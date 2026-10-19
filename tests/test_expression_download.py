"""Offline tests for download_expression.py (fake HPA server, fake rclone)."""
import contextlib
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cichlid_data.downloaders import expression as de
from helpers import coord_env
FAKE_RCLONE = ""        # unused: backups go through the coordinator's stand-in rclone (coord_env)


def zbytes(name, text):
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w") as z:
        z.writestr(name, text)
    return b.getvalue()


class R:
    def __init__(self, status, body=b""):
        self.status_code, self.body = status, body
    def iter_content(self, n):
        yield self.body
    def __enter__(self):
        return self
    def __exit__(self, *a):
        pass


class FakeHPA:
    """Files only exist at the second (older) location."""
    def __init__(self):
        self.urls = []
    def get(self, url, stream=False, timeout=None):
        self.urls.append(url)
        if "/download/tsv/" in url:
            return R(404)
        name = url.rsplit("/", 1)[-1]
        return R(200, zbytes(name.replace(".zip", ""), "Gene\tCell type\tnTPM\nENSG1\thepatocytes\t1\n"))


BGEE_HEAD = ("Gene ID\tGene name\tAnatomical entity ID\tAnatomical entity name\tDevelopmental stage ID\t"
             "Developmental stage name\tSex\tStrain\tExpression\tCall quality\tFDR\tExpression score\tExpression rank\n")
BGEE_ROWS = [
    ("ENSDARG1", "rho", "UBERON:0000966", "retina", "ZFS:0000044", "adult (Danio)", "female", "wild-type", "present", "gold quality", "0", "95.5", "10"),
    ("ENSDARG1", "rho", "UBERON:0000966", "retina", "ZFS:0000045", "post-juvenile adult stage", "male", "wild-type", "present", "silver quality", "0", "97.0", "8"),
    ("ENSDARG1", "rho", "UBERON:0000955", "brain", "ZFS:0000044", "adult (Danio)", "any", "wild-type", "absent", "gold quality", "0", "10", "99"),
    ("ENSDARG1", "rho", "UBERON:0002107", "liver", "ZFS:0000035", "larval stage", "any", "wild-type", "present", "gold quality", "0", "80", "20"),
    ("ENSDARG2", "vtg1", "UBERON:0002107", "liver", "ZFS:0000044", "adult (Danio)", "female", "wild-type", "present", "gold quality", "0", "99", "1"),
    # Bgee's actual naming of the zebrafish adult stage
    ("ENSDARG3", "myh7", "UBERON:0000948", "heart", "UBERON:0000113", "post-juvenile", "any", "wild-type", "present", "gold quality", "0", "91", "3"),
    ("ENSDARG3", "myh7", "UBERON:0002107", "liver", "ZFS:0000040", "Juvenile:Days 30-44 (Danio)", "any", "wild-type", "present", "gold quality", "0", "50", "60"),
]


def bgee_gz():
    import gzip as _g
    return _g.compress((BGEE_HEAD + "".join("\t".join(r) + "\n" for r in BGEE_ROWS)).encode())


class TestBgee(unittest.TestCase):
    def test_adult_stage_names(self):
        self.assertTrue(de.is_adult("post-juvenile"))
        self.assertTrue(de.is_adult("anything", "UBERON:0000113"))
        self.assertTrue(de.is_adult("adult (Danio)"))
        self.assertFalse(de.is_adult("Juvenile:Days 30-44 (Danio)"))
        self.assertFalse(de.is_adult("Larval:Day 5 (Danio)"))

    def test_adult_summary(self):
        import gzip as _g
        with tempfile.TemporaryDirectory() as d:
            src, dest = Path(d) / "calls.tsv.gz", Path(d) / "adult.tsv.gz"
            src.write_bytes(bgee_gz())
            genes, rows = de.summarize_bgee(src, dest)
            got = {(r.split("\t")[0], r.split("\t")[3]): r.rstrip("\n").split("\t") for r in
                   _g.decompress(dest.read_bytes()).decode().splitlines()[1:]}
            self.assertEqual((genes, rows), (3, 3))                      # absent, larval and juvenile calls dropped
            self.assertEqual(got[("ENSDARG3", "heart")][4], "91.00")     # 'post-juvenile' is the adult stage
            self.assertNotIn(("ENSDARG3", "liver"), got)
            r = got[("ENSDARG1", "retina")]
            self.assertEqual((r[4], r[5], r[6]), ("97.00", "gold", "female,male"))   # max score, best quality, sexes
            self.assertEqual(got[("ENSDARG2", "liver")][6], "female")

    def test_download_and_backup(self):
        class S:
            def get(self, url, stream=False, timeout=None):
                return R(200, bgee_gz())
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "cloud").mkdir()
            (d / "rclone").write_text(FAKE_RCLONE)
            (d / "rclone").chmod(0o755)
            env = coord_env(d)
            with mock.patch.dict(os.environ, env):
                de.main(["--source", "bgee"], session=S())
            out = d / "local" / "Inputs" / "Expression" / "Bgee"
            self.assertEqual(json.loads((out / "download.json").read_text())["adult_rows"], 3)
            class Never:
                def get(self, *a, **k):
                    raise AssertionError("should reuse the existing download")
            with mock.patch.dict(os.environ, env):
                de.main(["--source", "bgee", "--no-backup"], session=Never())        # summary only
            self.assertTrue((d / "cloud" / "Inputs/Expression/Bgee/bgee_adult.tsv.gz").is_file())


class TestExpressionDownload(unittest.TestCase):
    def test_fallback_record_and_backup(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "cloud").mkdir()
            (d / "rclone").write_text(FAKE_RCLONE)
            (d / "rclone").chmod(0o755)
            env = coord_env(d)
            hpa = FakeHPA()
            with mock.patch.dict(os.environ, env), mock.patch.object(de.time, "sleep"):
                de.main([], session=hpa)
            out = d / "local" / "Inputs" / "Expression" / "HPA"
            rec = json.loads((out / "download.json").read_text())
            self.assertEqual(sorted(rec["files"]), ["proteinatlas.tsv.zip", "rna_single_cell_type.tsv.zip"])
            self.assertTrue(rec["files"]["proteinatlas.tsv.zip"]["url"].endswith("/download/proteinatlas.tsv.zip"))
            self.assertTrue(zipfile.is_zipfile(out / "rna_single_cell_type.tsv.zip"))
            cloud = d / "cloud" / "Inputs/Expression/HPA"
            self.assertTrue((cloud / "proteinatlas.tsv.zip").is_file())
            self.assertIn("Karlsson", (cloud / "DATA_NOTICE.txt").read_text())

    def test_reports_when_nothing_works(self):
        class Down:
            def get(self, url, stream=False, timeout=None):
                return R(404)
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(IOError) as cm:
                de.fetch(Down(), ["https://a/x.zip", "https://b/x.zip"], Path(d) / "x.zip")
            self.assertIn("could not download", str(cm.exception))



class TestDaniocell(unittest.TestCase):
    FILES = {"Daniocell2023_SeuratV4.rds": b"rds", "cluster_annotations.csv": b"clust,tissue\n", "daniocell_load.R": b"# R"}

    class Session:
        def get(self, url, stream=False, timeout=None):
            return R(200, TestDaniocell.FILES[url.rsplit("/", 1)[-1]])

    def test_download_then_the_packaged_export_step(self):
        class Never:
            def get(self, *a, **k):
                raise AssertionError("should reuse the existing files")
        with tempfile.TemporaryDirectory() as d:
            env = coord_env(Path(d))
            err = io.StringIO()
            with mock.patch.dict(os.environ, env), contextlib.redirect_stderr(err):
                de.main(["--source", "daniocell"], session=self.Session())
            out = Path(d) / "local" / "Inputs" / "Expression" / "raw" / "daniocell"
            self.assertEqual(sorted(p.name for p in out.iterdir()), sorted([*self.FILES, "DATA_NOTICE.txt"]))
            self.assertIn("Next (in an environment with R and Seurat): Rscript " + str(de.EXPORT_SEURAT), err.getvalue())
            self.assertIn(f"{out / 'Daniocell2023_SeuratV4.rds'} {out / 'export'}", err.getvalue())
            err = io.StringIO()
            with mock.patch.dict(os.environ, env), contextlib.redirect_stderr(err):
                de.main(["--source", "daniocell"], session=Never())
            self.assertEqual(err.getvalue().count("already here"), 3)


class TestExportSeuratScript(unittest.TestCase):
    """export_seurat.R was moved verbatim from the pipeline; the Python reader depends on what it writes."""
    SCRIPT = de.EXPORT_SEURAT
    MD5 = "c6a5117bf2ab964d7af93c533215b119"          # of the script as moved; a deliberate change to it means changing this too

    def test_the_script_ships_with_the_package_and_is_unchanged(self):
        self.assertEqual(self.SCRIPT.parent.name, "downloaders")
        self.assertTrue(self.SCRIPT.is_file())
        self.assertEqual(hashlib.md5(self.SCRIPT.read_bytes()).hexdigest(), self.MD5)
        self.assertTrue(self.SCRIPT.read_text().startswith("#!/usr/bin/env Rscript\n"))

    def test_it_writes_exactly_the_files_summarize_h5ad_reads(self):
        written = set(re.findall(r'out\("([^"]+)"\)', self.SCRIPT.read_text()))
        from cichlid_data.downloaders import h5ad
        read = set(re.findall(r'd / "([a-z_]+\.(?:bin|txt|tsv))"', Path(h5ad.__file__).read_text()))
        self.assertEqual(written, {"counts_x.bin", "counts_i.bin", "counts_p.bin", "dims.txt", "genes.txt",
                                    "cells.txt", "meta.tsv", "source.txt"})
        self.assertEqual(written, read)

    def test_pyproject_ships_r_files(self):
        pyproject = self.SCRIPT.parents[3] / "pyproject.toml"
        if not pyproject.is_file():
            self.skipTest("pyproject.toml not present")
        self.assertIn('"cichlid_data.downloaders" = ["*.R"]', pyproject.read_text())


if __name__ == "__main__":
    unittest.main()
