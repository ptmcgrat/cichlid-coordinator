"""Offline tests for download_expression.py (fake HPA server, fake rclone)."""
import io
import json
import os
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


if __name__ == "__main__":
    unittest.main()


