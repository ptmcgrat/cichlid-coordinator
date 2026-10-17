"""Stage 1's tests (same fakes, same expectations), through the layout, plus the coordinator's additions."""
import hashlib
import io
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

from cichlid_data.downloaders import ncbi as s1
from cichlid_data import Layout
from cichlid_data.config import config_dir

ACC = "GCF_001858045.2"
GFF = b"##gff-version 3\nNC_031965.2\tGnomon\tgene\t1\t100\t.\t+\t.\tID=gene-a\n"
FAA = b">XP_1.1 protein A\nMKV\n>XP_2.1 protein B\nMKL\n"
SEQ = b'{"refseqAccession":"NC_027944.1","assignedMoleculeLocationType":"Mitochondrion"}\n'

REPORT = {"accession": ACC,
          "organism": {"organism_name": "Oreochromis niloticus", "tax_id": 8128},
          "assembly_info": {"assembly_name": "O_niloticus_UMD_NMBU", "assembly_status": "current",
                            "assembly_level": "Chromosome"},
          "assembly_stats": {"contig_n50": 2923640},
          "annotation_info": {"name": "GCF_001858045.2-RS_2025_06", "provider": "NCBI RefSeq",
                              "release_date": "2025-06-16",
                              "busco": {"busco_lineage": "actinopterygii_odb10",
                                        "complete": 0.994, "duplicated": 0.012},
                              "stats": {"gene_counts": {"protein_coding": 29562}}}}


def fake_zip(corrupt=False):
    files = {f"ncbi_dataset/data/{ACC}/genomic.gff": GFF,
             f"ncbi_dataset/data/{ACC}/protein.faa": FAA,
             f"ncbi_dataset/data/{ACC}/sequence_report.jsonl": SEQ,
             "ncbi_dataset/data/assembly_data_report.jsonl": b'{"accession":"%s"}\n' % ACC.encode()}
    catalog = {"apiVersion": "V2", "assemblies": [
        {"files": [{"filePath": "assembly_data_report.jsonl", "fileType": "DATA_REPORT"}]},
        {"accession": ACC, "files": [
            {"filePath": f"{ACC}/genomic.gff", "fileType": "GFF3"},
            {"filePath": f"{ACC}/protein.faa", "fileType": "PROTEIN_FASTA"},
            {"filePath": f"{ACC}/sequence_report.jsonl", "fileType": "SEQUENCE_REPORT"}]}]}
    md5 = "".join(f"{hashlib.md5(b).hexdigest() if not (corrupt and k.endswith('faa')) else '0'*32}  {k}\n"
                  for k, b in files.items())
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for k, b in files.items():
            z.writestr(k, b)
        z.writestr("ncbi_dataset/data/dataset_catalog.json", json.dumps(catalog))
        z.writestr("md5sum.txt", md5)
    return buf.getvalue()


class FakeResponse:
    def __init__(self, status=200, js=None, content=b""):
        self.status_code, self._js, self.content = status, js, content

    def json(self):
        return self._js

    def iter_content(self, n):
        yield self.content

    def raise_for_status(self):
        raise s1.requests.HTTPError(str(self.status_code))

    def __enter__(self):
        return self

    def __exit__(self, *a):
        pass


class FakeSession:
    def __init__(self, report=REPORT, zipbytes=None, fail_first=0):
        self.headers, self.report, self.zip = {}, report, zipbytes or fake_zip()
        self.fail_first, self.calls = fail_first, []

    def get(self, url, **kw):
        self.calls.append((url, kw.get("params")))
        if self.fail_first:
            self.fail_first -= 1
            return FakeResponse(503)
        if url.endswith("dataset_report"):
            return FakeResponse(js={"reports": [self.report] if self.report else []})
        return FakeResponse(content=self.zip)


def write_manifest(d, body):
    p = Path(d) / "species.tsv"
    p.write_text("# comment line\nprefix\tspecies\taccession\tsource\tannotation_format\n" + body)
    return p


class TestManifest(unittest.TestCase):
    def test_example_manifest_parses(self):
        import csv
        with open(config_dir() / "species.tsv") as fh:          # valid, whatever species it lists
            ncbi_rows = [r for r in csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t")
                         if r["source"] == "ncbi"]
        self.assertGreaterEqual(len(ncbi_rows), 9)
        with tempfile.TemporaryDirectory() as d:                # NCBI rows alone (local files are machine-specific)
            p = Path(d) / "species.tsv"
            p.write_text("prefix\tspecies\taccession\tsource\tannotation_format\n" +
                         "".join(f"{r['prefix']}\t{r['species']}\t{r['accession']}\tncbi\t{r.get('annotation_format') or 'gff3'}\n"
                                 for r in ncbi_rows))
            self.assertEqual(len(s1.read_manifest(p)), len(ncbi_rows))

    def test_all_errors_reported_together(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_manifest(d, "On_il\tX y\tGCF_1\tncbi\tgff3\n"
                                  "Mzeb\tM z\tGCF_041146795.1\tncbi\tbed\n"
                                  "mzeb\tM z\tGCF_041146795.1\tftp\tgff3\n")
            with self.assertRaises(s1.ManifestError) as cm:
                s1.read_manifest(p)
            msg = str(cm.exception)
            for needle in ("no underscores", "not like GCF", "gff3 or gtf", "duplicate prefix", "ncbi or local"):
                self.assertIn(needle, msg)

    def test_local_requires_existing_files(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_manifest(d, "Brk\tA b\t\tlocal\tgff3\n")
            with self.assertRaises(s1.ManifestError):
                s1.read_manifest(p)


class TestMetadata(unittest.TestCase):
    def test_summary_and_checks(self):
        s = s1.summarize_report(REPORT)
        self.assertEqual(s["protein_coding_genes"], 29562)
        row = {"species": "Oreochromis niloticus", "accession": ACC}
        self.assertEqual(s1.check_report(row, s), ([], []))

    def test_unannotated_is_error_and_mismatch_is_warning(self):
        rep = json.loads(json.dumps(REPORT))
        del rep["annotation_info"]
        rep["assembly_info"]["assembly_status"] = "suppressed"
        errs, warns = s1.check_report({"species": "Tilapia sp.", "accession": ACC},
                                      s1.summarize_report(rep))
        self.assertEqual(len(errs), 1)
        self.assertEqual(len(warns), 2)


class TestEndToEnd(unittest.TestCase):
    def run_main(self, d, session, *extra):
        """Stage 1's main, now download(): returns the downloads folder (the old --outdir)."""
        p = write_manifest(d, f"Onil\tOreochromis niloticus\t{ACC}\tncbi\tgff3\n")
        L = Layout(Path(d) / "data", "fake:")
        with mock.patch.object(s1.time, "sleep"):
            _, failed = s1.download(L, p, session=session, check_only="--check-only" in extra,
                                    force="--force" in extra, say=lambda *_: None)
        if failed:
            raise SystemExit(f"{len(failed)} species failed")
        return L.path("Inputs/Proteins/downloads").local

    def test_download_extract_and_skip(self):
        with tempfile.TemporaryDirectory() as d:
            sess = FakeSession(fail_first=1)  # one 503 then success: exercises retry
            out = self.run_main(d, sess)
            sp = out / "Onil"
            self.assertEqual((sp / "genomic.gff").read_bytes(), GFF)
            self.assertEqual((sp / "protein.faa").read_bytes(), FAA)
            self.assertTrue((sp / "sequence_report.jsonl").is_file())
            self.assertTrue((sp / "assembly_report.json").is_file())
            self.assertFalse(list(sp.glob("*.zip")))
            prov = json.loads((sp / "download.json").read_text())
            self.assertEqual(prov["protein_records"], 2)
            self.assertEqual(prov["ncbi_metadata"]["annotation_name"], "GCF_001858045.2-RS_2025_06")
            params = [c[1] for c in sess.calls if c[0].endswith("download")][0]
            self.assertIn(("include_annotation_type", "GENOME_GFF"), params)
            self.assertIn(("include_annotation_type", "PROT_FASTA"), params)
            self.assertIn("downloaded", (out / "download_summary.tsv").read_text())

            sess2 = FakeSession()
            self.run_main(d, sess2)
            self.assertFalse([c for c in sess2.calls if c[0].endswith("download")], "should skip")
            self.assertIn("up-to-date", (out / "download_summary.tsv").read_text())

            (sp / "protein.faa").write_text(">tampered\n")   # damaged file -> refetch
            sess3 = FakeSession()
            self.run_main(d, sess3)
            self.assertTrue([c for c in sess3.calls if c[0].endswith("download")])

    def test_check_only_reports_each_species(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_manifest(d, f"Onil\tOreochromis niloticus\t{ACC}\tncbi\tgff3\n")
            msgs = []
            with mock.patch.object(s1.time, "sleep"):
                s1.download(Layout(Path(d) / "data", "fake:"), p, session=FakeSession(), check_only=True, say=msgs.append)
            self.assertIn(f"[Onil] {ACC}: O_niloticus_UMD_NMBU (Chromosome), GCF_001858045.2-RS_2025_06, "
                          "29562 protein-coding genes", msgs)
            self.assertTrue(msgs[-1].startswith("1 species: 1 checked; wrote "))

    def test_check_only_downloads_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            sess = FakeSession()
            out = self.run_main(d, sess, "--check-only")
            self.assertFalse((out / "Onil").exists())
            self.assertIn("29562", (out / "check_summary.tsv").read_text())

    def test_md5_mismatch_fails_species(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                self.run_main(d, FakeSession(zipbytes=fake_zip(corrupt=True)))
            self.assertIn("md5 mismatch", (Path(d) / "data/Inputs/Proteins/downloads/download_summary.tsv").read_text())

    def test_missing_accession_fails_cleanly(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                self.run_main(d, FakeSession(report=None))
            self.assertIn("no record", (Path(d) / "data/Inputs/Proteins/downloads/download_summary.tsv").read_text())


if __name__ == "__main__":
    unittest.main()


class TestCoordinatorAdditions(unittest.TestCase):
    def test_pulls_from_cloud_instead_of_ncbi(self):
        from helpers import Env
        e = Env()
        try:
            man = write_manifest(e.a.local_root.parent, f"Onil\tOreochromis niloticus\t{ACC}\tncbi\tgff3\n")
            with mock.patch.object(s1.time, "sleep"):
                s1.download(e.a, man, session=FakeSession(), say=lambda *_: None)
            e.store(e.a).push("proteins/downloads")
            sess = FakeSession()                                  # machine b: in the cloud, not local
            out = []
            with mock.patch.object(s1.time, "sleep"):
                summary, failed = s1.download(e.b, man, session=sess, store=e.store(e.b), say=out.append)
            self.assertEqual((failed, summary[0]["status"]), ([], "up-to-date"))
            self.assertIn("[Onil] pulled from the cloud", out)
            self.assertFalse([c for c in sess.calls if c[0].endswith("download")], "nothing fetched from NCBI")
            self.assertEqual((e.b.download("Onil").local / "protein.faa").read_bytes(), FAA)
        finally:
            e.close()

    def test_summary_keeps_other_species_and_relative_local_paths(self):
        with tempfile.TemporaryDirectory() as d:
            L = Layout(Path(d) / "data", "fake:")
            src = L.local_source("Brk").local / "prepared"
            src.mkdir(parents=True)
            (src / "genomic.gff").write_bytes(GFF)
            (src / "protein.faa").write_bytes(FAA)
            p = Path(d) / "species.tsv"
            p.write_text("prefix\tspecies\taccession\tsource\tannotation_format\tgroup\tgff_path\tprotein_path\n"
                         f"Onil\tOreochromis niloticus\t{ACC}\tncbi\tgff3\t\t\t\n"
                         "Brk\tEtroplus suratensis\t\tlocal\tgff3\t\t"
                         "Inputs/Proteins/local_sources/Brk/prepared/genomic.gff\t"
                         "Inputs/Proteins/local_sources/Brk/prepared/protein.faa\n")
            with mock.patch.object(s1.time, "sleep"):
                s1.download(L, p, session=FakeSession(), say=lambda *_: None)
                s1.download(L, p, species=["Brk"], force=True, session=FakeSession(), say=lambda *_: None)
            self.assertEqual((L.download("Brk").local / "genomic.gff").read_bytes(), GFF)   # relative path resolved
            text = (L.path("Inputs/Proteins/downloads").local / "download_summary.tsv").read_text()
            self.assertIn("Onil", text)                                                     # kept
            self.assertIn("Brk", text)
            with self.assertRaises(s1.ManifestError):
                s1.download(L, p, species=["Nope"], session=FakeSession(), say=lambda *_: None)


class TestCommands(unittest.TestCase):
    def test_download_genomes_command(self):
        import io
        import os
        from contextlib import redirect_stdout
        from cichlid_data import Settings
        from cichlid_data.cli import main as cli
        from helpers import Env, releases_toml, write
        e = Env()
        try:
            cfg = e.a.local_root.parent / "cfg"
            write(cfg / "releases.toml", releases_toml({"OrthoFinder": 'species = ["Onil"]'}, deliver=()))
            write(cfg / "species.tsv", f"prefix\tspecies\taccession\tsource\tannotation_format\nOnil\tOreochromis niloticus\t{ACC}\tncbi\tgff3\n")
            write(cfg / "species_sets.toml", '[[organism]]\nprefix = "Onil"\nspecies = "Oreochromis niloticus"\nrole = "project"\n')
            os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
            s = Settings(local_root=e.a.local_root, cloud_remote="fake:")
            f = io.StringIO()
            with redirect_stdout(f), mock.patch.object(s1, "make_session", return_value=FakeSession()), \
                    mock.patch.object(s1.time, "sleep"):
                cli(["download", "genomes", "--push"], backend=e.store(e.a).backend, settings=s)
            out = f.getvalue()
            self.assertIn("Catalog updated: 1 species downloaded", out)
            self.assertIn("pushed proteins/downloads: uploaded", out)
            self.assertTrue((e.cloud / "Inputs/Proteins/downloads/Onil.tar").is_file())
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)
            e.close()

    def test_prepare_local_default_outdir(self):
        from cichlid_data import Settings
        from cichlid_data.cli import main as cli
        from cichlid_data.downloaders import local
        from helpers import Env
        e = Env()
        try:
            got = {}
            with mock.patch.object(local, "main", side_effect=lambda argv: got.setdefault("argv", argv)):
                cli(["prepare-local", "--prefix", "Esur", "--cds", "x.fa"], backend=e.store(e.a).backend,
                    settings=Settings(local_root=e.a.local_root, cloud_remote="fake:"))
            self.assertEqual(got["argv"][-2:], ["--outdir", str(e.a.local_source("Esur").local / "prepared")])
        finally:
            e.close()


