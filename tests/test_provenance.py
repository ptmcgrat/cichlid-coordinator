import json
import tempfile
import unittest
from pathlib import Path

from cichlid_data import Layout
from cichlid_data.config import config_dir
from cichlid_data.provenance import SourceRecord, inventory, load_sources, read_record, verify
from cichlid_data.provenance_cmd import check_records, write_records
from helpers import write


class TestProvenance(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_stage1_record(self):
        d = self.root / "Mzeb"
        write(d / "protein.faa", "P")
        (d / "download.json").write_text(json.dumps({
            "manifest_row": {"prefix": "Mzeb", "species": "Maylandia zebra", "accession": "GCF_041146795.1",
                             "source": "ncbi", "notes": "lab reference"},
            "ncbi_metadata": {"assembly_name": "Mzebra_GT3a"},
            "files": {"protein": {"file": "protein.faa", "md5": "abc", "bytes": 1}},
            "downloaded": "2026-10-07T12:00:00+00:00"}))
        r = read_record(d)
        self.assertEqual((r.name, r.type, r.accession, r.version, r.retrieved),
                         ("Maylandia zebra", "ncbi_refseq", "GCF_041146795.1", "Mzebra_GT3a", "2026-10-07T12:00:00+00:00"))
        self.assertEqual(r.files, {"protein.faa": {"bytes": 1, "md5": "abc"}})
        self.assertIn("GCF_041146795.1", r.url)
        write(d / "assembly_report.json", "{}")                                    # Stage 1 bookkeeping
        self.assertEqual(verify(d, r), [])

    def test_other_download_json_is_not_a_record(self):
        L = Layout(self.root / "data", "fake:")
        hpa = L.expression("HPA").local
        write(hpa / "rna_single_cell_type.tsv.zip", "z")
        (hpa / "download.json").write_text(json.dumps({"downloaded": "2026-10-09T10:00:00Z", "files": 2}))
        blumer = L.population("Blumer2025").local
        write(blumer / "chr10.vcf.gz", "v")
        (blumer / "download.json").write_text(json.dumps({"files": ["chr10.vcf.gz"]}))
        self.assertIsNone(read_record(hpa))                                       # not mistaken for Stage 1's
        self.assertIsNone(read_record(blumer))
        sources = self.root / "sources.toml"
        sources.write_text('[[source]]\npath = "Inputs/Expression/HPA"\nname = "HPA"\n'
                           '[[source]]\npath = "Inputs/Populations/Blumer2025"\nname = "Blumer"\n')
        self.assertEqual(write_records(L, sources_file=sources, say=lambda *_: None), 2)
        rec = read_record(hpa)
        self.assertEqual((rec.retrieved, list(rec.files)), ("2026-10-09T10:00:00Z", ["rna_single_cell_type.tsv.zip"]))
        self.assertIn("download.json", rec.notes)
        self.assertEqual(check_records(L, sources_file=sources, say=lambda *_: None), 0)

    def test_inventory_and_verify(self):
        d = self.root / "res"
        write(d / "small.txt", "aaaa")
        write(d / "sub/big.bin", "x" * 2000)
        inv = inventory(d, md5_limit=1000)
        self.assertIsNotNone(inv["small.txt"]["md5"])
        self.assertEqual(inv["sub/big.bin"], {"bytes": 2000, "md5": None})          # over the limit: size only
        rec = SourceRecord(name="r", files=inv)
        self.assertEqual(verify(d, rec), [])
        (d / "small.txt").write_text("bbbb")                                          # same size, new content
        self.assertEqual(verify(d, rec), [])
        self.assertEqual(verify(d, rec, deep=True), ["content changed: small.txt"])
        (d / "sub/big.bin").write_text("x")
        write(d / "extra.txt", "e")
        (d / "small.txt").unlink()
        self.assertEqual(verify(d, rec), ["missing: small.txt", "size changed: sub/big.bin",
                                          "new file (not in record): extra.txt"])

    def test_write_and_check(self):
        L = Layout(self.root / "data", "fake:")
        write(L.family_resource("TreeFam9").local / "hmm/TreeFam9", "hmm")
        write(L.annotations("GO").local / "go-basic.obo", "go")
        write(L.annotations("Mystery").local / "x.txt", "?")
        sources = self.root / "sources.toml"
        sources.write_text('[[source]]\npath = "Inputs/Families/TreeFam9"\nname = "TreeFam 9"\ntype = "treefam"\nversion = "9"\n'
                           '[[source]]\npath = "Inputs/Annotations/GO"\nname = "GO"\nlicence = "CC BY 4.0"\n'
                           '[[source]]\npath = "Inputs/Expression/HPA"\nname = "HPA"\n')
        out = []
        self.assertEqual(write_records(L, sources_file=sources, say=out.append), 2)
        text = "\n".join(out)
        self.assertIn("Inputs/Expression/HPA: not on this machine", text)
        self.assertIn("Inputs/Annotations/Mystery: no source description", text)
        rec = read_record(L.family_resource("TreeFam9").local)
        self.assertEqual((rec.name, rec.version, list(rec.files)), ("TreeFam 9", "9", ["hmm/TreeFam9"]))
        self.assertTrue(rec.recorded_by.startswith("cichlid-data"))
        out = []
        self.assertEqual(write_records(L, sources_file=sources, say=out.append), 0)    # already recorded
        self.assertEqual(check_records(L, sources_file=sources, say=lambda *_: None), 0)
        (L.annotations("GO").local / "go-basic.obo").write_text("changed!")
        msgs = []
        self.assertEqual(check_records(L, sources_file=sources, say=msgs.append), 1)
        self.assertTrue(any("size changed: go-basic.obo" in m for m in msgs))

    def test_repository_sources_file(self):
        entries = load_sources(config_dir() / "sources.toml")
        self.assertGreaterEqual(len(entries), 14)
        for e in entries:
            self.assertTrue(e["path"] and e["name"], e)
            SourceRecord.from_dict({k: v for k, v in e.items() if k != "path"})      # every field is known


if __name__ == "__main__":
    unittest.main()


