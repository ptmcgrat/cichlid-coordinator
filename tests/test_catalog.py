import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Layout, Settings
from cichlid_data.cli import main as cli
from cichlid_data.config import config_dir
from helpers import Env, write

SETS = '''
[methods.TreeFam2026]
queries = "project"
[methods.OrthoFinder2026]
species = ["Mzeb", "Acal", "Esur", "Drer"]
[[organism]]
prefix = "Mzeb"
species = "Maylandia zebra"
role = "project"
source = "ncbi_refseq"
[[organism]]
prefix = "Acal"
species = "Astatotilapia calliptera"
role = "project"
source = "ncbi_refseq"
[[organism]]
prefix = "Esur"
species = "Etroplus suratensis"
role = "project"
source = "authors"
[[organism]]
prefix = "Drer"
species = "Danio rerio"
common = "zebrafish"
role = "reference"
source = "ncbi_refseq"
[[organism]]
prefix = "Amca"
species = "Amia calva"
role = "candidate"
'''
MANIFEST = "prefix\tspecies\taccession\tsource\nMzeb\tMaylandia zebra\tGCF_041146795.1\tncbi\nAcal\tAstatotilapia calliptera\tGCF_964374335.1\tncbi\n"
PANEL = "# panel\nspecies\trole\tgroup\nhomo_sapiens\treference\tmammal\nesox_lucius\tpanel\tteleost\nsalmo_salar\texclude\tfish\n"


def setup(root):
    cfg = root / "config"
    write(cfg / "species_sets.toml", SETS)
    write(cfg / "species.tsv", MANIFEST)
    write(cfg / "compara_panel.tsv", PANEL)
    L = Layout(root / "data", "fake:")
    write(L.proteome("Mzeb").local, ">Mzeb_1\nM\n>Mzeb_2\nK\n")
    write(L.genes_table("Mzeb").local, "gene_id\tsymbol\nMzeb_1\ta\nMzeb_2\tb\n")
    write(L.download("Acal").local / "protein.faa", ">x\nM\n")
    (L.download("Acal").local / "download.json").write_text(json.dumps({
        "manifest_row": {"prefix": "Acal", "species": "Astatotilapia calliptera", "accession": "GCF_964374335.1",
                         "source": "ncbi"}, "ncbi_metadata": {"assembly_name": "fAstCal68.1"},
        "files": {"protein": {"file": "protein.faa", "md5": "x", "bytes": 5}}, "downloaded": "2026-10-07T00:00:00Z"}))
    write(L.download("Acal").local / "assembly_report.json", "{}")
    write(L.local_source("Esur").local / "Genes.fasta.gz", "g")
    return L, cfg


class TestCatalog(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L, self.cfg = setup(self.root)

    def tearDown(self):
        self._t.cleanup()

    def test_build_and_query(self):
        C = Catalog.build(self.L, self.cfg)
        m, a, e, d = C["Mzeb"], C["Acal"], C["Esur"], C["Drer"]
        self.assertEqual((m.status, m.counts, m.accession, m.in_manifest),
                         ("proteome ready", {"proteins": 2, "genes": 2}, "GCF_041146795.1", True))
        self.assertEqual(m.proteome, "Inputs/Proteins/proteomes/Mzeb.faa")
        self.assertEqual((a.status, a.version, a.retrieved), ("downloaded", "fAstCal68.1", "2026-10-07T00:00:00Z"))
        self.assertEqual((e.status, e.in_manifest), ("local source only", False))
        self.assertEqual(d.status, "not downloaded")
        self.assertEqual(C.methods["TreeFam2026"]["species"], ["Mzeb", "Acal", "Esur"])      # the project species
        self.assertEqual(m.used_by, ["TreeFam2026", "OrthoFinder2026"])
        self.assertEqual(d.used_by, ["OrthoFinder2026"])
        self.assertEqual([s.prefix for s in C.missing("OrthoFinder2026")], ["Acal", "Esur", "Drer"])
        self.assertEqual([s.prefix for s in C.species(role="reference")], ["Drer"])
        self.assertEqual([s["species"] for s in C.panel], ["homo_sapiens", "esox_lucius"])   # excluded rows dropped
        errors, warnings = C.validate()
        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["Esur: a project species not in species.tsv, so Stage 1 won't download it"])
        # write and load
        C.write(self.L)
        C2 = Catalog.load(self.L)
        self.assertEqual(C2["Mzeb"].to_dict(), m.to_dict())
        self.assertEqual(C2.methods, C.methods)
        self.assertIn("species_sets.toml", C2.meta["config"])

    def test_validation_errors(self):
        (self.cfg / "species_sets.toml").write_text(SETS.replace('species = ["Mzeb", "Acal", "Esur", "Drer"]',
                                                                   'species = ["Mzeb", "Amca", "Nope"]')
                                                    + '[[organism]]\nprefix = "Bad_1"\nspecies = "x"\nrole = "pet"\n')
        errors, _ = Catalog.build(self.L, self.cfg).validate()
        joined = "\n".join(errors)
        for bit in ("Amca is a candidate", "Nope has no organism entry", "Bad_1: a prefix must be letters",
                    "unknown role 'pet'"):
            self.assertIn(bit, joined)

    def test_repository_configuration(self):
        C = Catalog.build(self.L, config_dir())
        errors, warnings = C.validate()
        self.assertEqual(errors, [])
        self.assertEqual(len(C.methods["OrthoFinder2026"]["species"]), 35)
        self.assertEqual(len(C.species(role="project")), 12)
        self.assertEqual(len(C.methods["TreeFam2026"]["species"]), 12)
        self.assertGreater(len(C.panel), 50)


class TestCatalogCli(unittest.TestCase):
    def test_build_and_missing(self):
        e = Env()
        try:
            import os
            L, cfg = setup(e.a.local_root.parent / "x")
            os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
            s = Settings(local_root=L.local_root, cloud_remote="fake:")
            f = io.StringIO()
            with redirect_stdout(f):
                cli(["catalog", "build"], backend=e.store(e.a).backend, settings=s)
                cli(["catalog", "missing", "--method", "OrthoFinder2026"], backend=e.store(e.a).backend, settings=s)
            out = f.getvalue()
            self.assertIn("5 organisms (1 with proteomes ready)", out)
            self.assertIn("1 of 4 proteomes ready", out)
            self.assertIn("Drer", out)
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)
            e.close()


if __name__ == "__main__":
    unittest.main()


