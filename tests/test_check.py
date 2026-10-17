import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Settings, schemas
from cichlid_data.check import check
from cichlid_data.cli import main as cli
from cichlid_data.provenance import SourceRecord, inventory
from helpers import Env, write

SETS = '''[methods.TreeFam2026]
queries = "project"
[methods.OrthoFinder2026]
species = ["Mzeb", "Drer"]
[[organism]]
prefix = "Mzeb"
species = "Maylandia zebra"
role = "project"
[[organism]]
prefix = "Drer"
species = "Danio rerio"
role = "reference"
'''


def proteome(L, p, seq="MKLV"):
    write(L.proteome(p).local, f">{p}_1\n{seq}\n")
    schemas.GENES.write(L.genes_table(p).local, [{"gene_id": f"{p}_1", "seqid": "c1", "start": 1, "end": 9,
                                                   "mitochondrial": False, "chosen_protein": "XP_1"}])


def setup(root):
    cfg = Path(root) / "cfg"
    write(cfg / "species_sets.toml", SETS)
    write(cfg / "species.tsv", "prefix\tspecies\taccession\tsource\nMzeb\tMaylandia zebra\tGCF_1\tncbi\n")
    L = __import__("cichlid_data").Layout(Path(root) / "data", "fake:")
    tf = L.family_resource("TreeFam9").local
    write(tf / "hmm/TreeFam9", "hmm")
    write(tf / "trees/TF1.nhx", "(a,b);")
    write(tf / "protein_to_species.txt", "a\thuman\n")
    SourceRecord(name="TreeFam 9", files=inventory(tf)).write(tf)
    write(L.annotations("Legacy2018").local / "HumanEnsDescriptions.txt", "x")
    proteome(L, "Mzeb")
    return L, cfg


class TestCheck(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.L, self.cfg = setup(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def findings(self, **kw):
        return {(f.method, f.item): f for f in check(self.L, Catalog.build(self.L, self.cfg), **kw)}

    def test_complete_method_passes_and_problems_are_reported(self):
        f = self.findings()
        tf = {k: v for k, v in f.items() if k[0] == "TreeFam2026"}
        self.assertTrue(all(v.ok for v in tf.values()), [v for v in tf.values() if not v.ok])
        self.assertEqual(f[("TreeFam2026", "Inputs/Families/TreeFam9/hmm")].detail, "matches its source record")
        self.assertEqual(f[("TreeFam2026", "Inputs/Annotations/Legacy2018")].detail, "no source record")
        self.assertEqual((f[("OrthoFinder2026", "proteome Drer")].ok, f[("OrthoFinder2026", "proteome Drer")].detail),
                         (False, "missing"))
        # an invalid proteome, an empty folder, and a resource that no longer matches its record
        write(self.L.proteome("Mzeb").local, ">Mzeb_1\nMK*V\n")
        for t in (self.L.family_resource("TreeFam9").local / "trees").iterdir():
            t.unlink()
        (self.L.family_resource("TreeFam9").local / "protein_to_species.txt").write_text("changed contents")
        f = self.findings()
        self.assertIn("amino-acid alphabet", f[("TreeFam2026", "proteome Mzeb")].detail)
        self.assertEqual(f[("TreeFam2026", "Inputs/Families/TreeFam9/trees")].detail, "empty")
        self.assertIn("size changed: protein_to_species.txt",
                      f[("TreeFam2026", "Inputs/Families/TreeFam9/protein_to_species.txt")].detail)

    def test_deep_check_catches_same_size_changes(self):
        (self.L.family_resource("TreeFam9").local / "protein_to_species.txt").write_text("b\thuman\n")   # same size
        self.assertTrue(self.findings()[("TreeFam2026", "Inputs/Families/TreeFam9/protein_to_species.txt")].ok)
        self.assertFalse(self.findings(deep=True)[("TreeFam2026", "Inputs/Families/TreeFam9/protein_to_species.txt")].ok)

    def test_fetch_and_command(self):
        import os
        e = Env()
        try:
            L, cfg = setup(e.a.local_root.parent / "x")
            os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
            # Drer's proteome exists only in the cloud
            proteome(L, "Drer")
            from cichlid_data import LocalBackend, Store
            Store(L, LocalBackend(e.cloud)).push("proteins/proteomes")
            import shutil
            shutil.rmtree(L.proteomes().local)              # as on a machine that never had the proteomes
            shutil.rmtree(L.state_dir())
            s = Settings(local_root=L.local_root, cloud_remote="fake:")
            out = io.StringIO()
            with redirect_stdout(out):
                cli(["check", "--method", "OrthoFinder2026", "--fetch"], backend=LocalBackend(e.cloud), settings=s)
            self.assertIn("OrthoFinder2026: 2 of 2 inputs ready", out.getvalue())
            self.assertTrue(L.proteome("Drer").local.is_file())
            (L.family_resource("TreeFam9").local / "protein_to_species.txt").unlink()
            out = io.StringIO()
            with redirect_stdout(out), self.assertRaises(SystemExit):
                cli(["check"], backend=LocalBackend(e.cloud), settings=s)
            self.assertIn("Inputs/Families/TreeFam9/protein_to_species.txt: missing", out.getvalue())
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)
            e.close()


if __name__ == "__main__":
    unittest.main()


