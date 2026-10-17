import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Settings, schemas
from cichlid_data.check import check
from cichlid_data.cli import main as cli
from cichlid_data.provenance import SourceRecord, inventory
from helpers import Env, releases_toml, write

RELEASES = releases_toml({"TreeFam": 'queries = "project"', "OrthoFinder": 'species = ["Mzeb", "Drer"]',
                          "Homology": 'description = "calls"'})
SETS = '''[[organism]]
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
    write(cfg / "releases.toml", RELEASES)
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
        tf = {k: v for k, v in f.items() if k[0] == "v1/TreeFam"}
        self.assertTrue(all(v.ok for v in tf.values()), [v for v in tf.values() if not v.ok])
        self.assertEqual(f[("v1/TreeFam", "Inputs/Families/TreeFam9/hmm")].detail, "matches its source record")
        # the 2018 name tables are the Homology part's requirement, not TreeFam's
        self.assertEqual(f[("v1/Homology", "Inputs/Annotations/Legacy2018")].detail, "no source record")
        self.assertNotIn(("v1/TreeFam", "Inputs/Annotations/Legacy2018"), f)
        self.assertEqual((f[("v1/OrthoFinder", "proteome Drer")].ok, f[("v1/OrthoFinder", "proteome Drer")].detail),
                         (False, "missing"))
        # an invalid proteome, an empty folder, and a resource that no longer matches its record
        write(self.L.proteome("Mzeb").local, ">Mzeb_1\nMK*V\n")
        for t in (self.L.family_resource("TreeFam9").local / "trees").iterdir():
            t.unlink()
        (self.L.family_resource("TreeFam9").local / "protein_to_species.txt").write_text("changed contents")
        f = self.findings()
        self.assertIn("amino-acid alphabet", f[("v1/TreeFam", "proteome Mzeb")].detail)
        self.assertEqual(f[("v1/TreeFam", "Inputs/Families/TreeFam9/trees")].detail, "empty")
        self.assertIn("size changed: protein_to_species.txt",
                      f[("v1/TreeFam", "Inputs/Families/TreeFam9/protein_to_species.txt")].detail)

    def test_deep_check_catches_same_size_changes(self):
        (self.L.family_resource("TreeFam9").local / "protein_to_species.txt").write_text("b\thuman\n")   # same size
        self.assertTrue(self.findings()[("v1/TreeFam", "Inputs/Families/TreeFam9/protein_to_species.txt")].ok)
        self.assertFalse(self.findings(deep=True)[("v1/TreeFam", "Inputs/Families/TreeFam9/protein_to_species.txt")].ok)

    def test_release_and_part_selection(self):
        C = Catalog.build(self.L, self.cfg)
        self.assertEqual({f.method for f in check(self.L, C)}, {"v1/TreeFam", "v1/OrthoFinder", "v1/Homology"})
        self.assertEqual({f.method for f in check(self.L, C, release="v1", parts=["Homology"])}, {"v1/Homology"})
        self.assertEqual([(f.method, f.item, f.ok) for f in check(self.L, C, release="v1", parts=["Compara"])],
                         [("v1/Compara", "part", False)])
        with self.assertRaises(KeyError):
            check(self.L, C, release="v9")

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
                cli(["check", "--release", "v1", "--part", "OrthoFinder", "--fetch"], backend=LocalBackend(e.cloud), settings=s)
            self.assertIn("v1/OrthoFinder: 2 of 2 inputs ready", out.getvalue())
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


