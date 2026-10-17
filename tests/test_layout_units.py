import os
import tempfile
import unittest
from pathlib import Path

from cichlid_data import ConfigError, Layout, Registry, Settings
from cichlid_data.config import write_default_config
from helpers import write


class TestLayout(unittest.TestCase):
    def setUp(self):
        self.L = Layout("/data/root", "remote:/base/")

    def test_paths(self):
        L = self.L
        self.assertEqual(L.proteome("Mzeb").rel, "Inputs/Proteins/proteomes/Mzeb.faa")
        self.assertEqual(str(L.proteome("Mzeb").local), "/data/root/Inputs/Proteins/proteomes/Mzeb.faa")
        self.assertEqual(L.proteome("Mzeb").cloud, "remote:/base/Inputs/Proteins/proteomes/Mzeb.faa")
        self.assertEqual(os.fspath(L.genes_table("Mzeb")), "/data/root/Inputs/Proteins/proteomes/Mzeb.genes.tsv")
        self.assertEqual(L.download("Acal").rel, "Inputs/Proteins/downloads/Acal")
        self.assertEqual(L.family_resource("TreeFam9").rel, "Inputs/Families/TreeFam9")
        R = L.results("v1", "Compara")
        self.assertEqual(R.assign("Mzeb").rel, "Results/v1/Compara/assign/Mzeb.assign.tsv")
        self.assertEqual(R.orthology("Mzeb").rel, "Results/v1/Compara/orthology/Mzeb.orthologs.tsv")
        self.assertEqual(R.family("TF101001").rel, "Results/v1/Compara/families/TF101001")
        self.assertEqual(R.run_record().rel, "Results/v1/Compara/run.json")
        self.assertEqual(R.root().rel, "Results/v1/Compara")
        self.assertEqual((R.release, R.part), ("v1", "Compara"))
        self.assertEqual((L.analysis("Tissues") / "summary.tsv").rel, "Analyses/Tissues/summary.tsv")
        with self.assertRaises(KeyError):
            R.section("bogus")
        with self.assertRaises(KeyError):
            L.results("v1", "TreeFam2026")                               # a part, not a method-and-year
        with self.assertRaises(ValueError):
            L.results("2026", "TreeFam")                                 # a release is v1, v2, ...
        self.assertEqual(L.release("v1").manifest().rel, "Releases/v1/release.json")
        self.assertEqual(L.release("v2").proteome("Mzeb").rel, "Releases/v2/proteomes/Mzeb.faa")
        self.assertEqual(L.release("v2").genes_table("Mzeb").rel, "Releases/v2/genes/Mzeb.genes.tsv")
        self.assertEqual(L.release("v2").catalog().rel, "Releases/v2/catalog.json")
        self.assertEqual(L.release("v2").species_tree("t.nwk").rel, "Releases/v2/config/t.nwk")
        with self.assertRaises(ValueError):
            L.release("latest")


class TestRunsAndReleases(unittest.TestCase):
    def test_what_is_on_this_machine(self):
        with tempfile.TemporaryDirectory() as d:
            L = Layout(d, "fake:")
            self.assertEqual((L.runs(), L.releases()), ([], []))
            write(Path(d) / "Results/v1/TreeFam/assign/x.tsv")
            write(Path(d) / "Results/v1/OrthoFinder/run.json")
            write(Path(d) / "Results/v2/Compara/run.json")
            write(Path(d) / "Results/v1/Notes/x.txt")                       # not a part
            write(Path(d) / "Releases/v3/release.json")                      # frozen archive, no results yet
            write(Path(d) / "Releases/latest/release.json")                  # not a release name
            self.assertEqual(L.runs(), [("v1", "OrthoFinder"), ("v1", "TreeFam"), ("v2", "Compara")])
            self.assertEqual(L.releases(), ["v1", "v2", "v3"])


class TestSettings(unittest.TestCase):
    def test_load_and_errors(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config.toml"
            old = {k: os.environ.pop(k, None) for k in ("CICHLID_DATA_CONFIG", "CICHLID_DATA_ROOT", "CICHLID_DATA_REMOTE")}
            try:
                os.environ["CICHLID_DATA_CONFIG"] = str(cfg)
                with self.assertRaises(ConfigError) as e:
                    Settings.load()
                self.assertIn("cichlid-data init", str(e.exception))
                write_default_config(cfg, local_root=f"{d}/root")
                s = Settings.load()
                self.assertEqual((s.local_root, s.cloud_remote), (Path(f"{d}/root"), "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/TreeFamData"))
                os.environ["CICHLID_DATA_ROOT"] = f"{d}/other"            # environment wins
                self.assertEqual(Settings.load().local_root, Path(f"{d}/other"))
                with self.assertRaises(ConfigError):
                    write_default_config(cfg)                             # refuses to overwrite
            finally:
                for k, v in old.items():
                    os.environ.pop(k, None)
                    if v is not None:
                        os.environ[k] = v


class TestRegistry(unittest.TestCase):
    def test_units_and_paths(self):
        with tempfile.TemporaryDirectory() as d:
            L = Layout(d, "fake:")
            write(Path(d) / "Results/v1/Compara/families/TF101001/tree.nwk")
            write(Path(d) / "Analyses/Tissues/summary.tsv")
            write(Path(d) / "Releases/v1/release.json")
            write(Path(d) / "Results/v1/Notes/x.txt")                          # not a part: ignored
            R = Registry(L, index={"results/v2/TreeFam/orthology": {"fingerprint": "x"},
                                   "releases/v2": {"fingerprint": "y"}})
            names = {u.name for u in R.all()}
            self.assertIn("results/v1/Compara/families", names)           # found locally
            self.assertIn("results/v2/TreeFam/orthology", names)          # found in the cloud index
            self.assertIn("releases/v1", names)
            self.assertIn("releases/v2", names)
            self.assertFalse([n for n in names if "Notes" in n])
            self.assertEqual(R.get("releases/v1").policy, "immutable")
            self.assertEqual(R.get("releases/v1").packing, "files")
            self.assertEqual(R.get("results/v1/Compara/run").rel, "Results/v1/Compara/run.json")
            self.assertIn("analyses/Tissues", names)
            u, shard = R.for_path("Results/v1/Compara/families/TF101001/tree.nwk")
            self.assertEqual((u.name, shard), ("results/v1/Compara/families", "TF101"))
            self.assertEqual(R.for_path("Releases/v1/proteomes/Mzeb.faa")[0].name, "releases/v1")
            u, shard = R.for_path("Inputs/Proteins/downloads/Mzeb/protein.faa")
            self.assertEqual((u.name, shard), ("proteins/downloads", "Mzeb"))
            u, shard = R.for_path("Inputs/Proteins/proteomes/Mzeb.faa")
            self.assertEqual((u.name, shard, u.packing), ("proteins/proteomes", None, "tar"))
            u, _ = R.for_path("Inputs/Expression/raw/zebrahub/x.h5ad")
            self.assertEqual((u.name, u.policy), ("expression/raw", "on_demand"))  # the longer match wins
            self.assertEqual({u.name for u in R.matching("families/TreeFam9/*")},
                             {"families/TreeFam9/hmm", "families/TreeFam9/trees", "families/TreeFam9/protein_to_species",
                              "families/TreeFam9/source"})


if __name__ == "__main__":
    unittest.main()


