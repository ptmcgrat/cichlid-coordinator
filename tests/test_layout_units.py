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
        R = L.results("Compara2026")
        self.assertEqual(R.assign("Mzeb").rel, "Results/Compara2026/assign/Mzeb.assign.tsv")
        self.assertEqual(R.orthology("Mzeb").rel, "Results/Compara2026/orthology/Mzeb.orthologs.tsv")
        self.assertEqual(R.family("TF101001").rel, "Results/Compara2026/families/TF101001")
        self.assertEqual(R.run_record().rel, "Results/Compara2026/run.json")
        self.assertEqual((L.analysis("Tissues") / "summary.tsv").rel, "Analyses/Tissues/summary.tsv")
        with self.assertRaises(KeyError):
            R.part("bogus")


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
            write(Path(d) / "Results/Compara2026/families/TF101001/tree.nwk")
            write(Path(d) / "Analyses/Tissues/summary.tsv")
            R = Registry(L, index={"results/TreeFam2026/orthology": {"fingerprint": "x"}})
            names = {u.name for u in R.all()}
            self.assertIn("results/Compara2026/families", names)           # found locally
            self.assertIn("results/TreeFam2026/orthology", names)          # found in the cloud index
            self.assertIn("analyses/Tissues", names)
            u, shard = R.for_path("Results/Compara2026/families/TF101001/tree.nwk")
            self.assertEqual((u.name, shard), ("results/Compara2026/families", "TF101"))
            u, shard = R.for_path("Inputs/Proteins/downloads/Mzeb/protein.faa")
            self.assertEqual((u.name, shard), ("proteins/downloads", "Mzeb"))
            u, shard = R.for_path("Inputs/Proteins/proteomes/Mzeb.faa")
            self.assertEqual((u.name, shard, u.packing), ("proteins/proteomes", None, "tar"))
            u, _ = R.for_path("Inputs/Expression/raw/zebrahub/x.h5ad")
            self.assertEqual((u.name, u.policy), ("expression/raw", "on_demand"))  # the longer match wins
            self.assertEqual({u.name for u in R.matching("families/TreeFam9/*")},
                             {"families/TreeFam9/hmm", "families/TreeFam9/trees", "families/TreeFam9/protein_to_species"})


if __name__ == "__main__":
    unittest.main()
