"""OrthoFinder-style family IDs in the families unit, and Results/<release>/<part>/native/ for a part's own outputs."""
import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Layout, Registry, Settings
from cichlid_data.cli import main as cli
from cichlid_data.units import entry_shard, family_shard
from helpers import Env, touch_later, write

REL, PART = "v1", "OrthoFinder"
RUN = f"{REL}/{PART}"                      # as in unit names: results/v1/OrthoFinder/...


class TestFamilyShards(unittest.TestCase):
    def test_shard_names(self):
        for fid, shard in (("TF101001", "TF101"), ("C116_012345", "C116_012"), ("OG0001234", "OG0001"),
                           ("OG0000007", "OG0000"), ("N0.HOG0001234", "N0.HOG0001"), ("N12.HOG0034567", "N12.HOG0034")):
            self.assertEqual(family_shard(fid, True), shard, fid)
        self.assertEqual(family_shard("family_summary.tsv", False), "_top")
        self.assertEqual(family_shard("OG12", True), "_top")                           # too short to shard

    def test_about_a_thousand_families_per_shard(self):
        shards = {family_shard(f"N0.HOG{i:07d}", True) for i in range(5000)}
        self.assertEqual(shards, {f"N0.HOG000{i}" for i in range(5)})

    def test_for_path_resolves_family_ids_with_dots(self):
        """A folder named N0.HOG0001234 is a family even where it is not on this machine yet."""
        with tempfile.TemporaryDirectory() as d:
            R = Registry(Layout(d, "fake:"), index={f"results/{RUN}/families": {"fingerprint": "x"}})
            base = f"Results/{RUN}/families"
            for rel, shard in ((f"{base}/N0.HOG0001234", "N0.HOG0001"),
                               (f"{base}/N0.HOG0001234/tree.nwk", "N0.HOG0001"),
                               (f"{base}/N0.HOG0001234/", "N0.HOG0001"),
                               (f"{base}/OG0001234", "OG0001"),
                               (f"{base}/TF101001", "TF101"),
                               (f"{base}/family_summary.tsv", "_top"),
                               (f"{base}/summary.v2.tsv.gz", "_top"),
                               (f"{base}/ids.json", "_top")):
                u, got = R.for_path(rel)
                self.assertEqual((u.name, got), (f"results/{RUN}/families", shard), rel)

    def test_for_path_for_other_shard_units_is_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            R = Registry(Layout(d, "fake:"))
            self.assertEqual(R.for_path("Inputs/Proteins/downloads/Mzeb")[1], "Mzeb")
            self.assertEqual(R.for_path("Inputs/Proteins/downloads/Mzeb/protein.faa")[1], "Mzeb")
            self.assertEqual(R.for_path("Inputs/Proteins/downloads/notes.txt")[1], "_top")      # a loose file
            u, shard = R.for_path("Inputs/Proteins/downloads/download_summary.tsv")             # its own unit
            self.assertEqual((u.name, shard), ("proteins/download_summary", None))


class TestNativeUnit(unittest.TestCase):
    def test_layout_and_registry(self):
        with tempfile.TemporaryDirectory() as d:
            L = Layout(d, "fake:")
            R = L.results(REL, PART)
            self.assertEqual(R.native().rel, f"Results/{RUN}/native")
            self.assertEqual(R.section("native").rel, R.native().rel)
            self.assertIn("native", R.SECTIONS)
            write(R.native().local / "Results_Oct17" / "Orthogroups" / "Orthogroups.tsv", "x")
            reg = Registry(L)
            u = R.unit("native")
            self.assertEqual((u.name, u.policy, u.packing, u.rel), (f"results/{RUN}/native", "regenerable", "shards",
                                                                    f"Results/{RUN}/native"))
            self.assertEqual(reg.for_path(f"Results/{RUN}/native/Results_Oct17/Orthogroups/Orthogroups.tsv")[1],
                             "Results_Oct17")
            self.assertEqual(reg.for_path(f"Results/{RUN}/native/versions.txt")[1], "_top")
            self.assertEqual(u.cloud_tar("Results_Oct17"), f"Results/{RUN}/native/Results_Oct17.tar")
            self.assertEqual(entry_shard("own_species_tree", True), "own_species_tree")
            self.assertIn(f"results/{RUN}/native", {x.name for x in reg.matching(f"results/{RUN}/*")})
            for other in ("families", "assign", "orthology", "names", "residues", "site"):      # nothing else changed
                self.assertEqual(R.unit(other).name, f"results/{RUN}/{other}")

    def test_where_command(self):
        with tempfile.TemporaryDirectory() as d:
            s = Settings(local_root=Path(d), cloud_remote="fake:")
            f = io.StringIO()
            with redirect_stdout(f):
                cli(["where", "results", REL, PART, "native"], settings=s)
            self.assertEqual(f.getvalue().strip(), str(Path(d) / "Results" / REL / PART / "native"))

    def test_push_pull_and_only_changed_shards_upload(self):
        e = Env()
        try:
            a, b = e.a, e.b
            R = a.results(REL, PART)
            write(R.native().local / "Results_Oct17" / "Orthogroups" / "Orthogroups.tsv", "OG0000000\tMzeb_1\n")
            write(R.native().local / "Results_Oct17" / "Gene_Trees" / "OG0000000_tree.txt", "(a,b);")
            write(R.native().local / "own_species_tree" / "Results_Oct18" / "Species_Tree" / "t.txt", "(a,b);")
            write(R.native().local / "versions.txt", "OrthoFinder x\n")
            res = e.store(a).push(f"results/{RUN}/native")
            self.assertEqual(set(res), {f"results/{RUN}/native/{x}" for x in ("Results_Oct17", "own_species_tree", "_top")})
            self.assertTrue((e.cloud / "Results" / RUN / "native" / "Results_Oct17.tar").is_file())
            touch_later(R.native().local / "own_species_tree" / "Results_Oct18" / "Species_Tree" / "t.txt")
            res = e.store(a).push(f"results/{RUN}/native")
            self.assertEqual(res[f"results/{RUN}/native/own_species_tree"], "uploaded")
            self.assertEqual(res[f"results/{RUN}/native/Results_Oct17"], "up to date")
            e.store(b).pull(f"results/{RUN}/native")
            rb = b.results(REL, PART).native().local
            self.assertEqual((rb / "Results_Oct17" / "Orthogroups" / "Orthogroups.tsv").read_text(), "OG0000000\tMzeb_1\n")
            self.assertTrue((rb / "own_species_tree" / "Results_Oct18" / "Species_Tree" / "t.txt").read_text().endswith("+"))
            self.assertEqual((rb / "versions.txt").read_text(), "OrthoFinder x\n")
            self.assertEqual({st for _, st, _ in e.store(a).status(f"results/{RUN}/native")}, {"in sync"})
        finally:
            e.close()


class TestHierarchicalFamiliesThroughTheStore(unittest.TestCase):
    def test_push_pull_ensure_with_orthogroup_ids(self):
        e = Env()
        try:
            a, b = e.a, e.b
            R = a.results(REL, PART)
            ids = ("N0.HOG0000001", "N0.HOG0000002", "N0.HOG0001500", "N1.HOG0000003", "OG0000009")
            for fid in ids:
                write(R.family(fid).local / "tree.nwk", f"({fid},x);")
                write(R.family(fid).local / "done.json", "{}")
            write(R.family("N0.HOG0001500").local / "stage5" / "scratch.fa", "scratch")      # excluded
            res = e.store(a).push(f"results/{RUN}/families")
            self.assertEqual(set(res), {f"results/{RUN}/families/{s}"
                                        for s in ("N0.HOG0000", "N0.HOG0001", "N1.HOG0000", "OG0000")})
            self.assertTrue((e.cloud / "Results" / RUN / "families" / "N0.HOG0000.tar").is_file())
            # b fetches one family by the path of its folder: only the owning shard comes down
            got = e.store(b).ensure(b.results(REL, PART).family("N0.HOG0001500"))
            self.assertTrue((Path(got) / "tree.nwk").is_file())
            self.assertFalse((Path(got) / "stage5" / "scratch.fa").exists())
            self.assertFalse(b.results(REL, PART).family("N0.HOG0000001").local.exists())
            e.store(b).pull(f"results/{RUN}/families")
            for fid in ids:
                self.assertEqual((b.results(REL, PART).family(fid).local / "tree.nwk").read_text(), f"({fid},x);")
        finally:
            e.close()


if __name__ == "__main__":
    unittest.main()
