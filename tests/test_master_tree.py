"""The master species tree: it is reproducible, it keeps the agreed tree, and every ladder level is a clade."""
import random
import shutil
import tempfile
import unittest
from pathlib import Path

from cichlid_data import Catalog, Clades, Layout, master_tree as mt, newick
from cichlid_data.check import check
from cichlid_data.config import config_dir
from helpers import releases_toml, write


class TestCanonical(unittest.TestCase):
    def test_synonyms_spaces_and_case(self):
        self.assertEqual(mt.canonical("Canis_familiaris"), "canis_lupus_familiaris")
        self.assertEqual(mt.canonical("canis_lupus_familiaris"), "canis_lupus_familiaris")
        self.assertEqual(mt.canonical("Spermophilus_tridecemlineatus"), "ictidomys_tridecemlineatus")
        self.assertEqual(mt.canonical(" Anguilla anguilla "), "anguilla_anguilla")
        self.assertEqual(mt.canonical("schizosaccharomyces_pombe_972h-"), "schizosaccharomyces_pombe")
        self.assertEqual(len(mt.SYNONYMS), 8)                    # the eight spellings TreeFam's species tree uses


class TestTheCommittedTree(unittest.TestCase):
    def setUp(self):
        self.tree = mt.load_master_tree()
        self.CL = Clades.load()

    def test_it_is_what_build_gives(self):
        built = newick.write(mt.build(), internal_names=True) + "\n"
        self.assertEqual(mt.master_tree_path().read_text(), built)

    def test_no_problems(self):
        self.assertEqual(mt.check(self.tree), [])

    def test_the_species_it_holds(self):
        leaves = [l.name for l in self.tree.leaves()]
        self.assertEqual((len(leaves), len(set(leaves))), (154, 154))
        for sp in ("homo_sapiens", "danio_rerio", "maylandia_zebra", "callorhinchus_milii", "lepisosteus_oculatus",
                   "branchiostoma_floridae", "petromyzon_marinus", "tetraodon_nigroviridis", "ictidomys_tridecemlineatus"):
            self.assertIn(sp, leaves)
        for sp in ("canis_familiaris", "astatotilapia_burtoni", "metriaclima_zebra", "pholidichthys_leucotaenia", "amia_calva"):
            self.assertNotIn(sp, leaves)                          # synonyms, the owner's 2018 entries, candidates
        self.assertEqual(len(self.tree.children), 2)

    def test_the_agreed_tree_is_preserved_exactly(self):
        agreed, _, _, _, _ = mt.sources()
        keep = {l.name for l in agreed.leaves()}
        self.assertEqual(len(keep), 35)
        self.assertEqual(mt.clade_sets(self.tree, keep), mt.clade_sets(agreed, keep))

    def test_treefams_non_fish_clades_are_kept(self):
        """TreeFam's tree is a base for everything outside the fish: its clades are clades of the master."""
        _, tf, _, _, _ = mt.sources()
        fish = {"danio_rerio", "gadus_morhua", "tetraodon_nigroviridis", "takifugu_rubripes", "oreochromis_niloticus",
                "gasterosteus_aculeatus", "oryzias_latipes", "xiphophorus_maculatus"}
        keep = {l.name for l in tf.leaves()} - fish
        master_clades = mt.clade_sets(self.tree, keep)
        for c in mt.clade_sets(tf, keep):
            self.assertIn(c, master_clades, sorted(c)[:3])

    def test_every_ladder_level_is_a_clade(self):
        everything = frozenset(l.name for l in self.tree.leaves())
        sets = mt.clade_sets(self.tree)
        for level in self.CL.rows and ("cyclostome", "gnathostome", "teleost", "percomorph", "ovalentaria", "cichlid"):
            members = frozenset(mt.canonical(s) for s in self.CL.within(level)) & everything
            self.assertIn(members, sets, level)
        self.assertEqual(len(frozenset(mt.canonical(s) for s in self.CL.within("cichlid")) & everything), 12)

    def test_the_tree_and_the_ladder_date_a_node_the_same_way(self):
        """For any two species, the level of the clade where they meet in the tree is the ladder's common level."""
        leaves = [l.name for l in self.tree.leaves()]
        rng = random.Random(7)
        for _ in range(250):
            a, b = rng.sample(leaves, 2)
            node = mt._mrca(self.tree, [a, b])
            self.assertEqual(self.CL.common_level(mt._leafset(node)), self.CL.common_level([a, b]), (a, b))

    def test_known_relationships(self):
        t = self.tree
        names = lambda *sp: mt._leafset(mt._mrca(t, sp))
        self.assertEqual(self.CL.common_level(names("maylandia_zebra", "danio_rerio")), "teleost")
        self.assertEqual(self.CL.common_level(names("maylandia_zebra", "homo_sapiens")), "gnathostome")
        self.assertEqual(self.CL.common_level(names("maylandia_zebra", "branchiostoma_floridae")), "outside_vertebrates")
        self.assertEqual(self.CL.common_level(names("homo_sapiens", "petromyzon_marinus")), "cyclostome")
        # the gar is the sister of the teleosts, the reedfish of both, and the shark outside the bony fishes
        self.assertEqual(len(mt._leafset(mt._find(t, "Teleostei"))), 49)
        self.assertIn("lepisosteus_oculatus", mt._leafset(mt._find(t, "Actinopterygii")))
        self.assertNotIn("lepisosteus_oculatus", mt._leafset(mt._find(t, "Teleostei")))
        self.assertNotIn("callorhinchus_milii", mt._leafset(mt._find(t, "Euteleostomi")))
        # tilapia is Ovalentaria with medaka and platyfish (TreeFam's own tree separates them)
        self.assertEqual(self.CL.level("oreochromis_niloticus"), "cichlid")
        self.assertIn("oryzias_latipes", mt._leafset(mt._find(t, "Ovalentaria")))
        self.assertIn("oreochromis_niloticus", mt._leafset(mt._find(t, "Ovalentaria")))
        self.assertIn("tetraodon_nigroviridis", mt._leafset(mt._find(t, "Percomorpha")))


class TestCheckCatchesDamage(unittest.TestCase):
    def setUp(self):
        self.tree = mt.load_master_tree()

    def problems(self, tree):
        return mt.check(tree)

    def test_a_leaf_twice_missing_or_extra(self):
        t = mt.build()
        leaf = next(l for l in t.leaves() if l.name == "danio_rerio")
        leaf.parent.add(newick.Node("danio_rerio"))
        self.assertTrue(any("appear twice: danio_rerio" in p for p in self.problems(t)))
        t = mt.build()
        gone = next(l for l in t.leaves() if l.name == "mus_musculus")
        gone.parent.children.remove(gone)
        self.assertTrue(any("species missing from the tree: mus_musculus" in p for p in self.problems(t)))
        t = mt.build()
        next(l for l in t.leaves() if l.name == "mus_musculus").parent.add(newick.Node("unicorn"))
        probs = self.problems(t)
        self.assertTrue(any("leaves that no method has: unicorn" in p for p in probs))
        self.assertTrue(any("no row in species_clades.tsv: unicorn" in p for p in probs))

    def test_the_root_and_unary_nodes(self):
        t = mt.build()
        t.add(newick.Node("homo_sapiens"))                               # a third child at the root (and a repeat)
        self.assertTrue(any("must be rooted (two)" in p for p in self.problems(t)))
        t = mt.build()
        v = mt._find(t, "Vertebrata")
        wrapper = newick.Node("")
        mt._replace(v, wrapper)
        wrapper.add(v)
        self.assertIn("the tree has a node with a single child", self.problems(t))

    def test_the_agreed_tree_changed(self):
        t = mt.build()
        a = next(l for l in t.leaves() if l.name == "danio_rerio")
        b = next(l for l in t.leaves() if l.name == "gadus_morhua")
        pa, pb = a.parent, b.parent
        pa.children[pa.children.index(a)], pb.children[pb.children.index(b)] = b, a
        a.parent, b.parent = pb, pa
        self.assertIn("pruned to the 35 OrthoFinder species the tree is not the agreed tree", self.problems(t))

    def test_a_ladder_level_that_is_not_a_clade(self):
        t = mt.build()
        # move one cichlid out of the cichlid node, into the percomorph node: the cichlids no longer form a clade
        cichlid = next(l for l in t.leaves() if l.name == "pundamilia_nyererei")
        cichlid.parent.children.remove(cichlid)
        mt._find(t, "Percomorpha").add(cichlid)
        self.assertTrue(any("level cichlid and inside it do not form a clade" in p for p in self.problems(t)), self.problems(t))


class TestInTheCatalogAndCheck(unittest.TestCase):
    def test_catalog_validate_and_check_use_it(self):
        with tempfile.TemporaryDirectory() as d:
            C = Catalog.build(Layout(d, "fake:"))
            self.assertEqual(C.validate()[0], [])
            self.assertIn("species_tree_master.nwk", C.meta["config"])
            self.assertIn("species_tree_treefam9.phy", C.meta["config"])
            found = {f.item: f for f in check(Layout(d, "fake:"), C, release="v1", parts=["Homology"])}
            self.assertTrue(found["species_tree_master.nwk"].ok, found["species_tree_master.nwk"].detail)
            self.assertTrue(found["species_clades.tsv"].ok)

    def test_a_damaged_tree_is_an_error_in_a_copy_of_the_config(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config"
            shutil.copytree(config_dir(), cfg)
            write(cfg / "species_tree_master.nwk", "((a,b),c);\n")
            C = Catalog.build(Layout(Path(d) / "data", "fake:"), cfg)
            errors = C.validate()[0]
            self.assertTrue(any(e.startswith("species_tree_master.nwk: species missing from the tree") for e in errors), errors[:2])
            found = {f.item: f for f in check(Layout(Path(d) / "data", "fake:"), C, release="v1", parts=["Homology"])}
            self.assertFalse(found["species_tree_master.nwk"].ok)
            (cfg / "species_tree_master.nwk").unlink()
            found = {f.item: f for f in check(Layout(Path(d) / "data", "fake:"), Catalog.build(Layout(Path(d) / "data", "fake:"), cfg),
                                              release="v1", parts=["Homology"])}
            self.assertFalse(found["species_tree_master.nwk"].ok)                 # the Homology part needs it


if __name__ == "__main__":
    unittest.main()
