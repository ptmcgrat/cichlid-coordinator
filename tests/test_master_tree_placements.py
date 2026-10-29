"""Where the extra fish go in the master tree: the table, the four operations, the topology it gives, and the polytomies that
are left on purpose. A new polytomy, or a resolved one, must be a decision, so the set is pinned here."""
import shutil
import tempfile
import unittest
from pathlib import Path

from cichlid_data import Catalog, Clades, Layout, master_tree as mt, newick
from cichlid_data.check import check
from cichlid_data.config import config_dir
from helpers import write

PLACEMENTS = config_dir() / mt.PLACEMENTS_FILE
HEADER = "\t".join(mt.PLACEMENT_COLUMNS)


def row(species, op, anchor, group="", support="robust", source="a source", note=""):
    return "\t".join([species, op, anchor, group, support, source, note])


def toy():
    """((a,b),(c,d)) with the nodes named, as the root of a small tree."""
    return newick.parse("((a,b)AB,(c,d)CD)ROOT;")


def clade(t, *names):
    return mt._leafset(mt._mrca(t, names))


class TestTheTable(unittest.TestCase):
    def test_it_is_valid_and_every_row_is_cited(self):
        self.assertEqual(mt.placement_problems(), [])
        rows = mt.load_placements()
        self.assertEqual(len(rows), 22)
        for r in rows:
            self.assertTrue(r["source"], r["species"])
            self.assertIn(r["support"], mt.SUPPORT)
            self.assertIn(r["op"], mt.PLACEMENT_OPS)

    def test_it_places_every_missing_fish_that_treefam_does_not_already_put_beside_a_relative(self):
        agreed, tf, wanted, clades, _ = mt.sources()
        extras = {s for s in wanted - mt._leafset(agreed) if clades.level(s) in mt.TELEOST_LEVELS}
        in_treefam = extras & mt._leafset(tf)
        self.assertEqual(in_treefam, {"tetraodon_nigroviridis"})                        # the one fish the older rule places
        self.assertEqual({r["species"] for r in mt.load_placements()}, extras - in_treefam)
        # ... and it sits beside its TreeFam relative, as before
        T = mt.load_master_tree()
        self.assertEqual(mt._leafset(mt._mrca(T, ["takifugu_rubripes", "tetraodon_nigroviridis"])), {"takifugu_rubripes", "tetraodon_nigroviridis"})

    def test_the_four_nodes_project_3_asked_about_have_a_row_for_every_fish_in_them(self):
        placed = {r["species"] for r in mt.load_placements()}
        for sp in ("denticeps_clupeoides", "ictalurus_punctatus", "electrophorus_electricus", "myripristis_murdjan",   # teleost base
                   "sparus_aurata", "dicentrarchus_labrax", "labrus_bergylta", "seriola_dumerili", "cottoperca_gobio", "larimichthys_crocea",
                   "lates_calcarifer", "sander_lucioperca", "scophthalmus_maximus", "anabas_testudineus", "mastacembelus_armatus",
                   "hippocampus_comes", "cyclopterus_lumpus",                                                            # percomorphs
                   "acanthochromis_polyacanthus", "cyprinodon_variegatus", "kryptolebias_marmoratus",                   # Ovalentaria
                   "amphilophus_citrinellus", "pundamilia_nyererei"):                                                    # cichlids
            self.assertIn(sp, placed)

    def test_parse_errors(self):
        with self.assertRaises(ValueError):
            mt.parse_placements("species\top\n")
        with self.assertRaises(ValueError):
            mt.parse_placements("# only comments\n")
        with self.assertRaises(ValueError) as e:
            mt.parse_placements(HEADER + "\nx\tsister_of\n")
        self.assertIn("expected 7 columns, found 2", str(e.exception))
        self.assertEqual(mt.parse_placements("# c\n\n" + HEADER + "\n\n" + row("x", "child_of", "a") + "\n")[0]["species"], "x")


class TestTheOperations(unittest.TestCase):
    def apply(self, t, *rows, labels=None):
        return mt.apply_placements(t, labels or {}, [dict(zip(mt.PLACEMENT_COLUMNS, r.split("\t"))) for r in rows])

    def test_sister_of_makes_a_pair_in_the_place_of_the_anchor(self):
        t = toy()
        self.apply(t, row("x", "sister_of", "a", "XA"))
        node = mt._find(t, "XA")
        self.assertEqual(mt._leafset(node), {"a", "x"})
        self.assertIs(node.parent, mt._find(t, "AB"))
        self.assertEqual(len(node.parent.children), 2)               # the anchor's place is taken, not widened
        self.assertEqual(mt._leafset(mt._find(t, "AB")), {"a", "b", "x"})

    def test_sister_of_a_clade_keeps_the_clade_whole(self):
        t = toy()
        self.apply(t, row("x", "sister_of", "a,b"))
        self.assertEqual(clade(t, "a", "b"), {"a", "b"})
        self.assertEqual(clade(t, "a", "x"), {"a", "b", "x"})

    def test_join_group_makes_a_polytomy_that_holds_the_anchor(self):
        t = toy()
        self.apply(t, row("x", "join_group", "a", "G"), row("y", "join_group", "a", "G"), row("z", "join_group", "a", "G"))
        g = mt._find(t, "G")
        self.assertEqual({c.name for c in g.children}, {"a", "x", "y", "z"})
        self.assertEqual(mt._leafset(mt._find(t, "AB")), {"a", "b", "x", "y", "z"})
        self.assertEqual(len(mt._find(t, "AB").children), 2)         # (G, b)

    def test_join_group_of_a_clade_anchor_keeps_the_clade(self):
        t = toy()
        self.apply(t, row("x", "join_group", "a,b", "G"), row("y", "join_group", "a,b", "G"))
        g = mt._find(t, "G")
        self.assertEqual(len(g.children), 3)
        self.assertEqual(mt._leafset(g.children[0]), {"a", "b"})     # the agreed pair is intact inside the group

    def test_group_under_hangs_a_new_group_from_the_anchor(self):
        t = toy()
        self.apply(t, row("x", "group_under", "@AB", "G"), row("y", "group_under", "@AB", "G"), labels={"AB": mt._find(t, "AB")})
        g = mt._find(t, "G")
        self.assertEqual({c.name for c in g.children}, {"x", "y"})
        self.assertEqual({c.name for c in mt._find(t, "AB").children}, {"a", "b", "G"})

    def test_child_of_widens_the_anchor(self):
        t = toy()
        self.apply(t, row("x", "child_of", "a,b"))
        self.assertEqual({c.name for c in mt._find(t, "AB").children}, {"a", "b", "x"})

    def test_a_row_can_anchor_on_a_fish_an_earlier_row_placed(self):
        t = toy()
        self.apply(t, row("x", "sister_of", "a", "XA"), row("y", "sister_of", "a,x", "YXA"))
        self.assertEqual(mt._leafset(mt._find(t, "YXA")), {"a", "x", "y"})

    def test_errors(self):
        for rows, msg in (
                ((row("a", "child_of", "b"),), "already in the tree"),
                ((row("x", "sister_of", "nope"),), "anchor species not in the tree: nope"),
                ((row("x", "sister_of", "@Nope"),), "unknown anchor @Nope"),
                ((row("x", "sister_of", "a,b,c,d"),), "cannot be the sister of the root"),
                ((row("x", "join_group", "a,b,c,d", "G"),), "cannot group with the root"),
                ((row("x", "join_group", "a", "G"), row("y", "join_group", "c", "G")), "its anchor is not inside group G"),
                ((row("x", "wobble", "a"),), "unknown op"),
        ):
            with self.assertRaises(ValueError) as e:
                self.apply(toy(), *rows)
            self.assertIn(msg, str(e.exception), rows)
        t = toy()
        with self.assertRaises(ValueError) as e:
            self.apply(t, row("x", "group_under", "@AB", "G"), row("y", "group_under", "@CD", "G"),
                       labels={"AB": mt._find(t, "AB"), "CD": mt._find(t, "CD")})
        self.assertIn("hangs under a different node", str(e.exception))


class TestValidationOfTheTable(unittest.TestCase):
    def problems_with(self, body):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config"
            shutil.copytree(config_dir(), cfg)
            write(cfg / mt.PLACEMENTS_FILE, HEADER + "\n" + body)
            return mt.placement_problems(cfg)

    def test_each_kind_of_mistake(self):
        for body, msg in (
                (row("denticeps_clupeoides", "sister_of", "clupea_harengus", "X", "robust", "") + "\n", "no source"),
                (row("denticeps_clupeoides", "sister_of", "clupea_harengus", "X", "sure") + "\n", "support 'sure' must be one of"),
                (row("denticeps_clupeoides", "wobble", "clupea_harengus") + "\n", "op 'wobble' must be one of"),
                (row("denticeps_clupeoides", "sister_of", "") + "\n", "no anchor"),
                (row("denticeps_clupeoides", "sister_of", "not_a_fish") + "\n", "anchor not_a_fish is not in the agreed tree"),
                (row("denticeps_clupeoides", "sister_of", "@Nowhere") + "\n", "anchor @Nowhere must be one of"),
                (row("denticeps_clupeoides", "join_group", "clupea_harengus") + "\n", "join_group needs a group name"),
                (row("denticeps_clupeoides", "child_of", "clupea_harengus", "X") + "\n", "child_of makes no new node"),
                (row("clupea_harengus", "sister_of", "danio_rerio") + "\n", "not a fish that is missing from the agreed tree"),
                (row("homo_sapiens", "sister_of", "danio_rerio") + "\n", "not a fish that is missing from the agreed tree"),
                (row("denticeps_clupeoides", "child_of", "clupea_harengus") + "\n" + row("denticeps_clupeoides", "child_of", "clupea_harengus") + "\n",
                 "the species appears twice"),
                # an anchor must exist when the row is applied: the fish placed LATER cannot be an anchor
                (row("denticeps_clupeoides", "sister_of", "ictalurus_punctatus") + "\n" + row("ictalurus_punctatus", "child_of", "danio_rerio") + "\n",
                 "anchor ictalurus_punctatus is not in the agreed tree or placed by an earlier row"),
        ):
            probs = self.problems_with(body)
            self.assertTrue(any(msg in p for p in probs), (msg, probs))

    def test_a_bad_header_is_reported_not_raised(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config"
            shutil.copytree(config_dir(), cfg)
            write(cfg / mt.PLACEMENTS_FILE, "species\top\n")
            self.assertTrue(any("the header must be" in p for p in mt.placement_problems(cfg)))

    def test_an_invalid_table_stops_the_build_and_shows_in_check(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config"
            shutil.copytree(config_dir(), cfg)
            write(cfg / mt.PLACEMENTS_FILE, HEADER + "\n" + row("denticeps_clupeoides", "sister_of", "clupea_harengus", source="") + "\n")
            with self.assertRaises(ValueError) as e:
                mt.build(cfg)
            self.assertIn("the placement table is not valid", str(e.exception))
            self.assertTrue(any("no source" in p for p in mt.check(mt.load_master_tree(), cfg)))


class TestTheFallbackStillWorks(unittest.TestCase):
    def test_a_fish_without_a_row_is_attached_by_the_older_rule_and_the_tree_stays_valid(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config"
            shutil.copytree(config_dir(), cfg)
            text = PLACEMENTS.read_text().splitlines(True)
            kept = [ln for ln in text if not ln.startswith("hippocampus_comes\t")]
            self.assertEqual(len(kept), len(text) - 1)                                   # the row was there, and is gone
            write(cfg / mt.PLACEMENTS_FILE, "".join(kept))
            tree = mt.build(cfg)
            self.assertEqual(mt.check(tree, cfg), [])
            self.assertIn("hippocampus_comes", mt._leafset(mt._find(tree, "Percomorpha")))
            self.assertIs(mt._find(tree, "hippocampus_comes").parent, mt._find(tree, "Percomorpha"))   # a ladder-level child
            # and nothing else moved
            self.assertEqual(mt.clade_sets(tree, {r["species"] for r in mt.load_placements()} - {"hippocampus_comes"}),
                             mt.clade_sets(mt.load_master_tree(), {r["species"] for r in mt.load_placements()} - {"hippocampus_comes"}))

    def test_without_the_table_the_old_polytomies_come_back(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config"
            shutil.copytree(config_dir(), cfg)
            (cfg / mt.PLACEMENTS_FILE).unlink()
            tree = mt.build(cfg)
            self.assertEqual(mt.check(tree, cfg), [])
            self.assertEqual(len(mt._find(tree, "Percomorpha").children), 15)
            self.assertEqual(len(mt._find(tree, "Teleostei").children), 6)


class TestTheTopology(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.T = mt.load_master_tree()
        cls.CL = Clades.load()
        cls.node = {n.name: n for n in cls.T.iter() if n.name}

    def leaves(self, name):
        return mt._leafset(self.node[name])

    def test_the_teleost_base(self):
        tel = self.node["Teleostei"]
        self.assertEqual(len(tel.children), 2)                                        # was 6
        self.assertEqual({c.name for c in tel.children if not c.children}, {"anguilla_anguilla"})
        # herrings: Denticeps is sister to Clupea, and Clupeiformes is sister to Otophysi
        self.assertEqual(self.leaves("Clupeiformes"), {"clupea_harengus", "denticeps_clupeoides"})
        clup = self.node["Clupeiformes"]
        self.assertEqual(mt._leafset(clup.parent), {"clupea_harengus", "denticeps_clupeoides", "danio_rerio", "astyanax_mexicanus",
                                                    "ictalurus_punctatus", "electrophorus_electricus"})
        # Otophysi: zebrafish is basal; tetra, catfish and electric eel stay unresolved (three published hypotheses)
        chara = self.node["Characiphysi"]
        self.assertEqual({c.name for c in chara.children}, {"astyanax_mexicanus", "ictalurus_punctatus", "electrophorus_electricus"})
        self.assertEqual({c.name for c in chara.parent.children}, {"danio_rerio", "Characiphysi"})

    def test_the_soldierfish_is_sister_to_the_percomorphs_and_the_cod_is_outside_both(self):
        perc = self.node["Percomorpha"]
        self.assertNotIn("myripristis_murdjan", mt._leafset(perc))
        pair = perc.parent
        self.assertEqual({c.name or "X" for c in pair.children}, {"Percomorpha", "myripristis_murdjan"})
        outer = pair.parent
        self.assertEqual({c.name for c in outer.children if not c.children}, {"gadus_morhua"})
        self.assertEqual(self.CL.level("myripristis_murdjan"), "teleost")             # the ladder is unchanged: no acanthomorph level
        self.assertEqual(self.CL.level("gadus_morhua"), "teleost")

    def test_the_percomorphs(self):
        perc = self.node["Percomorpha"]
        self.assertEqual(len(perc.children), 4)                                        # was 15
        kids = {c.name or ("agreed:" + "+".join(sorted(mt._leafset(c))[:1])) for c in perc.children}
        self.assertEqual({c.name for c in perc.children if not c.children}, {"hippocampus_comes"})
        self.assertEqual(self.leaves("Carangaria"), {"seriola_dumerili", "scophthalmus_maximus", "lates_calcarifer"})
        self.assertIs(self.node["Carangaria"].parent, perc)
        euper = self.node["Eupercaria"]
        self.assertEqual(self.leaves("Eupercaria"), {"gasterosteus_aculeatus", "takifugu_rubripes", "tetraodon_nigroviridis", "sparus_aurata",
                                                      "dicentrarchus_labrax", "labrus_bergylta", "larimichthys_crocea", "sander_lucioperca",
                                                      "cottoperca_gobio", "cyclopterus_lumpus"})
        agreed_pair = [c for c in euper.children if c.children]
        self.assertEqual(len(agreed_pair), 1)                                          # the agreed stickleback and pufferfishes stay together
        self.assertEqual(mt._leafset(agreed_pair[0]), {"gasterosteus_aculeatus", "takifugu_rubripes", "tetraodon_nigroviridis"})
        self.assertIs(euper.parent, perc)
        del kids

    def test_anabantaria_sits_with_betta_and_is_sister_to_ovalentaria_as_in_the_agreed_tree(self):
        self.assertEqual(self.leaves("Anabantiformes"), {"betta_splendens", "anabas_testudineus"})
        self.assertEqual(self.leaves("Anabantaria"), {"betta_splendens", "anabas_testudineus", "mastacembelus_armatus"})
        self.assertEqual({c.name for c in self.node["Anabantaria"].parent.children}, {"Anabantaria", "Ovalentaria"})

    def test_ovalentaria_is_binary_and_complete(self):
        ov = self.node["Ovalentaria"]
        self.assertEqual(len(ov.children), 2)                                          # was 5
        self.assertEqual(self.leaves("Pomacentridae"), {"amphiprion_ocellaris", "acanthochromis_polyacanthus"})
        self.assertEqual(self.leaves("Aplocheiloidei"), {"nothobranchius_furzeri", "kryptolebias_marmoratus"})
        self.assertEqual(len(mt._leafset(ov)), 21)
        # Cyprinodon joins the killifish node of Fundulus and the poeciliids, unresolved
        node = mt._mrca(self.T, ["fundulus_heteroclitus", "poecilia_reticulata"])
        self.assertEqual(len(node.children), 3)
        self.assertEqual({c.name for c in node.children if not c.children}, {"fundulus_heteroclitus", "cyprinodon_variegatus"})

    def test_the_cichlids(self):
        ci = self.node["Cichlidae"]
        self.assertEqual(len(ci.children), 2)                                          # was 4
        self.assertEqual(len(mt._leafset(ci)), 12)
        self.assertEqual(self.leaves("Cichlinae"), {"archocentrus_centrarchus", "amphilophus_citrinellus"})
        # Pundamilia sits with Haplochromis, Maylandia and Astatotilapia calliptera, unresolved (Brawand 2014, Fig. 1)
        hap = mt._mrca(self.T, ["haplochromis_burtoni", "maylandia_zebra", "pundamilia_nyererei"])
        self.assertEqual(mt._leafset(hap), {"haplochromis_burtoni", "maylandia_zebra", "astatotilapia_calliptera", "pundamilia_nyererei"})
        self.assertEqual(len(hap.children), 3)
        self.assertEqual(self.CL.level("pundamilia_nyererei"), "cichlid")

    def test_where_the_cichlid_carriers_of_a_duplicate_meet_is_no_longer_the_whole_clade(self):
        """Project 3's B2: a duplicate shared by M. zebra, A. calliptera and Pundamilia used to meet at the whole cichlid clade (12)."""
        meet = mt._mrca(self.T, ["maylandia_zebra", "astatotilapia_calliptera", "pundamilia_nyererei"])
        self.assertEqual(len(mt._leafset(meet)), 4)
        meet = mt._mrca(self.T, ["maylandia_zebra", "amphilophus_citrinellus"])
        self.assertEqual(len(mt._leafset(meet)), 11)                                   # with the African cichlids and the Neotropical pair

    def test_the_polytomies_that_are_left_are_the_ones_we_decided_to_leave(self):
        fish = mt._leafset(self.node["Teleostei"])
        got = {(name or "(unnamed)", n, len(ls)) for name, n, ls in mt.polytomies(self.T) if ls <= fish}
        self.assertEqual(got, {("Percomorpha", 4, 38),          # the order of the percomorph series is not resolved in the sources used
                               ("Eupercaria", 8, 10),            # Betancur-R's unresolved bush
                               ("Carangaria", 3, 3),
                               ("Characiphysi", 3, 3),           # three published hypotheses disagree
                               ("(unnamed)", 3, 4),              # Cyprinodontoidei: Cyprinodon, Fundulus, the poeciliids
                               ("(unnamed)", 3, 4)})             # haplochromines: Haplochromis, Malawi, Pundamilia (Brawand 2014)

    def test_every_contested_row_has_a_polytomy_and_no_polytomy_comes_from_a_robust_row_alone(self):
        rows = mt.load_placements()
        self.assertTrue(any(r["support"] == "contested" for r in rows))
        for r in rows:
            if r["op"] == "sister_of":
                self.assertEqual(r["support"], "robust", r["species"])               # a pair is a resolution, so it must be supported

    def test_names_in_the_fish_clade_must_be_unique_but_treefams_own_duplicate_elsewhere_is_tolerated(self):
        """Project 3 looks nodes up by name inside the fish clade. TreeFam's tree has two nodes called Drosophila, which is not ours."""
        names = [n.name for n in self.T.iter() if n.children and n.name]
        self.assertGreater(names.count("Drosophila"), 1)                              # the duplicate is real ...
        self.assertEqual(mt.check(self.T), [])                                        # ... and tolerated, outside the fish
        t = mt.load_master_tree()
        {n.name: n for n in t.iter() if n.name}["Eupercaria"].name = "Percomorpha"     # a duplicate inside the fish clade
        self.assertTrue(any("internal node names in the fish clade that appear twice: Percomorpha" in p for p in mt.check(t)))

    def test_the_catalog_freezes_the_table(self):
        with tempfile.TemporaryDirectory() as d:
            C = Catalog.build(Layout(d, "fake:"))
            self.assertIn(mt.PLACEMENTS_FILE, C.meta["config"])
            found = {f.item: f for f in check(Layout(d, "fake:"), C, release="v1", parts=["Homology"])}
            self.assertTrue(found["species_tree_master.nwk"].ok, found["species_tree_master.nwk"].detail)


if __name__ == "__main__":
    unittest.main()
