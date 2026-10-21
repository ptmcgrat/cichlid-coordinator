"""Which species a leaf of a family tree belongs to, and the shadowed species."""
import tempfile
import tomllib
import unittest
from pathlib import Path

from cichlid_data import Catalog, Layout, LeafSpecies, master_tree as mt
from cichlid_data.check import check
from cichlid_data.config import config_dir
from cichlid_data.leaf_species import MAPPINGS, reference_species_in_mapping
from helpers import write
from test_species_tree import config

TREEFAM = ("222156|estExt_fgenesh1_pg.C_90040\tdrosophila_melanogaster\nBGIBMGA005661-TA\tbombyx_mori\n"
           "ENSGALP00000012291\tgallus_gallus\nISCW010711-PA\tixodes_scapularis\nENSMZEBP0001\tmaylandia_zebra\n"
           "CFAM0001\tcanis_familiaris\nSPTR0001\tSpermophilus_tridecemlineatus\nMzeb_decoy\tgallus_gallus\n")
COMPARA = ("ENSP00000354855\thomo_sapiens\nK08E5.2a.1\tcaenorhabditis_elegans\nENSLOCP00000001444\tlepisosteus_oculatus\n"
           "ENSCAFP00845023514\tcanis_lupus_familiaris\n")


class Base(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L = Layout(self.root / "data", "fake:")
        cfg = config(self.root)
        self.C = Catalog.build(self.L, cfg)
        write(self.L.family_resource("TreeFam9").local / "protein_to_species.txt", TREEFAM)
        write(self.L.family_resource("Compara116").local / "prepared" / "protein_to_species.txt", COMPARA)

    def tearDown(self):
        self._t.cleanup()

    def load(self, part, **kw):
        return LeafSpecies.load(self.L, part, catalog=self.C, **kw)


class TestSpecies(Base):
    def test_reference_leaves_of_either_method(self):
        tf, co = self.load("TreeFam"), self.load("Compara")
        self.assertEqual(tf.species("222156|estExt_fgenesh1_pg.C_90040"), "drosophila_melanogaster")
        self.assertEqual(tf.species("BGIBMGA005661-TA"), "bombyx_mori")
        self.assertEqual(tf.species("ENSGALP00000012291"), "gallus_gallus")
        self.assertEqual(co.species("K08E5.2a.1"), "caenorhabditis_elegans")
        self.assertEqual(co.species("ENSP00000354855"), "homo_sapiens")
        self.assertIsNone(co.species("ENSGALP00000012291"))                  # a TreeFam leaf is not in Compara's file

    def test_names_are_the_master_trees(self):
        tf = self.load("TreeFam")
        self.assertEqual(tf.species("CFAM0001"), "canis_lupus_familiaris")             # canis_familiaris folded
        self.assertEqual(tf.species("SPTR0001"), "ictidomys_tridecemlineatus")         # old genus name folded
        self.assertEqual(self.load("Compara").species("ENSCAFP00845023514"), "canis_lupus_familiaris")

    def test_our_genes_by_prefix_for_every_part(self):
        for part in MAPPINGS:
            ls = self.load(part)
            self.assertEqual(ls.species("Mzeb_101467881"), "maylandia_zebra", part)
            self.assertEqual(ls.species("Drer_5"), "danio_rerio", part)           # a catalog organism, not one of the project's
            self.assertTrue(ls.is_project_gene("Mzeb_101467881"))

    def test_what_is_not_a_leaf_of_either_kind(self):
        ls = self.load("TreeFam")
        for leaf in ("no such leaf", "", "Mzeb", "Mzeb_", "_12", "mzeb_1", "Unknown_5"):         # prefixes are exact-case
            self.assertIsNone(ls.species(leaf), repr(leaf))
            self.assertFalse(ls.is_project_gene(leaf), repr(leaf))

    def test_the_reference_mapping_wins_over_the_prefix_rule(self):
        ls = self.load("TreeFam")
        self.assertEqual(ls.species("Mzeb_decoy"), "gallus_gallus")             # listed as a reference leaf
        self.assertFalse(ls.is_project_gene("Mzeb_decoy"))
        self.assertEqual(self.load("Compara").species("Mzeb_decoy"), "maylandia_zebra")   # not in Compara's file: ours

    def test_of_leaves(self):
        ls = self.load("TreeFam")
        self.assertEqual(ls.of_leaves(["ENSGALP00000012291", "Mzeb_1"]),
                         {"ENSGALP00000012291": "gallus_gallus", "Mzeb_1": "maylandia_zebra"})
        with self.assertRaises(KeyError) as e:
            ls.of_leaves(["Mzeb_1"] + [f"x{i}" for i in range(7)])
        self.assertIn("no species for 7 leaves: x0, x1, x2, x3, x4 (and 2 more)", str(e.exception))

    def test_only_keeps_what_is_asked_for(self):
        ls = self.load("TreeFam", only={"ENSGALP00000012291", "BGIBMGA005661-TA"})
        self.assertEqual(sorted(ls.reference), ["BGIBMGA005661-TA", "ENSGALP00000012291"])
        self.assertEqual(ls.species("ENSGALP00000012291"), "gallus_gallus")
        self.assertIsNone(ls.species("ISCW010711-PA"))                          # not kept, not one of ours


class TestShadowed(Base):
    def test_reference_species_that_are_project_species(self):
        self.assertEqual(self.load("TreeFam").shadowed(), {"maylandia_zebra"})
        self.assertEqual(self.load("Compara").shadowed(), set())
        self.assertEqual(self.load("OrthoFinder").shadowed(), set())              # every species once: our proteome
        self.assertEqual(self.load("TreeFam").project,
                         frozenset({"maylandia_zebra", "astatotilapia_calliptera", "etroplus_suratensis"}))

    def test_reference_species_in_the_mapping(self):
        self.assertEqual(reference_species_in_mapping(self.L, "Compara"),
                         {"homo_sapiens", "caenorhabditis_elegans", "lepisosteus_oculatus", "canis_lupus_familiaris"})
        self.assertIn("ictidomys_tridecemlineatus", reference_species_in_mapping(self.L, "TreeFam"))


class TestErrors(Base):
    def test_a_missing_mapping_file(self):
        (self.L.family_resource("TreeFam9").local / "protein_to_species.txt").unlink()
        with self.assertRaises(FileNotFoundError) as e:
            self.load("TreeFam")
        self.assertIn("check --fetch", str(e.exception))
        self.load("OrthoFinder")                                               # needs no file

    def test_malformed_and_conflicting_lines(self):
        f = self.L.family_resource("TreeFam9").local / "protein_to_species.txt"
        write(f, "a\tx_y\nbroken line\n")
        with self.assertRaises(ValueError) as e:
            self.load("TreeFam")
        self.assertIn("line 2: expected protein<TAB>species", str(e.exception))
        write(f, "a\tx_y\na\tz_w\n")
        with self.assertRaises(ValueError) as e:
            self.load("TreeFam")
        self.assertIn("line 2: a is given two species", str(e.exception))
        write(f, "a\tx_y\n\na\tX_Y\r\nb\tz_w\r\n")                               # blanks, CRLF, a repeat of the same species
        self.assertEqual(self.load("TreeFam").species("b"), "z_w")

    def test_unknown_part(self):
        with self.assertRaises(KeyError):
            LeafSpecies.load(self.L, "Homology", catalog=self.C)


class TestRepositoryConfiguration(unittest.TestCase):
    """The shadowed species of release 1, from the files that define the trees."""

    def project(self):
        sets = tomllib.loads((config_dir() / "species_sets.toml").read_text())
        return {mt.canonical(o["species"]) for o in sets["organism"] if o.get("role") == "project"}

    def test_compara_shadows_betta_and_stickleback_only(self):
        """The five cichlids with an Ensembl copy have the `validation` role: removed from the trees when prepared."""
        roles = {}
        for ln in (config_dir() / "compara_panel.tsv").read_text().splitlines():
            f = ln.split("\t")
            if not ln.startswith("#") and len(f) >= 3 and f[0] != "species":
                roles[mt.canonical(f[0])] = f[1]
        same = {sp: roles[sp] for sp in self.project() if sp in roles}
        self.assertEqual({sp for sp, r in same.items() if r == "validation"},
                         {"oreochromis_niloticus", "neolamprologus_brichardi", "haplochromis_burtoni",
                          "maylandia_zebra", "astatotilapia_calliptera"})
        self.assertEqual({sp for sp, r in same.items() if r != "validation"}, {"betta_splendens", "gasterosteus_aculeatus"})

    def test_treefam_shadows_tilapia_and_stickleback(self):
        tf = {l.name for l in mt.sources()[1].leaves()}
        self.assertEqual(tf & self.project(), {"oreochromis_niloticus", "gasterosteus_aculeatus"})

    def test_pundamilia_and_amphilophus_are_reference_species_not_project_species(self):
        self.assertFalse({"pundamilia_nyererei", "amphilophus_citrinellus"} & self.project())

    def test_every_synonym_resolves_to_a_leaf_of_the_master_tree(self):
        leaves = {l.name for l in mt.load_master_tree().leaves()}
        for alt, target in mt.SYNONYMS.items():
            self.assertIn(target, leaves, alt)


class TestCheckFinding(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def run_check(self, treefam, compara):
        L = Layout(self.root / "d", "fake:")
        if treefam is not None:
            write(L.family_resource("TreeFam9").local / "protein_to_species.txt", treefam)
        if compara is not None:
            write(L.family_resource("Compara116").local / "prepared" / "protein_to_species.txt", compara)
        C = Catalog.build(L)                                                    # the repository's own configuration
        return {f.item: f for f in check(L, C, release="v1", parts=["Homology"])}

    def test_all_species_resolve(self):
        f = self.run_check("ENSGALP00000012291\tgallus_gallus\nX\tSpermophilus_tridecemlineatus\nY\tAstatotilapia_burtoni\n",
                           "ENSP1\thomo_sapiens\n")
        r = f["reference leaf species in the master tree"]
        self.assertTrue(r.ok, r.detail)
        self.assertEqual(r.detail, "all in the master tree (TreeFam 3 species, Compara 1 species)")
        self.assertTrue(f["Inputs/Families/TreeFam9/protein_to_species.txt"].ok)
        self.assertTrue(f["Inputs/Families/Compara116/prepared/protein_to_species.txt"].ok)

    def test_a_species_that_is_not_in_the_tree_is_named(self):
        f = self.run_check("a\tunicorn_equus\n", "b\thomo_sapiens\nc\tdragon_draco\n")
        r = f["reference leaf species in the master tree"]
        self.assertFalse(r.ok)
        self.assertIn("TreeFam: unicorn_equus", r.detail)
        self.assertIn("Compara: dragon_draco", r.detail)

    def test_missing_files_fail(self):
        f = self.run_check(None, None)
        self.assertFalse(f["reference leaf species in the master tree"].ok)
        self.assertFalse(f["Inputs/Families/TreeFam9/protein_to_species.txt"].ok)


if __name__ == "__main__":
    unittest.main()
