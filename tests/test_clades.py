"""The species-to-clade table: its API, its validation, and its agreement with the species tree and the Compara panel."""
import io
import os
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Clades, Layout, Settings, newick
from cichlid_data.check import check
from cichlid_data.cli import main as cli
from cichlid_data.clades import COLUMNS, DISPUTED, LADDER
from cichlid_data.config import config_dir
from helpers import releases_toml, write
from test_species_tree import SETS, config

HEADER = "\t".join(COLUMNS) + "\n"


def table(*rows):
    return "# a comment\n" + HEADER + "".join("\t".join(r) + "\n" for r in rows)


GOODROWS = [("maylandia_zebra", "Mzeb", "cichlid", "", "x", ""), ("astatotilapia_calliptera", "Acal", "cichlid", "", "x", ""),
            ("etroplus_suratensis", "Esur", "cichlid", "", "x", ""), ("danio_rerio", "Drer", "teleost", "", "x", ""),
            ("homo_sapiens", "", "gnathostome", "sarcopterygian", "x", ""), ("ciona_intestinalis", "", "outside_vertebrates", "tunicate", "x", "")]


class TestLadderAndLookups(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.cfg = Path(self._t.name)
        write(self.cfg / "species_clades.tsv", table(*GOODROWS))
        self.CL = Clades.load(self.cfg)

    def tearDown(self):
        self._t.cleanup()

    def test_the_ladder(self):
        self.assertEqual(LADDER, ("outside_vertebrates", "cyclostome", "gnathostome", "teleost", "percomorph", "ovalentaria", "cichlid"))
        self.assertEqual([self.CL.rank(x) for x in LADDER], list(range(7)))
        self.assertEqual(DISPUTED, ("cyclostome",))
        with self.assertRaises(KeyError):
            self.CL.rank("fish")

    def test_level_by_prefix_or_name_in_any_case(self):
        for name in ("Mzeb", "mzeb", "maylandia_zebra", "Maylandia_Zebra", "MAYLANDIA_ZEBRA"):
            self.assertEqual(self.CL.level(name), "cichlid", name)
        self.assertEqual(self.CL.level("homo_sapiens"), "gnathostome")
        self.assertEqual((self.CL.lineage("ciona_intestinalis"), self.CL.lineage("Mzeb")), ("tunicate", ""))
        self.assertIn("Drer", self.CL)
        self.assertNotIn("unicorn", self.CL)
        with self.assertRaises(KeyError) as e:
            self.CL.level("unicorn")
        self.assertIn("'unicorn' is not in species_clades.tsv", str(e.exception))

    def test_common_level_is_the_outermost_of_the_levels(self):
        c = self.CL.common_level
        self.assertEqual(c(["Mzeb", "Acal"]), "cichlid")                  # two cichlids meet inside the cichlids
        self.assertEqual(c(["Mzeb", "Drer"]), "teleost")                  # a cichlid and a zebrafish
        self.assertEqual(c(["Mzeb", "homo_sapiens"]), "gnathostome")      # a cichlid and a human
        self.assertEqual(c(["Mzeb", "Drer", "homo_sapiens", "ciona_intestinalis"]), "outside_vertebrates")
        self.assertEqual(c(["Mzeb"]), "cichlid")
        with self.assertRaises(ValueError):
            c([])

    def test_within(self):
        self.assertEqual(self.CL.within("cichlid"), ["maylandia_zebra", "astatotilapia_calliptera", "etroplus_suratensis"])
        self.assertIn("danio_rerio", self.CL.within("teleost"))
        self.assertNotIn("homo_sapiens", self.CL.within("teleost"))
        self.assertEqual(len(self.CL.within("outside_vertebrates")), 6)

    def test_load_errors(self):
        with self.assertRaises(FileNotFoundError):
            Clades.load(self.cfg / "nowhere")
        write(self.cfg / "species_clades.tsv", "species\tclade\nx\tcichlid\n")
        with self.assertRaises(ValueError) as e:
            Clades.load(self.cfg)
        self.assertIn("the header must be species, prefix, clade, lineage, basis, notes", str(e.exception))
        write(self.cfg / "species_clades.tsv", HEADER + "a\tA\tcichlid\t\tb\tn\textra\n")
        with self.assertRaises(ValueError):
            Clades.load(self.cfg)
        write(self.cfg / "species_clades.tsv", HEADER + "a\tA\tcichlid\n")          # short rows are padded
        self.assertEqual(Clades.load(self.cfg).rows[0]["notes"], "")


class TestValidation(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def problems(self, rows, **kw):
        write(self.root / "species_clades.tsv", table(*rows))
        return Clades.load(self.root).validate(**kw)

    def test_good(self):
        self.assertEqual(self.problems(GOODROWS), [])

    def test_row_problems(self):
        bad = self.problems([("a_b", "", "fish", "", "", "")])
        self.assertEqual(bad, ["a_b: clade 'fish' must be one of " + ", ".join(LADDER)])
        self.assertIn("must be lower case with underscores", self.problems([("Homo sapiens", "", "gnathostome", "", "", "")])[0])
        self.assertIn("must be lower case with underscores", self.problems([("Homo_sapiens", "", "gnathostome", "", "", "")])[0])
        dup = self.problems([("a_b", "Ab", "cichlid", "", "", ""), ("A_b", "", "cichlid", "", "", "")])
        self.assertTrue(any("species 'a_b' appears twice (rows 1 and 2)" in p for p in dup), dup)
        dup = self.problems([("a_b", "Ab", "cichlid", "", "", ""), ("c_d", "AB", "cichlid", "", "", "")])
        self.assertEqual(dup, ["prefix 'ab' appears twice (rows 1 and 2); case does not count"])

    def test_every_catalog_organism_needs_a_row(self):
        cfg = config(self.root)
        write(cfg / "species_clades.tsv", table(*GOODROWS[:2]))                      # Esur and Drer missing
        C = Catalog.build(Layout(self.root / "data", "fake:"), cfg)
        probs = Clades.load(cfg).validate(C)
        self.assertEqual(probs, ["Esur (Etroplus suratensis) is in the catalog but has no row",
                                 "Drer (Danio rerio) is in the catalog but has no row"])

    def test_compara_panel_leaves_need_a_row_but_not_excluded_or_validation_species(self):
        write(self.root / "compara_panel.tsv", "# c\nspecies\trole\tgroup\n"
              "homo_sapiens\treference\treference\nfugu_x\tpanel\tpercomorph\nan_ensembl_copy\tvalidation\tcichlid\n"
              "carassius_y\texclude\tfish\n")
        probs = self.problems(GOODROWS, cfg=self.root)
        self.assertEqual(probs, ["fugu_x (Compara panel, panel) has no row"])

    def test_catalog_validate_and_check_use_the_table_when_it_is_there(self):
        cfg = config(self.root, releases=releases_toml({"TreeFam": 'queries = "project"', "Homology": 'description = "calls"'}))
        L = Layout(self.root / "data", "fake:")
        self.assertEqual(Catalog.build(L, cfg).validate()[0], [])                    # no table: nothing to check
        write(cfg / "species_clades.tsv", table(*GOODROWS))
        C = Catalog.build(L, cfg)
        self.assertEqual(C.validate()[0], [])
        self.assertIn("species_clades.tsv", C.meta["config"])
        found = {f.item: f for f in check(L, C, release="v1", parts=["Homology"])}
        self.assertTrue(found["species_clades.tsv"].ok, found["species_clades.tsv"].detail)
        write(cfg / "species_clades.tsv", table(*GOODROWS[:3]))                      # Drer missing now
        C = Catalog.build(L, cfg)
        self.assertEqual(C.validate()[0], ["species_clades.tsv: Drer (Danio rerio) is in the catalog but has no row"])
        found = {f.item: f for f in check(L, C, release="v1", parts=["Homology"])}
        self.assertFalse(found["species_clades.tsv"].ok)
        write(cfg / "species_clades.tsv", "nonsense\n")
        self.assertIn("the header must be", Catalog.build(L, cfg).validate()[0][0])
        (cfg / "species_clades.tsv").unlink()
        found = {f.item: f for f in check(L, Catalog.build(L, cfg), release="v1", parts=["Homology"])}
        self.assertFalse(found["species_clades.tsv"].ok)                             # the Homology part needs it
        self.assertIn("not found", found["species_clades.tsv"].detail)

    def test_catalog_validate_command_reports_it(self):
        cfg = config(self.root)
        write(cfg / "species_clades.tsv", table(*GOODROWS[:2]))
        s = Settings(local_root=self.root / "data", cloud_remote="fake:")
        os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
        try:
            out = io.StringIO()
            with redirect_stdout(out), self.assertRaises(SystemExit):
                cli(["catalog", "validate"], settings=s)
            self.assertIn("species_clades.tsv: Esur", out.getvalue())
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)


class TestTheRepositoryTable(unittest.TestCase):
    """config/species_clades.tsv against the sources it was derived from."""

    def setUp(self):
        self.CL = Clades.load()

    def test_valid_against_the_catalog_and_the_panel(self):
        with tempfile.TemporaryDirectory() as t:
            C = Catalog.build(Layout(t, "fake:"))
            self.assertEqual(self.CL.validate(C, config_dir()), [])
            self.assertEqual(C.validate()[0], [])

    def test_the_orthofinder_species_have_the_levels_the_species_tree_gives(self):
        tree = newick.parse((config_dir() / "species_tree_orthofinder.nwk").read_text())

        def clade(*landmarks):
            best = None
            for node in tree.iter():
                ls = {l.name for l in node.leaves()}
                if set(landmarks) <= ls and (best is None or len(ls) < len(best)):
                    best = ls
            return best
        inner = [("cichlid", clade("Esur", "Acal")), ("ovalentaria", clade("Mzeb", "Xmac")), ("percomorph", clade("Gacu", "Mzeb")),
                 ("teleost", clade("Aang", "Mzeb")), ("gnathostome", clade("Cmil", "Hsap"))]
        self.assertEqual([len(c) for _, c in inner], [10, 16, 19, 26, 34])
        for leaf in (l.name for l in tree.leaves()):
            want = next((name for name, c in inner if leaf in c), "outside_vertebrates")
            self.assertEqual(self.CL.level(leaf), want, leaf)
        self.assertEqual(self.CL.level("Bflo"), "outside_vertebrates")                # amphioxus is the only one outside

    def test_the_compara_panel_labels_agree(self):
        gmap = {"cichlid": "cichlid", "ovalentaria": "ovalentaria", "percomorph": "percomorph", "teleost": "teleost",
                "non-teleost fish": "gnathostome", "outgroup": "outside_vertebrates", "vertebrate": "gnathostome"}
        checked = 0
        for ln in (config_dir() / "compara_panel.tsv").read_text().splitlines():
            f = ln.split("\t")
            if ln.startswith("#") or len(f) < 3 or f[0] == "species" or f[1] in ("exclude", "validation") or f[2] == "reference":
                continue
            want = "cyclostome" if f[0] == "petromyzon_marinus" else gmap[f[2]]
            self.assertEqual(self.CL.level(f[0]), want, f[0])
            checked += 1
        self.assertGreaterEqual(checked, 50)

    def test_known_placements_and_the_disputed_level(self):
        for name, level in (("homo_sapiens", "gnathostome"), ("mus_musculus", "gnathostome"), ("danio_rerio", "teleost"),
                            ("lepisosteus_oculatus", "gnathostome"), ("callorhinchus_milii", "gnathostome"),
                            ("branchiostoma_floridae", "outside_vertebrates"), ("ciona_intestinalis", "outside_vertebrates"),
                            ("petromyzon_marinus", "cyclostome"), ("gasterosteus_aculeatus", "percomorph"),
                            ("oryzias_latipes", "ovalentaria"), ("Metriaclima_zebra", "cichlid"), ("drosophila_melanogaster", "outside_vertebrates")):
            self.assertEqual(self.CL.level(name), level, name)
        self.assertEqual(self.CL.lineage("callorhinchus_milii"), "chondrichthyan")
        self.assertEqual(self.CL.lineage("lepisosteus_oculatus"), "non-teleost actinopterygian")
        self.assertIn("DISPUTED", (config_dir() / "species_clades.tsv").read_text())
        self.assertEqual(self.CL.common_level(["Mzeb", "danio_rerio", "lepisosteus_oculatus"]), "gnathostome")

    def test_all_treefam_cichlid_entries_of_2018_are_cichlids(self):
        for name in ("Astatotilapia_burtoni", "Metriaclima_zebra", "Metriaclima_zebra_b", "Neolamprologus_brichardi", "Pundamilia_nyererei"):
            self.assertEqual(self.CL.level(name), "cichlid")


if __name__ == "__main__":
    unittest.main()
