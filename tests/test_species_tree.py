"""The species tree a method names in species_sets.toml: validated against the catalog's species."""
import io
import re
import tempfile
import tomllib
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Layout, Settings, newick
from cichlid_data.catalog import check_species_tree, species_tree_path
from cichlid_data.check import check
from cichlid_data.cli import main as cli
from cichlid_data.config import config_dir
from helpers import write

SETS = '''
[methods.TreeFam2026]
queries = "project"
[methods.OrthoFinder2026]
species = ["Mzeb", "Acal", "Esur", "Drer"]
species_tree = "tree.nwk"
[[organism]]
prefix = "Mzeb"
species = "Maylandia zebra"
role = "project"
[[organism]]
prefix = "Acal"
species = "Astatotilapia calliptera"
role = "project"
[[organism]]
prefix = "Esur"
species = "Etroplus suratensis"
role = "project"
[[organism]]
prefix = "Drer"
species = "Danio rerio"
role = "reference"
'''
GOOD = "(Drer,(Esur,(Mzeb,Acal)));\n"


def config(root, tree=GOOD, sets=SETS):
    cfg = Path(root) / "config"
    write(cfg / "species_sets.toml", sets)
    write(cfg / "species.tsv", "prefix\tspecies\taccession\tsource\n")
    write(cfg / "compara_panel.tsv", "species\trole\tgroup\n")
    if tree is not None:
        write(cfg / "tree.nwk", tree)
    return cfg


class TestCheckSpeciesTree(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)
        self.sp = ["Mzeb", "Acal", "Esur", "Drer"]

    def tearDown(self):
        self._t.cleanup()

    def problems(self, text, species=None):
        return check_species_tree(write(self.d / "t.nwk", text), species or self.sp)

    def test_good(self):
        self.assertEqual(self.problems(GOOD), [])
        self.assertEqual(self.problems("((Drer,Esur),(Mzeb,Acal));"), [])        # leaf order does not matter
        self.assertEqual(self.problems("(Drer:1,(Esur:2,(Mzeb:3,Acal:4)90:5):6);"), [])   # lengths and labels allowed

    def test_missing_file_and_unparseable(self):
        self.assertEqual(check_species_tree(self.d / "none.nwk", self.sp), [f"none.nwk not found in {self.d}"])
        self.assertIn("does not parse", self.problems("(Drer,(Esur,(Mzeb,Acal))")[0])
        self.assertIn("does not parse", self.problems("")[0])

    def test_leaves_must_be_exactly_the_species(self):
        probs = self.problems("(Drer,(Esur,(Mzeb,Hsap)));")
        self.assertEqual(probs, ["leaves that are not species of the method: Hsap",
                                 "species of the method missing from the tree: Acal"])
        self.assertEqual(self.problems("(Drer,(Esur,(Mzeb,(Acal,Drer))));")[0], "duplicate leaves: Drer")
        self.assertIn("a leaf has no name", self.problems("(Drer,(Esur,(Mzeb,(Acal,))));")[0])

    def test_must_be_rooted_and_fully_resolved(self):
        probs = self.problems("(Drer,Esur,(Mzeb,Acal));")                          # trifurcating root
        self.assertEqual(len(probs), 1)
        self.assertIn("the root has 3 children", probs[0])
        probs = self.problems("(Drer,(Esur,Mzeb,Acal));")                          # a polytomy inside
        self.assertEqual(len(probs), 1)
        self.assertIn("a node (Esur, Mzeb, Acal, ...) has 3 children", probs[0])
        probs = self.problems("((((Drer,Esur,Mzeb,Acal))));")                      # unary nodes
        self.assertTrue(any("has 1 children" in p for p in probs))
        many = "(a,b,c,d,(e,f,g,(h,i,j,(k,l,m,(n,o,p)))));"
        probs = self.problems(many, list("abcdefghijklmnop"))
        self.assertEqual(len(probs), 4)
        self.assertEqual(probs[-1], "... and 2 more unresolved nodes")


class TestCatalogSpeciesTree(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L = Layout(self.root / "data", "fake:")

    def tearDown(self):
        self._t.cleanup()

    def test_valid_tree_passes_and_is_reachable(self):
        cfg = config(self.root)
        C = Catalog.build(self.L, cfg)
        self.assertEqual(C.validate()[0], [])
        self.assertEqual(C.species_tree("OrthoFinder2026"), cfg / "tree.nwk")
        self.assertEqual(species_tree_path("OrthoFinder2026", cfg), cfg / "tree.nwk")
        self.assertEqual(C.methods["OrthoFinder2026"]["species_tree"], "tree.nwk")
        self.assertIn("tree.nwk", C.meta["config"])                                # its checksum is recorded
        self.assertEqual(len(C.meta["config"]["tree.nwk"]), 32)

    def test_species_list_and_tree_must_agree(self):
        cfg = config(self.root, sets=SETS.replace('"Drer"]', '"Drer", "Acen"]') +
                     '[[organism]]\nprefix = "Acen"\nspecies = "Archocentrus centrarchus"\nrole = "project"\n')
        errors, _ = Catalog.build(self.L, cfg).validate()
        self.assertEqual(errors, ["OrthoFinder2026: species tree tree.nwk: species of the method missing from the "
                                  "tree: Acen"])
        cfg = config(self.root, tree="(Drer,(Esur,(Mzeb,(Acal,Hsap))));\n")
        errors, _ = Catalog.build(self.L, cfg).validate()
        self.assertEqual(errors, ["OrthoFinder2026: species tree tree.nwk: leaves that are not species of the "
                                  "method: Hsap"])

    def test_missing_tree_file_is_an_error(self):
        cfg = config(self.root, tree=None)
        errors, _ = Catalog.build(self.L, cfg).validate()
        self.assertEqual(len(errors), 1)
        self.assertIn("tree.nwk not found", errors[0])

    def test_methods_without_a_tree_are_unaffected(self):
        cfg = config(self.root)
        C = Catalog.build(self.L, cfg)
        with self.assertRaises(KeyError) as e:
            C.species_tree("TreeFam2026")
        self.assertIn("no species_tree", str(e.exception))
        with self.assertRaises(KeyError):
            C.species_tree("Nope")
        with self.assertRaises(KeyError):
            species_tree_path("TreeFam2026", cfg)
        cfg = config(self.root, sets=SETS.replace('species_tree = "tree.nwk"\n', ""), tree=None)
        self.assertEqual(Catalog.build(self.L, cfg).validate()[0], [])

    def test_catalog_build_refuses_to_write_an_invalid_tree(self):
        cfg = config(self.root, tree="(Drer,Esur,Mzeb,Acal);\n")
        s = Settings(local_root=self.L.local_root, cloud_remote="fake:")
        import os
        os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
        try:
            with self.assertRaises(SystemExit) as e:
                with redirect_stdout(io.StringIO()):
                    cli(["catalog", "build"], settings=s)
            self.assertIn("catalog not written", str(e.exception))
            self.assertIn("species tree tree.nwk", str(e.exception))
            self.assertFalse(self.L.catalog().local.exists())
            write(cfg / "tree.nwk", GOOD)
            with redirect_stdout(io.StringIO()):
                cli(["catalog", "build"], settings=s)
            self.assertTrue(self.L.catalog().local.is_file())
            C = Catalog.load(self.L)                                                # loaded: tree found via config_dir
            self.assertEqual(C.species_tree("OrthoFinder2026"), cfg / "tree.nwk")
            f = io.StringIO()
            with redirect_stdout(f):
                cli(["catalog", "tree", "--method", "OrthoFinder2026"], settings=s)
            self.assertEqual(f.getvalue().strip(), str(cfg / "tree.nwk"))
            with self.assertRaises(SystemExit) as e:
                cli(["catalog", "tree", "--method", "TreeFam2026"], settings=s)
            self.assertIn("no species_tree", str(e.exception))
            with self.assertRaises(SystemExit) as e:
                cli(["catalog", "tree"], settings=s)
            self.assertIn("needs --method", str(e.exception))
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)


class TestCheckFinding(unittest.TestCase):
    def test_check_reports_the_species_tree(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            L = Layout(root / "data", "fake:")
            for p in ("Mzeb", "Acal", "Esur", "Drer"):
                write(L.proteome(p).local, f">{p}_1\nMK\n")
            cfg = config(root)
            C = Catalog.build(L, cfg)
            found = [f for f in check(L, C, ["OrthoFinder2026"]) if f.item.startswith("species tree")]
            self.assertEqual([(f.item, f.ok, f.detail) for f in found], [("species tree tree.nwk", True, "valid")])
            write(cfg / "tree.nwk", "(Drer,(Esur,(Mzeb,Hsap)));\n")
            found = [f for f in check(L, C, ["OrthoFinder2026"]) if f.item.startswith("species tree")]
            self.assertFalse(found[0].ok)
            self.assertIn("not species of the method: Hsap", found[0].detail)
            self.assertEqual([f.item for f in check(L, C, ["TreeFam2026"]) if f.item.startswith("species tree")], [])


class TestRepositoryFiles(unittest.TestCase):
    """The repository's own configuration."""

    def test_orthofinder_species_tree_matches_the_species_list(self):
        sets = tomllib.loads((config_dir() / "species_sets.toml").read_text())
        m = sets["methods"]["OrthoFinder2026"]
        self.assertEqual(m["species_tree"], "species_tree_orthofinder.nwk")
        self.assertEqual(check_species_tree(config_dir() / m["species_tree"], m["species"]), [])
        root = newick.parse((config_dir() / m["species_tree"]).read_text())
        self.assertEqual(len(root.leaves()), 35)
        self.assertEqual(len(root.children), 2)

    def test_species_tree_agrees_with_species_final(self):
        doc = config_dir().parent / "docs" / "species_final.md"
        if not doc.is_file():
            self.skipTest("docs/species_final.md not present")
        block = re.search(r"```\n(\(Bflo.*?;)\n```", doc.read_text(), re.S)
        self.assertIsNotNone(block, "species_final.md has no species tree block")
        self.assertEqual(block.group(1), (config_dir() / "species_tree_orthofinder.nwk").read_text().strip())


if __name__ == "__main__":
    unittest.main()
