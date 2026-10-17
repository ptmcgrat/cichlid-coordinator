"""The family contract: what any method's Results/<run>/families/<family>/ must hold."""
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Settings, schemas
from cichlid_data.cli import main as cli
from helpers import write

TREE = "((Mzeb_1,Acal_1),(Hsap_9,ENSDARP00000012345));"
ROOTED = "(((Mzeb_1,Acal_1),Hsap_9),ENSDARP00000012345);"
UNROOTED = "(Mzeb_1,Acal_1,(Hsap_9,ENSDARP00000012345));"


def family(parent, fid, tree=TREE, done=None, alignment=None):
    d = Path(parent) / fid
    write(d / "tree.nwk", tree)
    write(d / "done.json", json.dumps({} if done is None else done))
    if alignment is not None:
        write(d / "alignment.fa", alignment)
    return d


class TestFamilyIds(unittest.TestCase):
    def test_valid_and_invalid(self):
        for ok in ("TF101001", "C116_012345", "OG0001234", "N0.HOG0001234", "N3.HOG0000007", "a-b", "7"):
            self.assertTrue(schemas.valid_family_id(ok), ok)
        for bad in ("", ".hidden", "_top", "-x", "a b", "a/b", "OG:1", "N0.HOG 1", "fam\u00e9"):
            self.assertFalse(schemas.valid_family_id(bad), bad)

    def test_family_folder_name_is_checked(self):
        with tempfile.TemporaryDirectory() as t:
            good = family(t, "N0.HOG0001234")
            bad = family(t, "N0.HOG 1234")
            self.assertEqual(schemas.validate(good), [])
            self.assertEqual(len(schemas.validate(bad)), 1)
            self.assertIn("family ID 'N0.HOG 1234'", schemas.validate(bad)[0])

    def test_assign_for_orthofinder(self):
        """best_family is the orthogroup; the score columns may be empty."""
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "Mzeb.assign.tsv"
            rows = [{"gene_id": "Mzeb_1", "protein_length": 400, "best_family": "N0.HOG0001234"},
                    {"gene_id": "Mzeb_2", "protein_length": 90, "best_family": "OG0000007"},
                    {"gene_id": "Mzeb_3", "protein_length": 50}]                    # in no family
            schemas.ASSIGN.write(p, rows, prefix="Mzeb")
            self.assertEqual(schemas.validate(p), [])
            self.assertEqual(schemas.ASSIGN.read(p)[0]["score"], "")
            with self.assertRaises(schemas.SchemaError):
                schemas.ASSIGN.write(p, [{"gene_id": "Mzeb_1", "best_family": "OG 12"}])
            rows = [{"gene_id": "Mzeb_1", "best_family": "TF101001", "second_family": "bad family"}]
            self.assertIn("second_family 'bad family' is not a family ID",
                          "\n".join(schemas.ASSIGN.check_rows(rows)))
            for fam in ("TF101001", "C116_012345"):                                    # the existing methods' IDs
                self.assertEqual(schemas.ASSIGN.check_rows([{"gene_id": "Mzeb_1", "best_family": fam}]), [])


class TestLeaves(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_good_tree(self):
        self.assertEqual(schemas.validate(family(self.d, "OG0000001")), [])
        self.assertEqual(schemas.validate(family(self.d, "OG0000002"), kind="family", prefix=None), [])

    def test_orthofinder_internal_names_refused(self):
        f = family(self.d, "OG0000001", "((0_123,0_124),(3_5,Hsap_9));")
        probs = schemas.validate(f)
        self.assertEqual(len(probs), 3)
        self.assertIn("'0_123' looks like an OrthoFinder internal sequence name", probs[0])

    def test_duplicate_empty_and_whitespace_leaves(self):
        self.assertEqual(schemas.validate(family(self.d, "OG1000001", "((a,a),b);")), ["tree.nwk: duplicate leaf 'a'"])
        self.assertEqual(schemas.validate(family(self.d, "OG1000002", "((a,),b);")), ["tree.nwk has a leaf with no name"])
        self.assertIn("contains whitespace", schemas.validate(family(self.d, "OG1000003", "((a b,c),d);"))[0])

    def test_known_prefix_must_be_a_gene_id(self):
        f = family(self.d, "OG0000001", "((Mzeb_XP_1,Acal_1),(Hsap_9,ENSDARP00000012345));")
        self.assertEqual(schemas.validate(f), [])                                      # could be a reference ID
        probs = schemas.validate(f, prefixes={"Mzeb", "Acal", "Hsap"})
        self.assertEqual(len(probs), 1)
        self.assertIn("'Mzeb_XP_1' starts with a species prefix but is not a gene ID", probs[0])
        self.assertEqual(schemas.validate(family(self.d, "OG0000002"), prefixes={"Mzeb", "Acal", "Hsap"}), [])

    def test_tree_that_does_not_parse(self):
        probs = schemas.validate(family(self.d, "OG0000003", "((a:x,b),c);"))
        self.assertEqual(len(probs), 1)
        self.assertIn("tree.nwk does not parse", probs[0])
        # the original messages are unchanged
        self.assertEqual(schemas.validate(family(self.d, "OG0000004", "((a,b),c;")), ["tree.nwk has unbalanced parentheses"])
        self.assertEqual(schemas.validate(family(self.d, "OG0000005", "((a,b),c)")), ["tree.nwk does not end with ';'"])

    def test_trees_with_lengths_supports_and_nhx(self):
        t = "((Mzeb_1:0.1,Acal_1:0.2)0.95:0.05,(Hsap_9:0.3,ENSDARP1:0.4)[&&NHX:D=Y]:0.1);"
        self.assertEqual(schemas.validate(family(self.d, "OG0000006", t)), [])


class TestRooting(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_rooted_true_needs_how_and_a_rooted_tree(self):
        ok = family(self.d, "OG0000001", ROOTED, {"rooted": True, "rooting": "orthofinder species tree"})
        self.assertEqual(schemas.validate(ok), [])
        self.assertEqual(schemas.family_rooting(ok), (True, "orthofinder species tree"))
        no_how = family(self.d, "OG0000002", ROOTED, {"rooted": True})
        self.assertIn('needs a "rooting" text', schemas.validate(no_how)[0])
        blank = family(self.d, "OG0000003", ROOTED, {"rooted": True, "rooting": "  "})
        self.assertIn('needs a "rooting" text', schemas.validate(blank)[0])
        liar = family(self.d, "OG0000004", UNROOTED, {"rooted": True, "rooting": "midpoint"})
        self.assertEqual(schemas.validate(liar),
                         ["done.json says rooted but the root of tree.nwk has 3 children (a rooted tree has two)"])

    def test_rooted_false_or_absent_is_unrooted(self):
        for i, done in enumerate(({}, {"rooted": False}, {"rooted": False, "rooting": "x"}, {"inputs": "abc"})):
            f = family(self.d, f"OG000010{i}", UNROOTED, done)
            self.assertEqual(schemas.validate(f), [], done)
            self.assertEqual(schemas.family_rooting(f), (False, ""))
        # even a tree that happens to have a two-child root is unrooted unless done.json says so
        self.assertEqual(schemas.family_rooting(family(self.d, "OG0000200", ROOTED, {})), (False, ""))

    def test_rooted_must_be_a_boolean(self):
        for bad in ("yes", 1, None, "true"):
            probs = schemas.validate(family(self.d, "OG0000300", ROOTED, {"rooted": bad}))
            self.assertEqual(len(probs), 1, bad)
            self.assertIn("rooted must be true or false", probs[0])
            # a malformed value is never taken as "rooted": consumers then root the tree themselves
            (self.d / "OG0000300" / "done.json").write_text(json.dumps({"rooted": bad, "rooting": "x"}))
            self.assertEqual(schemas.family_rooting(self.d / "OG0000300"), (False, ""), bad)
            (self.d / "OG0000300" / "done.json").unlink()
            (self.d / "OG0000300" / "tree.nwk").unlink()
            (self.d / "OG0000300").rmdir()

    def test_done_json_must_be_an_object(self):
        f = family(self.d, "OG0000400")
        (f / "done.json").write_text("[1, 2]")
        self.assertEqual(schemas.validate(f), ["done.json must be a JSON object"])
        self.assertEqual(schemas.family_rooting(f), (False, ""))
        self.assertIsNone(schemas.read_done(f))
        (f / "done.json").write_text("{not json")
        self.assertEqual(schemas.validate(f), ["done.json is not valid JSON"])
        self.assertEqual(schemas.family_rooting(f), (False, ""))
        self.assertEqual(schemas.family_rooting(self.d / "none"), (False, ""))


class TestAlignment(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    ALN = ">Mzeb_1 x\nMK-L\n>Acal_1\nMKQL\n>Hsap_9\nMRQL\n>ENSDARP00000012345\nMKQL\n"

    def test_good_and_optional(self):
        self.assertEqual(schemas.validate(family(self.d, "OG0000001", alignment=self.ALN)), [])
        self.assertEqual(schemas.validate(family(self.d, "OG0000002")), [])            # no alignment file: fine

    def test_problems(self):
        f = family(self.d, "OG0000003", alignment=">Mzeb_1\nMKL\n>Acal_1\nMKQL\n>Hsap_9\nMRQL\n")
        probs = schemas.validate(f)
        self.assertEqual(len(probs), 2)
        self.assertIn("rows have different lengths (3 to 4)", probs[0])
        self.assertIn("tree leaves without a row: ENSDARP00000012345", probs[1])
        f = family(self.d, "OG0000004", alignment="")
        self.assertEqual(schemas.validate(f), ["alignment.fa is empty"])
        f = family(self.d, "OG0000005", alignment=self.ALN + ">Acal_1\nMKQL\n")
        self.assertEqual(schemas.validate(f), ["alignment.fa: duplicate row 'Acal_1'"])
        f = family(self.d, "OG0000006", alignment=">0_1\nMKQL\n>0_2\nMKQL\n")           # names not mapped back
        self.assertIn("tree leaves without a row", schemas.validate(f)[0])

    def test_extra_rows_are_allowed_and_check_can_be_skipped(self):
        f = family(self.d, "OG0000007", alignment=self.ALN + ">Extra_1\nMKQL\n")
        self.assertEqual(schemas.validate_family(f), [])
        f = family(self.d, "OG0000008", alignment=">only\nMK\n")
        self.assertNotEqual(schemas.validate_family(f), [])
        self.assertEqual(schemas.validate_family(f, check_alignment=False), [])


class TestFamiliesFolder(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)
        self.fams = self.d / "Results" / "OrthoFinder2026" / "families"

    def tearDown(self):
        self._t.cleanup()

    def test_kind_is_inferred_from_the_folder_name(self):
        family(self.fams, "N0.HOG0000001")
        self.assertEqual(schemas.kind_of(self.fams), "families")
        self.assertEqual(schemas.kind_of(self.fams / "N0.HOG0000001"), "family")

    def test_problems_name_the_family(self):
        family(self.fams, "N0.HOG0000001")
        family(self.fams, "N0.HOG0000002", "((0_1,0_2),a);")
        family(self.fams, "N0.HOG0000003", done={"rooted": "yes"})
        write(self.fams / "family_summary.tsv", "x")                                   # loose files are ignored
        probs = schemas.validate(self.fams)
        self.assertEqual(len(probs), 3)
        self.assertTrue(all(p.startswith(("N0.HOG0000002: ", "N0.HOG0000003: ")) for p in probs))
        self.assertEqual(schemas.validate(self.fams, kind="families", prefixes={"Mzeb"}), probs)

    def test_good_empty_and_missing(self):
        family(self.fams, "OG0000001")
        family(self.fams, "OG0000002")
        self.assertEqual(schemas.validate(self.fams), [])
        empty = self.d / "empty" / "families"
        empty.mkdir(parents=True)
        self.assertIn("has no family folders", schemas.validate(empty)[0])
        self.assertIn("is not a folder", schemas.validate_families(self.d / "nothing")[0])

    def test_stops_after_many_problems(self):
        for i in range(40):
            family(self.fams, f"OG{i:07d}", "((0_1,0_2),a);")
        probs = schemas.validate(self.fams)
        self.assertEqual(len(probs), schemas.MAX_PROBLEMS + 1)
        self.assertIn("stopping after too many problems", probs[-1])


class TestCli(unittest.TestCase):
    def test_validate_families_with_and_without_catalog(self):
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            s = Settings(local_root=t / "data", cloud_remote="fake:")
            fams = s.local_root / "Results" / "OrthoFinder2026" / "families"
            family(fams, "N0.HOG0000001")
            family(fams, "N0.HOG0000002", "((Mzeb_XP_1,Acal_1),a);")
            f = io.StringIO()
            with redirect_stdout(f):
                cli(["validate", str(fams)], settings=s)                                # no catalog: prefix check off
            self.assertIn("no catalog.json here", f.getvalue())
            self.assertIn("families: ok", f.getvalue())
            write(s.local_root / "Inputs/Proteins/catalog.json", json.dumps(
                {"methods": {}, "compara_panel": [], "organisms": {
                    p: {"prefix": p, "species": p, "role": "project"} for p in ("Mzeb", "Acal")}}))
            f = io.StringIO()
            with redirect_stdout(f), self.assertRaises(SystemExit):
                cli(["validate", str(fams)], settings=s)
            self.assertIn("1 problem(s)", f.getvalue())
            self.assertIn("N0.HOG0000002: tree.nwk: leaf 'Mzeb_XP_1' starts with a species prefix", f.getvalue())
            self.assertNotIn("no catalog.json", f.getvalue())
            f = io.StringIO()
            with redirect_stdout(f):                                                    # a single good family
                cli(["validate", str(fams / "N0.HOG0000001"), "--kind", "family"], settings=s)
            self.assertIn("N0.HOG0000001: ok", f.getvalue())


if __name__ == "__main__":
    unittest.main()
