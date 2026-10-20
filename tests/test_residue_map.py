"""residue_map.tsv (the placement parts' per-protein map to alignment columns), tree_mode, and proteomes_md5."""
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Layout, Settings, schemas
from cichlid_data.cli import main as cli
from helpers import write
from test_catalog import MANIFEST, SETS, setup as catalog_setup

HEADER = "gene_id\tprotein_length\taln_length\tresidues_not_in_alignment\tcolumns\n"


def rmap(*rows):
    return HEADER + "".join("\t".join(str(x) for x in r) + "\n" for r in rows)


GOOD = rmap(("Acal_1", 4, 10, 1, "2,3,-,7"), ("Mzeb_5", 3, 10, 0, "1,5,10"), ("Mzeb_9", 2, 10, 2, "-,-"))


class TestValidateResidueMap(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def problems(self, text, aln_length=None):
        return schemas.validate_residue_map(write(self.d / "residue_map.tsv", text), aln_length)

    def test_good(self):
        self.assertEqual(self.problems(GOOD), [])
        self.assertEqual(self.problems(GOOD, aln_length=10), [])
        self.assertEqual(self.problems(HEADER), [])                              # no rows: nothing to complain about
        self.assertEqual(self.problems(GOOD.replace("\n", "\r\n")), [])            # Windows line ends
        self.assertEqual(self.problems(GOOD + "\n\n"), [])                       # blank lines are ignored

    def test_kind_is_inferred_from_the_name(self):
        f = write(self.d / "residue_map.tsv", GOOD)
        self.assertEqual(schemas.kind_of(f), "residue_map")
        self.assertEqual(schemas.validate(f), [])
        write(f, "nonsense")
        self.assertEqual(schemas.validate(f), ["header must be gene_id, protein_length, aln_length, "
                                               "residues_not_in_alignment, columns"])
        self.assertEqual(schemas.validate(self.d / "missing" / "residue_map.tsv")[0][:11], "[Errno 2] N")

    def test_every_kind_of_row_problem(self):
        cases = [
            (rmap(("Acal_1", 3, 10, 0, "1,2")), "columns has 2 entries for a protein of length 3"),
            (rmap(("Acal_1", 2, 10, 0, "1,2,3")), "columns has 3 entries for a protein of length 2"),
            (rmap(("Acal_1", 3, 10, 1, "1,2,3")), "residues_not_in_alignment is 1 but columns has 0 '-' entries"),
            (rmap(("Acal_1", 3, 10, 0, "1,-,3")), "residues_not_in_alignment is 0 but columns has 1 '-' entries"),
            (rmap(("Acal_1", 3, 10, 0, "1,3,3")), "columns must increase along the sequence (3 after 3)"),
            (rmap(("Acal_1", 3, 10, 0, "1,5,4")), "columns must increase along the sequence (4 after 5)"),
            (rmap(("Acal_1", 3, 10, 0, "1,2,11")), "column '11' must be '-' or a whole number from 1 to aln_length (10)"),
            (rmap(("Acal_1", 3, 10, 0, "0,2,3")), "column '0' must be"),
            (rmap(("Acal_1", 3, 10, 0, "1,x,3")), "column 'x' must be"),
            (rmap(("Acal_1", 3, 10, 0, "1,2.5,3")), "column '2.5' must be"),
            (rmap(("Acal_1", "three", 10, 0, "1,2,3")), "must be whole numbers"),
            (rmap(("Acal_1", 0, 10, 0, "")), "must be at least 1"),
            (rmap(("Acal 1", 3, 10, 0, "1,2,3")), "is not a gene ID"),
            (rmap(("ENSDARP0000001", 3, 10, 0, "1,2,3")), "is not a gene ID"),
            (HEADER + "Acal_1\t3\t10\n", "line 2: expected 5 columns, found 3"),
        ]
        for text, want in cases:
            probs = self.problems(text)
            self.assertEqual(len(probs), 1, (text, probs))
            self.assertIn(want, probs[0])
            self.assertTrue(probs[0].startswith("line 2: "))

    def test_rows_are_sorted_and_unique(self):
        probs = self.problems(rmap(("Mzeb_5", 1, 10, 0, "1"), ("Acal_1", 1, 10, 0, "1")))
        self.assertEqual(probs, ["line 3: Acal_1 is not after Mzeb_5 (rows are sorted by gene_id, with no duplicates)"])
        probs = self.problems(rmap(("Acal_1", 1, 10, 0, "1"), ("Acal_1", 1, 10, 0, "2")))
        self.assertIn("is not after Acal_1", probs[0])

    def test_aln_length_must_agree(self):
        probs = self.problems(rmap(("Acal_1", 1, 10, 0, "1"), ("Mzeb_5", 1, 12, 0, "1")))
        self.assertEqual(probs, ["aln_length differs between rows (10 to 12)"])
        self.assertEqual(self.problems(GOOD, aln_length=11), ["aln_length is 10 but combined.aln is 11 columns long"])

    def test_header_must_be_exact(self):
        for bad in (HEADER.replace("columns", "cols"), "gene_id\tprotein_length\n", HEADER.replace("\t", ",")):
            self.assertIn("header must be", self.problems(bad + "Acal_1\t1\t10\t0\t1\n")[0])

    def test_a_giant_protein_exceeds_csv_defaults(self):
        n = 200_000                                                              # titin-sized; csv's limit is 131,072
        cols = ",".join(str(i) for i in range(1, n + 1))
        self.assertEqual(self.problems(rmap(("Acal_1", n, n, 0, cols)), aln_length=n), [])

    def test_stops_after_many_problems(self):
        probs = self.problems(rmap(*[(f"Acal_{i:03d}", 3, 10, 0, "1,2") for i in range(50)]))
        self.assertEqual(len(probs), schemas.MAX_PROBLEMS + 1)
        self.assertIn("stopping after too many problems", probs[-1])


def family(parent, name="TF101001", done=None, aln=">a\nMKQLVV\n>b\nMKQLVV\n>c\nMKQLVV\n"):
    d = Path(parent) / name
    write(d / "tree.nwk", "((a,b),c);")
    write(d / "combined.aln", aln)
    write(d / "done.json", json.dumps({} if done is None else done))
    return d


class TestInAFamilyFolder(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_checked_against_combined_aln_when_present(self):
        f = family(self.d)
        self.assertEqual(schemas.validate(f), [])                                # no map: fine, it is optional
        write(f / "residue_map.tsv", rmap(("Mzeb_1", 3, 6, 0, "1,2,6")))
        self.assertEqual(schemas.validate(f), [])
        write(f / "residue_map.tsv", rmap(("Mzeb_1", 3, 7, 0, "1,2,6")))
        self.assertEqual(schemas.validate(f), ["residue_map.tsv: aln_length is 7 but combined.aln is 6 columns long"])
        self.assertEqual(schemas.validate_family(f, check_alignment=False), [])      # alignment check skipped
        write(f / "residue_map.tsv", rmap(("Mzeb_1", 3, 6, 0, "1,2")))
        self.assertIn("residue_map.tsv: line 2: Mzeb_1: columns has 2 entries", schemas.validate(f)[0])
        write(f / "residue_map.tsv", rmap(("Mzeb_1", 3, 6, 0, "1,2,6")))
        write(f / "combined.aln", ">a\nMKQLVV\n>b\nMKQ\n>c\nMKQLVV\n")        # a ragged alignment: no length to compare
        probs = schemas.validate(f)
        self.assertEqual(len(probs), 1)
        self.assertIn("rows have different lengths", probs[0])

    def test_tree_mode(self):
        for mode in ("place", "denovo"):
            self.assertEqual(schemas.validate(family(self.d, f"TF{mode}", {"tree_mode": mode})), [])
        for bad in ("de novo", "Place", "orthofinder", "", None, 1):
            f = family(self.d, "TF999999", {"tree_mode": bad})
            self.assertEqual(schemas.validate(f), [f"done.json: tree_mode {bad!r} must be one of place, denovo"], bad)
            for x in f.iterdir():
                x.unlink()
            f.rmdir()
        # the checks of a done.json are independent of each other
        probs = schemas.validate(family(self.d, "TF888888", {"tree_mode": "x", "rooted": "yes", "release": "latest"}))
        self.assertEqual(len(probs), 3)

    def test_cli_validates_a_residue_map(self):
        s = Settings(local_root=self.d / "data", cloud_remote="fake:")
        f = write(self.d / "x" / "residue_map.tsv", GOOD)
        out = io.StringIO()
        with redirect_stdout(out):
            cli(["validate", str(f)], settings=s)
        self.assertIn("residue_map.tsv: ok", out.getvalue())
        write(f, GOOD.replace("2,3,-,7", "2,3,-,2"))
        out = io.StringIO()
        with redirect_stdout(out), self.assertRaises(SystemExit):
            cli(["validate", str(f), "--kind", "residue_map"], settings=s)
        self.assertIn("columns must increase along the sequence (2 after 3)", out.getvalue())


class TestProteomesMd5(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L, self.cfg = catalog_setup(self.root)

    def tearDown(self):
        self._t.cleanup()

    def test_stable_across_rebuilds_and_independent_of_order(self):
        import time
        C = Catalog.build(self.L, self.cfg)
        C.write(self.L)
        m1, f1 = C.proteomes_md5(["Mzeb"]), C.md5()
        time.sleep(1.1)
        C2 = Catalog.build(self.L, self.cfg)
        C2.write(self.L)
        self.assertNotEqual(C2.md5(), f1)                                        # the file checksum moves on every rebuild
        self.assertEqual(C2.proteomes_md5(["Mzeb"]), m1)                         # the proteomes' checksum does not
        self.assertRegex(m1, r"^[0-9a-f]{32}$")
        self.assertEqual(C2.proteomes_md5(), m1)                                 # only Mzeb has a proteome in this fixture

    def test_it_follows_the_proteomes_and_the_species_asked_for(self):
        from test_check import proteome
        proteome(self.L, "Acal", seq="MKLVAAA")
        C = Catalog.build(self.L, self.cfg)
        both = C.proteomes_md5(["Mzeb", "Acal"])
        self.assertEqual(both, C.proteomes_md5(["Acal", "Mzeb"]))
        self.assertNotEqual(both, C.proteomes_md5(["Mzeb"]))
        self.assertNotEqual(C.proteomes_md5(["Acal"]), C.proteomes_md5(["Mzeb"]))
        write(self.L.proteome("Acal").local, ">Acal_1\nMKLVAAAAAAAA\n")            # a better annotation
        self.assertNotEqual(Catalog.build(self.L, self.cfg).proteomes_md5(["Acal"]), C.proteomes_md5(["Acal"]))
        self.assertEqual(Catalog.build(self.L, self.cfg).proteomes_md5(["Mzeb"]), C.proteomes_md5(["Mzeb"]))

    def test_a_species_without_a_proteome_is_refused(self):
        C = Catalog.build(self.L, self.cfg)
        with self.assertRaises(KeyError) as e:
            C.proteomes_md5(["Mzeb", "Drer", "Nope"])
        self.assertIn("Drer, Nope", str(e.exception))

    def test_the_recipe_is_documented_and_checkable_by_hand(self):
        import hashlib
        C = Catalog.build(self.L, self.cfg)
        want = hashlib.md5(f"Mzeb {C['Mzeb'].proteome_md5}\n".encode()).hexdigest()
        self.assertEqual(C.proteomes_md5(["Mzeb"]), want)

    def test_run_record_accepts_and_checks_the_field(self):
        run = self.root / "Results" / "v1" / "TreeFam" / "run.json"
        rec = dict(release="v1", method="TreeFam", method_version="x", parameters={}, started="a", finished="b",
                   coordinator_version="0.10.0", catalog_md5="m", species=["Mzeb"])
        schemas.write_run(run, **rec)
        self.assertEqual(schemas.validate(run), [])                              # optional
        schemas.write_run(run, proteomes_md5=Catalog.build(self.L, self.cfg).proteomes_md5(["Mzeb"]), **rec)
        self.assertEqual(schemas.validate(run), [])
        schemas.write_run(run, proteomes_md5="not a checksum", **rec)
        self.assertEqual(schemas.validate(run), ["proteomes_md5 must be a 32-character MD5 (Catalog.proteomes_md5)"])


if __name__ == "__main__":
    unittest.main()
