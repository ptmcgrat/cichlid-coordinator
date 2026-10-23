"""The names of the reference proteins of TreeFam and Compara: the reader, the validator, the coverage, and the check finding."""
import io
import os
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Layout, ReferenceNames, Settings, schemas
from cichlid_data.check import RECORDED, REQUIREMENTS, check
from cichlid_data.cli import main as cli
from cichlid_data.config import config_dir
from cichlid_data.reference_names import (HAS_GENE, HEADER, PARTS, RELIABLE, SPECIES, STATUSES, Record, coverage,
                                          table_path, validate_reference_names)
from helpers import write

HEAD = "\t".join(HEADER) + "\n"


def row(protein, species, status, route="", mapped="", transcript="", gene="", symbol="", source="", official="", description=""):
    return "\t".join([protein, species, status, route, mapped, transcript, gene, symbol, source, official, description]) + "\n"


# rows in the format of the real files (the first two human rows are from the TreeFam table, the zebrafish ones are made up)
TREEFAM = (HEAD
           + row("ENSP00000409362", "homo_sapiens", "legacy_id", "accession", "", "", "ENSG00000237613", "FAM138A", "HGNC", "HGNC:32334", "family with sequence similarity 138 member A")
           + row("ENSP00000329982", "homo_sapiens", "current", "ensembl", "ENSP00000329982", "ENST00000332831", "ENSG00000284662", "OR4F16", "HGNC", "HGNC:15079", "olfactory receptor family 4 subfamily F member 16")
           + row("ENSP00000111111", "homo_sapiens", "grch37", "grch37", "ENSP00000222222", "ENST00000333333", "ENSG00000444444", "ABC1", "HGNC", "HGNC:1")
           + row("ENSP00000555555", "homo_sapiens", "unmapped")
           + row("ENSDARP00000000001", "danio_rerio", "current", "ensembl", "ENSDARP00000000001", "ENSDART00000000001", "ENSDARG00000000001", "fixa", "ZFIN", "ZDB-GENE-000000-1", "a fixture gene")
           + row("ENSDARP00000000002", "danio_rerio", "replaced", "archive", "ENSDARP00000000009", "ENSDART00000000009", "ENSDARG00000000009", "fixb", "Ensembl", "ZDB-GENE-000000-2")
           + row("ENSDARP00000000003", "danio_rerio", "legacy_symbol", "symbol", "", "", "", "fixc", "legacy"))
def without(text, protein):
    """The table without the line of `protein`; fails loudly if there was no such line, so a test cannot pass vacuously."""
    lines = text.splitlines(True)
    kept = [ln for ln in lines if not ln.startswith(protein + "\t")]
    assert len(kept) == len(lines) - 1, f"no single line for {protein}"
    return "".join(kept)


MAPPING = ("ENSP00000409362\thomo_sapiens\nENSP00000329982\thomo_sapiens\nENSP00000111111\thomo_sapiens\nENSP00000555555\thomo_sapiens\n"
           "ENSDARP00000000001\tdanio_rerio\nENSDARP00000000002\tdanio_rerio\nENSDARP00000000003\tdanio_rerio\nENSMUSP00000000001\tmus_musculus\n"
           "ENSGALP00000000001\tgallus_gallus\n")


class Base(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L = Layout(self.root / "data", "fake:")
        write(table_path(self.L, "TreeFam"), TREEFAM)
        write(table_path(self.L, "Compara"), HEAD + row("ENSP00000000009", "homo_sapiens", "current", "ensembl", "ENSP00000000009", "T", "G", "X1", "HGNC", "HGNC:9"))

    def tearDown(self):
        self._t.cleanup()


class TestTheConstants(unittest.TestCase):
    def test_the_vocabulary(self):
        self.assertEqual(STATUSES, ("current", "grch37", "legacy_id", "replaced", "legacy_symbol", "unmapped"))
        self.assertEqual(RELIABLE, ("current", "grch37", "legacy_id"))
        self.assertEqual(SPECIES, ("homo_sapiens", "mus_musculus", "danio_rerio"))
        self.assertEqual(PARTS, ("TreeFam", "Compara"))
        self.assertEqual(len(HEADER), 11)
        self.assertTrue(set(RELIABLE) <= set(HAS_GENE) <= set(STATUSES))

    def test_the_folders_are_recorded_and_required_for_the_homology_part(self):
        for part in PARTS:
            folder = f"Inputs/Annotations/ReferenceNames/{part}"
            self.assertIn(folder, RECORDED)
            self.assertEqual(REQUIREMENTS["Homology"][folder], "dir")
        self.assertNotIn("Inputs/Annotations/ReferenceNames/TreeFam", REQUIREMENTS["TreeFam"])

    def test_the_source_records_exist_and_say_what_is_not_known(self):
        import tomllib
        d = {x["path"]: x for x in tomllib.loads((config_dir() / "sources.toml").read_text())["source"]}
        for part in PARTS:
            e = d[f"Inputs/Annotations/ReferenceNames/{part}"]
            self.assertEqual(e["type"], "derived")
            self.assertTrue(e["licence"].startswith("to confirm"))
            self.assertTrue(e["citation"].startswith("to confirm"))
            self.assertIn("stage6a_names.py", e["notes"])


class TestReader(Base):
    def test_get_and_the_fields(self):
        RN = ReferenceNames.load(self.L, "TreeFam")
        self.assertEqual(len(RN), 7)
        r = RN.get("ENSP00000409362")
        self.assertEqual((r.symbol, r.official_id, r.status, r.route, r.gene_id, r.symbol_source),
                         ("FAM138A", "HGNC:32334", "legacy_id", "accession", "ENSG00000237613", "HGNC"))
        self.assertEqual(r.mapped_protein, "")
        self.assertEqual(r.description, "family with sequence similarity 138 member A")
        self.assertIsInstance(r, Record)
        self.assertIn("ENSP00000329982", RN)
        self.assertIsNone(RN.get("ENSP_NOT_THERE"))
        self.assertEqual(RN.get("ENSDARP00000000001").official_id, "ZDB-GENE-000000-1")

    def test_reliable_and_named(self):
        RN = ReferenceNames.load(self.L, "TreeFam")
        got = {p: RN.get(p).reliable for p in RN.proteins()}
        self.assertEqual(got, {"ENSP00000409362": True, "ENSP00000329982": True, "ENSP00000111111": True, "ENSP00000555555": False,
                               "ENSDARP00000000001": True, "ENSDARP00000000002": False, "ENSDARP00000000003": False})
        self.assertTrue(RN.get("ENSP00000329982").named)
        self.assertFalse(RN.get("ENSP00000555555").named)

    def test_counts_and_filters(self):
        RN = ReferenceNames.load(self.L, "TreeFam")
        self.assertEqual(RN.status_counts("homo_sapiens"), {"legacy_id": 1, "current": 1, "grch37": 1, "unmapped": 1})
        self.assertEqual(RN.status_counts("danio_rerio"), {"current": 1, "replaced": 1, "legacy_symbol": 1})
        self.assertEqual(sum(RN.status_counts().values()), 7)
        self.assertEqual(sorted(RN.proteins("danio_rerio")), ["ENSDARP00000000001", "ENSDARP00000000002", "ENSDARP00000000003"])
        only = ReferenceNames.load(self.L, "TreeFam", only={"ENSP00000329982", "ENSDARP00000000001"})
        self.assertEqual(sorted(only.proteins()), ["ENSDARP00000000001", "ENSP00000329982"])

    def test_the_two_parts_are_separate(self):
        self.assertEqual(len(ReferenceNames.load(self.L, "Compara")), 1)
        self.assertIsNone(ReferenceNames.load(self.L, "Compara").get("ENSP00000409362"))

    def test_errors(self):
        with self.assertRaises(KeyError) as e:
            ReferenceNames.load(self.L, "OrthoFinder")
        self.assertIn("genes tables", str(e.exception))
        table_path(self.L, "Compara").unlink()
        with self.assertRaises(FileNotFoundError) as e:
            ReferenceNames.load(self.L, "Compara")
        self.assertIn("check --fetch", str(e.exception))
        write(table_path(self.L, "Compara"), "nonsense\n")
        with self.assertRaises(ValueError):
            ReferenceNames.load(self.L, "Compara")
        write(table_path(self.L, "Compara"), HEAD + row("A", "homo_sapiens", "current") + row("A", "homo_sapiens", "current"))
        with self.assertRaises(ValueError) as e:
            ReferenceNames.load(self.L, "Compara")
        self.assertIn("line 3: A appears twice", str(e.exception))
        write(table_path(self.L, "Compara"), HEAD + "A\thomo_sapiens\n")
        with self.assertRaises(ValueError) as e:
            ReferenceNames.load(self.L, "Compara")
        self.assertIn("expected 11 columns, found 2", str(e.exception))


class TestValidation(Base):
    def problems(self, text):
        return validate_reference_names(write(self.root / "reference_names.tsv", text))

    def test_good_and_empty_and_crlf(self):
        self.assertEqual(self.problems(TREEFAM), [])
        self.assertEqual(self.problems(HEAD), [])
        self.assertEqual(self.problems(TREEFAM.replace("\n", "\r\n")), [])
        self.assertEqual(self.problems(TREEFAM + "\n\n"), [])

    def test_every_kind_of_problem(self):
        cases = [
            (HEAD.replace("status", "state"), "header must be"),
            (HEAD + "A\thomo_sapiens\n", "line 2: expected 11 columns, found 2"),
            (HEAD + row("", "homo_sapiens", "current", gene="G"), "line 2: protein_id is empty"),
            (HEAD + row("A", "homo_sapiens", "current", gene="G") + row("A", "homo_sapiens", "current", gene="G"), "line 3: A appears twice"),
            (HEAD + row("A", "gallus_gallus", "current", gene="G"), "species 'gallus_gallus' must be one of"),
            (HEAD + row("A", "homo_sapiens", "fine", gene="G"), "status 'fine' must be one of"),
            (HEAD + row("A", "homo_sapiens", "current", gene="G", official="ZDB-GENE-000000-1"), "'ZDB-GENE-000000-1' is not an official ID of homo_sapiens"),
            (HEAD + row("A", "danio_rerio", "current", gene="G", official="HGNC:5"), "'HGNC:5' is not an official ID of danio_rerio"),
            (HEAD + row("A", "mus_musculus", "current", gene="G", official="MGI:abc"), "is not an official ID of mus_musculus"),
            (HEAD + row("A", "homo_sapiens", "current"), "A has status current but no gene_id"),
            (HEAD + row("A", "homo_sapiens", "replaced"), "has status replaced but no gene_id"),
        ]
        for text, want in cases:
            probs = self.problems(text)
            self.assertTrue(any(want in p for p in probs), (want, probs))

    def test_what_is_allowed(self):
        # unmapped and legacy_symbol rows need no gene; an empty official ID is fine; mouse and every ZFIN prefix are accepted
        ok = (HEAD + row("A", "homo_sapiens", "unmapped") + row("B", "homo_sapiens", "legacy_symbol", symbol="X")
              + row("C", "mus_musculus", "current", gene="G", official="MGI:12345")
              + row("D", "danio_rerio", "current", gene="G", official="ZDB-LINCRNAG-000000-12"))
        self.assertEqual(self.problems(ok), [])

    def test_stops_after_many_problems(self):
        probs = self.problems(HEAD + "".join(row(f"P{i}", "homo_sapiens", "fine") for i in range(60)))
        self.assertEqual(len(probs), schemas.MAX_PROBLEMS + 1)
        self.assertIn("stopping after too many problems", probs[-1])

    def test_through_validate_and_the_command_line(self):
        f = write(self.root / "x" / "reference_names.tsv", TREEFAM)
        self.assertEqual(schemas.kind_of(f), "reference_names")
        self.assertEqual(schemas.validate(f), [])
        s = Settings(local_root=self.root / "data", cloud_remote="fake:")
        out = io.StringIO()
        with redirect_stdout(out):
            cli(["validate", str(f)], settings=s)
        self.assertIn("reference_names.tsv: ok", out.getvalue())
        write(f, TREEFAM.replace("HGNC:15079", "HGNC:x"))
        out = io.StringIO()
        with redirect_stdout(out), self.assertRaises(SystemExit):
            cli(["validate", str(f), "--kind", "reference_names"], settings=s)
        self.assertIn("is not an official ID of homo_sapiens", out.getvalue())


class TestCoverage(Base):
    def setUp(self):
        super().setUp()
        write(self.L.family_resource("TreeFam9").local / "protein_to_species.txt", MAPPING)

    def test_counts_per_species(self):
        self.assertEqual(coverage(self.L, "TreeFam"),
                         {"homo_sapiens": (4, 4), "mus_musculus": (1, 0), "danio_rerio": (3, 3)})

    def test_a_protein_without_a_row_shows(self):
        write(table_path(self.L, "TreeFam"), without(TREEFAM, "ENSP00000555555"))
        self.assertEqual(coverage(self.L, "TreeFam")["homo_sapiens"], (4, 3))


class TestCheckFinding(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.L = Layout(Path(self._t.name) / "d", "fake:")
        write(self.L.family_resource("TreeFam9").local / "protein_to_species.txt", MAPPING)
        write(self.L.family_resource("Compara116").local / "prepared" / "protein_to_species.txt", "ENSP00000000009\thomo_sapiens\n")
        write(table_path(self.L, "TreeFam"), TREEFAM)
        write(table_path(self.L, "Compara"), HEAD + row("ENSP00000000009", "homo_sapiens", "current", gene="G", official="HGNC:9"))

    def tearDown(self):
        self._t.cleanup()

    def run_check(self):
        C = Catalog.build(self.L)
        return {f.item: f for f in check(self.L, C, release="v1", parts=["Homology"])}

    def test_all_covered(self):
        f = self.run_check()
        r = f["reference_names.tsv of TreeFam and Compara"]
        self.assertTrue(r.ok, r.detail)
        self.assertIn("TreeFam: homo 4 of 4, danio 3 of 3", r.detail)
        self.assertIn("Compara: homo 1 of 1, danio 0 of 0", r.detail)
        self.assertTrue(f["Inputs/Annotations/ReferenceNames/TreeFam"].ok)
        self.assertTrue(f["Inputs/Annotations/ReferenceNames/Compara"].ok)

    def test_a_missing_row_fails_and_says_how_many(self):
        write(table_path(self.L, "TreeFam"), without(TREEFAM, "ENSDARP00000000003"))
        r = self.run_check()["reference_names.tsv of TreeFam and Compara"]
        self.assertFalse(r.ok)
        self.assertIn("danio 2 of 3", r.detail)

    def test_an_invalid_table_fails_with_the_problem(self):
        write(table_path(self.L, "Compara"), HEAD + row("ENSP00000000009", "homo_sapiens", "current", gene="G", official="ZDB-GENE-000000-1"))
        r = self.run_check()["reference_names.tsv of TreeFam and Compara"]
        self.assertFalse(r.ok)
        self.assertIn("is not an official ID of homo_sapiens", r.detail)

    def test_missing_tables_fail(self):
        table_path(self.L, "TreeFam").unlink()
        f = self.run_check()
        self.assertFalse(f["reference_names.tsv of TreeFam and Compara"].ok)
        self.assertFalse(f["Inputs/Annotations/ReferenceNames/TreeFam"].ok)
        self.assertTrue(f["Inputs/Annotations/ReferenceNames/Compara"].ok)


if __name__ == "__main__":
    unittest.main()
