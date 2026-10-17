import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Settings, schemas
from cichlid_data.cli import main as cli
from helpers import Env, write


def gene(gid, **kw):
    r = {"gene_id": gid, "seqid": "NC_1", "start": 10, "end": 90, "strand": "+", "mitochondrial": False,
         "n_transcripts": 1, "n_proteins": 1, "chosen_protein": "XP_1", "protein_length": 5,
         "partial": False, "low_quality": False, "residues_replaced": 0, "internal_stops": 0}
    r.update(kw)
    return r


class TestTables(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_write_read_round_trip(self):
        p = self.d / "Mzeb.genes.tsv"
        schemas.GENES.write(p, [gene("Mzeb_1"), gene("Mzeb_2.2", strand="-")], prefix="Mzeb")
        rows = schemas.GENES.read(p)
        self.assertEqual((rows[0]["mitochondrial"], rows[1]["gene_id"], rows[1]["strand"]), ("False", "Mzeb_2.2", "-"))
        self.assertEqual(schemas.validate(p), [])

    def test_row_problems(self):
        rows = [gene("Mzeb_1", start="ten"), gene("Mzeb_1"), gene("Mzeb_bad_id"), gene("Acal_3"),
                gene("Mzeb_4", partial="yes", strand="up"), gene("Mzeb_5", seqid="")]
        probs = "\n".join(schemas.GENES.check_rows(rows, prefix="Mzeb"))
        for bit in ("start is not an integer", "duplicate gene_id 'Mzeb_1'", "'Mzeb_bad_id' is not a gene ID",
                    "'Acal_3' does not start with the file's prefix Mzeb_", "partial is not True/False",
                    "strand 'up' is not one of", "seqid is empty"):
            self.assertIn(bit, probs)
        with self.assertRaises(schemas.SchemaError):
            schemas.GENES.write(self.d / "x.genes.tsv", rows)
        self.assertFalse((self.d / "x.genes.tsv").exists())                     # nothing written

    def test_header_problems(self):
        p = self.d / "Mzeb.assign.tsv"
        cols = schemas.ASSIGN.names
        p.write_text("\t".join(cols[:3] + ["extra"]) + "\n")
        self.assertIn("missing columns", schemas.validate(p)[0])
        p.write_text("\t".join([cols[1], cols[0]] + cols[2:]) + "\n")
        self.assertIn("wrong order", schemas.validate(p)[0])
        with self.assertRaises(schemas.SchemaError):
            schemas.ASSIGN.read(p)

    def test_orthologs_any_reference_set(self):
        S = schemas.orthologs_schema(("human", "zebrafish"))
        row = {"gene_id": "Mzeb_1", "human_relationship": "Ortholog", "human_ortholog_type": "OneToTwo",
               "zebrafish_relationship": "NoHomolog", "confidence": "high", "placement_lwr": "0.93"}
        p = self.d / "Mzeb.orthologs.tsv"
        S.write(p, [row])
        self.assertEqual(schemas.refs_in_header(p.read_text().split("\n")[0].split("\t")), ["human", "zebrafish"])
        self.assertEqual(schemas.validate(p), [])
        row["human_relationship"] = "Sibling"
        row["confidence"] = "certain"
        probs = S.check_rows([row])
        self.assertTrue(any("human_relationship 'Sibling'" in x for x in probs))
        self.assertTrue(any("confidence 'certain'" in x for x in probs))
        self.assertEqual(schemas.ORTHOLOGS.names[15:17], ["human_relationship", "human_ortholog_type"])


class TestOtherKinds(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_proteome(self):
        p = write(self.d / "Mzeb.faa", ">Mzeb_1\nMKLV\n>Mzeb_2\nMK*L\n>Acal_3\nMK\n>Mzeb_1\nM\n>Mzeb_4\n\n")
        probs = "\n".join(schemas.validate(p))
        for bit in ("Mzeb_2: characters outside the amino-acid alphabet: *", "'Acal_3' does not start with Mzeb_",
                    "duplicate protein 'Mzeb_1'", "Mzeb_4: empty sequence"):
            self.assertIn(bit, probs)
        # cross-check with the genes table beside it
        q = write(self.d / "Hbur.faa", ">Hbur_1\nMK\n>Hbur_9\nMK\n")
        schemas.GENES.write(self.d / "Hbur.genes.tsv", [gene("Hbur_1"), gene("Hbur_2"), gene("Hbur_3", chosen_protein="")])
        probs = schemas.validate(q)
        self.assertEqual(probs, ["Hbur_2: has a chosen protein in the genes table but none in the proteome",
                                 "Hbur_9: in the proteome but has no chosen protein in the genes table"])

    def test_run_and_family(self):
        run = self.d / "run.json"
        with self.assertRaises(schemas.SchemaError):
            schemas.write_run(run, method="OrthoFinder")
        schemas.write_run(run, method="OrthoFinder", method_version="3.0", parameters={"-M": "msa"}, started="t0",
                          finished="t1", coordinator_version="0.4.0", catalog_md5="abc", species=["Mzeb", "Hsap"])
        self.assertEqual(schemas.validate(run), [])
        run.write_text(json.dumps({"method": "x", "species": ["bad_prefix"]}))
        self.assertIn("species must be a list of prefixes", schemas.validate(run))
        fam = self.d / "families" / "TF101001"
        fam.mkdir(parents=True)
        self.assertEqual(schemas.validate(fam), ["tree.nwk is missing", "done.json is missing (the family is not finished)",
                                                 "combined.aln is missing"])
        (fam / "tree.nwk").write_text("((a,b),c);")
        (fam / "done.json").write_text("{}")
        self.assertEqual(schemas.validate(fam), ["combined.aln is missing"])
        (fam / "combined.aln").write_text(">a\nMK\n>b\nMK\n>c\nMK\n")
        self.assertEqual(schemas.validate(fam), [])
        (fam / "tree.nwk").write_text("((a,b),c;")
        self.assertEqual(schemas.validate(fam), ["tree.nwk has unbalanced parentheses"])

    def test_cli(self):
        e = Env()
        try:
            good = self.d / "Mzeb.genes.tsv"
            schemas.GENES.write(good, [gene("Mzeb_1")])
            bad = write(self.d / "Acal.faa", ">Acal_1\nMK*\n")
            s = Settings(local_root=e.a.local_root, cloud_remote="fake:")
            f = io.StringIO()
            with redirect_stdout(f), self.assertRaises(SystemExit):
                cli(["validate", str(good), str(bad)], backend=e.store(e.a).backend, settings=s)
            out = f.getvalue()
            self.assertIn("Mzeb.genes.tsv: ok", out)
            self.assertIn("Acal.faa: 1 problem(s)", out)
        finally:
            e.close()


if __name__ == "__main__":
    unittest.main()


