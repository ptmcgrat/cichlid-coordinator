"""Offline tests for prepare_local_species.py, checked against stage 2's own parser."""
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

from cichlid_data.downloaders import local as pl
from cichlid_data import proteomes as s2

GENOME = "A" * 9 + "ATGAAACCCGGGTTTAAACCCGGGTTTAAA" + "C" * 100 + "CCCGGGTTTAAACCCGGGTTTAAACCCTAA" + "G" * 230


def setup(d):
    cds1 = GENOME[9:39] + GENOME[139:169]                     # two exons, 60 nt
    with gzip.open(d / "Genes.fasta.gz", "wt") as fh:
        fh.write(f">evm.model.Chr_1.1 F-box only protein 48 [Archocentrus centrarchus]\n{cds1}\n")
        fh.write(">evm.model.Chr_1.2 something\nATGTGGTGGTAA\n")              # will not align
        fh.write(">evm.model.Chr_1.3 broken\nATGAAGTAGCTATCCAAGGATGGATACAAAGCCATAGGCATTAATGACGTACTTTAGACA\n")  # excluded
    with gzip.open(d / "Genome.fa.gz", "wt") as fh:
        fh.write(f">Chr_1\n{GENOME}\n>Scaffold_9\n{'N' * 50}\n")
    (d / "func.gff").write_text(
        "evm.model.Chr_1.1\tBlast2GO\tCDS\t1\t60\t.\t.\t.\tID=evm.model.Chr_1.1;Description=PREDICTED: F-box only protein 48 [Archocentrus centrarchus];Gene=FBXO48\n"
        "evm.model.Chr_1.2\tBlast2GO\tCDS\t1\t12\t.\t.\t.\tID=evm.model.Chr_1.2;Description=---NA---;Gene=W5NGD6\n")
    fake = d / "minimap2"
    fake.write_text("#!/bin/sh\nprintf 'Chr1-1\\t60\\t0\\t60\\t+\\tChr_1\\t399\\t9\\t169\\t60\\t60\\t60\\ttp:A:P\\tcg:Z:30M100N30M\\n'\n")
    fake.chmod(0o755)
    return fake


class TestPrepareLocal(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(pl.clean_id("evm.model.Chr_12.345"), "Chr12-345")
        self.assertEqual(pl.translate("ATGAAATAA"), "MK")
        self.assertEqual(pl.exons_from_cg(9, "30M100N30M"), [(10, 39), (140, 169)])
        self.assertEqual(pl.exons_from_cg(0, "10M2I5M3D10M"), [(1, 28)])            # insertion has no span

    def test_protein_keys(self):
        self.assertEqual(pl.protein_keys("protein|evm.model.Chr_1.1243", "ID=evm.model.Chr_1.1243|Parent=evm.TU.Chr_1.1243"),
                         ["protein|evm.model.Chr_1.1243", "evm.model.Chr_1.1243", "evm.model.Chr_1.1243"])

    def test_protein_repair(self):
        body = "GCT" * 30                                                              # 30 alanines
        self.assertEqual(pl.protein_for("ATG" + body + "TAA"), ("M" + "A" * 30, "clean"))
        # partial model starting mid-codon: frame 1 has stops, frame 2 is clean
        self.assertEqual(pl.protein_for("C" + "ATG" + "AAT" * 20), ("M" + "N" * 20, "reframed (frame 2)"))
        # runs on past the true stop (late stop, then non-coding sequence)
        p, how = pl.protein_for("ATG" + body + "TAG" + "CCCTGA")
        self.assertEqual((p, how), ("M" + "A" * 30, "cut at first stop"))
        # early stop and no clean frame: left out
        bad = "ATGAAGTAGCTATCCAAGGATGGATACAAAGCCATAGGCATTAATGACGTACTTTAGACA"   # stops in all six frames
        self.assertEqual(pl.protein_for(bad), (None, "first stop too early, no clean frame"))

    def test_prepare_and_stage2(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            fake = setup(d)
            out = d / "prepared"
            pl.main(["--prefix", "Esur", "--cds", str(d / "Genes.fasta.gz"), "--genome", str(d / "Genome.fa.gz"),
                     "--functional", str(d / "func.gff"), "--outdir", str(out), "--minimap2", str(fake)])
            summ = json.loads((out / "prepare_summary.json").read_text())
            self.assertEqual((summ["genes"], summ["aligned"], summ["unplaced"], summ["symbols"], summ["excluded"]), (2, 1, 1, 1, 1))
            self.assertIn("Chr1-3\tevm.model.Chr_1.3\tfirst stop too early", (out / "excluded.tsv").read_text())
            self.assertNotIn("Chr1-3", (out / "protein.faa").read_text())
            rep = [json.loads(l) for l in (out / "sequence_report.jsonl").read_text().splitlines()]
            roles = {r["sequenceName"]: (r["role"], r.get("chrName")) for r in rep}
            self.assertEqual(roles["Chr_1"], ("assembled-molecule", "1"))
            self.assertEqual(roles["Scaffold_9"][0], "unplaced-scaffold")
            self.assertEqual(roles["unplaced_genes"][0], "unplaced-scaffold")
            # stage 2 reads the result like any annotation
            ann = s2.Annotation()
            ann.load_gff3(out / "genomic.gff")
            fasta = s2.read_fasta(out / "protein.faa")
            res, stats = s2.select_proteins(ann, fasta)
            got = {s2.output_id("Esur", g): g for g in res}
            self.assertEqual(sorted(got), ["Esur_Chr1-1", "Esur_Chr1-2"])
            g1 = got["Esur_Chr1-1"]
            self.assertEqual((g1["seqid"], g1["start"], g1["end"], g1["symbol"]), ("Chr_1", 10, 169, "FBXO48"))
            self.assertEqual(g1["description"], "F-box only protein 48")
            self.assertEqual(g1["chosen"]["seq"], "MKPGFKPGFKPGFKPGFKP")                 # exons joined, stop removed
            g2 = got["Esur_Chr1-2"]
            self.assertEqual((g2["seqid"], g2["symbol"]), ("unplaced_genes", "Chr1-2"))  # UniProt accession not a symbol


class TestSuppliedProteins(unittest.TestCase):
    def test_supplied_proteins(self):
        import io, contextlib
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            fake = setup(d)
            # Chr_1.1 under its original ID (trailing '*'); Chr_1.3 under the cleaned ID, with an internal
            # stop; Chr_1.2 missing
            # header in the authors' style: 'protein|' prefix and an ID= field
            (d / "proteins.faa").write_text(">protein|evm.model.Chr_1.1 ID=evm.model.Chr_1.1|Parent=evm.TU.Chr_1.1|Name=EVM prediction\n"
                                           "MKPGFKPGFK\nPGFKPGFKP*\n>Chr1-3\nMKX*Q\n")
            out = d / "prepared"
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                pl.main(["--prefix", "Esur", "--cds", str(d / "Genes.fasta.gz"), "--genome", str(d / "Genome.fa.gz"),
                         "--functional", str(d / "func.gff"), "--outdir", str(out), "--minimap2", str(fake),
                         "--proteins", str(d / "proteins.faa")])
            prots = {n: s for n, _, s in pl.read_fasta(out / "protein.faa")}
            self.assertEqual(prots, {"Chr1-1": "MKPGFKPGFKPGFKPGFKP"})
            ex = (out / "excluded.tsv").read_text()
            self.assertIn("Chr1-2\tevm.model.Chr_1.2\tnot in the protein file", ex)
            self.assertIn("Chr1-3\tevm.model.Chr_1.3\tinternal stop in the supplied protein", ex)
            summ = json.loads((out / "prepare_summary.json").read_text())
            self.assertEqual(summ["comparison_with_cds_translation"], {"ours clean: identical": 1})
            self.assertIn("Read 2 proteins", err.getvalue())


class TestStage1Local(unittest.TestCase):
    def test_copies_sequence_report(self):
        from cichlid_data.downloaders import ncbi as s1
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "src").mkdir()
            (d / "src" / "genomic.gff").write_text("##gff-version 3\n")
            (d / "src" / "protein.faa").write_text(">x\nM\n")
            (d / "src" / "sequence_report.jsonl").write_text('{"sequenceName": "Chr_1", "role": "assembled-molecule"}\n')
            (d / "out").mkdir()
            files = s1.process_local({"annotation_format": "gff3", "gff_path": str(d / "src" / "genomic.gff"),
                                      "protein_path": str(d / "src" / "protein.faa")}, d / "out")
            self.assertIn("sequence_report", files)
            self.assertTrue((d / "out" / "sequence_report.jsonl").is_file())


if __name__ == "__main__":
    unittest.main()
