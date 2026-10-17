"""Offline tests for the Compara preparation (Stage 3c's tests stay in the pipeline)."""
import csv
import gzip
import io
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from cichlid_data.downloaders import compara as cp
from cichlid_data import newick


def _text(path):
    return io.StringIO(Path(path).read_text())


PANEL = """# test panel
species\trole\tgroup\tnote
homo_sapiens\treference\treference\t
mus_musculus\treference\treference\t
danio_rerio\treference\treference\t
oryzias_latipes\tpanel\tovalentaria\t
lepisosteus_oculatus\tpanel\tnon-teleost fish\t
maylandia_zebra\tvalidation\tcichlid\t
salmo_salar\texclude\tfish\t
"""

T1 = ("((((ENSORLP1:0.1[&&NHX:D=N:G=ENSORLG1:T=8090],ENSMZEP1:0.1[&&NHX:D=N:G=ENSMZEG1:T=106582])"
      ":0.1[&&NHX:D=N:B=90:T=1],ENSSSAP1:0.3[&&NHX:D=N:G=ENSSSAG1:T=8030]):0.1[&&NHX:D=Y:B=68:T=1],"
      "ENSDARP1:0.2[&&NHX:D=N:G=ENSDARG1:T=7955]):0.1,(ENSLOCP1:0.2[&&NHX:D=N:G=ENSLOCG1:T=7918],"
      "(ENSP1:0.1[&&NHX:D=N:G=ENSG1:T=9606],ENSMUSP1:0.1[&&NHX:D=N:G=ENSMUSG1:T=10090]):0.1):0.1);")
SEQ = {"ENSORLP1": "MKV-LLAA", "ENSMZEP1": "MKVALLAA", "ENSSSAP1": "MKVQLLAAW", "ENSDARP1": "MKV-LLAG",
       "ENSLOCP1": "MRV-LLAA", "ENSP1": "MKV-LLSA", "ENSMUSP1": "MKV-LLSA"}
SPECIES = {"ENSORLP1": "oryzias_latipes", "ENSMZEP1": "maylandia_zebra", "ENSSSAP1": "salmo_salar",
           "ENSDARP1": "danio_rerio", "ENSLOCP1": "lepisosteus_oculatus", "ENSP1": "homo_sapiens",
           "ENSMUSP1": "mus_musculus"}


def seq_line(p, sp, label="NULL"):
    return f"SEQ {sp} {p} 1 100 200 1 {p.replace('P', 'G', 1)} {label}\n"


def make_dumps(d):
    nhx, aln = [], []
    # family 1: normal; salmon column 9 is the only residue there -> dropped after pruning
    for p in SEQ:
        nhx.append(seq_line(p, SPECIES[p], "BMP2" if p == "ENSP1" else "NULL"))
    nhx += ["DATA\n", T1 + "\n", "//\n"]
    for p, s in SEQ.items():
        aln.append(f">{p}\n{s.ljust(9, '-')}\n")
    aln.append("//\n")
    # family 2: only excluded species
    nhx += [seq_line("ENSSSAP2", "salmo_salar"), seq_line("ENSSSAP3", "salmo_salar"), "DATA\n",
            "(ENSSSAP2:0.1[&&NHX:D=N:T=8030],ENSSSAP3:0.1[&&NHX:D=N:T=8030]);\n", "//\n"]
    aln += [">ENSSSAP2\nMK\n>ENSSSAP3\nMK\n//\n"]
    # family 3: tree and alignment disagree
    nhx += [seq_line("ENSP3", "homo_sapiens"), seq_line("ENSDARP3", "danio_rerio"), "DATA\n",
            "(ENSP3:0.1,ENSDARP3:0.1);\n", "//\n"]
    aln += [">ENSP3\nMK\n>ENSDARP9\nMK\n//\n"]
    # family 4: tiny (2 kept leaves)
    nhx += [seq_line("ENSP4", "homo_sapiens"), seq_line("ENSDARP4", "danio_rerio"), "DATA\n",
            "(ENSP4:0.1[&&NHX:T=9606],ENSDARP4:0.1[&&NHX:T=7955]);\n", "//\n"]
    aln += [">ENSP4\nMKLLW\n>ENSDARP4\nMKLLW\n//\n"]
    # family 5: only Ensembl's M. zebra copies (validation) -> cichlid_only
    nhx += [seq_line("ENSMZEP5", "maylandia_zebra"), seq_line("ENSMZEP6", "maylandia_zebra"), "DATA\n",
            "(ENSMZEP5:0.1[&&NHX:T=106582],ENSMZEP6:0.1[&&NHX:T=106582]);\n", "//\n"]
    aln += [">ENSMZEP5\nMKQ\n>ENSMZEP6\nMKR\n//\n"]
    with gzip.open(d / "Compara.116.protein_default.nhx.emf.gz", "wt") as fh:
        fh.writelines(nhx)
    with gzip.open(d / "Compara.116.protein_default.aa.fasta.gz", "wt") as fh:
        fh.writelines(aln)
    (d / "panel.tsv").write_text(PANEL)
    return d / "Compara.116.protein_default.nhx.emf.gz", d / "Compara.116.protein_default.aa.fasta.gz"


class TestPrepare(unittest.TestCase):
    def run_prepare(self, d, *extra):
        nhx, aln = make_dumps(d) if not (d / "panel.tsv").exists() else (
            d / "Compara.116.protein_default.nhx.emf.gz", d / "Compara.116.protein_default.aa.fasta.gz")
        cp.main(["--nhx", str(nhx), "--aln", str(aln), "--panel", str(d / "panel.tsv"),
                 "--outdir", str(d / "prep"), *extra])
        return d / "prep"

    def test_split_prune_and_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            out = self.run_prepare(Path(d))
            audit = {r["family"]: r for r in csv.DictReader(_text(out / "family_audit.tsv"), delimiter="\t")}
            self.assertEqual({k: v["status"] for k, v in audit.items()},
                             {"C116_000001": "ok", "C116_000002": "no_panel_members",
                              "C116_000003": "mismatch", "C116_000004": "tiny_tree",
                              "C116_000005": "cichlid_only"})
            self.assertEqual((out / "trees" / "C116_000005.aln").read_text(), "")   # de novo marker
            self.assertIn(">ENSMZEP5\nMKQ\n", (out / "search_db.fa").read_text())
            self.assertEqual(audit["C116_000001"]["members_kept"], "5")
            self.assertEqual(audit["C116_000001"]["validation_members"], "1")

            tree = newick.parse((out / "trees" / "C116_000001.nhx").read_text())
            leaves = sorted(l.name for l in tree.leaves())
            self.assertEqual(leaves, ["ENSDARP1", "ENSLOCP1", "ENSMUSP1", "ENSORLP1", "ENSP1"])
            self.assertNotIn("NHX", (out / "trees" / "C116_000001.nhx").read_text())
            aln = {}
            for line in (out / "trees" / "C116_000001.aln").read_text().split(">")[1:]:
                n, s = line.split("\n", 1)
                aln[n] = s.strip()
            self.assertEqual(sorted(aln), leaves)
            # columns with residues only in removed proteins (M. zebra copy at 4, salmon at 9) dropped
            self.assertEqual({len(s) for s in aln.values()}, {7})
            self.assertEqual(aln["ENSORLP1"], "MKVLLAA")

            members = {r["protein"]: r for r in csv.DictReader(_text(out / "members.tsv"), delimiter="\t")}
            self.assertEqual(members["ENSMZEP1"]["role"], "validation")
            self.assertNotIn("ENSSSAP1", members)
            self.assertEqual((members["ENSP1"]["label"], members["ENSP1"]["gene"], members["ENSP1"]["family"]),
                             ("BMP2", "ENSG1", "C116_000001"))
            db = (out / "search_db.fa").read_text()
            self.assertIn(">ENSMZEP1\nMKVALLAA\n", db)                   # validation searchable, ungapped
            self.assertNotIn("ENSSSAP1", db)
            self.assertIn("ENSP1\thomo_sapiens", (out / "protein_to_species.txt").read_text())

    def test_skip_and_redo_on_panel_change(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            out = self.run_prepare(d)
            stamp = (out / "trees" / "C116_000001.aln").stat().st_mtime_ns
            self.run_prepare(d)                                          # unchanged: skipped
            self.assertEqual((out / "trees" / "C116_000001.aln").stat().st_mtime_ns, stamp)
            (d / "panel.tsv").write_text(PANEL.replace("salmo_salar\texclude", "salmo_salar\tpanel"))
            self.run_prepare(d)
            audit = {r["family"]: r for r in csv.DictReader(_text(out / "family_audit.tsv"), delimiter="\t")}
            self.assertEqual(audit["C116_000002"]["status"], "tiny_tree")    # salmon now kept

    def test_bad_panel_role(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "panel.tsv"
            p.write_text("species\trole\ngallus_gallus\tmaybe\n")
            with self.assertRaises(ValueError):
                cp.read_panel(p)


if __name__ == "__main__":
    unittest.main()


