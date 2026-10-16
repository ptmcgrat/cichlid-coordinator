"""Proteomes: Stage 2's tests (same fixtures, same expectations) through the layout, plus the
coordinator's additions."""
import csv
import io
import json
import os
import tempfile
import time
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Layout, Settings, schemas
from cichlid_data import proteomes as pr
from cichlid_data.cli import main as cli
from helpers import Env, write

# RefSeq-shaped GFF3 covering the cases that bit the 2018 code
GFF = """##gff-version 3
NC_1.1\tRefSeq\tregion\t1\t100000\t.\t+\t.\tID=NC_1.1:1..100000;Dbxref=taxon:8128
NC_1.1\tGnomon\tgene\t100\t5000\t.\t+\t.\tID=gene-bmp2;Dbxref=GeneID:111;Name=bmp2;gene=bmp2;gene_biotype=protein_coding;description=bone morphogenetic protein 2
NC_1.1\tGnomon\tmRNA\t100\t5000\t.\t+\t.\tID=rna-XM_1.1;Parent=gene-bmp2;Dbxref=GeneID:111
NC_1.1\tGnomon\texon\t100\t400\t.\t+\t.\tID=exon-XM_1.1-1;Parent=rna-XM_1.1
NC_1.1\tGnomon\tCDS\t150\t400\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1;protein_id=XP_1.1;product=bmp2 isoform X1
NC_1.1\tGnomon\tCDS\t600\t900\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1;protein_id=XP_1.1;product=bmp2 isoform X1
NC_1.1\tGnomon\tmRNA\t100\t5000\t.\t+\t.\tID=rna-XM_2.1;Parent=gene-bmp2;partial=true;start_range=.,100
NC_1.1\tGnomon\tCDS\t100\t900\t.\t+\t0\tID=cds-XP_2.1;Parent=rna-XM_2.1;protein_id=XP_2.1;partial=true
NC_1.1\tGnomon\tmRNA\t100\t5000\t.\t+\t.\tID=rna-XM_3.1;Parent=gene-bmp2
NC_1.1\tGnomon\tCDS\t150\t900\t.\t+\t0\tID=cds-XP_3.1;Parent=rna-XM_3.1;protein_id=XP_3.1
NC_1.1\tGnomon\tgene\t6000\t7000\t.\t-\t.\tID=gene-fraglike;Dbxref=GeneID:222;gene=fraglike;gene_biotype=protein_coding;description=only partial%2C sadly
NC_1.1\tGnomon\tmRNA\t6000\t7000\t.\t-\t.\tID=rna-XM_4.1;Parent=gene-fraglike;partial=true
NC_1.1\tGnomon\tCDS\t6000\t7000\t.\t-\t0\tID=cds-XP_4.1;Parent=rna-XM_4.1;protein_id=XP_4.1;partial=true
NC_1.1\tGnomon\tmRNA\t6000\t7000\t.\t-\t.\tID=rna-XM_5.1;Parent=gene-fraglike;partial=true
NC_1.1\tGnomon\tCDS\t6000\t6500\t.\t-\t0\tID=cds-XP_5.1;Parent=rna-XM_5.1;protein_id=XP_5.1;partial=true
NC_1.1\tGnomon\tgene\t8000\t9000\t.\t+\t.\tID=gene-sel;Dbxref=GeneID:333;gene=selenop;gene_biotype=protein_coding
NC_1.1\tGnomon\tmRNA\t8000\t9000\t.\t+\t.\tID=rna-XM_6.1;Parent=gene-sel
NC_1.1\tGnomon\tCDS\t8000\t9000\t.\t+\t0\tID=cds-XP_6.1;Parent=rna-XM_6.1;protein_id=XP_6.1
NC_1.1\tGnomon\tgene\t9500\t9900\t.\t+\t.\tID=gene-lnc;Dbxref=GeneID:444;gene_biotype=lncRNA
NC_1.1\tGnomon\tlnc_RNA\t9500\t9900\t.\t+\t.\tID=rna-XR_1.1;Parent=gene-lnc
NC_1.1\tGnomon\tpseudogene\t10000\t11000\t.\t+\t.\tID=gene-ps;Dbxref=GeneID:555;gene_biotype=pseudogene
NC_1.1\tGnomon\tgene\t12000\t13000\t.\t+\t.\tID=gene-ighv;Dbxref=GeneID:666;gene_biotype=V_segment
NC_1.1\tGnomon\tV_gene_segment\t12000\t13000\t.\t+\t.\tID=id-ighv;Parent=gene-ighv
NC_1.1\tGnomon\tCDS\t12000\t13000\t.\t+\t0\tID=cds-ighv;Parent=id-ighv;protein_id=XP_9.1
NC_1.1\tGnomon\tgene\t14000\t15000\t.\t+\t.\tID=gene-lost;Dbxref=GeneID:777;gene=lost;gene_biotype=protein_coding
NC_1.1\tGnomon\tmRNA\t14000\t15000\t.\t+\t.\tID=rna-XM_7.1;Parent=gene-lost
NC_1.1\tGnomon\tCDS\t14000\t15000\t.\t+\t0\tID=cds-XP_7.1;Parent=rna-XM_7.1;protein_id=XP_7.1
NC_MT.1\tRefSeq\tgene\t1\t1500\t.\t+\t.\tID=gene-COX1;Dbxref=GeneID:888;gene=COX1;gene_biotype=protein_coding
NC_MT.1\tRefSeq\tCDS\t1\t1500\t.\t+\t0\tID=cds-YP_1.1;Parent=gene-COX1;protein_id=YP_1.1;product=cytochrome c oxidase subunit I
"""

FAA = """>XP_1.1 bone morphogenetic protein 2 isoform X1 [O. niloticus]
MAAAAAAAAA
>XP_2.1 bone morphogenetic protein 2 isoform X2, partial [O. niloticus]
MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>XP_3.1 LOW QUALITY PROTEIN: bone morphogenetic protein 2 isoform X3 [O. niloticus]
MAAAAAAAAA
>XP_4.1 fraglike, partial [O. niloticus]
MKKKKKKKKK
>XP_5.1 fraglike, partial [O. niloticus]
MKK
>XP_6.1 selenoprotein P [O. niloticus]
MUUAC*
>XP_9.1 ig heavy chain V [O. niloticus]
MQQQ
>YP_1.1 cytochrome c oxidase subunit I (mitochondrion) [O. niloticus]
MFITRWLF*
>XP_99.1 orphan protein not in the GFF
MSSS
"""

SEQREPORT = '{"refseqAccession":"NC_MT.1","assignedMoleculeLocationType":"Mitochondrion"}\n' \
            '{"refseqAccession":"NC_1.1","assignedMoleculeLocationType":"Chromosome"}\n'

# Braker-style GTF: no gene lines, no protein_id, FASTA keyed by transcript ID
GTF = """chr1\tAUGUSTUS\tCDS\t100\t400\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
chr1\tAUGUSTUS\tCDS\t600\t900\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
chr1\tAUGUSTUS\tCDS\t100\t950\t.\t+\t0\ttranscript_id "g1.t2"; gene_id "g1";
chr2\tAUGUSTUS\tCDS\t50\t80\t.\t-\t0\ttranscript_id "g2.t1"; gene_id "g2";
"""
GTF_FAA = ">g1.t1\nMAAA\n>g1.t2\nMAAAAAA\n>g2.t1\nMK\n"


MANIFEST = ("prefix\tspecies\taccession\tsource\tannotation_format\n"
            "Onil\tOreochromis niloticus\tGCF_001858045.2\tncbi\tgff3\n"
            "Brk\tEtroplus suratensis\t\tlocal\tgtf\n")


def setup(root, species=("Onil",)):
    L = Layout(Path(root) / "data", "fake:")
    if "Onil" in species:
        dl = L.download("Onil").local
        write(dl / "genomic.gff", GFF)
        write(dl / "protein.faa", FAA)
        write(dl / "sequence_report.jsonl", SEQREPORT)
        write(dl / "download.json", json.dumps({"ncbi_metadata": {"protein_coding_genes": 5}}))
    if "Brk" in species:
        dl = L.download("Brk").local
        write(dl / "genomic.gtf", GTF)
        write(dl / "protein.faa", GTF_FAA)
    man = write(Path(root) / "species.tsv", MANIFEST)
    return L, man


def read_genes(path):
    with open(path) as fh:
        return {r["symbol"] or r["annotation_gene_key"]: r for r in csv.DictReader(fh, delimiter="\t")}


class TestStage2Rules(unittest.TestCase):
    """Stage 2's own tests, unchanged in substance."""
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def run_species(self, mode="complete-longest"):
        L, man = setup(self.root)
        pr.build(L, man, species=["Onil"], mode=mode, say=lambda *_: None)
        return L.proteomes().local

    def test_selection_rules(self):
        out = self.run_species()
        g = read_genes(out / "Onil.genes.tsv")
        self.assertEqual(g["bmp2"]["chosen_protein"], "XP_1.1")
        self.assertEqual(g["bmp2"]["n_proteins"], "3")
        self.assertEqual(g["fraglike"]["chosen_protein"], "XP_4.1")
        self.assertEqual(g["fraglike"]["partial"], "True")
        self.assertEqual(g["fraglike"]["description"], "only partial, sadly")

    def test_longest_mode_ignores_completeness(self):
        g = read_genes(self.run_species("longest") / "Onil.genes.tsv")
        self.assertEqual(g["bmp2"]["chosen_protein"], "XP_2.1")

    def test_mito_selenocysteine_and_biotypes(self):
        out = self.run_species()
        g = read_genes(out / "Onil.genes.tsv")
        self.assertEqual(g["COX1"]["chosen_protein"], "YP_1.1")
        self.assertEqual(g["COX1"]["mitochondrial"], "True")
        self.assertEqual(g["selenop"]["residues_replaced"], "2")
        for skipped in ("gene-lnc", "gene-ps", "gene-ighv"):
            self.assertNotIn(skipped, g)
        self.assertEqual(g["lost"]["chosen_protein"], "")
        seqs = pr.read_fasta(out / "Onil.faa")
        self.assertEqual(seqs["Onil_333"][1], "MXXAC")
        self.assertEqual(seqs["Onil_888"][1], "MFITRWLF")
        self.assertEqual(len(seqs), 4)
        self.assertTrue(all(k.startswith("Onil_") and "_" not in k[5:] for k in seqs))
        with open(out / "stage2_summary.tsv") as fh:
            s = next(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(s["protein_coding_genes"], "5")
        self.assertEqual(s["coding_genes_without_protein"], "1")
        self.assertEqual(s["fasta_proteins_not_in_annotation"], "1")
        self.assertEqual(s["cds_protein_not_in_fasta"], "1")
        self.assertEqual(s["mitochondrial_proteins"], "1")
        self.assertIn("lncRNA=1", s["skipped_biotypes"])
        self.assertIn("V_segment=1", s["skipped_biotypes"])

    def test_braker_style_gtf(self):
        L, man = setup(self.root, species=("Brk",))
        summaries, failed = pr.build(L, man, species=["Brk"], say=lambda *_: None)
        self.assertEqual((failed, summaries[0]["proteins_written"]), ([], 2))
        g = read_genes(L.genes_table("Brk").local)
        self.assertEqual(g["g1"]["chosen_protein"], "g1.t2")
        self.assertEqual((g["g1"]["start"], g["g1"]["end"]), ("100", "950"))
        self.assertIn("Brk_g1", pr.read_fasta(L.proteome("Brk").local))


class TestPrimaryAssemblyOnly(unittest.TestCase):
    def test_alternate_loci_and_patches_are_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            L, man = setup(d)
            dl = L.download("Onil").local
            # a second copy of bmp2 (same GeneID) on an alternate locus, and a gene on a patch
            alt = ("NT_ALT.1\tGnomon\tgene\t100\t5000\t.\t+\t.\tID=gene-bmp2-alt;Dbxref=GeneID:111;gene=bmp2;gene_biotype=protein_coding\n"
                   "NT_ALT.1\tGnomon\tmRNA\t100\t5000\t.\t+\t.\tID=rna-XM_alt.1;Parent=gene-bmp2-alt\n"
                   "NT_ALT.1\tGnomon\tCDS\t150\t900\t.\t+\t0\tID=cds-XP_1.1b;Parent=rna-XM_alt.1;protein_id=XP_1.1\n"
                   "NW_PATCH.1\tGnomon\tgene\t1\t500\t.\t+\t.\tID=gene-p;Dbxref=GeneID:999;gene=patchy;gene_biotype=protein_coding\n")
            (dl / "genomic.gff").write_text(GFF + alt)
            # the primary unit is named after the strain, as in mouse GRCm39: roles decide, not names
            (dl / "sequence_report.jsonl").write_text(
                SEQREPORT.replace('"Mitochondrion"}', '"Mitochondrion","assemblyUnit":"non-nuclear","role":"assembled-molecule"}')
                         .replace('"Chromosome"}', '"Chromosome","assemblyUnit":"C57BL/6J","role":"assembled-molecule"}')
                + '{"refseqAccession":"NT_ALT.1","assemblyUnit":"ALT_REF_LOCI_1","role":"alt-scaffold"}\n'
                + '{"refseqAccession":"NW_PATCH.1","assemblyUnit":"PATCHES","role":"fix-patch"}\n')
            pr.build(L, man, species=["Onil"], say=lambda *_: None)
            ids = [r["gene_id"] for r in schemas.GENES.read(L.genes_table("Onil").local)]
            self.assertNotIn("Onil_111.2", ids)                                 # no false duplicate
            self.assertNotIn("Onil_999", ids)
            self.assertIn("Onil_888", ids)                                      # mitochondrion kept
            with open(L.proteomes().local / "stage2_summary.tsv") as fh:
                s = next(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual((s["non_primary_genes_skipped"], s["protein_coding_genes"],
                              s["coding_genes_without_protein"], s["duplicate_gene_ids"]), ("2", "5", "1", "0"))


    def test_refuses_when_most_genes_would_be_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            L, man = setup(d)
            dl = L.download("Onil").local
            (dl / "sequence_report.jsonl").write_text(
                '{"refseqAccession":"NC_1.1","role":"alt-scaffold"}\n{"refseqAccession":"NC_MT.1","role":"assembled-molecule",'
                '"assignedMoleculeLocationType":"Mitochondrion"}\n')
            out = []
            summaries, failed = pr.build(L, man, species=["Onil"], say=out.append)
            self.assertEqual(failed, ["Onil"])
            self.assertIn("would be skipped as alternate loci or patches; refusing", out[-1])
            self.assertFalse(L.proteome("Onil").local.exists())


class TestCoordinatorAdditions(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    def test_outputs_validate(self):
        L, man = setup(self.root, species=("Onil", "Brk"))
        pr.build(L, man, say=lambda *_: None)
        for p in ("Onil", "Brk"):
            self.assertEqual(schemas.validate(L.proteome(p).local), [])
            self.assertEqual(schemas.validate(L.genes_table(p).local), [])

    def test_skip_force_and_selection(self):
        L, man = setup(self.root, species=("Onil", "Brk"))
        out = []
        pr.build(L, man, species=["onil", "Nope"], say=out.append)        # prefixes are case-insensitive
        self.assertTrue(L.proteome("Onil").local.is_file())
        self.assertFalse(L.proteome("Brk").local.exists())
        self.assertIn("[nope] not in", out[0])
        out = []
        pr.build(L, man, species=["Onil"], say=out.append)
        self.assertEqual(out, ["[Onil] up to date, skipping"])
        out = []
        pr.build(L, man, species=["Onil"], force=True, say=out.append)
        self.assertIn("4 proteins from 5 coding genes", out[0])
        # a newer download makes it out of date again
        t = time.time() + 10
        os.utime(L.download("Onil").local / "protein.faa", (t, t))
        self.assertFalse(pr.up_to_date(L, "Onil"))

    def test_invalid_output_is_not_written(self):
        L, man = setup(self.root)
        orig = pr.output_id
        pr.output_id = lambda prefix, g: f"{prefix}_bad_id"                # underscores: not a valid gene ID
        try:
            summaries, failed = pr.build(L, man, species=["Onil"], say=lambda *_: None)
        finally:
            pr.output_id = orig
        self.assertEqual((summaries, failed), ([], ["Onil"]))
        self.assertFalse(L.genes_table("Onil").local.exists())
        self.assertFalse(L.proteome("Onil").local.exists())

    def test_summary_keeps_other_species(self):
        L, man = setup(self.root, species=("Onil", "Brk"))
        pr.build(L, man, say=lambda *_: None)
        pr.build(L, man, species=["Onil"], force=True, say=lambda *_: None)
        with open(L.proteomes().local / "stage2_summary.tsv") as fh:
            self.assertEqual([r["prefix"] for r in csv.DictReader(fh, delimiter="\t")], ["Brk", "Onil"])

    def test_fetches_missing_download_and_updates_catalog(self):
        e = Env()
        try:
            L, man = setup(e.a.local_root.parent / "seed", species=("Onil",))
            # machine a has the download only in the cloud
            import shutil
            dst = e.a.download("Onil").local
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(L.download("Onil").local, dst)
            e.store(e.a).push("proteins/downloads")
            shutil.rmtree(e.a.download("Onil").local)
            cfg = e.a.local_root.parent / "cfg"
            write(cfg / "species.tsv", MANIFEST)
            write(cfg / "species_sets.toml", '[methods.OrthoFinder2026]\nspecies = ["Onil"]\n[[organism]]\nprefix = "Onil"\n'
                                             'species = "Oreochromis niloticus"\nrole = "project"\n[[organism]]\n'
                                             'prefix = "Brk"\nspecies = "Etroplus suratensis"\nrole = "project"\n')
            os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
            s = Settings(local_root=e.a.local_root, cloud_remote="fake:")
            f = io.StringIO()
            with redirect_stdout(f):
                try:
                    cli(["proteomes", "--species", "Onil", "--push"], backend=e.store(e.a).backend, settings=s)
                except SystemExit:
                    pass
            out = f.getvalue()
            self.assertTrue(e.a.proteome("Onil").local.is_file())             # download fetched, proteome built
            self.assertIn("Catalog updated: 1 proteomes ready", out)
            self.assertIn("pushed proteins/proteomes: uploaded", out)
            cat = json.loads(e.a.catalog().local.read_text())
            self.assertEqual(cat["organisms"]["Onil"]["status"], "proteome ready")
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)
            e.close()


if __name__ == "__main__":
    unittest.main()
