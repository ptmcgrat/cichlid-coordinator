"""Offline tests for blumer_samples.py (fake VCF, fake sheet, fake rclone)."""
import csv
import gzip
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cichlid_data.downloaders import blumer_samples as bs
from helpers import coord_env
FAKE_RCLONE = ""        # unused: backups go through the coordinator's stand-in rclone (coord_env)

SAMPLES = [{"SampleID": "cichlid001", "Species": "Maylandia zebra", "Ecogroup": "Mbuna", "Subgroup": "",
            "Sex": "m", "DoB": ""},
           {"SampleID": "cichlid002", "Species": "Diplotaxodon limnothrissa", "Ecogroup": "Deep pelagic",
            "Subgroup": "", "Sex": "f", "DoB": ""},
           {"SampleID": "lab999", "Species": "Maylandia zebra", "Ecogroup": "Mbuna", "Subgroup": "", "Sex": "f"}]
READS = [{"SampleID": "cichlid001", "ProjectID": "Blumer2025"}, {"SampleID": "cichlid001", "ProjectID": "Malinsky2018"},
         {"SampleID": "cichlid002", "ProjectID": "Blumer2025"}]


class TestSamples(unittest.TestCase):
    def test_snapshot(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            data = d / "local" / "Inputs" / "Populations" / "Blumer2025"
            data.mkdir(parents=True)
            with gzip.open(data / "chr1.vcf.gz", "wt") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                         "cichlid001\tcichlid002\tcichlid003\nchr1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n")
            (d / "cloud").mkdir()
            (d / "rclone").write_text(FAKE_RCLONE)
            (d / "rclone").chmod(0o755)
            creds = d / "creds.json"
            creds.write_text("{}")
            env = coord_env(d)
            with mock.patch.dict(os.environ, env):
                bs.main(["--creds", str(creds)], sheet_reader=lambda c: (SAMPLES, READS))   # no Table S1: sheet
            rows = {r["SampleID"]: r for r in csv.DictReader(io.StringIO((data / "samples.tsv").read_text()),
                                                              delimiter="\t")}
            self.assertEqual(list(rows), ["cichlid001", "cichlid002", "cichlid003"])        # VCF order
            self.assertEqual((rows["cichlid001"]["Ecogroup"], rows["cichlid001"]["ProjectIDs"]),
                             ("Mbuna", "Blumer2025,Malinsky2018"))
            self.assertEqual((rows["cichlid003"]["in_database"], rows["cichlid003"]["Species"]), ("False", ""))
            self.assertNotIn("lab999", rows)                                                # not in the VCF
            report = (data / "samples_report.txt").read_text()
            self.assertIn("missing: 1", report)
            self.assertIn("cichlid003", report)
            cloud = d / "cloud" / "Inputs/Populations/Blumer2025"
            self.assertTrue((cloud / "samples.tsv").is_file())
            self.assertFalse(any("creds" in p.name for p in (d / "cloud").rglob("*")))       # key never uploaded

    def test_table_s1_first_then_sheet(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            data = d / "Blumer2025"
            data.mkdir()
            with gzip.open(data / "chr1.vcf.gz", "wt") as fh:
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcichlid7020305\tILBCDS5438962\tcichlid001\tlost1\n")
            (data / "Table_S1.tsv").write_text(
                "sequence_id\tspecimen_id\tseq_depth\tgenus\tspecies\tsex\tclade\tlocation\tsublocation\tBiosample\tBioproject\tInversion10\r\n"
                "cichlid7020305\tD19-E05\t16.35\tAstatotilapia\tcalliptera\tM\tAstCal\tLake_Malombe\tChimwala\tSAMEA117806043\tPRJEB21343\t\r\n"
                "ILBCDS5438962\tX1\t10\tAlticorpus\tgeoffreyi\tU\tDeep\tNkhata_Bay\tNA\tSAMEA1\tPRJEB1254\tHomozygous\r\n")
            calls = []
            def sheet(c):
                calls.append(c)
                return SAMPLES, READS
            bs.main(["--dir", str(data), "--no-backup", "--creds", "x"], sheet_reader=sheet)
            self.assertEqual(calls, [])                                    # Table S1 alone: sheet not read
            rows = {r["SampleID"]: r for r in csv.DictReader(io.StringIO((data / "samples.tsv").read_text()), delimiter="\t")}
            a = rows["cichlid7020305"]
            self.assertEqual((a["Species"], a["Ecogroup"], a["Location"], a["source"]),
                             ("Astatotilapia calliptera", "AstCal", "Lake_Malombe", "Table S1"))
            self.assertEqual((rows["ILBCDS5438962"]["Inversion10"], rows["ILBCDS5438962"]["Sublocation"]), ("Homozygous", ""))
            self.assertEqual(rows["cichlid001"]["in_database"], "False")
            bs.main(["--dir", str(data), "--no-backup", "--creds", "x", "--sheet"], sheet_reader=sheet)
            rows = {r["SampleID"]: r for r in csv.DictReader(io.StringIO((data / "samples.tsv").read_text()), delimiter="\t")}
            self.assertEqual((rows["cichlid001"]["Ecogroup"], rows["cichlid001"]["source"]), ("Mbuna", "lab sheet"))
            self.assertIn("missing: 1", (data / "samples_report.txt").read_text())       # lost1 only

    def test_vcf_without_header(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.vcf.gz"
            with gzip.open(p, "wt") as fh:
                fh.write("chr1\t5\t.\tA\tG\n")
            with self.assertRaises(ValueError):
                bs.vcf_samples(p)


if __name__ == "__main__":
    unittest.main()


