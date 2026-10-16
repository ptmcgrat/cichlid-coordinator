"""Offline tests for summarize_tissue_rnaseq.py."""
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

from cichlid_data.downloaders import rnaseq as st


def layout(root):
    def tpm(tissue, header, rows):
        p = root / tissue / "output" / "star_salmon"
        p.mkdir(parents=True)
        (p / "salmon.merged.gene_tpm.tsv").write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    tpm("heart", ["gene_id", "gene_name", "Mz_heart_Rpool"], [["avpr1aa", "avpr1aa", 12.5], ["LOC999", "LOC999", 0]])
    tpm("retina", ["gene_id", "gene_name", "Mz_retina_R1", "Mz_retina_R2"], [["avpr1aa", "avpr1aa", 2, 4]])
    tpm("all_tissues_pool", ["gene_id", "gene_name", "Mz_all"], [["avpr1aa", "avpr1aa", 1]])
    (root / "kidney" / "output").mkdir(parents=True)                     # still uploading: no table yet


class TestTissueRNAseq(unittest.TestCase):
    def test_combine(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            layout(d / "Brawand_data")
            st.main(["--local", str(d / "Brawand_data"), "--dataset", "Brawand", "--species", "Mzeb",
                     "--outdir", str(d / "out"), "--no-backup"])
            rows = [l.split("\t") for l in gzip.decompress((d / "out" / "tissue_tpm.tsv.gz").read_bytes()).decode().splitlines()[1:]]
            got = {(r[0], r[2]): (r[3], r[4], r[5]) for r in rows}
            self.assertEqual(got[("avpr1aa", "heart")], ("12.5", "1", "12.5"))
            self.assertEqual(got[("avpr1aa", "retina")], ("3", "2", "2,4"))         # mean of two samples
            self.assertNotIn(("avpr1aa", "all_tissues_pool"), got)                  # pools skipped
            info = json.loads((d / "out" / "summary.json").read_text())
            self.assertEqual(info["tissues"], ["heart", "retina"])
            self.assertTrue(any("kidney" in s for s in info["skipped"]))
            self.assertTrue(any("pool" in s for s in info["skipped"]))
            self.assertIn("Mz_retina_R2", (d / "out" / "samples.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
