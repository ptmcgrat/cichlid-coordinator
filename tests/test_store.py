import json
import unittest
from pathlib import Path

from cichlid_data import TransferError
from helpers import Env, touch_later, write


class TestStore(unittest.TestCase):
    def setUp(self):
        self.e = Env()

    def tearDown(self):
        self.e.close()

    def test_tar_round_trip_and_conflicts(self):
        a, b = self.e.a, self.e.b
        write(a.proteome("Mzeb").local, ">Mzeb_1\nM\n")
        Sa = self.e.store(a)
        self.assertEqual(Sa.push("proteins/proteomes"), {"proteins/proteomes": "uploaded"})
        self.assertTrue((self.e.cloud / "Inputs/Proteins/proteomes.tar").is_file())
        Sb = self.e.store(b)
        self.assertEqual(Sb.pull("proteins/proteomes"), {"proteins/proteomes": "downloaded"})
        self.assertEqual(b.proteome("Mzeb").local.read_text(), ">Mzeb_1\nM\n")
        self.assertEqual(self.e.store(b).pull("proteins/proteomes"), {"proteins/proteomes": "up to date"})
        # both machines change it: b pushes first, then a's push is refused, and a's pull too
        touch_later(b.proteome("Mzeb").local)
        self.e.store(b).push("proteins/proteomes")
        touch_later(a.proteome("Mzeb").local)
        touch_later(a.proteome("Mzeb").local)                                   # a different edit from b's
        with self.assertRaises(TransferError) as err:
            self.e.store(a).push("proteins/proteomes")
        self.assertIn("CONFLICT", str(err.exception))
        self.assertIn("CONFLICT", self.e.store(a).pull("proteins/proteomes")["proteins/proteomes"])
        self.assertEqual(self.e.store(a, force=True).pull("proteins/proteomes"), {"proteins/proteomes": "downloaded"})

    def test_shards_only_upload_what_changed(self):
        a = self.e.a
        R = a.results("Compara2026")
        write(R.family("TF101001").local / "tree.nwk", "(a,b);")
        write(R.family("TF101002").local / "tree.nwk", "(c,d);")
        write(R.family("TF314000").local / "tree.nwk", "(e,f);")
        write(R.family("TF314000").local / "stage5/x.fa", "scratch")                  # excluded
        write(a.local_root / "Results/Compara2026/families/family_summary.tsv", "s")
        res = self.e.store(a).push("results/Compara2026/families")
        self.assertEqual(set(res), {"results/Compara2026/families/TF101", "results/Compara2026/families/TF314",
                                    "results/Compara2026/families/_top"})
        touch_later(R.family("TF314000").local / "tree.nwk")
        res = self.e.store(a).push("results/Compara2026/families")
        self.assertEqual(res["results/Compara2026/families/TF314"], "uploaded")
        self.assertEqual(res["results/Compara2026/families/TF101"], "up to date")
        b = self.e.b
        self.e.store(b).pull("results/Compara2026/families")
        self.assertTrue((b.results("Compara2026").family("TF101002").local / "tree.nwk").is_file())
        self.assertFalse((b.results("Compara2026").family("TF314000").local / "stage5/x.fa").exists())  # excluded

    def test_immutable_files_and_locks(self):
        a, b = self.e.a, self.e.b
        write(a.family_resource("TreeFam9").local / "hmm/TreeFam9", "hmm v1")
        self.e.store(a).push("families/TreeFam9/hmm")
        self.assertEqual((self.e.cloud / "Inputs/Families/TreeFam9/hmm/TreeFam9").read_text(), "hmm v1")
        touch_later(a.family_resource("TreeFam9").local / "hmm/TreeFam9")
        with self.assertRaises(TransferError) as err:
            self.e.store(a).push("families/TreeFam9/hmm")
        self.assertIn("REFUSED (", str(err.exception))
        self.assertIn("immutable", str(err.exception))
        self.assertEqual((self.e.cloud / "Inputs/Families/TreeFam9/hmm/TreeFam9").read_text(), "hmm v1")
        # adding a file to an immutable file-by-file unit is allowed (e.g. a source.json)
        (a.family_resource("TreeFam9").local / "hmm/TreeFam9").write_text("hmm v1")       # undo the edit
        write(a.family_resource("TreeFam9").local / "hmm/README.txt", "notes")
        self.assertEqual(self.e.store(a).push("families/TreeFam9/hmm"), {"families/TreeFam9/hmm": "uploaded"})
        self.assertEqual((self.e.cloud / "Inputs/Families/TreeFam9/hmm/README.txt").read_text(), "notes")
        # tar-packed immutable units stay strict
        write(a.download("Mzeb").local / "protein.faa", "P")
        self.e.store(a).push("proteins/downloads")
        write(a.download("Mzeb").local / "extra.txt", "x")
        with self.assertRaises(TransferError):
            self.e.store(a).push("proteins/downloads")
        # a stale-looking but fresh lock blocks a push
        write(b.proteome("Hbur").local, ">Hbur_1\nM\n")
        lock = self.e.cloud / ".cichlid-data/locks/proteins__proteomes.lock"
        write(lock, json.dumps({"user": "someone", "machine": "m", "time": "now", "epoch": 9e12}))
        with self.assertRaises(TransferError):
            self.e.store(b).push("proteins/proteomes")
        lock.unlink()
        self.e.store(b).push("proteins/proteomes")
        self.assertFalse(lock.exists())                                           # released after the push

    def test_ensure_and_status(self):
        a, b = self.e.a, self.e.b
        write(a.download("Mzeb").local / "protein.faa", "P")
        write(a.download("Acal").local / "protein.faa", "Q")
        write(a.population("Blumer2025").local / "chr10.vcf.gz", "VCF")
        Sa = self.e.store(a)
        Sa.push(["proteins/downloads", "populations"])
        Sb = self.e.store(b)
        # ensure fetches just the owning shard
        p = Sb.ensure(b.download("Mzeb") / "protein.faa")
        self.assertEqual(p.read_text(), "P")
        self.assertFalse(b.download("Acal").local.exists())
        # on-demand units are skipped by bulk pulls but fetched by ensure, one file at a time
        Sb.pull("*")
        self.assertFalse(b.population("Blumer2025").local.exists())
        self.assertEqual(Sb.ensure(b.population("Blumer2025") / "chr10.vcf.gz").read_text(), "VCF")
        st = dict((k, s) for k, s, _ in self.e.store(b).status("proteins/*"))
        self.assertEqual(st["proteins/downloads/Mzeb"], "in sync")
        self.assertEqual(st["proteins/downloads/Acal"], "in sync")              # pulled by the bulk pull
        log = (b.state_dir() / "transfers.log").read_text()
        self.assertIn("proteins/downloads/Mzeb", log)

    def test_download_reports_and_retire(self):
        a = self.e.a
        write(a.download("Mzeb").local / "protein.faa", "P")
        write(a.download("Mzeb").local.parent / "download_summary.tsv", "v1")
        res = self.e.store(a).push(["proteins/downloads", "proteins/download_summary"])
        self.assertEqual(set(res), {"proteins/downloads/Mzeb", "proteins/download_summary"})   # no _top shard
        write(a.download("Mzeb").local.parent / "download_summary.tsv", "v2, rewritten by the next run")
        res = self.e.store(a).push(["proteins/downloads", "proteins/download_summary"])
        self.assertEqual(res["proteins/download_summary"], "uploaded")                         # not immutable
        self.assertEqual((self.e.cloud / "Inputs/Proteins/downloads/download_summary.tsv").read_text(),
                         "v2, rewritten by the next run")
        # an old _top shard, as left by earlier versions, can be retired
        S = self.e.store(a)
        S.index()["proteins/downloads/_top"] = {"fingerprint": "old"}
        write(self.e.cloud / "Inputs/Proteins/downloads/_top.tar", "old tar")
        S._save()
        self.assertEqual(self.e.store(a).retire("proteins/downloads/_top"), "deleted")
        self.assertFalse((self.e.cloud / "Inputs/Proteins/downloads/_top.tar").exists())
        self.assertNotIn("proteins/downloads/_top", self.e.store(a).index())
        self.assertEqual(self.e.store(a).retire("proteins/downloads/_top"), "not in the cloud index")

    def test_dry_run_changes_nothing(self):
        a = self.e.a
        write(a.proteome("Mzeb").local, ">Mzeb_1\nM\n")
        self.assertEqual(self.e.store(a, dry_run=True).push("proteins/proteomes"), {"proteins/proteomes": "would upload"})
        self.assertEqual(list(self.e.cloud.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
