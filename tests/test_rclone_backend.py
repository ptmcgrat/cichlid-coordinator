"""The store tests again, through RcloneBackend and a stand-in rclone."""
import os
import sys
import unittest
from pathlib import Path

from cichlid_data import RcloneBackend, Store
import test_store

FAKE = str(Path(__file__).resolve().parent / "fake_rclone.py")


class RcloneEnv(test_store.Env):
    def store(self, layout, **kw):
        os.environ["FAKE_CLOUD"] = str(self.cloud)
        return Store(layout, RcloneBackend("fake:", rclone=FAKE, retries=1), **kw)


class TestStoreThroughRclone(test_store.TestStore):
    def setUp(self):
        self.e = RcloneEnv()

    def test_backend_basics(self):
        os.environ["FAKE_CLOUD"] = str(self.e.cloud)
        B = RcloneBackend("fake:", rclone=FAKE, retries=1)
        self.assertEqual(B.list("nowhere"), [])
        B.write_text("a/b.txt", "hello")
        self.assertEqual((B.read_text("a/b.txt"), B.list("a"), B.exists("a/b.txt")), ("hello", ["b.txt"], True))
        B.move("a/b.txt", "c/d.txt")
        self.assertEqual((B.exists("a/b.txt"), B.read_text("c/d.txt")), (False, "hello"))
        B.delete("c/d.txt")
        self.assertIsNone(B.read_text("c/d.txt"))

    def test_move_that_errs_but_completes(self):
        os.environ["FAKE_CLOUD"] = str(self.e.cloud)
        B = RcloneBackend("fake:", rclone=FAKE, retries=1)
        B.write_text("Data/Big/a.txt", "a")
        os.environ["FAKE_MOVE_ERRS_BUT_COMPLETES"] = "1"
        try:
            B.move("Data/Big", "Archive/Big")                                  # no exception
        finally:
            os.environ.pop("FAKE_MOVE_ERRS_BUT_COMPLETES")
        self.assertEqual(B.read_text("Archive/Big/a.txt"), "a")

    def test_large_folder_move_falls_back_to_files(self):
        os.environ["FAKE_CLOUD"] = str(self.e.cloud)
        B = RcloneBackend("fake:", rclone=FAKE, retries=1)
        for i in range(5):
            B.write_text(f"Data/NewTreeFiles/sub/TF{i}.nex", f"tree {i}")
        os.environ["FAKE_FAIL_MOVEDIR"] = "1"
        try:
            B.move("Data/NewTreeFiles", "Archive/Project2018/NewTreeFiles")
        finally:
            os.environ.pop("FAKE_FAIL_MOVEDIR")
        self.assertEqual(B.read_text("Archive/Project2018/NewTreeFiles/sub/TF3.nex"), "tree 3")
        self.assertEqual(len(B.list("Archive/Project2018/NewTreeFiles/sub")), 5)
        self.assertNotIn("NewTreeFiles/", B.list("Data"))                       # source folder removed


if __name__ == "__main__":
    unittest.main()
