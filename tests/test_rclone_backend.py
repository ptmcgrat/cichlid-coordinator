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


if __name__ == "__main__":
    unittest.main()
