"""Offline tests for download_blumer.py with a fake DSpace repository and fake rclone."""
import hashlib
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cichlid_data.downloaders import blumer as db
from helpers import coord_env
FAKE_RCLONE = ""        # unused: backups go through the coordinator's stand-in rclone (coord_env)

DATA = {"chr1.vcf.gz": b"A" * 5000, "chr2.vcf.gz": b"B" * 3000,
        "fAstCal1.2.ancestral.fa.gz": b"C" * 2000, "fDipLim2.1.fa.gz": b"D" * 1000}


class Resp:
    def __init__(self, status=200, js=None, body=b"", fail_after=None):
        self.status_code, self._js, self.body, self.fail_after = status, js, body, fail_after
    def json(self):
        return self._js
    def raise_for_status(self):
        raise db.requests.HTTPError(str(self.status_code))
    def iter_content(self, n):
        sent = 0
        for i in range(0, len(self.body), 1000):
            if self.fail_after is not None and sent >= self.fail_after:
                raise db.requests.ConnectionError("dropped")
            yield self.body[i:i + 1000]
            sent += 1000
    def __enter__(self):
        return self
    def __exit__(self, *a):
        pass


class FakeRepo:
    """Two pages of bitstreams; chr1 drops once mid-download; md5 of chr2 can be wrong."""
    def __init__(self, bad_md5=False):
        self.dropped, self.bad_md5, self.ranges = False, bad_md5, []

    def get(self, url, params=None, headers=None, stream=False, timeout=None):
        if url.endswith("/bundles") and url.startswith("https://www."):   # only the API host answers
            return Resp(status=404)
        if url.endswith("/bundles"):
            return Resp(js={"_embedded": {"bundles": [
                {"name": "LICENSE", "_links": {"bitstreams": {"href": "X/lic"}}},
                {"name": "ORIGINAL", "_links": {"bitstreams": {"href": "X/orig"}}}]}})
        if url == "X/orig":
            names = sorted(DATA)
            page = (params or {}).get("page", 0)
            chunk = names[page * 2:(page + 1) * 2]
            bs = [{"name": n, "sizeBytes": len(DATA[n]),
                   "checkSum": {"checkSumAlgorithm": "MD5",
                                "value": "0" * 32 if (self.bad_md5 and n == "chr2.vcf.gz") else hashlib.md5(DATA[n]).hexdigest()},
                   "_links": {"content": {"href": "F/" + n}}} for n in chunk]
            return Resp(js={"_embedded": {"bitstreams": bs}, "page": {"totalPages": 2, "number": page}})
        if url.startswith("F/"):
            name = url[2:]
            rng = (headers or {}).get("Range")
            self.ranges.append((name, rng))
            start = int(rng.split("=")[1].rstrip("-")) if rng else 0
            if name == "chr1.vcf.gz" and not self.dropped:
                self.dropped = True
                return Resp(status=200, body=DATA[name], fail_after=2000)
            return Resp(status=206 if rng else 200, body=DATA[name][start:])
        return Resp(status=404)


class TestBlumer(unittest.TestCase):
    def run_main(self, d, repo, *argv):
        env = coord_env(Path(d))
        with mock.patch.dict(os.environ, env), mock.patch.object(db.time, "sleep"):
            db.main(list(argv), session=repo)
        return Path(d) / "local" / "Inputs" / "Populations" / "Blumer2025"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        (d / "cloud").mkdir()
        (d / "rclone").write_text(FAKE_RCLONE)
        (d / "rclone").chmod(0o755)

    def tearDown(self):
        self.tmp.cleanup()

    def test_listing_selection_resume_and_backup(self):
        repo = FakeRepo()
        out = self.run_main(self.tmp.name, repo)
        names = sorted(p.name for p in out.iterdir())
        self.assertIn("DATA_NOTICE.txt", names)
        self.assertIn("files.tsv", names)
        self.assertNotIn("fDipLim2.1.fa.gz", names)                  # assemblies not needed by default
        for n in ("chr1.vcf.gz", "chr2.vcf.gz", "fAstCal1.2.ancestral.fa.gz"):
            self.assertEqual((out / n).read_bytes(), DATA[n])
        self.assertIn(("chr1.vcf.gz", "bytes=2000-"), repo.ranges)   # resumed, not restarted
        self.assertFalse(list(out.glob("*.part")))
        cloud = Path(self.tmp.name) / "cloud" / "Inputs/Populations/Blumer2025"
        self.assertEqual((cloud / "chr1.vcf.gz").read_bytes(), DATA["chr1.vcf.gz"])
        self.assertIn("Benefit Sharing agreement", (cloud / "DATA_NOTICE.txt").read_text())
        # rerun: everything verified, nothing downloaded again
        repo2 = FakeRepo()
        self.run_main(self.tmp.name, repo2)
        self.assertEqual(repo2.ranges, [])

    def test_checksum_failure_stops_before_backup(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(self.tmp.name, FakeRepo(bad_md5=True))
        self.assertIn("chr2.vcf.gz", str(cm.exception))
        self.assertFalse((Path(self.tmp.name) / "cloud" / "COS").exists())

    def test_falls_back_and_reports(self):
        files = db.list_files(FakeRepo(), apis=["https://www.example/server/api", "https://api.example/server/api"])
        self.assertEqual(len(files), 4)
        with self.assertRaises(IOError) as cm:
            db.list_files(FakeRepo(), apis=["https://www.example/server/api"])
        self.assertIn("could not reach the repository API", str(cm.exception))

    def test_list_and_only(self):
        repo = FakeRepo()
        self.run_main(self.tmp.name, repo, "--list")
        self.assertEqual(repo.ranges, [])
        out = self.run_main(self.tmp.name, FakeRepo(), "--only", "chr2", "--no-backup")
        self.assertEqual(sorted(p.name for p in out.glob("*.gz")), ["chr2.vcf.gz"])


if __name__ == "__main__":
    unittest.main()
