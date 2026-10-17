"""The annotations downloader (GO and PANTHER), with a fake web: paths, release discovery, reuse, backup."""
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from cichlid_data.downloaders import annotations as an
from helpers import coord_env

LISTING = '<a href="PTHR18.0_human">x</a> <a href="PTHR19.0_human">x</a> <a href="PTHR19.0_zebrafish">x</a>'


class Resp:
    def __init__(self, body, status=200):
        self.body, self.status_code = body.encode(), status
        self.text = body

    def iter_content(self, n):
        yield self.body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class FakeWeb:
    def __init__(self, fail=()):
        self.calls, self.fail = [], set(fail)

    def get(self, url, **kw):
        self.calls.append(url)
        if any(f in url for f in self.fail):
            return Resp("", 404)
        if url.endswith("PANTHER_Sequence_Classification_files/"):
            return Resp(LISTING)
        return Resp(f"content of {url.rsplit('/', 1)[-1]}")


class TestAnnotations(unittest.TestCase):
    def test_download_reuse_and_backup(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "cloud").mkdir()
            with mock.patch.dict(os.environ, coord_env(d)):
                web = FakeWeb(fail=("purl.obolibrary.org",))            # first GO mirror down: falls back
                an.main([], session=web)
                out = d / "local" / "Inputs" / "Annotations"
                self.assertEqual((out / "GO" / "go-basic.obo").read_text(), "content of go-basic.obo")
                self.assertEqual((out / "PANTHER" / "human.tsv").read_text(), "content of PTHR19.0_human")  # newest release
                rec = json.loads((out / "download.json").read_text())
                self.assertEqual(rec["files"]["GO/go-basic.obo"], "https://current.geneontology.org/ontology/go-basic.obo")
                self.assertTrue((d / "cloud" / "Inputs/Annotations/PANTHER/zebrafish.tsv").is_file())   # backed up
                again = FakeWeb()
                an.main(["--no-backup"], session=again)
                self.assertEqual(again.calls, [])                       # everything already here


class TestPassthroughCommands(unittest.TestCase):
    def test_routing_and_compara_defaults(self):
        from cichlid_data import Settings
        from cichlid_data.cli import main as cli
        from cichlid_data.downloaders import compara, expression
        with tempfile.TemporaryDirectory() as d:
            s = Settings(local_root=Path(d), cloud_remote="fake:")
            got = {}
            with mock.patch.object(expression, "main", side_effect=lambda argv: got.setdefault("expr", argv)):
                cli(["download", "expression", "--source", "bgee"], settings=s)
            self.assertEqual(got["expr"], ["--source", "bgee"])                  # options passed through untouched
            with mock.patch.object(compara, "main", side_effect=lambda argv: got.setdefault("cmp", argv)):
                cli(["prepare", "compara", "--outdir", "/elsewhere"], settings=s)
            argv = got["cmp"]
            opt = dict(zip(argv[::2], argv[1::2]))
            self.assertEqual(opt["--outdir"], "/elsewhere")                      # explicit option wins
            self.assertEqual(opt["--nhx"], str(Path(d) / "Inputs/Families/Compara116/dumps/Compara.116.protein_default.nhx.emf.gz"))
            self.assertTrue(opt["--panel"].endswith("config/compara_panel.tsv"))


