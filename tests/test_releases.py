"""Releases: freezing one, its archive Releases/<release>/, drift from it, and the `release` command."""
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from cichlid_data import Catalog, LocalBackend, Settings, Store, TransferError, __version__, releases, schemas
from cichlid_data.check import check
from cichlid_data.cli import main as cli
from cichlid_data.releases import ReleaseError, drift, freeze, load_manifest, verify_archive
from helpers import Env, releases_toml, touch_later, write
from test_check import proteome

SETS = '''
[[organism]]
prefix = "Mzeb"
species = "Maylandia zebra"
role = "project"
[[organism]]
prefix = "Acal"
species = "Astatotilapia calliptera"
role = "project"
[[organism]]
prefix = "Drer"
species = "Danio rerio"
role = "reference"
'''
TREE = "(Drer,(Mzeb,Acal));\n"
PARTS = {"TreeFam": 'queries = "project"',
         "OrthoFinder": 'species = ["Mzeb", "Acal", "Drer"]\nspecies_tree = "tree.nwk"',
         "Homology": 'description = "calls"'}


def make_config(root, status="frozen", parts=PARTS, tree=TREE, deliver=("Mzeb",)):
    cfg = Path(root) / "config"
    write(cfg / "species_sets.toml", SETS)
    write(cfg / "releases.toml", releases_toml(parts, status=status, deliver=deliver,
                                               inputs={"treefam": "TreeFam 9", "compara": "Ensembl Compara 116"}))
    write(cfg / "species.tsv", "prefix\tspecies\taccession\tsource\nMzeb\tMaylandia zebra\tGCF_1\tncbi\n")
    write(cfg / "compara_panel.tsv", "species\trole\tgroup\n")
    if tree is not None:
        write(cfg / "tree.nwk", tree)
    return cfg


class Base(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        from cichlid_data import Layout
        self.L = Layout(self.root / "data", "fake:")
        self.cfg = make_config(self.root)
        for p in ("Mzeb", "Acal", "Drer"):
            proteome(self.L, p, seq="MKLV" + "A" * len(p))

    def tearDown(self):
        self._t.cleanup()

    def catalog(self):
        return Catalog.build(self.L, self.cfg)

    def frozen(self):
        freeze(self.L, self.catalog(), "v1", say=lambda m: None)
        return self.L.release("v1").root().local


class TestFreeze(Base):
    def test_the_archive(self):
        said = []
        root = freeze(self.L, self.catalog(), "v1", say=said.append)
        self.assertEqual(root, self.root / "data" / "Releases" / "v1")
        self.assertEqual(sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()),
                         ["catalog.json", "config/tree.nwk",
                          "genes/Acal.genes.tsv", "genes/Drer.genes.tsv", "genes/Mzeb.genes.tsv",
                          "proteomes/Acal.faa", "proteomes/Drer.faa", "proteomes/Mzeb.faa", "release.json"])
        self.assertIn("cichlid-data push releases/v1", said[-1])
        self.assertFalse((root.parent / ".v1.partial").exists())
        m = load_manifest(self.L, "v1")
        self.assertEqual((m["release"], m["status"], m["deliver"], m["coordinator_version"]),
                         ("v1", "frozen", ["Mzeb"], __version__))
        self.assertEqual(m["inputs"], {"treefam": "TreeFam 9", "compara": "Ensembl Compara 116"})
        self.assertEqual(sorted(m["species"]), ["Acal", "Drer", "Mzeb"])
        self.assertEqual(m["parts"]["OrthoFinder"]["species_tree"], "tree.nwk")
        self.assertEqual(m["parts"]["OrthoFinder"]["species"], ["Mzeb", "Acal", "Drer"])
        self.assertEqual(m["species"]["Mzeb"]["proteome_md5"], self.catalog()["Mzeb"].proteome_md5)
        self.assertEqual(m["catalog_md5"], m["files"]["catalog.json"]["md5"])
        self.assertEqual(m["files"]["proteomes/Mzeb.faa"]["size"], (root / "proteomes/Mzeb.faa").stat().st_size)
        self.assertNotIn("release.json", m["files"])
        # the archive is a faithful copy, and carries the catalog and tree as they were
        self.assertEqual((root / "proteomes/Mzeb.faa").read_bytes(), self.L.proteome("Mzeb").local.read_bytes())
        self.assertEqual((root / "config/tree.nwk").read_text(), TREE)
        self.assertIn("releases", json.loads((root / "catalog.json").read_text()))
        self.assertEqual(verify_archive(self.L, "v1"), [])
        self.assertEqual(drift(self.L, "v1", self.cfg), [])
        # the layout API reaches the archive's files
        self.assertEqual(self.L.release("v1").proteome("Drer").local, root / "proteomes/Drer.faa")

    def test_published_may_be_frozen_too(self):
        make_config(self.root, status="published")
        self.assertTrue(self.frozen().is_dir())

    def test_only_the_owner_decision_allows_it(self):
        for status in ("planned", "building", "superseded"):
            make_config(self.root, status=status)
            with self.assertRaises(ReleaseError) as e:
                freeze(self.L, self.catalog(), "v1", say=lambda m: None)
            self.assertIn(f"status {status!r}", str(e.exception) if status != "superseded" else f"status {status!r}")
            self.assertFalse((self.root / "data" / "Releases").exists())

    def test_refuses_what_is_not_valid_or_not_ready(self):
        with self.assertRaises(KeyError):
            freeze(self.L, self.catalog(), "v9", say=lambda m: None)
        make_config(self.root, tree="(Drer,Mzeb,Acal);\n")                       # a polytomy: invalid
        with self.assertRaises(ReleaseError) as e:
            freeze(self.L, self.catalog(), "v1", say=lambda m: None)
        self.assertIn("is not valid", str(e.exception))
        make_config(self.root)
        self.L.proteome("Acal").local.unlink()
        with self.assertRaises(ReleaseError) as e:
            freeze(self.L, self.catalog(), "v1", say=lambda m: None)
        self.assertIn("no proteome yet for Acal", str(e.exception))
        self.assertFalse((self.root / "data" / "Releases").exists())

    def test_a_frozen_release_is_never_frozen_again(self):
        root = self.frozen()
        before = (root / "release.json").read_text()
        with self.assertRaises(ReleaseError) as e:
            freeze(self.L, self.catalog(), "v1", say=lambda m: None)
        self.assertIn("already frozen", str(e.exception))
        self.assertEqual((root / "release.json").read_text(), before)

    def test_the_archive_appears_whole_or_not_at_all(self):
        real = shutil.copy2
        calls = []

        def flaky(src, dst, *a, **kw):
            calls.append(dst)
            if len(calls) == 3:
                raise OSError("disk full")
            return real(src, dst, *a, **kw)

        with mock.patch.object(releases.shutil, "copy2", flaky), self.assertRaises(OSError):
            freeze(self.L, self.catalog(), "v1", say=lambda m: None)
        self.assertEqual(list((self.root / "data" / "Releases").iterdir()), [])        # nothing, not even a partial
        self.assertIsNotNone(freeze(self.L, self.catalog(), "v1", say=lambda m: None))   # and a retry works

    def test_the_catalog_in_the_archive_is_the_one_frozen(self):
        root = self.frozen()
        write(self.L.proteome("Mzeb").local, ">Mzeb_1\nMKLVAAAAAAAAAA\n")               # the annotation improves later
        self.assertEqual((root / "proteomes/Mzeb.faa").read_text(), ">Mzeb_1\nMKLVAAAA\n")
        self.assertNotEqual(self.catalog()["Mzeb"].proteome_md5, load_manifest(self.L, "v1")["species"]["Mzeb"]["proteome_md5"])


class TestVerifyAndDrift(Base):
    def test_verify_archive_finds_every_kind_of_change(self):
        root = self.frozen()
        self.assertEqual(verify_archive(self.L, "v1"), [])
        f = root / "proteomes/Mzeb.faa"
        text = f.read_text()
        f.write_text(text.replace("A", "C", 1))                                          # same size, other content
        self.assertEqual(verify_archive(self.L, "v1", deep=False), [])
        self.assertEqual(verify_archive(self.L, "v1"), ["proteomes/Mzeb.faa: content differs from when it was frozen"])
        f.write_text(text + "extra")
        self.assertIn("size", verify_archive(self.L, "v1", deep=False)[0])
        f.write_text(text)
        (root / "genes/Acal.genes.tsv").unlink()
        self.assertEqual(verify_archive(self.L, "v1"), ["genes/Acal.genes.tsv: missing"])
        shutil.copy2(self.L.genes_table("Acal").local, root / "genes/Acal.genes.tsv")
        write(root / "proteomes/New.faa", ">x\nM\n")
        self.assertEqual(verify_archive(self.L, "v1"), ["proteomes/New.faa: not in the manifest (a frozen release does not gain files)"])

    def test_not_frozen_here(self):
        self.assertIn("not frozen on this machine", verify_archive(self.L, "v1")[0])
        self.assertIsNone(load_manifest(self.L, "v1"))
        self.assertEqual(drift(self.L, "v1"), [])

    def test_drift_reports_newer_proteomes_and_configuration(self):
        self.frozen()
        self.assertEqual(drift(self.L, "v1", self.cfg), [])
        write(self.L.proteome("Acal").local, ">Acal_1\nMKLVAAAAAAAAAAAAAA\n")
        write(self.cfg / "tree.nwk", "(Drer,(Acal,Mzeb));\n")
        out = drift(self.L, "v1", self.cfg)
        self.assertEqual(len(out), 2)
        self.assertIn("proteome Acal: differs from the frozen v1 (a new annotation or assembly: make a new release)", out)
        self.assertIn("config/tree.nwk: differs from the frozen v1", out)
        self.L.proteome("Drer").local.unlink()
        self.assertTrue(any(x.startswith("proteome Drer: not on this machine") for x in drift(self.L, "v1", self.cfg)))


class TestCheckOfAFrozenRelease(Base):
    def findings(self, **kw):
        return {(f.method, f.item): f for f in check(self.L, self.catalog(), **kw)}

    def test_a_frozen_release_without_its_archive_fails(self):
        f = self.findings()
        self.assertFalse(f[("v1", "frozen archive")].ok)
        self.assertIn("not on this machine", f[("v1", "frozen archive")].detail)

    def test_intact_archive_and_no_drift_pass_then_drift_fails(self):
        self.frozen()
        f = self.findings()
        self.assertEqual((f[("v1", "frozen archive")].ok, f[("v1", "frozen archive")].detail), (True, "intact"))
        self.assertEqual((f[("v1", "inputs against the frozen archive")].ok,
                          f[("v1", "inputs against the frozen archive")].detail), (True, "the same as when frozen"))
        write(self.L.proteome("Mzeb").local, ">Mzeb_1\nMKLVAAAAAAAAAAAA\n")
        f = self.findings()
        self.assertFalse(f[("v1", "inputs against the frozen archive")].ok)
        self.assertIn("proteome Mzeb: differs from the frozen v1", f[("v1", "inputs against the frozen archive")].detail)
        (self.L.release("v1").root().local / "genes/Drer.genes.tsv").unlink()
        self.assertFalse(self.findings()[("v1", "frozen archive")].ok)

    def test_a_release_that_is_building_has_no_archive_findings(self):
        make_config(self.root, status="building")
        self.assertEqual({k for k in self.findings() if k[0] == "v1"}, set())
        self.assertEqual({f.method for f in check(self.L, self.catalog())},
                         {"v1/TreeFam", "v1/OrthoFinder", "v1/Homology"})
        make_config(self.root, status="superseded")
        self.assertEqual(check(self.L, self.catalog()), [])                          # superseded: not checked by default
        self.assertTrue(check(self.L, self.catalog(), release="v1"))                  # unless asked for

    def test_fetch_brings_the_archive_from_the_cloud(self):
        e = Env()
        try:
            cloud_L = e.a
            for p in ("Mzeb", "Acal", "Drer"):
                proteome(cloud_L, p, seq="MKLV" + "A" * len(p))
            cat = Catalog.build(cloud_L, self.cfg)
            freeze(cloud_L, cat, "v1", say=lambda m: None)
            e.store(cloud_L).push("releases/v1")
            shutil.rmtree(cloud_L.release("v1").root().local)                          # a machine without the archive
            shutil.rmtree(cloud_L.state_dir())
            f = {(x.method, x.item): x for x in check(cloud_L, Catalog.build(cloud_L, self.cfg))}
            self.assertFalse(f[("v1", "frozen archive")].ok)
            f = {(x.method, x.item): x for x in check(cloud_L, Catalog.build(cloud_L, self.cfg), fetch=True,
                                                      store=e.store(cloud_L))}
            self.assertTrue(f[("v1", "frozen archive")].ok, f[("v1", "frozen archive")].detail)
            self.assertTrue((cloud_L.release("v1").proteome("Mzeb").local).is_file())
        finally:
            e.close()


class TestArchiveIsImmutableInTheStore(unittest.TestCase):
    def test_push_pull_and_refusal(self):
        e = Env()
        try:
            cfg = make_config(e.a.local_root.parent / "x")
            for p in ("Mzeb", "Acal", "Drer"):
                proteome(e.a, p, seq="MKLV" + "A" * len(p))
            freeze(e.a, Catalog.build(e.a, cfg), "v1", say=lambda m: None)
            self.assertEqual(e.store(e.a).push("releases/v1"), {"releases/v1": "uploaded"})
            self.assertTrue((e.cloud / "Releases/v1/release.json").is_file())
            e.store(e.b).pull("releases/v1")                                           # another machine gets all of it
            self.assertEqual(verify_archive(e.b, "v1"), [])
            self.assertEqual(load_manifest(e.b, "v1")["release"], "v1")
            touch_later(e.a.release("v1").proteome("Mzeb").local)                      # an edit to a frozen archive
            with self.assertRaises(TransferError) as err:
                e.store(e.a).push("releases/v1")
            self.assertIn("REFUSED", str(err.exception))
            self.assertIn("immutable", str(err.exception))
        finally:
            e.close()


class TestRunRecordsAgainstTheArchive(Base):
    """run.json says which proteomes the run used; a frozen release can tell whether that is what it froze."""

    def write_run(self, part, species, **extra):
        rec = dict(release="v1", method=part, method_version="x", parameters={}, started="a", finished="b",
                   coordinator_version=__version__, catalog_md5="m", species=species)
        rec.update(extra)
        schemas.write_run(self.L.results("v1", part).run_record().local, **rec)

    def test_manifest_records_the_proteomes_checksum(self):
        self.frozen()
        m = load_manifest(self.L, "v1")
        self.assertEqual(m["proteomes_md5"], self.catalog().proteomes_md5(["Mzeb", "Acal", "Drer"]))

    def test_matching_mismatching_missing_and_foreign_runs(self):
        from cichlid_data.releases import run_records
        self.assertEqual(run_records(self.L, "v1"), ([], []))                    # not frozen: nothing to compare
        self.frozen()
        C = self.catalog()
        self.write_run("TreeFam", ["Mzeb", "Acal"], proteomes_md5=C.proteomes_md5(["Mzeb", "Acal"]))
        self.assertEqual(run_records(self.L, "v1"), ([], []))
        self.write_run("OrthoFinder", ["Mzeb", "Acal", "Drer"], proteomes_md5="0" * 32)
        problems, notes = run_records(self.L, "v1")
        self.assertEqual(problems, ["OrthoFinder: its proteomes differ from the ones frozen in v1"])
        self.write_run("OrthoFinder", ["Mzeb", "Acal", "Drer"], proteomes_md5=C.proteomes_md5(["Mzeb", "Acal", "Drer"]))
        self.write_run("Compara", ["Mzeb"])                                      # written before proteomes_md5 existed
        problems, notes = run_records(self.L, "v1")
        self.assertEqual(problems, [])
        self.assertEqual(notes, ["Compara: run.json has no proteomes_md5, so its proteomes cannot be checked against "
                                 "the frozen v1"])
        self.write_run("Homology", ["Mzeb", "Hsap"], proteomes_md5="0" * 32)
        self.assertEqual(run_records(self.L, "v1")[0],
                         ["Homology: the run used species the frozen v1 does not have: Hsap"])
        write(self.L.results("v1", "Analysis").run_record().local, "{not json")
        self.assertIn("Analysis/run.json is not valid JSON", run_records(self.L, "v1")[0])

    def test_a_newer_annotation_does_not_disturb_a_run_of_the_frozen_release(self):
        from cichlid_data.releases import run_records
        self.frozen()
        self.write_run("TreeFam", ["Mzeb"], proteomes_md5=self.catalog().proteomes_md5(["Mzeb"]))
        write(self.L.proteome("Mzeb").local, ">Mzeb_1\nMKLVAAAAAAAAAAAAAA\n")            # the annotation improves later
        self.assertEqual(run_records(self.L, "v1"), ([], []))                    # the run still matches what was frozen


class TestReleaseCommand(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.cfg = make_config(self.root / "x", status="building")
        os.environ["CICHLID_DATA_CONFIG_DIR"] = str(self.cfg)
        self.s = Settings(local_root=self.root / "data", cloud_remote="fake:")
        from cichlid_data import Layout
        self.L = Layout(self.root / "data", "fake:")
        for p in ("Mzeb", "Acal", "Drer"):
            proteome(self.L, p, seq="MKLV" + "A" * len(p))

    def tearDown(self):
        os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)
        self._t.cleanup()

    def run_cli(self, *args, fail=False):
        f = io.StringIO()
        with redirect_stdout(f):
            if fail:
                with self.assertRaises(SystemExit) as e:
                    cli(list(args), settings=self.s)
                return f.getvalue() + "\n" + str(e.exception)
            cli(list(args), settings=self.s)
        return f.getvalue()

    def set_status(self, status):
        make_config(self.root / "x", status=status)

    def test_show(self):
        out = self.run_cli("release", "show")
        self.assertIn("v1 (building): test release", out)
        self.assertIn("delivers: Mzeb", out)
        self.assertIn("input compara: Ensembl Compara 116", out)
        self.assertIn("TreeFam: 2 species", out)                                  # the project species
        self.assertIn("OrthoFinder: 3 species, species tree tree.nwk", out)
        self.assertIn("Homology: 0 species", out)
        self.assertIn("archive: not frozen on this machine", out)
        self.assertIn("v1 (building)", self.run_cli("release", "show", "v1"))
        self.assertIn("unknown release 'v7'", self.run_cli("release", "show", "v7", fail=True))

    def test_freeze_validate_and_show_the_archive(self):
        self.assertIn("status 'building'", self.run_cli("release", "freeze", "v1", fail=True))
        self.assertIn("needs a release", self.run_cli("release", "freeze", fail=True))
        self.set_status("frozen")
        out = self.run_cli("release", "freeze", "v1")
        self.assertIn("Froze v1: 3 species", out)
        self.assertIn("cichlid-data push releases/v1", out)
        self.assertIn("already frozen", self.run_cli("release", "freeze", "v1", fail=True))
        self.assertIn("v1 (frozen): valid", self.run_cli("release", "validate"))
        self.assertIn("archive: frozen ", self.run_cli("release", "show", "v1"))
        # the archive is damaged: validate fails, and says what
        write(self.L.release("v1").proteome("Mzeb").local, ">Mzeb_1\nXXXXXXXXXXXX\n")
        out = self.run_cli("release", "validate", "v1", fail=True)
        self.assertIn("v1 (frozen): 1 problem(s)", out)
        self.assertIn("proteomes/Mzeb.faa", out)
        # a newer annotation is a note, not an error
        shutil.copy2(self.L.proteome("Mzeb").local, self.L.release("v1").proteome("Mzeb").local)
        write(self.L.proteome("Acal").local, ">Acal_1\nMKLVAAAAAAAAAAAAAAAAAAAA\n")
        out = self.run_cli("release", "validate", "v1")
        self.assertIn("v1 (frozen): valid", out)
        self.assertIn("note: proteome Acal: differs from the frozen v1", out)
        self.assertIn("--fast", self.run_cli("release", "--help", fail=True))

    def test_release_validate_checks_the_runs_against_the_archive(self):
        self.set_status("frozen")
        self.run_cli("release", "freeze", "v1")
        from cichlid_data import Catalog as _C
        good = _C.build(self.L, self.cfg_dir()).proteomes_md5(["Mzeb", "Acal"])
        rec = dict(release="v1", method="TreeFam", method_version="x", parameters={}, started="a", finished="b",
                   coordinator_version=__version__, catalog_md5="m", species=["Mzeb", "Acal"])
        run = self.L.results("v1", "TreeFam").run_record().local
        schemas.write_run(run, **rec)                                           # no proteomes_md5: a note only
        out = self.run_cli("release", "validate", "v1")
        self.assertIn("v1 (frozen): valid", out)
        self.assertIn("note: TreeFam: run.json has no proteomes_md5", out)
        schemas.write_run(run, proteomes_md5=good, **rec)
        out = self.run_cli("release", "validate", "v1")
        self.assertIn("v1 (frozen): valid", out)
        self.assertNotIn("proteomes_md5", out)
        schemas.write_run(run, proteomes_md5="1" * 32, **rec)                   # a run on other proteomes: a problem
        out = self.run_cli("release", "validate", "v1", fail=True)
        self.assertIn("v1 (frozen): 1 problem(s)", out)
        self.assertIn("TreeFam: its proteomes differ from the ones frozen in v1", out)

    def cfg_dir(self):
        return self.root / "x" / "config"

    def test_status_and_archive_must_agree(self):
        self.set_status("frozen")
        self.assertIn("frozen archive", self.run_cli("release", "validate", fail=True) + "frozen archive")
        out = self.run_cli("release", "validate", fail=True)
        self.assertIn("not frozen on this machine", out)                       # frozen status, no archive
        self.run_cli("release", "freeze", "v1")
        self.set_status("building")
        out = self.run_cli("release", "validate", fail=True)
        self.assertIn("has an archive but its status is 'building'", out)

    def test_where_and_check_commands(self):
        self.assertEqual(self.run_cli("where", "release", "v1").strip(), str(self.root / "data" / "Releases" / "v1"))
        self.assertEqual(self.run_cli("where", "results", "v1", "OrthoFinder", "families").strip(),
                         str(self.root / "data" / "Results" / "v1" / "OrthoFinder" / "families"))
        out = self.run_cli("check", "--release", "v1", "--part", "OrthoFinder")
        self.assertIn("v1/OrthoFinder: 4 of 4 inputs ready", out)                 # 3 proteomes and the species tree
        self.assertIn("unknown release 'v7'", self.run_cli("check", "--release", "v7", fail=True))


class TestReleaseInRecords(unittest.TestCase):
    """run.json and done.json say which release they belong to, and it must match where they are."""

    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()

    RUN = dict(method="TreeFam", method_version="x", parameters={}, started="a", finished="b",
               coordinator_version="0.9.0", catalog_md5="m", species=["Mzeb"])

    def test_run_record(self):
        run = self.d / "Results" / "v1" / "TreeFam" / "run.json"
        schemas.write_run(run, release="v1", **self.RUN)
        self.assertEqual(schemas.validate(run), [])
        with self.assertRaises(schemas.SchemaError):
            schemas.write_run(self.d / "x.json", **self.RUN)                          # release is required
        schemas.write_run(run, release="v2", **self.RUN)
        self.assertEqual(schemas.validate(run), ["release is v2 but the file is in the results of v1"])
        schemas.write_run(run, release="2026", **self.RUN)
        self.assertEqual(schemas.validate(run), ["release '2026' must be like v1"])
        elsewhere = self.d / "run.json"
        schemas.write_run(elsewhere, release="v2", **self.RUN)
        self.assertEqual(schemas.validate(elsewhere), [])                              # not in a results folder: no cross-check
        wrong_part = self.d / "Results" / "v1" / "Stuff" / "run.json"
        schemas.write_run(wrong_part, release="v2", **self.RUN)
        self.assertEqual(schemas.validate(wrong_part), [])

    def family(self, release, done):
        d = self.d / "Results" / release / "OrthoFinder" / "families" / "OG0000001"
        write(d / "tree.nwk", "((a,b),c);")
        write(d / "combined.aln", ">a\nMK\n>b\nMK\n>c\nMK\n")
        write(d / "done.json", json.dumps(done))
        return d

    def test_done_json(self):
        self.assertEqual(schemas.validate(self.family("v1", {"release": "v1"})), [])
        self.assertEqual(schemas.validate(self.family("v1", {})), [])                  # still optional
        self.assertEqual(schemas.validate(self.family("v1", {"release": "v2"})),
                         ["done.json says release v2 but the family is in the results of v1"])
        self.assertEqual(schemas.validate(self.family("v1", {"release": "latest"})),
                         ["done.json: release 'latest' must be like v1"])
        # the release check and the rooting check are independent
        probs = schemas.validate(self.family("v1", {"release": "v2", "rooted": "yes"}))
        self.assertEqual(len(probs), 2)
        self.assertEqual(schemas.release_of_location(self.family("v3", {}), "family"), "v3")
        self.assertIsNone(schemas.release_of_location(self.d, "family"))


if __name__ == "__main__":
    unittest.main()
