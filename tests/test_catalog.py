import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Catalog, Layout, Settings
from cichlid_data.cli import main as cli
from cichlid_data.config import config_dir
from helpers import Env, releases_toml, write

SETS = '''
[[organism]]
prefix = "Mzeb"
species = "Maylandia zebra"
role = "project"
source = "ncbi_refseq"
[[organism]]
prefix = "Acal"
species = "Astatotilapia calliptera"
role = "project"
source = "ncbi_refseq"
[[organism]]
prefix = "Esur"
species = "Etroplus suratensis"
role = "project"
source = "authors"
[[organism]]
prefix = "Drer"
species = "Danio rerio"
common = "zebrafish"
role = "reference"
source = "ncbi_refseq"
[[organism]]
prefix = "Amca"
species = "Amia calva"
role = "candidate"
'''
RELEASES = releases_toml({"TreeFam": 'queries = "project"', "OrthoFinder": 'species = ["Mzeb", "Acal", "Esur", "Drer"]'})
MANIFEST = "prefix\tspecies\taccession\tsource\nMzeb\tMaylandia zebra\tGCF_041146795.1\tncbi\nAcal\tAstatotilapia calliptera\tGCF_964374335.1\tncbi\n"
PANEL = "# panel\nspecies\trole\tgroup\nhomo_sapiens\treference\tmammal\nesox_lucius\tpanel\tteleost\nsalmo_salar\texclude\tfish\n"


def setup(root):
    cfg = root / "config"
    write(cfg / "species_sets.toml", SETS)
    write(cfg / "releases.toml", RELEASES)
    write(cfg / "species.tsv", MANIFEST)
    write(cfg / "compara_panel.tsv", PANEL)
    L = Layout(root / "data", "fake:")
    write(L.proteome("Mzeb").local, ">Mzeb_1\nM\n>Mzeb_2\nK\n")
    write(L.genes_table("Mzeb").local, "gene_id\tsymbol\nMzeb_1\ta\nMzeb_2\tb\n")
    write(L.download("Acal").local / "protein.faa", ">x\nM\n")
    (L.download("Acal").local / "download.json").write_text(json.dumps({
        "manifest_row": {"prefix": "Acal", "species": "Astatotilapia calliptera", "accession": "GCF_964374335.1",
                         "source": "ncbi"}, "ncbi_metadata": {"assembly_name": "fAstCal68.1"},
        "files": {"protein": {"file": "protein.faa", "md5": "x", "bytes": 5}}, "downloaded": "2026-10-07T00:00:00Z"}))
    write(L.download("Acal").local / "assembly_report.json", "{}")
    write(L.local_source("Esur").local / "Genes.fasta.gz", "g")
    return L, cfg


class TestCatalog(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L, self.cfg = setup(self.root)

    def tearDown(self):
        self._t.cleanup()

    def test_build_and_query(self):
        C = Catalog.build(self.L, self.cfg)
        m, a, e, d = C["Mzeb"], C["Acal"], C["Esur"], C["Drer"]
        self.assertEqual((m.status, m.counts, m.accession, m.in_manifest),
                         ("proteome ready", {"proteins": 2, "genes": 2}, "GCF_041146795.1", True))
        self.assertEqual(m.proteome, "Inputs/Proteins/proteomes/Mzeb.faa")
        self.assertEqual((a.status, a.version, a.retrieved), ("downloaded", "fAstCal68.1", "2026-10-07T00:00:00Z"))
        self.assertEqual((e.status, e.in_manifest), ("local source only", False))
        self.assertEqual(d.status, "not downloaded")
        self.assertEqual(C.releases["v1"]["parts"]["TreeFam"]["species"], ["Mzeb", "Acal", "Esur"])   # the project species
        self.assertEqual(m.used_by, ["v1/TreeFam", "v1/OrthoFinder"])
        self.assertEqual(d.used_by, ["v1/OrthoFinder"])
        self.assertEqual([s.prefix for s in C.missing("v1", "OrthoFinder")], ["Acal", "Esur", "Drer"])
        self.assertEqual([s.prefix for s in C.missing("v1")], ["Acal", "Esur", "Drer"])           # the whole release
        self.assertEqual([s.prefix for s in C.species(role="reference")], ["Drer"])
        self.assertEqual([s["species"] for s in C.panel], ["homo_sapiens", "esox_lucius"])   # excluded rows dropped
        errors, warnings = C.validate()
        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["Esur: a project species not in species.tsv, so Stage 1 won't download it"])
        # write and load
        C.write(self.L)
        C2 = Catalog.load(self.L)
        self.assertEqual(C2["Mzeb"].to_dict(), m.to_dict())
        self.assertEqual(C2.releases, C.releases)
        self.assertIn("species_sets.toml", C2.meta["config"])
        self.assertIn("releases.toml", C2.meta["config"])

    def test_validation_errors(self):
        (self.cfg / "releases.toml").write_text(RELEASES.replace('species = ["Mzeb", "Acal", "Esur", "Drer"]',
                                                                   'species = ["Mzeb", "Amca", "Nope"]'))
        (self.cfg / "species_sets.toml").write_text(SETS + '[[organism]]\nprefix = "Bad_1"\nspecies = "x"\nrole = "pet"\n')
        errors, _ = Catalog.build(self.L, self.cfg).validate()
        joined = "\n".join(errors)
        for bit in ("Amca is a candidate", "Nope has no organism entry", "Bad_1: a prefix must be letters",
                    "unknown role 'pet'"):
            self.assertIn(bit, joined)

    def test_repository_configuration(self):
        C = Catalog.build(self.L, config_dir())
        errors, warnings = C.validate()
        self.assertEqual(errors, [])
        r = C.releases["v1"]
        self.assertEqual(len(r["parts"]["OrthoFinder"]["species"]), 35)
        self.assertEqual(len(C.species(role="project")), 12)
        self.assertEqual(len(r["parts"]["TreeFam"]["species"]), 12)
        self.assertEqual(len(r["parts"]["Compara"]["species"]), 12)
        self.assertEqual(r["deliver"], ["Mzeb", "Acal"])
        self.assertIn(r["status"], ("planned", "building", "frozen", "published", "superseded"))
        self.assertGreater(len(C.panel), 50)


class TestReleases(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L, self.cfg = setup(self.root)

    def tearDown(self):
        self._t.cleanup()

    def errors(self, releases=None, sets=None):
        if releases is not None:
            (self.cfg / "releases.toml").write_text(releases)
        if sets is not None:
            (self.cfg / "species_sets.toml").write_text(sets)
        return Catalog.build(self.L, self.cfg).validate()[0]

    def test_a_valid_release(self):
        C = Catalog.build(self.L, self.cfg)
        r = C.releases["v1"]
        self.assertEqual((r["status"], r["deliver"], r["inputs"]), ("building", ["Mzeb"], {"treefam": "TreeFam 9"}))
        self.assertEqual(list(r["parts"]), ["TreeFam", "OrthoFinder"])
        self.assertEqual(C.validate(release="v1")[0], [])

    def test_status_must_be_one_of_the_five(self):
        for status in ("planned", "building", "frozen", "published", "superseded"):
            self.assertEqual(self.errors(RELEASES.replace('"building"', f'"{status}"')), [], status)
        self.assertEqual(self.errors(RELEASES.replace('"building"', '"done"')),
                         ["v1: status 'done' must be one of planned, building, frozen, published, superseded"])
        self.assertEqual(self.errors(RELEASES.replace('status = "building"\n', "")),
                         ["v1: status '' must be one of planned, building, frozen, published, superseded"])

    def test_release_and_part_names(self):
        errs = self.errors(RELEASES.replace("releases.v1", "releases.first"))
        self.assertTrue(any("a release must be 'v' and a number" in e for e in errs), errs)
        errs = self.errors(RELEASES + '[releases.v1.parts.Browser]\ndescription = "x"\n')
        self.assertEqual(errs, ["v1/Browser: unknown part (one of TreeFam, Compara, OrthoFinder, Homology, Analysis)"])
        self.assertEqual(self.errors(releases_toml({}, deliver=())), ["v1: no parts"])
        self.assertEqual(self.errors(releases_toml({"TreeFam": 'species = []'}, deliver=())), ["v1/TreeFam: no species"])

    def test_parts_without_species_are_fine_for_homology_and_analysis(self):
        text = RELEASES + '[releases.v1.parts.Homology]\ndescription = "calls"\n[releases.v1.parts.Analysis]\n'
        self.assertEqual(self.errors(text), [])
        C = Catalog.build(self.L, self.cfg)
        self.assertEqual(C.releases["v1"]["parts"]["Homology"]["species"], [])
        self.assertEqual(len(C.species(release="v1")), 4)                         # the union over its parts

    def test_delivered_species_must_be_project_species_in_every_tree_part(self):
        errs = self.errors(RELEASES.replace('deliver = ["Mzeb"]', 'deliver = ["Mzeb", "Drer", "Nope"]'))
        joined = "\n".join(errs)
        self.assertIn("delivered species Nope has no organism entry", joined)
        self.assertIn("delivered species Drer is not a project species", joined)
        self.assertIn("delivered species Drer is not in v1/TreeFam", joined)
        errs = self.errors(RELEASES.replace('species = ["Mzeb", "Acal", "Esur", "Drer"]', 'species = ["Acal", "Esur", "Drer"]'))
        self.assertEqual(errs, ["v1: delivered species Mzeb is not in v1/OrthoFinder"])

    def test_two_releases_are_validated_independently(self):
        two = RELEASES + releases_toml({"TreeFam": 'queries = "project"'}, release="v2", status="planned",
                                       deliver=("Mzeb", "Acal"))
        bad = two.replace('status = "planned"', 'status = "later"')
        C = Catalog.build(self.L, self.cfg)
        (self.cfg / "releases.toml").write_text(bad)
        C = Catalog.build(self.L, self.cfg)
        self.assertEqual(C.validate(release="v1")[0], [])
        self.assertEqual(len(C.validate(release="v2")[0]), 1)
        self.assertEqual(len(C.validate()[0]), 1)
        self.assertEqual(C.validate(release="v3")[0], ["unknown release 'v3' (one of v1, v2)"])

    def test_release_name_and_part_lookups(self):
        C = Catalog.build(self.L, self.cfg)
        self.assertEqual(C.release_name(), "v1")                                   # the only live release
        self.assertEqual(C.release_name("v1"), "v1")
        with self.assertRaises(KeyError):
            C.release_name("v9")
        with self.assertRaises(KeyError):
            C.part("v1", "Browser")
        self.assertEqual([s.prefix for s in C.species(part="OrthoFinder")], ["Mzeb", "Acal", "Esur", "Drer"])
        (self.cfg / "releases.toml").write_text(RELEASES + releases_toml({"TreeFam": 'queries = "project"'}, release="v2"))
        C = Catalog.build(self.L, self.cfg)
        with self.assertRaises(KeyError) as e:
            C.release_name()                                                       # two live releases: be explicit
        self.assertIn("give a release", str(e.exception))
        (self.cfg / "releases.toml").write_text(RELEASES.replace('"building"', '"superseded"') +
                                                releases_toml({"TreeFam": 'queries = "project"'}, release="v2"))
        self.assertEqual(Catalog.build(self.L, self.cfg).release_name(), "v2")      # superseded ones do not count

    def test_missing_releases_toml_and_old_catalog_files(self):
        (self.cfg / "releases.toml").unlink()
        with self.assertRaises(FileNotFoundError):
            Catalog.build(self.L, self.cfg)
        write(self.L.catalog().local, json.dumps({"catalog_version": 1, "methods": {}, "organisms": {}}))
        with self.assertRaises(ValueError) as e:
            Catalog.load(self.L)
        self.assertIn("catalog build", str(e.exception))


class TestCatalogCli(unittest.TestCase):
    def test_build_and_missing(self):
        e = Env()
        try:
            import os
            L, cfg = setup(e.a.local_root.parent / "x")
            os.environ["CICHLID_DATA_CONFIG_DIR"] = str(cfg)
            s = Settings(local_root=L.local_root, cloud_remote="fake:")
            f = io.StringIO()
            with redirect_stdout(f):
                cli(["catalog", "build"], backend=e.store(e.a).backend, settings=s)
                cli(["catalog", "missing", "--release", "v1", "--part", "OrthoFinder"], backend=e.store(e.a).backend, settings=s)
            out = f.getvalue()
            self.assertIn("5 organisms (1 with proteomes ready)", out)
            self.assertIn("releases: v1 (building)", out)
            self.assertIn("v1/OrthoFinder: 1 of 4 proteomes ready", out)
            f = io.StringIO()
            with redirect_stdout(f):
                cli(["catalog", "missing"], backend=e.store(e.a).backend, settings=s)          # the only live release
                cli(["catalog", "show", "--release", "v1", "--part", "TreeFam", "--role", "project"],
                    backend=e.store(e.a).backend, settings=s)
            out = f.getvalue()
            self.assertIn("v1: 1 of 4 proteomes ready; missing:", out)
            self.assertIn("v1/TreeFam,v1/OrthoFinder", out)                                       # used_by, for Mzeb
            with self.assertRaises(SystemExit) as err, redirect_stdout(io.StringIO()):
                cli(["catalog", "missing", "--release", "v9"], backend=e.store(e.a).backend, settings=s)
            self.assertIn("unknown release 'v9'", str(err.exception))
            with self.assertRaises(SystemExit) as err, redirect_stdout(io.StringIO()):
                cli(["catalog", "show", "--part", "Browser"], backend=e.store(e.a).backend, settings=s)
            self.assertIn("no part 'Browser'", str(err.exception))
            self.assertIn("Drer", out)
        finally:
            os.environ.pop("CICHLID_DATA_CONFIG_DIR", None)
            e.close()


if __name__ == "__main__":
    unittest.main()


