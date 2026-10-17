"""Catalog.md5(): the catalog_md5 of a run record."""
import hashlib
import tempfile
import unittest
from pathlib import Path

from cichlid_data import Catalog, Layout, __version__, schemas
from helpers import write
from test_species_tree import config


class TestCatalogMd5(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.root = Path(self._t.name)
        self.L = Layout(self.root / "data", "fake:")
        self.cfg = config(self.root)

    def tearDown(self):
        self._t.cleanup()

    def file_md5(self):
        return hashlib.md5(self.L.catalog().local.read_bytes()).hexdigest()

    def test_a_catalog_that_was_only_built_has_no_file_to_checksum(self):
        C = Catalog.build(self.L, self.cfg)
        with self.assertRaises(FileNotFoundError) as e:
            C.md5()
        self.assertIn("catalog build", str(e.exception))
        with self.assertRaises(FileNotFoundError):
            C.md5(self.L)                                                       # a layout with no catalog.json either

    def test_md5_of_the_file_after_write_and_after_load(self):
        C = Catalog.build(self.L, self.cfg)
        C.write(self.L)
        self.assertEqual(C.path, self.L.catalog().local)
        self.assertEqual(C.md5(), self.file_md5())
        self.assertRegex(C.md5(), r"^[0-9a-f]{32}$")
        loaded = Catalog.load(self.L)
        self.assertEqual(loaded.path, self.L.catalog().local)
        self.assertEqual(loaded.md5(), C.md5())

    def test_built_catalog_can_name_the_layout(self):
        Catalog.build(self.L, self.cfg).write(self.L)
        fresh = Catalog.build(self.L, self.cfg)                                  # built again, not written
        self.assertEqual(fresh.md5(self.L), self.file_md5())

    def test_the_checksum_follows_the_file(self):
        C = Catalog.build(self.L, self.cfg)
        C.write(self.L)
        before = C.md5()
        write(self.cfg / "tree.nwk", "(Drer,(Mzeb,(Esur,Acal)));\n")             # a different species tree
        C2 = Catalog.build(self.L, self.cfg)
        C2.write(self.L)
        self.assertNotEqual(C2.md5(), before)                                    # its checksum is in the catalog
        self.assertEqual(C2.md5(), self.file_md5())
        self.assertEqual(C.md5(), C2.md5())                                      # read from the file at call time

    def test_a_run_record_built_from_the_helper_validates(self):
        Catalog.build(self.L, self.cfg).write(self.L)
        C = Catalog.load(self.L)
        run = self.L.results("v1", "OrthoFinder").run_record().local
        schemas.write_run(run, release="v1", method="OrthoFinder", method_version="x", parameters={"-M": "msa"},
                          started="t0", finished="t1", coordinator_version=__version__, catalog_md5=C.md5(),
                          species=[s.prefix for s in C.species(release="v1", part="OrthoFinder")])
        self.assertEqual(schemas.validate(run), [])
        self.assertIn(C.md5(), run.read_text())


if __name__ == "__main__":
    unittest.main()
