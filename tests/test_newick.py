"""Five tests of newick (prepare_reference, parse, midpoint_root), from the cichlid-trees repository.

They were written for the pipeline's newick.py, which is now cichlid_data.newick, and test it as stage 5 uses it
(prepare_reference: prune a TreeFam tree to given rows, unroot, resolve polytomies; duplicates and too-small trees are
refused; a 5,000-leaf caterpillar round-trips) and as the calling code uses it (midpoint_root). Unchanged except the
import. Source: tests/test_newick_contract.py at cichlid-trees 0.5.1; originally TestNewick in test_stage5.py and
TestMidpoint in test_stage6.py of the legacy-2026-10 tag.
"""
import unittest

from cichlid_data import newick

NHX = ("((ENSP1:0.1[&&NHX:S=homo_sapiens:B=100],ENSDARP2:0.2)90:0.05,"
       "((ENSMUSP3:0.1,FBpp4:0.4,CE5:0.3):0.1,ENSGALP6:0.2):0.1);")


class TestNewick(unittest.TestCase):
    def test_prepare_reference(self):
        txt, rep = newick.prepare_reference(NHX, {"ENSP1", "ENSDARP2", "ENSMUSP3", "FBpp4", "CE5"})
        root = newick.parse(txt)
        self.assertEqual(sorted(l.name for l in root.leaves()),
                         ["CE5", "ENSDARP2", "ENSMUSP3", "ENSP1", "FBpp4"])
        self.assertEqual(len(root.children), 3)                          # unrooted trifurcation
        for n in root.iter():
            if n.children and n is not root:
                self.assertEqual(len(n.children), 2)                     # strictly binary
            if n is not root:
                self.assertIsNotNone(n.length)
        self.assertEqual(rep["pruned"], ["ENSGALP6"])
        self.assertNotIn("90", txt)                                      # support labels removed

    def test_duplicates_and_too_small(self):
        with self.assertRaises(ValueError):
            newick.prepare_reference("((A,B),(A,C));", {"A", "B", "C"})
        with self.assertRaises(ValueError):
            newick.prepare_reference(NHX, {"ENSP1", "ENSDARP2", "CE5"})

    def test_roundtrip_and_deep_tree(self):
        t = "((A:1,B:2)x:0.5,(C:1,'D d':1):0.3);"
        self.assertEqual(newick.write(newick.parse(t)), "((A:1,B:2):0.5,(C:1,D d:1):0.3);")
        s = "L0:1"
        for i in range(1, 5000):
            s = f"({s},L{i}:1):1"
        self.assertEqual(len(newick.parse(s + ";").leaves()), 5000)


class TestMidpoint(unittest.TestCase):
    def test_root_on_longest_path_midpoint(self):
        r = newick.midpoint_root(newick.parse("(A:1,B:1,(C:1,(D:1,E:9):1):1);"))
        d = newick.depths(r)
        leaf = {l.name: d[l] for l in r.leaves()}
        self.assertAlmostEqual(leaf["E"], 6.0)
        self.assertAlmostEqual(leaf["A"], 6.0)
        self.assertEqual(len(r.children), 2)
        self.assertEqual(sorted(leaf), list("ABCDE"))

    def test_zero_lengths_fall_back_to_balanced_split(self):
        r = newick.midpoint_root(newick.parse("((A,B),(C,D),(E,F));"))
        self.assertEqual(len(r.children), 2)
        self.assertEqual(len(r.leaves()), 6)


if __name__ == "__main__":
    unittest.main()
