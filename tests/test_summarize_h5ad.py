"""Offline tests for summarize_h5ad.py (core maths on synthetic data; script on a stand-in atlas)."""
import gzip
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cichlid_data.downloaders import h5ad as sh
from helpers import coord_env
FAKE_RCLONE = ""        # unused: backups go through the coordinator's stand-in rclone (coord_env)

# 5 cells x 3 genes; cell types: A A B B (unlabelled)
COUNTS = np.array([[10, 0, 0],      # A: lib 10 -> gene0 10000 cp10k
                   [5, 5, 0],       # A: lib 10 -> 5000, 5000
                   [0, 0, 4],       # B: lib 4  -> gene2 10000
                   [0, 2, 2],       # B: lib 4  -> 5000, 5000
                   [9, 9, 9]], float)
TYPES = ["A", "A", "B", "B", None]
STAGES = ["10hpf"] * 2 + ["2dpf"] * 2 + ["5dpf"]


def chunks(M, size):
    M = sparse.csr_matrix(M)
    return [(s, M[s:s + size]) for s in range(0, M.shape[0], size)]


class TestCore(unittest.TestCase):
    def test_means_and_fractions(self):
        types, mean, pct, n, st = sh.aggregate(chunks(COUNTS, 5), TYPES, STAGES, min_stage_cells=1)
        self.assertEqual(types, ["A", "B"])
        np.testing.assert_allclose(mean[0], [7500, 2500, 0])
        np.testing.assert_allclose(mean[1], [0, 2500, 7500])
        np.testing.assert_allclose(pct[0], [1.0, 0.5, 0.0])
        np.testing.assert_allclose(pct[1], [0.0, 0.5, 1.0])
        self.assertEqual(list(n), [2, 2])                                   # unlabelled cell skipped
        self.assertEqual(st, [["10hpf"], ["2dpf"]])

    def test_chunking_and_log_input(self):
        a = sh.aggregate(chunks(COUNTS, 5), TYPES)
        b = sh.aggregate(chunks(COUNTS, 2), TYPES)                           # chunks split mid type
        np.testing.assert_allclose(a[1], b[1])
        c = sh.aggregate(chunks(np.log1p(COUNTS), 3), TYPES, already_log=True)
        np.testing.assert_allclose(a[1], c[1], rtol=1e-9)

    def test_gene_columns(self):
        ids = [f"ENSDARG{i:011d}" for i in range(150)]
        var = pd.DataFrame({"gene_symbols": [f"g{i}" for i in range(150)]}, index=ids)
        self.assertEqual(sh.gene_columns(var.index, var)[1][:2], ["g0", "g1"])
        var2 = pd.DataFrame({"gene_ids": ids}, index=[f"g{i}" for i in range(150)])
        got = sh.gene_columns(var2.index, var2)
        self.assertEqual((got[0][0], got[1][0]), ("ENSDARG00000000000", "g0"))


class FakeAnnData:
    def __init__(self):
        self.X = sparse.csr_matrix(COUNTS)
        self.n_obs, self.n_vars = COUNTS.shape
        self.obs = pd.DataFrame({"cell_type": TYPES, "timepoint": STAGES})
        self.var = pd.DataFrame({"gene_ids": ["ENSDARG00000000001", "ENSDARG00000000002", "ENSDARG00000000003"]},
                                index=["avpr1aa", "bmp2b", "irx1a"])
        self.var_names = self.var.index
        self.layers = {}


OBO = """format-version: 1.2

[Term]
id: ZFA:0000001
name: anatomical structure

[Term]
id: ZFA:0000396
name: nervous system
is_a: ZFA:0000001

[Term]
id: ZFA:0001127
name: visual system
is_a: ZFA:0000001

[Term]
id: ZFA:0000152
name: retina
relationship: part_of ZFA:0001127

[Term]
id: ZFA:0009248
name: retinal neuron
is_a: ZFA:0009999
relationship: part_of ZFA:0000152

[Term]
id: ZFA:0009999
name: neuron
relationship: part_of ZFA:0000396

[Term]
id: ZFA:0000434
name: skin
"""


class TestSystemsAndStages(unittest.TestCase):
    def test_system_assignment(self):
        names, parents = sh.parse_obo(OBO)
        self.assertEqual(sh.system_of("ZFA:0009248", names, parents), "Eye")       # nearest: retina, not neuron
        self.assertEqual(sh.system_of("ZFA:0009999", names, parents), "Nervous system")
        self.assertEqual(sh.system_of("ZFA:0000434", names, parents), "Skin")
        self.assertEqual(sh.system_of("ZFA:0000001", names, parents), "Other")

    def test_stage_hours(self):
        self.assertEqual([sh.stage_hours(s) for s in ("10hpf", "2dpf", "10dpf", "5 somites")], [10, 48, 240, None])


def seurat_export(d):
    """A tiny export_seurat.R output: 3 genes x 4 cells, two clusters, two stage ranges."""
    from scipy import sparse as sp
    counts = np.array([[10, 0, 4, 0],      # rho
                       [0, 5, 0, 2],       # krt4
                       [0, 5, 4, 2]], float)   # actb
    m = sp.csc_matrix(counts)
    d.mkdir(parents=True)
    m.data.astype(np.float64).tofile(d / "counts_x.bin")
    m.indices.astype(np.int32).tofile(d / "counts_i.bin")
    m.indptr.astype(np.int32).tofile(d / "counts_p.bin")
    (d / "dims.txt").write_text(f"3\n4\n{m.nnz}\n")
    (d / "genes.txt").write_text("rho\nkrt4\nactb\n")
    (d / "cells.txt").write_text("c1\nc2\nc3\nc4\n")
    (d / "meta.tsv").write_text("cell\tclust\tstage.nice\n"
                                "c1\teye.1\t24-34hpf\nc2\tepi.3\t24-34hpf\nc3\teye.1\t3-4hpf\nc4\tepi.3\t3-4hpf\n")
    ann = d.parent / "cluster_annotations.csv"
    ann.write_text("clust,identity.super,tissue.name\neye.1,photoreceptor,eye\nepi.3,periderm,epidermis\n")
    return ann


class TestSeuratExport(unittest.TestCase):
    def test_daniocell_style_export(self):
        with tempfile.TemporaryDirectory() as d:
            exp = Path(d) / "export"
            ann = seurat_export(exp)
            out = Path(d) / "out"
            sh.main([str(exp), "--dataset", "Daniocell", "--celltype-col", "ann_identity.super",
                     "--stage-col", "stage.nice", "--by-stage", "--system-col", "ann_tissue.name",
                     "--annotations", str(ann), "--join", "clust=clust", "--outdir", str(out), "--no-backup"])
            with gzip.open(out / "expression_by_stage.tsv.gz", "rt") as fh:
                rows = {(r.split("\t")[1], r.split("\t")[2], r.split("\t")[3]): r.rstrip("\n").split("\t")[4:]
                        for r in list(fh)[1:]}
            # c1: rho 10 of 10 counts -> 10000 cp10k; c3 (3-4hpf): rho 4 of 8 -> 5000
            self.assertEqual(rows[("rho", "photoreceptor", "24-34hpf")], ["1e+04", "100"])
            self.assertEqual(rows[("rho", "photoreceptor", "3-4hpf")], ["5000", "100"])
            self.assertEqual(rows[("krt4", "periderm", "3-4hpf")], ["5000", "100"])
            types = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in (out / "cell_types.tsv").read_text().splitlines()[1:]}
            self.assertEqual((types["photoreceptor"][4], types["periderm"][4]), ("eye", "epidermis"))
            with gzip.open(out / "expression.tsv.gz", "rt") as fh:
                self.assertIn("\trho\tphotoreceptor\t7500\t100", fh.read())    # gene_id empty: symbols only
            buf = io.StringIO()
            with mock.patch("sys.stdout", buf):
                sh.main([str(exp), "--inspect", "--annotations", str(ann), "--join", "clust=clust"])
            self.assertIn("ann_tissue.name", buf.getvalue())

    def test_combined_celltype_labels(self):
        import pandas as pd
        obs = pd.DataFrame({"s": ["notochord", "notochord", "neural", ""], "u": ["sheath cells", "", "neural", "x"]})
        self.assertEqual(sh.celltype_values(obs, "s+u"), ["notochord (sheath cells)", "notochord", "neural", ""])
        self.assertEqual(sh.celltype_values(obs, "s"), ["notochord", "notochord", "neural", ""])

    def test_normalized_export(self):
        # an export of the 'data' slot (log1p of counts per 10,000) gives the same summary as raw counts
        with tempfile.TemporaryDirectory() as d:
            exp = Path(d) / "export"
            ann = seurat_export(exp)
            from scipy import sparse as sp
            raw = sp.csc_matrix(np.array([[10, 0, 4, 0], [0, 5, 0, 2], [0, 5, 4, 2]], float))
            lib = np.asarray(raw.sum(axis=0)).ravel()
            norm = sp.csc_matrix(raw.multiply(1e4 / lib))
            norm.data = np.log1p(norm.data)
            norm.data.astype(np.float64).tofile(exp / "counts_x.bin")
            (exp / "source.txt").write_text("RNA\ndata\n")
            out = Path(d) / "out"
            sh.main([str(exp), "--dataset", "Daniocell", "--celltype-col", "ann_identity.super",
                     "--annotations", str(ann), "--join", "clust=clust", "--outdir", str(out), "--no-backup"])
            with gzip.open(out / "expression.tsv.gz", "rt") as fh:
                self.assertIn("\trho\tphotoreceptor\t7500\t100", fh.read())
            self.assertTrue(json.loads((out / "summary.json").read_text())["log_normalized_input"])

    def test_stage_ranges_order(self):
        self.assertEqual([sh.stage_hours(s) for s in ("3-4hpf", "24-34hpf", "5dpf")], [3, 24, 120])
        # Daniocell stage groups: bare hour ranges with padding
        self.assertEqual([sh.stage_label(s) for s in ("  3-4", "10-12", "120", "5 somites", "2dpf")],
                         ["3-4 hpf", "10-12 hpf", "120 hpf", "5 somites", "2dpf"])
        self.assertEqual([sh.stage_hours(sh.stage_label(s)) for s in ("  3-4", "10-12", "120")], [3, 10, 120])


class FakeStage:
    def __init__(self, counts, ids, syms, types, zfa=None):
        self.X = sparse.csr_matrix(np.array(counts, float))
        self.n_obs, self.n_vars = self.X.shape
        self.obs = pd.DataFrame({"cell_type": types, **({"zfa": zfa} if zfa else {})})
        self.var = pd.DataFrame({"gene_ids": ids}, index=syms)
        self.var_names = self.var.index
        self.layers = {}


class TestScript(unittest.TestCase):
    def test_per_stage_files_merge(self):
        files = {
            # 2dpf: genes (avpr1aa, bmp2b); one A cell: lib 10 -> avpr1aa 10000
            "zf_atlas_2dpf_v1_release.h5ad": FakeStage([[10, 0]], ["ENSDARG00000000001", "ENSDARG00000000002"],
                                                       ["avpr1aa", "bmp2b"], ["A"]),
            # 5dpf: genes in another order plus irx1a; one A cell (bmp2b 5000, avpr1aa 5000), one B cell (irx1a)
            "zf_atlas_5dpf_v1_release.h5ad": FakeStage([[5, 5, 0], [0, 0, 3]],
                                                       ["ENSDARG00000000002", "ENSDARG00000000001", "ENSDARG00000000003"],
                                                       ["bmp2b", "avpr1aa", "irx1a"], ["A", "B"]),
        }
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            with mock.patch.object(sh, "open_h5ad", side_effect=lambda p, *a: files[Path(p).name]):
                sh.main(list(files) + ["--dataset", "Zebrahub", "--celltype-col", "cell_type", "--outdir", str(out),
                                       "--no-backup"])
            with gzip.open(out / "expression.tsv.gz", "rt") as fh:
                rows = {(r.split("\t")[1], r.split("\t")[2]): r.rstrip("\n").split("\t")[3:] for r in list(fh)[1:]}
            self.assertEqual(rows[("avpr1aa", "A")], ["7500", "100"])     # (10000 + 5000) / 2, both cells
            self.assertEqual(rows[("bmp2b", "A")], ["2500", "50"])
            self.assertEqual(rows[("irx1a", "B")], ["1e+04", "100"])
            types = (out / "cell_types.tsv").read_text()
            self.assertIn("A\t2\t", types)
            info = json.loads((out / "summary.json").read_text())
            self.assertEqual((info["genes"], info["stage_col"]), (3, "from file names"))
        self.assertEqual(sh.stage_from_name("zf_atlas_10hpf_v1_release.h5ad"), "10hpf")
        self.assertEqual(sh.stage_from_name("zf_atlas_2dpf_v1_release.h5ad"), "2dpf")

    def test_by_stage_and_systems(self):
        files = {
            "zf_atlas_2dpf_v1_release.h5ad": FakeStage([[10, 0], [0, 5]], ["ENSDARG00000000001", "ENSDARG00000000002"],
                                                       ["rho", "skin1"], ["retinal neuron", "skin cell"],
                                                       ["ZFA:0009248", "ZFA:0000434"]),
            "zf_atlas_10hpf_v1_release.h5ad": FakeStage([[4, 4]], ["ENSDARG00000000001", "ENSDARG00000000002"],
                                                        ["rho", "skin1"], ["retinal neuron"], ["ZFA:0009248"]),
        }
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            obo = Path(d) / "zfa.obo"
            obo.write_text(OBO)
            with mock.patch.object(sh, "open_h5ad", side_effect=lambda p, *a: files[Path(p).name]):
                sh.main(list(files) + ["--dataset", "Zebrahub", "--celltype-col", "cell_type", "--by-stage",
                                       "--zfa-col", "zfa", "--zfa-obo", str(obo), "--outdir", str(out), "--no-backup"])
            with gzip.open(out / "expression_by_stage.tsv.gz", "rt") as fh:
                rows = {(r.split("\t")[1], r.split("\t")[2], r.split("\t")[3]): r.rstrip("\n").split("\t")[4:]
                        for r in list(fh)[1:]}
            self.assertEqual(rows[("rho", "retinal neuron", "2dpf")], ["1e+04", "100"])
            self.assertEqual(rows[("rho", "retinal neuron", "10hpf")], ["5000", "100"])
            self.assertNotIn(("skin1", "retinal neuron", "2dpf"), rows)                # not expressed there
            types = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in (out / "cell_types.tsv").read_text().splitlines()[1:]}
            self.assertEqual(types["retinal neuron"][3:], ["ZFA:0009248", "Eye"])
            self.assertEqual(types["skin cell"][4], "Skin")
            per = (out / "cell_types_by_stage.tsv").read_text()
            self.assertIn("retinal neuron\t10hpf\t1", per)

    def test_summarize_and_backup(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "cloud").mkdir()
            (d / "rclone").write_text(FAKE_RCLONE)
            (d / "rclone").chmod(0o755)
            env = coord_env(d)
            with mock.patch.dict(os.environ, env), mock.patch.object(sh, "open_h5ad", return_value=FakeAnnData()):
                sh.main(["atlas.h5ad", "--dataset", "Zebrahub", "--celltype-col", "cell_type",
                         "--stage-col", "timepoint"])
                buf = io.StringIO()
                with mock.patch("sys.stdout", buf):
                    sh.main(["atlas.h5ad", "--inspect"])
            self.assertIn("timepoint", buf.getvalue())
            out = d / "local" / "Inputs" / "Expression" / "Zebrahub"
            with gzip.open(out / "expression.tsv.gz", "rt") as fh:
                rows = [l.rstrip("\n").split("\t") for l in fh][1:]
            self.assertIn(["ENSDARG00000000001", "avpr1aa", "A", "7500", "100"], rows)
            self.assertNotIn("irx1a", [r[1] for r in rows if r[2] == "A"])   # not expressed in A: left out
            info = json.loads((out / "summary.json").read_text())
            self.assertEqual((info["cell_types"], info["log_normalized_input"]), (2, False))
            cloud = d / "cloud" / "Inputs/Expression/Zebrahub"
            self.assertTrue((cloud / "cell_types.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
