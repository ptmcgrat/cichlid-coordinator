import io
import json
import shutil
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from cichlid_data import Settings
from cichlid_data.cli import main as cli
from cichlid_data.migrate import Migration
from cichlid_data.store import fingerprint, make_tar
from helpers import Env, write


def old_layout(root, cloud):
    """A miniature of the current layout, locally and in the cloud as the old sync left it."""
    files = ["Data/hmm_lib/TreeFam9", "Data/OriginalTreeFiles/TF1.aln", "Data/protein_to_species.txt",
             "Data/Compara116/Compara.116.protein_default.aa.fasta.gz", "Data/Compara116/prepared/trees/C1.nhx",
             "Data/Annotations/GO/go-basic.obo", "Data/Expression/HPA/x.tsv", "Data/Blumer2025/chr10.vcf.gz",
             "Data/Blumer2025/coding/coding_snps.tsv.gz", "Data/Blumer2025/inversions/x.tsv",
             "Rebuild2026/genomes/Mzeb/protein.faa", "Rebuild2026/proteomes/Mzeb.faa",
             "Rebuild2026/treefam/Mzeb.treefam.tsv", "Rebuild2026/families/TF101001/tree.nwk",
             "Rebuild2026/orthology/Mzeb.orthologs.tsv", "Compara2026/treefam/Mzeb.treefam.tsv",
             "Compara2026/families/C116_012345/tree.nwk", "Combined2026/site/site.sqlite",
             "Combined2026/analysis/tissues/summary.tsv"]
    for f in files:
        write(root / f, f)
    # cloud: inputs copied as they are, the TreeFam trees as a tar
    for d in ("Data/hmm_lib", "Data/Annotations", "Data/Expression", "Data/Blumer2025"):
        shutil.copytree(root / d, cloud / d)
    write(cloud / "Data/protein_to_species.txt", "p2s")
    write(cloud / "Data/Compara116/Compara.116.protein_default.aa.fasta.gz", "dump")
    make_tar(cloud / "Data/Compara116/prepared.tar", root / "Data/Compara116", ["prepared"])
    make_tar(cloud / "Data/OriginalTreeFiles.tar", root / "Data", ["OriginalTreeFiles"])
    # results: one tar per unit, paths relative to the old run folder; index and sync state
    for old, units in (("Rebuild2026", ["genomes", "proteomes", "treefam", "orthology"]),
                       ("Compara2026", ["treefam"]), ("Combined2026", ["site"])):
        base, index = root / old, {}
        for u in units:
            make_tar(cloud / old / f"{u}.tar", base, [u])
            index[u] = {"fingerprint": fingerprint(base, [u])}
        fams = [p.name for p in (base / "families").iterdir()] if (base / "families").is_dir() else []
        for fam in fams:
            shard = fam[:-3]
            (cloud / old / "families").mkdir(parents=True, exist_ok=True)
            make_tar(cloud / old / "families" / f"{shard}.tar", base, [f"families/{fam}"])
            index[f"families/{shard}"] = {"fingerprint": fingerprint(base, [f"families/{fam}"])}
        write(cloud / old / "index.json", json.dumps(index))
        write(base / ".sync_state.json", json.dumps({k: v["fingerprint"] for k, v in index.items()}))


class TestMigration(unittest.TestCase):
    def setUp(self):
        self.e = Env()
        old_layout(self.e.a.local_root, self.e.cloud)

    def tearDown(self):
        self.e.close()

    def test_dry_run_changes_nothing(self):
        before = sorted(str(p) for p in self.e.cloud.rglob("*"))
        out = []
        steps, repush = Migration(self.e.a, self.e.store(self.e.a).backend).run(lambda: None, apply=False, say=out.append)
        self.assertTrue(any("cloud move  Rebuild2026/families/TF101.tar -> Results/TreeFam2026/families/TF101.tar"
                            in s.replace("  ", " ").replace("would ", "") or "TF101.tar" in s for s in out))
        self.assertEqual(sorted(str(p) for p in self.e.cloud.rglob("*")), before)
        self.assertTrue((self.e.a.local_root / "Rebuild2026/genomes").is_dir())

    def test_apply(self):
        a, cloud = self.e.a, self.e.cloud
        backend = self.e.store(a).backend
        Migration(a, backend).run(lambda: self.e.store(a), apply=True, cleanup=True, say=lambda *_: None)
        r = a.local_root
        for rel in ("Inputs/Families/TreeFam9/hmm/TreeFam9", "Inputs/Families/TreeFam9/trees/TF1.aln",
                    "Inputs/Families/TreeFam9/protein_to_species.txt",
                    "Inputs/Families/Compara116/dumps/Compara.116.protein_default.aa.fasta.gz",
                    "Inputs/Families/Compara116/prepared/trees/C1.nhx", "Inputs/Annotations/GO/go-basic.obo",
                    "Inputs/Expression/HPA/x.tsv", "Inputs/Populations/Blumer2025/chr10.vcf.gz",
                    "Analyses/Variants/coding/coding_snps.tsv.gz", "Inputs/Proteins/downloads/Mzeb/protein.faa",
                    "Inputs/Proteins/proteomes/Mzeb.faa", "Results/TreeFam2026/assign/Mzeb.assign.tsv",
                    "Results/Compara2026/assign/Mzeb.assign.tsv", "Results/TreeFam2026/families/TF101001/tree.nwk",
                    "Results/Combined2026/site/site.sqlite", "Analyses/tissues/summary.tsv"):
            self.assertTrue((r / rel).is_file(), rel)
        self.assertFalse((r / "Data/Blumer2025/coding").exists())
        # cloud: moved server-side, or repacked and pushed
        for rel in ("Inputs/Families/TreeFam9/hmm/TreeFam9", "Inputs/Proteins/proteomes.tar",
                    "Inputs/Families/Compara116/prepared.tar", "Results/TreeFam2026/families/TF101.tar",
                    "Results/TreeFam2026/orthology.tar", "Results/Combined2026/site.tar",
                    "Inputs/Populations/Blumer2025/chr10.vcf.gz", "Analyses/Variants/coding/coding_snps.tsv.gz",
                    "Inputs/Proteins/downloads/Mzeb.tar", "Results/TreeFam2026/assign.tar",
                    "Results/Compara2026/assign.tar", "Inputs/Families/TreeFam9/trees.tar", "Analyses/tissues.tar"):
            self.assertTrue((cloud / rel).is_file(), rel)
        for gone in ("Rebuild2026/genomes.tar", "Rebuild2026/treefam.tar", "Data/OriginalTreeFiles.tar"):
            self.assertFalse((cloud / gone).exists(), gone)                         # --cleanup
        # moved units are in sync straight away, with their old fingerprints
        st = {k: s for k, s, _ in self.e.store(a).status("*")}
        for k in ("proteins/proteomes", "results/TreeFam2026/families/TF101", "results/TreeFam2026/orthology",
                  "results/Combined2026/site", "proteins/downloads/Mzeb", "results/TreeFam2026/assign"):
            self.assertEqual(st.get(k), "in sync", k)
        # a second machine can pull the new layout
        b = self.e.b
        self.e.store(b).pull("results/TreeFam2026/*")
        self.assertTrue((b.local_root / "Results/TreeFam2026/families/TF101001/tree.nwk").is_file())
        self.assertTrue((b.local_root / "Results/TreeFam2026/assign/Mzeb.assign.tsv").is_file())
        # running it again finds nothing to do
        out = []
        Migration(a, backend).run(lambda: self.e.store(a), apply=True, say=out.append)
        self.assertIn("Nothing to migrate", out[0])


class TestCli(unittest.TestCase):
    def setUp(self):
        self.e = Env()
        self.s = Settings(local_root=self.e.a.local_root, cloud_remote="fake:")

    def tearDown(self):
        self.e.close()

    def run_cli(self, *args):
        f = io.StringIO()
        with redirect_stdout(f):
            cli(list(args), backend=self.e.store(self.e.a).backend, settings=self.s)
        return f.getvalue()

    def test_where_env_units_push_status(self):
        self.assertEqual(self.run_cli("where", "proteome", "Mzeb").strip(),
                         str(self.e.a.local_root / "Inputs/Proteins/proteomes/Mzeb.faa"))
        self.assertEqual(self.run_cli("where", "results", "Compara2026", "orthology", "--cloud").strip(),
                         "fake:/Results/Compara2026/orthology")
        self.assertIn("export CICHLID_PROTEOMES=", self.run_cli("env"))
        self.assertIn("families/TreeFam9/hmm", self.run_cli("units", "families/*"))
        write(self.e.a.proteome("Mzeb").local, ">Mzeb_1\nM\n")
        self.assertIn("proteins/proteomes: would upload", self.run_cli("push", "proteins/proteomes", "--dry-run"))
        self.assertIn("proteins/proteomes: uploaded", self.run_cli("push", "proteins/proteomes"))
        self.assertIn("in sync", self.run_cli("status", "proteins/*"))


if __name__ == "__main__":
    unittest.main()
