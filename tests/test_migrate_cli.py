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
             "Combined2026/analysis/tissues/summary.tsv", "Combined2026/analysis/teleost_copies.tsv",
             "Combined2026/analysis/teleost_copies.png", "Data/Species/legacy2018.txt", "Data/Mystery/x.txt"]
    for f in files:
        write(root / f, f)
    # cloud: inputs copied as they are, the TreeFam trees as a tar
    for d in ("Data/hmm_lib", "Data/Annotations", "Data/Expression", "Data/Blumer2025"):
        shutil.copytree(root / d, cloud / d)
    write(cloud / "Data/protein_to_species.txt", "p2s")
    write(cloud / "Data/Compara116/Compara.116.protein_default.aa.fasta.gz", "dump")
    make_tar(cloud / "Data/Compara116/prepared.tar", root / "Data/Compara116", ["prepared"])
    make_tar(cloud / "Data/OriginalTreeFiles.tar", root / "Data", ["OriginalTreeFiles"])
    # the 2018 resource's files, in the cloud only
    for f in ("Data/Proteins/Metriaclima_zebra.BROADMZ2.pep.all.fa", "Data/Species/homo_sapiens.fa",
              "Data/AllCichlidProteins.fa", "Data/AllCichlidProteins.fa.fai", "Data/AllCichlidProteins.treefam",
              "Data/NewTreeFiles/TF101001.added.nex"):
        write(cloud / f, "2018")
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

    def test_dry_run_changes_nothing_and_reports(self):
        before = sorted(str(p) for p in self.e.cloud.rglob("*"))
        out = []
        Migration(self.e.a, self.e.store(self.e.a).backend).run(lambda: None, apply=False, say=out.append)
        import re
        text = re.sub(r" +", " ", "\n".join(out))
        self.assertEqual(sorted(str(p) for p in self.e.cloud.rglob("*")), before)
        self.assertTrue((self.e.a.local_root / "Rebuild2026/genomes").is_dir())
        self.assertIn("cloud move Rebuild2026/proteomes.tar -> Inputs/Proteins/proteomes.tar", text)
        design = text.split("Left in place by design")[1].split("Not moved and not expected")[0]
        unexpected = text.split("Not moved and not expected")[1]
        for item in ("local Rebuild2026/", "cloud Rebuild2026/", "local Compara2026/", "local Combined2026/",
                     "local Data/Blumer2025/coding/", "cloud Data/Blumer2025/inversions/"):
            self.assertIn(item, design, item)
        self.assertIn("local Data/Mystery/", unexpected)                     # unknown folders are reported
        self.assertNotIn("Data/Species/", unexpected)                         # the 2018 files are archived
        self.assertIn("cloud move Data/NewTreeFiles -> Archive/Project2018/NewTreeFiles", text)
        for absent in ("Results/", "Analyses/", "families/TF101.tar ->", "Rebuild2026/treefam", "Compara2026/treefam"):
            self.assertNotIn(absent, text.split("Left in place")[0], absent)

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
                    "Inputs/Proteins/downloads/Mzeb/protein.faa", "Inputs/Proteins/proteomes/Mzeb.faa"):
            self.assertTrue((r / rel).is_file(), rel)
        # all results left exactly where they were, locally and in the cloud
        for rel in ("Rebuild2026/treefam/Mzeb.treefam.tsv", "Rebuild2026/families/TF101001/tree.nwk",
                    "Rebuild2026/orthology/Mzeb.orthologs.tsv", "Compara2026/treefam/Mzeb.treefam.tsv",
                    "Combined2026/site/site.sqlite", "Combined2026/analysis/tissues/summary.tsv",
                    "Data/Blumer2025/coding/coding_snps.tsv.gz", "Data/Blumer2025/inversions/x.tsv"):
            self.assertTrue((r / rel).is_file(), rel)
        for rel in ("Rebuild2026/families/TF101.tar", "Rebuild2026/treefam.tar", "Rebuild2026/orthology.tar",
                    "Compara2026/treefam.tar", "Combined2026/site.tar", "Rebuild2026/index.json"):
            self.assertTrue((cloud / rel).is_file(), rel)
        self.assertFalse((r / "Results").exists())
        self.assertFalse((r / "Analyses").exists())
        # inputs in the cloud: moved server-side or repacked; old repacked tars removed by --cleanup
        for rel in ("Inputs/Families/TreeFam9/hmm/TreeFam9", "Inputs/Proteins/proteomes.tar",
                    "Inputs/Families/Compara116/prepared.tar", "Inputs/Populations/Blumer2025/chr10.vcf.gz",
                    "Inputs/Proteins/downloads/Mzeb.tar", "Inputs/Families/TreeFam9/trees.tar"):
            self.assertTrue((cloud / rel).is_file(), rel)
        for rel in ("Archive/Project2018/Proteins/Metriaclima_zebra.BROADMZ2.pep.all.fa", "Archive/Project2018/Species/homo_sapiens.fa",
                    "Archive/Project2018/AllCichlidProteins.fa", "Archive/Project2018/AllCichlidProteins.fa.fai",
                    "Archive/Project2018/AllCichlidProteins.treefam", "Archive/Project2018/NewTreeFiles/TF101001.added.nex"):
            self.assertTrue((cloud / rel).is_file(), rel)
        self.assertEqual(sorted(p.name for p in (cloud / "Data").iterdir()), ["Blumer2025"])   # only the variant analyses remain
        self.assertEqual(sorted(p.name for p in (cloud / "Data/Blumer2025").iterdir()), ["coding", "inversions"])
        self.assertEqual(sorted(p.name for p in (r / "Data").iterdir()), ["Blumer2025", "Mystery"])  # unknown left alone
        self.assertTrue((r / "Archive/Project2018/Species/legacy2018.txt").is_file())
        for gone in ("Rebuild2026/genomes.tar", "Data/OriginalTreeFiles.tar", "Rebuild2026/proteomes.tar"):
            self.assertFalse((cloud / gone).exists(), gone)
        st = {k: s for k, s, _ in self.e.store(a).status("proteins/*")}
        self.assertEqual((st["proteins/proteomes"], st["proteins/downloads/Mzeb"]), ("in sync", "in sync"))
        b = self.e.b
        self.e.store(b).pull("proteins/*")
        self.assertTrue(b.proteome("Mzeb").local.is_file())
        # an empty folder left in Data/ (as after an interrupted run): not "unexpected", removed by a rerun
        (r / "Data/Compara116").mkdir(parents=True, exist_ok=True)
        (cloud / "Data/Compara116").mkdir(parents=True, exist_ok=True)
        dry = []
        Migration(a, backend).run(lambda: None, apply=False, say=dry.append)
        self.assertFalse(any("Data/Compara116" in x for x in dry))
        out = []
        Migration(a, backend).run(lambda: self.e.store(a), apply=True, say=out.append)
        self.assertIn("Nothing to migrate", out[-1])
        self.assertFalse((r / "Data/Compara116").exists())
        self.assertFalse((cloud / "Data/Compara116").exists())


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


