"""One-time move from the old layout (Data/, Rebuild2026/, Compara2026/, Combined2026/) to the
new one (Inputs/, Results/, Analyses/), locally and in the cloud.

    cichlid-data migrate              # dry run: print every step
    cichlid-data migrate --apply      # do it
    cichlid-data migrate --apply --cleanup   # and delete the old cloud tars that were repacked

Cloud copies whose packing is unchanged are moved server-side (no re-upload). Units whose
internal paths change are repacked from the local copy and pushed: the genomes (now one shard
per species), the TreeFam tree files, and each run's assignment tables (treefam/ -> assign/,
<prefix>.treefam.tsv -> <prefix>.assign.tsv). Old cloud indexes and local sync state are carried
over under the new unit names, so moved units are recognized as in sync.
"""
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

OLD_RUNS = {"Rebuild2026": "TreeFam2026", "Compara2026": "Compara2026", "Combined2026": "Combined2026"}
PARTS_MOVED = ("names", "orthology", "residues", "site")


@dataclass
class Step:
    where: str          # local | cloud
    action: str         # move | repack | rename
    src: str
    dst: str
    note: str = ""

    def __str__(self):
        return f"{self.where:<5} {self.action:<6} {self.src} -> {self.dst}" + (f"   ({self.note})" if self.note else "")


class Migration:
    def __init__(self, layout, backend):
        self.L, self.backend = layout, backend
        self.root = Path(layout.local_root)

    # ---- what exists
    def _local(self, rel):
        return (self.root / rel).exists()

    def _cloud(self, rel):
        parent, _, name = rel.rpartition("/")
        return any(x.rstrip("/") == name for x in self.backend.list(parent))

    def plan(self):
        steps, repush = [], []

        def both(src, dst, note=""):
            if self._local(src):
                steps.append(Step("local", "move", src, dst, note))
            if self._cloud(src):
                steps.append(Step("cloud", "move", src, dst, note))

        tf, cm = "Inputs/Families/TreeFam9", "Inputs/Families/Compara116"
        both("Data/hmm_lib", f"{tf}/hmm")
        if self._local("Data/OriginalTreeFiles"):
            steps.append(Step("local", "move", "Data/OriginalTreeFiles", f"{tf}/trees"))
            if self._cloud("Data/OriginalTreeFiles.tar"):
                repush.append(("families/TreeFam9/trees", ["Data/OriginalTreeFiles.tar"]))
        both("Data/protein_to_species.txt", f"{tf}/protein_to_species.txt")
        names = set()
        if self._local("Data/Compara116"):
            names |= {p.name for p in (self.root / "Data/Compara116").iterdir()}
        names |= {x.rstrip("/") for x in self.backend.list("Data/Compara116")}
        for n in sorted(names):
            src = f"Data/Compara116/{n}"
            dst = f"{cm}/{n}" if n.startswith("prepared") else f"{cm}/dumps/{n}"
            both(src, dst)
        both("Data/Annotations", "Inputs/Annotations")
        both("Data/Expression", "Inputs/Expression")
        for sub in ("coding", "inversions"):
            both(f"Data/Blumer2025/{sub}", f"Analyses/Variants/{sub}", "results, not inputs")
        both("Data/Blumer2025", "Inputs/Populations/Blumer2025")
        both("Data/Esur", "Inputs/Proteins/local_sources/Esur")
        # runs
        for old, run in OLD_RUNS.items():
            R = f"Results/{run}"
            if old == "Rebuild2026":
                if self._local(f"{old}/genomes"):
                    steps.append(Step("local", "move", f"{old}/genomes", "Inputs/Proteins/downloads"))
                    if self._cloud(f"{old}/genomes.tar"):
                        repush.append(("proteins/downloads", [f"{old}/genomes.tar"]))
                if self._local(f"{old}/proteomes"):
                    steps.append(Step("local", "move", f"{old}/proteomes", "Inputs/Proteins/proteomes"))
                if self._cloud(f"{old}/proteomes.tar"):
                    steps.append(Step("cloud", "move", f"{old}/proteomes.tar", "Inputs/Proteins/proteomes.tar"))
            if self._local(f"{old}/treefam"):
                steps.append(Step("local", "move", f"{old}/treefam", f"{R}/assign"))
                steps.append(Step("local", "rename", f"{R}/assign", f"{R}/assign", "*.treefam*.tsv -> *.assign*.tsv"))
                if self._cloud(f"{old}/treefam.tar"):
                    repush.append((f"results/{run}/assign", [f"{old}/treefam.tar"]))
            if self._local(f"{old}/families"):
                steps.append(Step("local", "move", f"{old}/families", f"{R}/families"))
            for shard in self.backend.list(f"{old}/families"):
                if shard.endswith(".tar"):
                    steps.append(Step("cloud", "move", f"{old}/families/{shard}", f"{R}/families/{shard}"))
            for part in PARTS_MOVED:
                if self._local(f"{old}/{part}"):
                    steps.append(Step("local", "move", f"{old}/{part}", f"{R}/{part}"))
                if self._cloud(f"{old}/{part}.tar"):
                    steps.append(Step("cloud", "move", f"{old}/{part}.tar", f"{R}/{part}.tar"))
            if self._local(f"{old}/analysis"):
                for a in sorted(p.name for p in (self.root / old / "analysis").iterdir() if p.is_dir()):
                    steps.append(Step("local", "move", f"{old}/analysis/{a}", f"Analyses/{a}"))
                    repush.append((f"analyses/{a}", []))
        return steps, repush

    # ---- index and state carried over
    def _remap(self, old, key):
        run = OLD_RUNS[old]
        if key == "proteomes" and old == "Rebuild2026":
            return "proteins/proteomes"
        if key in PARTS_MOVED:
            return f"results/{run}/{key}"
        if key.startswith("families/"):
            return f"results/{run}/{key}"
        return None                              # genomes, treefam: repacked, recomputed on push

    def carry_over(self, apply):
        new_index = json.loads(self.backend.read_text(".cichlid-data/index.json") or "{}")
        state_path = self.L.state_dir() / "state.json"
        new_state = json.loads(state_path.read_text()) if state_path.is_file() else {}
        moved = 0
        for old in OLD_RUNS:
            idx = json.loads(self.backend.read_text(f"{old}/index.json") or "{}")
            for k, v in idx.items():
                nk = self._remap(old, k)
                if nk and nk not in new_index:
                    new_index[nk] = v
                    moved += 1
            sp = self.root / old / ".sync_state.json"
            if sp.is_file():
                for k, v in json.loads(sp.read_text()).items():
                    nk = self._remap(old, k)
                    if nk and nk not in new_state:
                        new_state[nk] = v
        if apply:
            self.backend.write_text(".cichlid-data/index.json", json.dumps(new_index, indent=2, sort_keys=True))
            self.L.state_dir().mkdir(parents=True, exist_ok=True)
            state_path.write_text(json.dumps(new_state, indent=2, sort_keys=True))
        return moved

    # ---- apply
    def run(self, store_factory, apply=False, cleanup=False, say=print):
        steps, repush = self.plan()
        if not steps and not repush:
            say("Nothing to migrate: no old-layout folders found locally or in the cloud.")
            return steps, repush
        for s in steps:
            say(("  " if apply else "  would ") + str(s))
        for unit, old in repush:
            say(f"  {'' if apply else 'would '}repack and push {unit}" + (f" (replacing {', '.join(old)})" if old else ""))
        if not apply:
            say(f"Index entries to carry over: {self.carry_over(False)}. Rerun with --apply to migrate.")
            return steps, repush
        for s in steps:
            src, dst = self.root / s.src, self.root / s.dst
            if s.where == "local" and s.action == "move":
                if dst.exists():
                    if dst.is_dir() and src.is_dir():              # merge folder contents
                        for e in src.iterdir():
                            if not (dst / e.name).exists():
                                shutil.move(str(e), str(dst / e.name))
                        if not any(src.iterdir()):
                            src.rmdir()
                    continue
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
            elif s.where == "local" and s.action == "rename":
                for f in sorted(dst.iterdir()):
                    if ".treefam" in f.name:
                        f.rename(f.with_name(f.name.replace(".treefam", ".assign")))
            elif s.where == "cloud":
                if not self._cloud(s.dst):
                    self.backend.move(s.src, s.dst)
        n = self.carry_over(True)
        say(f"Carried over {n} index entries.")
        store = store_factory()
        for unit, old in repush:
            res = store.push(unit)
            say(f"  pushed {unit}: " + ", ".join(sorted(set(res.values()))))
            if cleanup:
                for o in old:
                    self.backend.delete(o)
                    say(f"  deleted old {o}")
        return steps, repush
