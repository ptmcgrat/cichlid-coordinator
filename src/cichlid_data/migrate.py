"""One-time move of the project's INPUTS from the old layout (Data/, Rebuild2026/genomes and
proteomes) to the new one (Inputs/), locally and in the cloud.

    cichlid-data migrate              # dry run: print every step
    cichlid-data migrate --apply      # do it
    cichlid-data migrate --apply --cleanup   # and delete the old cloud tars that were repacked

Moved: downloads (genomes and annotations), proteomes, TreeFam 9, Compara 116, annotations,
expression, population data (VCFs, tables, reference) and local species sources.
Left in place by design: all results (Rebuild2026/, Compara2026/, Combined2026/: phylogenies,
calls, browser site, analyses) and the variant analyses in Data/Blumer2025/coding and inversions.
The TreeFam and Compara runs are redone with the new species set, so their old results stay
where they are, for comparison.

Cloud copies whose packing is unchanged are moved server-side (no re-upload). Two units are
repacked from the local copy and pushed: the downloads (now one shard per species) and the
TreeFam tree files (folder renamed). Old cloud index entries and local sync state for the
proteomes are carried over, so they are in sync immediately.
"""
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

OLD = ("Data", "Rebuild2026", "Compara2026", "Combined2026")
BLUMER_ANALYSES = ("coding", "inversions")


@dataclass
class Step:
    where: str          # local | cloud
    action: str         # move
    src: str
    dst: str

    def __str__(self):
        return f"{self.where:<5} {self.action:<6} {self.src} -> {self.dst}"


class Migration:
    def __init__(self, layout, backend):
        self.L, self.backend = layout, backend
        self.root = Path(layout.local_root)

    def _local(self, rel):
        return (self.root / rel).exists()

    def _cloud(self, rel):
        parent, _, name = rel.rpartition("/")
        return any(x.rstrip("/") == name for x in self.backend.list(parent))

    def _names(self, rel):
        names = {e.name for e in (self.root / rel).iterdir()} if self._local(rel) else set()
        return names | {x.rstrip("/") for x in self.backend.list(rel)}

    def plan(self):
        steps, repush = [], []

        def both(src, dst):
            if self._local(src):
                steps.append(Step("local", "move", src, dst))
            if self._cloud(src):
                steps.append(Step("cloud", "move", src, dst))

        tf, cm = "Inputs/Families/TreeFam9", "Inputs/Families/Compara116"
        both("Data/hmm_lib", f"{tf}/hmm")
        if self._local("Data/OriginalTreeFiles"):
            steps.append(Step("local", "move", "Data/OriginalTreeFiles", f"{tf}/trees"))
            if self._cloud("Data/OriginalTreeFiles.tar"):
                repush.append(("families/TreeFam9/trees", ["Data/OriginalTreeFiles.tar"]))
        both("Data/protein_to_species.txt", f"{tf}/protein_to_species.txt")
        for n in sorted(self._names("Data/Compara116")):
            both(f"Data/Compara116/{n}", f"{cm}/{n}" if n.startswith("prepared") else f"{cm}/dumps/{n}")
        both("Data/Annotations", "Inputs/Annotations")
        both("Data/Expression", "Inputs/Expression")
        for n in sorted(self._names("Data/Blumer2025") - set(BLUMER_ANALYSES)):
            both(f"Data/Blumer2025/{n}", f"Inputs/Populations/Blumer2025/{n}")
        both("Data/Esur", "Inputs/Proteins/local_sources/Esur")
        if self._local("Rebuild2026/genomes"):
            steps.append(Step("local", "move", "Rebuild2026/genomes", "Inputs/Proteins/downloads"))
            if self._cloud("Rebuild2026/genomes.tar"):
                repush.append(("proteins/downloads", ["Rebuild2026/genomes.tar"]))
        if self._local("Rebuild2026/proteomes"):
            steps.append(Step("local", "move", "Rebuild2026/proteomes", "Inputs/Proteins/proteomes"))
        if self._cloud("Rebuild2026/proteomes.tar"):
            steps.append(Step("cloud", "move", "Rebuild2026/proteomes.tar", "Inputs/Proteins/proteomes.tar"))
        return steps, repush

    def leftovers(self, steps):
        """-> (left in place by design, unexpected)"""
        moved = {(s.where, s.src) for s in steps}
        handled = {"Data/OriginalTreeFiles", "Data/OriginalTreeFiles.tar", "Rebuild2026/genomes", "Rebuild2026/genomes.tar"}
        design, unexpected = [], []
        for where in ("local", "cloud"):
            for top in ("Rebuild2026", "Compara2026", "Combined2026"):
                if (where == "local" and self._local(top)) or (where == "cloud" and self._cloud(top)):
                    design.append((where, f"{top}/ (results: phylogenies, calls, site, analyses)"))
            for sub in BLUMER_ANALYSES:
                rel = f"Data/Blumer2025/{sub}"
                if (where == "local" and self._local(rel)) or (where == "cloud" and self._cloud(rel)):
                    design.append((where, f"{rel}/ (variant analyses)"))
            names = ([e.name + ("/" if e.is_dir() else "") for e in (self.root / "Data").iterdir()]
                     if where == "local" and self._local("Data") else self.backend.list("Data") if where == "cloud" else [])
            for n in sorted(names):
                rel = f"Data/{n.rstrip('/')}"
                if n.startswith(".") or rel in handled or (where, rel) in moved or \
                        any(w == where and src.startswith(rel + "/") for w, src in moved) or rel == "Data/Blumer2025":
                    continue
                unexpected.append((where, rel + ("/" if n.endswith("/") else "")))
        return design, unexpected

    def carry_over(self, apply):
        """The proteomes' cloud index entry and sync state, under the new unit name."""
        new_index = json.loads(self.backend.read_text(".cichlid-data/index.json") or "{}")
        state_path = self.L.state_dir() / "state.json"
        new_state = json.loads(state_path.read_text()) if state_path.is_file() else {}
        old_index = json.loads(self.backend.read_text("Rebuild2026/index.json") or "{}")
        n = 0
        if "proteomes" in old_index and "proteins/proteomes" not in new_index:
            new_index["proteins/proteomes"] = old_index["proteomes"]
            n += 1
        sp = self.root / "Rebuild2026/.sync_state.json"
        if sp.is_file():
            old_state = json.loads(sp.read_text())
            if "proteomes" in old_state:
                new_state.setdefault("proteins/proteomes", old_state["proteomes"])
        if apply:
            self.backend.write_text(".cichlid-data/index.json", json.dumps(new_index, indent=2, sort_keys=True))
            self.L.state_dir().mkdir(parents=True, exist_ok=True)
            state_path.write_text(json.dumps(new_state, indent=2, sort_keys=True))
        return n

    def run(self, store_factory, apply=False, cleanup=False, say=print):
        steps, repush = self.plan()
        if not steps and not repush:
            say("Nothing to migrate: no old-layout inputs found locally or in the cloud.")
            return steps, repush
        for s in steps:
            say(("  " if apply else "  would ") + str(s))
        for unit, old in repush:
            say(f"  {'' if apply else 'would '}repack and push {unit} (replacing {', '.join(old)})")
        design, unexpected = self.leftovers(steps)
        if design:
            say("Left in place by design:")
            for where, rel in design:
                say(f"  {where:<5} {rel}")
        if unexpected:
            say("Not moved and not expected (decide about these; they stay where they are):")
            for where, rel in unexpected:
                say(f"  {where:<5} {rel}")
        if not apply:
            say("Rerun with --apply to migrate.")
            return steps, repush
        for s in steps:
            src, dst = self.root / s.src, self.root / s.dst
            if s.where == "local":
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
            elif not self._cloud(s.dst):
                self.backend.move(s.src, s.dst)
        say(f"Carried over {self.carry_over(True)} index entries.")
        store = store_factory()
        for unit, old in repush:
            res = store.push(unit)
            say(f"  pushed {unit}: " + ", ".join(sorted(set(res.values()))))
            if cleanup:
                for o in old:
                    self.backend.delete(o)
                    say(f"  deleted old {o}")
        return steps, repush
