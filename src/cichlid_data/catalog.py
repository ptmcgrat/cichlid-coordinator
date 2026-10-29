"""The species catalog: every organism, its role, source, files, counts, and the releases using it.

    C = Catalog.build(get_layout())      # from config/ and what is on disk; C.write() saves catalog.json
    C = Catalog.load(get_layout())
    C.species(release="v1", part="OrthoFinder")   # -> [Species]
    C.missing("v1", "OrthoFinder")       # organisms a part needs whose proteome is not ready
    C["Mzeb"].proteome                   # relative path of its proteome
    C.species_tree("v1", "OrthoFinder")  # path of the rooted species tree its run must use (-s)
    C.releases["v1"]["deliver"]          # the species whose calls the release delivers
    C.md5()                              # checksum of catalog.json, for run.json's catalog_md5

A release (v1, v2, ...) is one state of the whole resource, shared by every project: its status
(planned, building, frozen, published, superseded), the species it delivers, the versions of its inputs,
and its parts (TreeFam, Compara, OrthoFinder, Homology, Analysis), each with the species it runs on.
A part may name a species tree (species_tree = "species_tree_orthofinder.nwk", a file in config/). It must
be rooted, fully resolved, and have exactly the part's species as leaves; `validate` checks this, so the
tree and the species list cannot drift apart.

Built from config/species_sets.toml (organisms and roles), config/releases.toml (releases and their parts),
config/species.tsv (what Stage 1 downloads), config/compara_panel.tsv (Compara's reference panel), and the
data: each organism's source record (Stage 1's download.json, or source.json), and whether its download
and proteome exist.
"""
import csv
import hashlib
import json
import re
import time
import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path

from . import __version__, newick
from .config import config_dir
from .layout import RELEASE_RE, RESULT_PARTS
from .provenance import md5sum, read_record

ROLES = ("project", "relative", "reference", "outgroup", "candidate")
PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
CATALOG_VERSION = 2
CLADES_FILE = "species_clades.tsv"
STATUSES = ("planned", "building", "frozen", "published", "superseded")
FROZEN_STATUSES = ("frozen", "published")        # these have an archive in Releases/<release>/
TREE_PARTS = ("TreeFam", "Compara", "OrthoFinder")   # the parts that run on a list of species


@dataclass
class Species:
    prefix: str
    species: str
    role: str
    lineage: str = ""
    common: str = ""
    source: str = ""
    accession: str = ""
    version: str = ""
    retrieved: str = ""
    in_manifest: bool = False             # listed in species.tsv, so Stage 1 downloads it
    status: str = "not downloaded"        # proteome ready | downloaded | local source only | not downloaded
    proteome: str = ""                    # relative paths, when present
    genes_table: str = ""
    proteome_md5: str = ""
    counts: dict = field(default_factory=dict)
    used_by: list = field(default_factory=list)
    notes: str = ""

    def to_dict(self):
        return asdict(self)


def _count_fasta(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def _count_rows(path):
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def read_manifest(path):
    if not Path(path).is_file():
        return {}
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t")]
    return {r["prefix"]: r for r in rows if r.get("prefix")}


def read_panel(path):
    if not Path(path).is_file():
        return []
    with open(path, newline="") as fh:
        rows = list(csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t"))
    return [{"species": r["species"], "role": r["role"], "group": r.get("group", "")} for r in rows
            if r.get("role") in ("reference", "panel")]


def check_species_tree(path, species):
    """-> problems with a species tree file for a method using exactly `species` (empty if fine).

    A species tree must parse, be rooted and fully resolved (the root and every other internal node
    have exactly two children), and have each of `species` as a leaf once, and no other leaves."""
    path = Path(path)
    if not path.is_file():
        return [f"{path.name} not found in {path.parent}"]
    try:
        root = newick.parse(path.read_text())
    except ValueError as e:
        return [f"{path.name} does not parse: {e}"]
    problems, leaves = [], [n.name for n in root.leaves()]
    dup = sorted({x for x in leaves if leaves.count(x) > 1})
    if "" in leaves:
        problems.append("a leaf has no name")
    if dup:
        problems.append(f"duplicate leaves: {', '.join(dup)}")
    extra, absent = sorted(set(leaves) - set(species) - {""}), sorted(set(species) - set(leaves))
    if extra:
        problems.append(f"leaves that are not species of the method: {', '.join(extra)}")
    if absent:
        problems.append(f"species of the method missing from the tree: {', '.join(absent)}")
    unresolved = [n for n in root.iter() if not n.is_leaf() and len(n.children) != 2]
    for n in unresolved[:3]:
        where = "the root" if n is root else "a node (" + ", ".join(x.name for x in n.leaves()[:4]) + ", ...)"
        problems.append(f"{where} has {len(n.children)} children; the tree must be rooted and fully resolved "
                        "(two children at every node)")
    if len(unresolved) > 3:
        problems.append(f"... and {len(unresolved) - 3} more unresolved nodes")
    return problems


def proteomes_digest(checksums):
    """MD5 of a set of proteomes: one line `<prefix> <proteome md5>` per species, sorted by prefix.

    `checksums` maps prefix -> the proteome's MD5. A run records this for its species (run.json `proteomes_md5`), and a
    frozen release can recompute it from its manifest, so a run can be tied to exactly the proteomes it used."""
    lines = [f"{p} {checksums[p]}" for p in sorted(checksums)]
    return hashlib.md5(("\n".join(lines) + "\n").encode()).hexdigest()


def species_tree_path(release, part, cfg=None):
    """Path of the species tree file a part names in releases.toml; KeyError if it names none."""
    cfg = Path(cfg or config_dir())
    rels = tomllib.loads((cfg / "releases.toml").read_text()).get("releases", {})
    if release is None:
        live = [r for r, v in rels.items() if v.get("status") != "superseded"]
        if len(live) != 1:
            raise KeyError("give a release (one of " + ", ".join(rels) + ")")
        release = live[0]
    name = rels.get(release, {}).get("parts", {}).get(part, {}).get("species_tree")
    if not name:
        raise KeyError(f"{release}/{part} has no species_tree in releases.toml")
    return cfg / name


def read_releases(cfg, organisms):
    """config/releases.toml -> release -> {status, description, deliver, inputs, parts: part -> {species, ...}}."""
    path = Path(cfg) / "releases.toml"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found (releases and their parts are defined there)")
    out = {}
    for rel, r in tomllib.loads(path.read_text()).get("releases", {}).items():
        parts = {}
        for part, m in r.get("parts", {}).items():
            members = ([p for p, s in organisms.items() if s.role == "project"] if m.get("queries") == "project"
                       else list(m.get("species", [])))
            parts[part] = {k: v for k, v in m.items() if k != "species"}
            parts[part]["species"] = members
        out[rel] = {"status": r.get("status", ""), "description": r.get("description", ""),
                    "deliver": list(r.get("deliver", [])), "inputs": dict(r.get("inputs", {})), "parts": parts}
    return out


class Catalog:
    def __init__(self, organisms, releases, panel, meta=None, cfg=None, path=None):
        self.path = Path(path) if path else None    # the catalog.json it was loaded from or written to
        self.organisms = organisms              # prefix -> Species
        self.releases = releases                # release -> {"status", "deliver", "inputs", "parts": {part: {"species": [...], ...}}}
        self.panel = panel
        self.meta = meta or {}
        self.cfg = Path(cfg) if cfg else None   # config folder it was built from (None when loaded)

    # ---- building
    @classmethod
    def build(cls, layout, cfg=None):
        cfg = Path(cfg or config_dir())
        sets = tomllib.loads((cfg / "species_sets.toml").read_text())
        manifest = read_manifest(cfg / "species.tsv")
        panel = read_panel(cfg / "compara_panel.tsv")
        organisms = {}
        for o in sets.get("organism", []):
            organisms[o["prefix"]] = Species(prefix=o["prefix"], species=o["species"], role=o["role"],
                                             lineage=o.get("lineage", ""), common=o.get("common", ""),
                                             source=o.get("source", ""), accession=o.get("accession", ""),
                                             notes=o.get("notes", ""))
        for p, row in manifest.items():                       # manifest species missing from species_sets
            if p not in organisms:
                organisms[p] = Species(prefix=p, species=row.get("species", ""), role="project",
                                       notes="in species.tsv but not in species_sets.toml")
        releases = read_releases(cfg, organisms)
        root = Path(layout.local_root)
        for p, s in organisms.items():
            row = manifest.get(p)
            s.in_manifest = row is not None
            if row and row.get("accession") and not s.accession:
                s.accession = row["accession"]
            rec = read_record(layout.download(p).local)
            if rec:
                s.accession = s.accession or rec.accession
                s.version, s.retrieved = rec.version, rec.retrieved
            prot, genes = layout.proteome(p), layout.genes_table(p)
            if prot.local.is_file():
                s.status, s.proteome = "proteome ready", prot.rel
                s.proteome_md5 = md5sum(prot.local)
                s.counts["proteins"] = _count_fasta(prot.local)
                if genes.local.is_file():
                    s.genes_table = genes.rel
                    s.counts["genes"] = _count_rows(genes.local)
            elif layout.download(p).local.is_dir():
                s.status = "downloaded"
            elif layout.local_source(p).local.is_dir():
                s.status = "local source only"
            s.used_by = [f"{rel}/{part}" for rel, r in releases.items() for part, m in r["parts"].items()
                         if p in m["species"]]
        files = ["species_sets.toml", "releases.toml", "species.tsv", "compara_panel.tsv", "species_clades.tsv",
                 "species_tree_treefam9.phy", "species_tree_master.nwk", "species_tree_placements.tsv"] + sorted(
            {m["species_tree"] for r in releases.values() for m in r["parts"].values() if m.get("species_tree")})
        meta = {"catalog_version": CATALOG_VERSION, "built": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "built_by": f"cichlid-data {__version__}",
                "config": {f: md5sum(cfg / f) for f in files if (cfg / f).is_file()}}
        return cls(organisms, releases, panel, meta, cfg)

    # ---- files
    def to_dict(self):
        return {**self.meta, "releases": self.releases, "compara_panel": self.panel,
                "organisms": {p: s.to_dict() for p, s in self.organisms.items()}}

    def write(self, layout):
        path = layout.catalog().local
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2) + "\n")
        self.path = path
        return path

    @classmethod
    def load(cls, layout):
        path = layout.catalog().local
        if not path.is_file():
            raise FileNotFoundError(f"no catalog at {path}; run `cichlid-data catalog build`")
        d = json.loads(path.read_text())
        if "releases" not in d:
            raise ValueError(f"{path} was written by an older version (it has methods, not releases); "
                             "run `cichlid-data catalog build`")
        organisms = {p: Species(**v) for p, v in d.pop("organisms").items()}
        releases, panel = d.pop("releases"), d.pop("compara_panel", [])
        return cls(organisms, releases, panel, d, path=path)

    def proteomes_md5(self, prefixes=None):
        """MD5 of the proteomes of `prefixes` (default: every organism with a proteome): the `proteomes_md5` of a run
        record, which stays the same however often the catalog is rebuilt (see proteomes_digest). KeyError if a
        prefix has no proteome checksum."""
        names = list(prefixes) if prefixes is not None else [p for p, s in self.organisms.items() if s.proteome_md5]
        missing = [p for p in names if p not in self.organisms or not self.organisms[p].proteome_md5]
        if missing:
            raise KeyError(f"no proteome checksum in the catalog for {', '.join(missing)}")
        return proteomes_digest({p: self.organisms[p].proteome_md5 for p in names})

    def md5(self, layout=None):
        """MD5 of the catalog.json file this catalog was loaded from or written to. It identifies that one file: the
        file records when it was built, so the value changes on every rebuild even if nothing else did. To tie a run to
        the proteomes it used, use proteomes_md5. A catalog that was only built has no file yet (write() it first); `layout` names the
        file to hash instead."""
        path = self.path or (layout.catalog().local if layout is not None else None)
        if path is None or not Path(path).is_file():
            raise FileNotFoundError("this catalog has no catalog.json to checksum; run `cichlid-data catalog build`, "
                                    "then Catalog.load(layout)")
        return md5sum(path)

    # ---- queries
    def __getitem__(self, prefix):
        return self.organisms[prefix]

    def release_name(self, release=None):
        """`release`, checked; or, if none is given and exactly one release is not superseded, that one."""
        if release:
            if release not in self.releases:
                raise KeyError(f"unknown release {release!r} (one of {', '.join(self.releases) or 'none'})")
            return release
        live = [r for r, v in self.releases.items() if v["status"] != "superseded"]
        if len(live) != 1:
            raise KeyError("give a release (one of " + ", ".join(self.releases) + ")")
        return live[0]

    def part(self, release, part):
        rel = self.release_name(release)
        if part not in self.releases[rel]["parts"]:
            raise KeyError(f"{rel} has no part {part!r} (one of {', '.join(self.releases[rel]['parts'])})")
        return self.releases[rel]["parts"][part]

    def species(self, release=None, part=None, role=None, status=None):
        """Organisms: all; or those of a release (the union of its parts); or of one part of a release."""
        out = list(self.organisms.values())
        if release or part:
            rel = self.release_name(release)
            if part:
                wanted = self.part(rel, part)["species"]
            else:
                wanted = list(dict.fromkeys(p for m in self.releases[rel]["parts"].values() for p in m["species"]))
            out = [self.organisms[p] for p in wanted if p in self.organisms]
        if role:
            out = [s for s in out if s.role == role]
        if status:
            out = [s for s in out if s.status == status]
        return out

    def species_tree(self, release, part):
        """Path of the species tree file the part must use (config/<species_tree>)."""
        name = self.part(release, part).get("species_tree")
        if not name:
            raise KeyError(f"{self.release_name(release)}/{part} has no species_tree in releases.toml")
        return (self.cfg or config_dir()) / name

    def missing(self, release, part=None):
        """Organisms a release (or one of its parts) needs whose proteome is not ready."""
        return [s for s in self.species(release=release, part=part) if s.status != "proteome ready"]

    def validate(self, release=None):
        """-> (errors, warnings). The organisms always; the releases (all, or just `release`)."""
        errors, warnings = [], []
        for p, s in self.organisms.items():
            if not PREFIX_RE.match(p):
                errors.append(f"{p}: a prefix must be letters and digits, starting with a letter")
            if s.role not in ROLES:
                errors.append(f"{p}: unknown role {s.role!r} (one of {', '.join(ROLES)})")
            if s.notes == "in species.tsv but not in species_sets.toml":
                warnings.append(f"{p}: in species.tsv but not in species_sets.toml; add it there with its role")
            if s.role == "project" and not s.in_manifest:
                warnings.append(f"{p}: a project species not in species.tsv, so Stage 1 won't download it")
        lower = {}
        for p in self.organisms:
            if p.lower() in lower:
                errors.append(f"{p} and {lower[p.lower()]}: prefixes differ only in case")
            lower[p.lower()] = p
        if self.cfg is not None and (self.cfg / CLADES_FILE).is_file():            # optional: required by the Homology part
            from .clades import Clades
            try:
                errors += [f"{CLADES_FILE}: {pr}" for pr in Clades.load(self.cfg).validate(self, self.cfg)]
            except ValueError as e:
                errors.append(str(e))
        if self.cfg is not None and (self.cfg / "species_tree_master.nwk").is_file():      # optional: required by the Homology part
            from . import master_tree
            try:
                errors += [f"species_tree_master.nwk: {pr}" for pr in
                           master_tree.check(master_tree.load_master_tree(self.cfg), self.cfg)]
            except (OSError, ValueError, KeyError) as e:
                errors.append(f"species_tree_master.nwk: {e}")
        names = [release] if release else list(self.releases)
        for rel in names:
            if rel not in self.releases:
                errors.append(f"unknown release {rel!r} (one of {', '.join(self.releases) or 'none'})")
                continue
            errors += self._validate_release(rel)
        return errors, warnings

    def _validate_release(self, rel):
        r, errors = self.releases[rel], []
        if not RELEASE_RE.match(rel):
            errors.append(f"{rel}: a release must be 'v' and a number (v1, v2, ...)")
        if r["status"] not in STATUSES:
            errors.append(f"{rel}: status {r['status']!r} must be one of {', '.join(STATUSES)}")
        if not r["parts"]:
            errors.append(f"{rel}: no parts")
        for part, m in r["parts"].items():
            name = f"{rel}/{part}"
            if part not in RESULT_PARTS:
                errors.append(f"{name}: unknown part (one of {', '.join(RESULT_PARTS)})")
            if part in TREE_PARTS and not m["species"]:
                errors.append(f"{name}: no species")
            for p in m["species"]:
                if p not in self.organisms:
                    errors.append(f"{name}: {p} has no organism entry")
                elif self.organisms[p].role == "candidate":
                    errors.append(f"{name}: {p} is a candidate; change its role before using it")
            if m.get("species_tree"):
                errors += [f"{name}: species tree {m['species_tree']}: {pr}"
                           for pr in check_species_tree(self.species_tree(rel, part), m["species"])]
        for p in r["deliver"]:
            if p not in self.organisms:
                errors.append(f"{rel}: delivered species {p} has no organism entry")
                continue
            if self.organisms[p].role != "project":
                errors.append(f"{rel}: delivered species {p} is not a project species")
            for part, m in r["parts"].items():
                if part in TREE_PARTS and p not in m["species"]:
                    errors.append(f"{rel}: delivered species {p} is not in {rel}/{part}")
        return errors
