"""The species catalog: every organism, its role, source, files, counts, and the methods using it.

    C = Catalog.build(get_layout())      # from config/ and what is on disk; C.write() saves catalog.json
    C = Catalog.load(get_layout())
    C.species(method="OrthoFinder2026")  # -> [Species]
    C.missing("OrthoFinder2026")         # organisms a method needs whose proteome is not ready
    C["Mzeb"].proteome                   # relative path of its proteome

Built from config/species_sets.toml (organisms, roles, method membership), config/species.tsv
(what Stage 1 downloads), config/compara_panel.tsv (Compara's reference panel), and the data:
each organism's source record (Stage 1's download.json, or source.json), and whether its download
and proteome exist.
"""
import csv
import json
import re
import time
import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path

from . import __version__
from .config import config_dir
from .provenance import md5sum, read_record

ROLES = ("project", "relative", "reference", "outgroup", "candidate")
PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
CATALOG_VERSION = 1


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


class Catalog:
    def __init__(self, organisms, methods, panel, meta=None):
        self.organisms = organisms              # prefix -> Species
        self.methods = methods                  # name -> {"description", "species": [prefixes], ...}
        self.panel = panel
        self.meta = meta or {}

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
        methods = {}
        for name, m in sets.get("methods", {}).items():
            members = ([p for p, s in organisms.items() if s.role == "project"] if m.get("queries") == "project"
                       else list(m.get("species", [])))
            methods[name] = {k: v for k, v in m.items() if k != "species"}
            methods[name]["species"] = members
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
            s.used_by = [name for name, m in methods.items() if p in m["species"]]
        meta = {"catalog_version": CATALOG_VERSION, "built": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "built_by": f"cichlid-data {__version__}",
                "config": {f: md5sum(cfg / f) for f in ("species_sets.toml", "species.tsv", "compara_panel.tsv")
                           if (cfg / f).is_file()}}
        return cls(organisms, methods, panel, meta)

    # ---- files
    def to_dict(self):
        return {**self.meta, "methods": self.methods, "compara_panel": self.panel,
                "organisms": {p: s.to_dict() for p, s in self.organisms.items()}}

    def write(self, layout):
        path = layout.catalog().local
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2) + "\n")
        return path

    @classmethod
    def load(cls, layout):
        path = layout.catalog().local
        if not path.is_file():
            raise FileNotFoundError(f"no catalog at {path}; run `cichlid-data catalog build`")
        d = json.loads(path.read_text())
        organisms = {p: Species(**v) for p, v in d.pop("organisms").items()}
        methods, panel = d.pop("methods"), d.pop("compara_panel", [])
        return cls(organisms, methods, panel, d)

    # ---- queries
    def __getitem__(self, prefix):
        return self.organisms[prefix]

    def species(self, method=None, role=None, status=None):
        out = list(self.organisms.values())
        if method:
            if method not in self.methods:
                raise KeyError(f"unknown method {method!r} (one of {', '.join(self.methods)})")
            wanted = self.methods[method]["species"]
            out = [self.organisms[p] for p in wanted if p in self.organisms]
        if role:
            out = [s for s in out if s.role == role]
        if status:
            out = [s for s in out if s.status == status]
        return out

    def missing(self, method):
        """Organisms a method needs whose proteome is not ready."""
        return [s for s in self.species(method=method) if s.status != "proteome ready"]

    def validate(self):
        """-> (errors, warnings)"""
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
        for name, m in self.methods.items():
            for p in m["species"]:
                if p not in self.organisms:
                    errors.append(f"{name}: {p} has no organism entry")
                elif self.organisms[p].role == "candidate":
                    errors.append(f"{name}: {p} is a candidate; change its role before using it")
        return errors, warnings
