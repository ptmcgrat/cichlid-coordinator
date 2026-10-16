"""cichlid-data check [--method M] [--fetch] [--deep]: is everything a method needs here and sound?

For each method, its required inputs (below) and the proteomes of its species (from the catalog):
    proteomes   present, and valid against the schemas (including the genes-table cross-check)
    resources   present (files, or non-empty folders), and matching their source record if they have
                one (sizes; --deep also checksums)
With --fetch, anything missing locally is pulled from the cloud first. Exit status 1 if anything
is missing or invalid, so run scripts can call it before starting.
"""
from dataclasses import dataclass
from pathlib import Path

from . import schemas
from .provenance import read_record, verify

# relative path -> what it must be: "file" or "dir" (non-empty)
TREEFAM = {"Inputs/Families/TreeFam9/hmm": "dir", "Inputs/Families/TreeFam9/trees": "dir",
           "Inputs/Families/TreeFam9/protein_to_species.txt": "file"}
COMPARA = {"Inputs/Families/Compara116/prepared/trees": "dir",
           "Inputs/Families/Compara116/prepared/search_db.fa": "file",
           "Inputs/Families/Compara116/prepared/members.tsv": "file",
           "Inputs/Families/Compara116/prepared/protein_to_species.txt": "file",
           "Inputs/Families/Compara116/prepared/family_audit.tsv": "file"}
NAMES = {"Inputs/Annotations/Legacy2018": "dir"}
REQUIREMENTS = {"TreeFam2026": {**TREEFAM, **NAMES}, "Compara2026": {**COMPARA, **NAMES}, "OrthoFinder2026": {}}
# resources whose source record is verified (the folder holding source.json)
RECORDED = {"Inputs/Families/TreeFam9", "Inputs/Families/Compara116", "Inputs/Annotations/Legacy2018"}


@dataclass
class Finding:
    method: str
    item: str
    ok: bool
    detail: str = ""


def _present(path, kind):
    if kind == "file":
        return path.is_file() and path.stat().st_size > 0
    return path.is_dir() and any(path.iterdir())


def check(layout, catalog, methods=None, fetch=False, deep=False, store=None, say=None):
    """-> list of Findings. `catalog` is a Catalog (its methods list each one's species)."""
    say = say or (lambda m: None)
    root = Path(layout.local_root)
    methods = methods or [m for m in REQUIREMENTS if m in catalog.methods]
    findings, verified = [], {}
    for m in methods:
        if m not in catalog.methods:
            findings.append(Finding(m, "method", False, "not in species_sets.toml"))
            continue
        for rel, kind in REQUIREMENTS.get(m, {}).items():
            path = root / rel
            if not _present(path, kind) and fetch and store is not None:
                try:
                    store.ensure(layout.path(rel))
                except Exception as e:                    # report below as missing
                    say(f"  could not fetch {rel}: {e}")
            if not _present(path, kind):
                findings.append(Finding(m, rel, False, "missing" if not path.exists() else "empty"))
                continue
            rec_dir = next((r for r in RECORDED if rel == r or rel.startswith(r + "/")), None)
            detail = ""
            if rec_dir:
                if rec_dir not in verified:
                    rec = read_record(root / rec_dir)
                    verified[rec_dir] = ([] if rec is None else
                                         [p for p in verify(root / rec_dir, rec, deep=deep) if not p.startswith("new file")])
                    if rec is None:
                        verified[rec_dir] = None
                probs = verified[rec_dir]
                if probs:
                    findings.append(Finding(m, rel, False, f"differs from its source record: {'; '.join(probs[:3])}"))
                    continue
                detail = "no source record" if probs is None else "matches its source record"
            findings.append(Finding(m, rel, True, detail))
        for s in catalog.species(method=m):
            faa = layout.proteome(s.prefix)
            if not faa.local.is_file() and fetch and store is not None:
                try:
                    store.ensure(faa)
                except Exception as e:
                    say(f"  could not fetch {faa.rel}: {e}")
            if not faa.local.is_file():
                findings.append(Finding(m, f"proteome {s.prefix}", False, "missing"))
                continue
            probs = schemas.validate(faa.local)
            findings.append(Finding(m, f"proteome {s.prefix}", not probs,
                                    "valid" if not probs else "; ".join(probs[:3])))
    return findings
