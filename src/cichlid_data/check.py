"""cichlid-data check [--release R] [--part P] [--fetch] [--deep]: is everything a release's part needs here and sound?

For each part of a release, its required inputs (below) and the proteomes of its species (from the catalog):
    proteomes   present, and valid against the schemas (including the genes-table cross-check)
    resources   present (files, or non-empty folders), and matching their source record if they have
                one (sizes; --deep also checksums)
species tree  for a part that names one (releases.toml): rooted, fully resolved, and over exactly the
                part's species
A frozen (or published) release is also checked against its archive Releases/<release>/: the archive is
there and intact, and the current proteomes and configuration have not drifted from it.
With --fetch, anything missing locally is pulled from the cloud first. Exit status 1 if anything
is missing or invalid, so run scripts can call it before starting.
"""
from dataclasses import dataclass
from pathlib import Path

from . import schemas
from .catalog import FROZEN_STATUSES, check_species_tree
from .provenance import read_record, verify
from .releases import drift, load_manifest, verify_archive

# relative path -> what it must be: "file" or "dir" (non-empty)
TREEFAM = {"Inputs/Families/TreeFam9/hmm": "dir", "Inputs/Families/TreeFam9/trees": "dir",
           "Inputs/Families/TreeFam9/protein_to_species.txt": "file"}
COMPARA = {"Inputs/Families/Compara116/prepared/trees": "dir",
           "Inputs/Families/Compara116/prepared/search_db.fa": "file",
           "Inputs/Families/Compara116/prepared/members.tsv": "file",
           "Inputs/Families/Compara116/prepared/protein_to_species.txt": "file",
           "Inputs/Families/Compara116/prepared/family_audit.tsv": "file"}
NAMES = {"Inputs/Annotations/Legacy2018": "dir"}
# by part: what a part needs besides the proteomes of its species. The 2018 name tables are for the
# orthology calls (reference gene names), which the Homology part makes.
REQUIREMENTS = {"TreeFam": TREEFAM, "Compara": COMPARA, "OrthoFinder": {}, "Homology": NAMES, "Analysis": {}}
# resources whose source record is verified (the folder holding source.json)
RECORDED = {"Inputs/Families/TreeFam9", "Inputs/Families/Compara116", "Inputs/Annotations/Legacy2018"}


@dataclass
class Finding:
    method: str                 # "<release>/<part>", or "<release>" for a finding about the release itself
    item: str
    ok: bool
    detail: str = ""


def _present(path, kind):
    if kind == "file":
        return path.is_file() and path.stat().st_size > 0
    return path.is_dir() and any(path.iterdir())


def _frozen_findings(layout, rel, fetch, store, say):
    """Findings about a frozen release's archive and drift."""
    manifest_path = layout.release(rel).manifest()
    if not manifest_path.local.is_file() and fetch and store is not None:
        try:
            store.ensure(layout.release(rel).root())
        except Exception as e:
            say(f"  could not fetch the archive of {rel}: {e}")
    if load_manifest(layout, rel) is None:
        return [Finding(rel, "frozen archive", False, f"Releases/{rel}/ not on this machine "
                        "(status is frozen: run `cichlid-data release freeze`, or --fetch)")]
    probs = verify_archive(layout, rel, deep=False)
    out = [Finding(rel, "frozen archive", not probs, "intact" if not probs else "; ".join(probs[:3]))]
    drifted = drift(layout, rel)
    out.append(Finding(rel, "inputs against the frozen archive", not drifted,
                       "the same as when frozen" if not drifted else "; ".join(drifted[:3])))
    return out


def check(layout, catalog, release=None, parts=None, fetch=False, deep=False, store=None, say=None):
    """-> list of Findings. `catalog` is a Catalog (its releases list each part's species). Without `release`,
    every release that is not superseded is checked; without `parts`, every part of it."""
    say = say or (lambda m: None)
    root = Path(layout.local_root)
    releases = ([catalog.release_name(release)] if release else
                [r for r, v in catalog.releases.items() if v["status"] != "superseded"])
    findings, verified = [], {}
    for rel in releases:
        if catalog.releases[rel]["status"] in FROZEN_STATUSES:
            findings += _frozen_findings(layout, rel, fetch, store, say)
        for part in (parts or list(catalog.releases[rel]["parts"])):
            m = f"{rel}/{part}"
            if part not in catalog.releases[rel]["parts"]:
                findings.append(Finding(m, "part", False, "not in releases.toml"))
                continue
            findings += _check_part(layout, catalog, root, rel, part, m, fetch, deep, store, say, verified)
    return findings


def _check_part(layout, catalog, root, rel, part, m, fetch, deep, store, say, verified):
    findings = []
    for res, kind in REQUIREMENTS.get(part, {}).items():
        path = root / res
        if not _present(path, kind) and fetch and store is not None:
            try:
                store.ensure(layout.path(res))
            except Exception as e:                    # report below as missing
                say(f"  could not fetch {res}: {e}")
        if not _present(path, kind):
            findings.append(Finding(m, res, False, "missing" if not path.exists() else "empty"))
            continue
        rec_dir = next((r for r in RECORDED if res == r or res.startswith(r + "/")), None)
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
                findings.append(Finding(m, res, False, f"differs from its source record: {'; '.join(probs[:3])}"))
                continue
            detail = "no source record" if probs is None else "matches its source record"
        findings.append(Finding(m, res, True, detail))
    if catalog.part(rel, part).get("species_tree"):
        name = catalog.part(rel, part)["species_tree"]
        probs = check_species_tree(catalog.species_tree(rel, part),
                                   [s.prefix for s in catalog.species(release=rel, part=part)])
        findings.append(Finding(m, f"species tree {name}", not probs, "valid" if not probs else "; ".join(probs[:3])))
    for s in catalog.species(release=rel, part=part):
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
