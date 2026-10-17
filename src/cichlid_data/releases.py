"""Releases: freezing one, verifying its archive, and checking that the inputs have not drifted from it.

A release (v1, v2, ...) is defined in config/releases.toml. When its owner sets its status to frozen (or
published), `cichlid-data release freeze v1` writes the archive Releases/v1/:

    release.json        what the release is, with the size and MD5 of every file below
    catalog.json        the catalog as it was
    proteomes/<prefix>.faa, genes/<prefix>.genes.tsv     the proteomes the release was built from
    config/             its species trees and reference panel

The archive never changes (its unit is immutable). A better annotation of a species does not change a
frozen release: it makes a new one. Runs of a frozen release can read the proteomes from the archive
(`layout.release("v1").proteome("Mzeb")`).
"""
import json
import os
import shutil
import time
from pathlib import Path

from . import __version__
from .catalog import FROZEN_STATUSES
from .config import config_dir
from .provenance import md5sum

MANIFEST = "release.json"
CONFIG_KEYS = ("species_tree", "reference_panel")      # a part's keys that name a file in config/


class ReleaseError(Exception):
    pass


def load_manifest(layout, release):
    """The release.json of a frozen release on this machine, or None."""
    path = layout.release(release).manifest().local
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def _config_files(catalog, release):
    cfg = catalog.cfg or config_dir()
    names = sorted({m[k] for m in catalog.releases[release]["parts"].values() for k in CONFIG_KEYS if m.get(k)})
    return cfg, names


def _entry(path):
    return {"size": Path(path).stat().st_size, "md5": md5sum(path)}


def freeze(layout, catalog, release, say=print):
    """Write Releases/<release>/ from a freshly built `catalog`. -> the archive's path.

    Refuses unless the release exists and is valid, its status is frozen or published, every species it
    uses has a proteome, and it has not been frozen already."""
    rel = catalog.release_name(release)
    errors, _ = catalog.validate(release=rel)
    if errors:
        raise ReleaseError(f"{rel} is not valid:\n  " + "\n  ".join(errors))
    status = catalog.releases[rel]["status"]
    if status not in FROZEN_STATUSES:
        raise ReleaseError(f"{rel} has status {status!r}: set status = \"frozen\" in config/releases.toml "
                           "first (freezing is the owner's decision), then freeze")
    target = layout.release(rel).root().local
    if target.exists():
        raise ReleaseError(f"{target} exists: {rel} is already frozen, and a frozen release never changes. "
                           "Make a new release instead.")
    species = catalog.species(release=rel)
    missing = [s.prefix for s in species if s.status != "proteome ready"]
    if missing:
        raise ReleaseError(f"{rel}: no proteome yet for {', '.join(missing)}")
    cfg, config_names = _config_files(catalog, rel)
    for name in config_names:
        if not (cfg / name).is_file():
            raise ReleaseError(f"{rel}: config/{name} not found")

    partial = target.parent / f".{rel}.partial"
    if partial.exists():
        shutil.rmtree(partial)
    (partial / "proteomes").mkdir(parents=True)
    (partial / "genes").mkdir()
    (partial / "config").mkdir()
    try:
        (partial / "catalog.json").write_text(json.dumps(catalog.to_dict(), indent=2) + "\n")
        for s in species:
            shutil.copy2(layout.proteome(s.prefix).local, partial / "proteomes" / f"{s.prefix}.faa")
            genes = layout.genes_table(s.prefix).local
            if genes.is_file():
                shutil.copy2(genes, partial / "genes" / f"{s.prefix}.genes.tsv")
        for name in config_names:
            shutil.copy2(cfg / name, partial / "config" / name)
        files = {str(f.relative_to(partial)): _entry(f) for f in sorted(partial.rglob("*")) if f.is_file()}
        r = catalog.releases[rel]
        manifest = {
            "release": rel, "status": status, "frozen": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "frozen_by": f"cichlid-data {__version__}", "coordinator_version": __version__,
            "description": r["description"], "deliver": r["deliver"], "inputs": r["inputs"], "parts": r["parts"],
            "catalog_md5": files["catalog.json"]["md5"],
            "species": {s.prefix: {"species": s.species, "role": s.role, "accession": s.accession,
                                   "version": s.version, "proteome_md5": s.proteome_md5, "counts": s.counts}
                        for s in species},
            "files": files,
        }
        (partial / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
        os.replace(partial, target)                      # the archive appears whole or not at all
    except BaseException:
        shutil.rmtree(partial, ignore_errors=True)
        raise
    say(f"Froze {rel}: {len(species)} species, {len(files)} files in {target}")
    say(f"Now push it:  cichlid-data push releases/{rel}")
    return target


def verify_archive(layout, release, deep=True):
    """-> problems with the archive on this machine (empty if it is as it was frozen)."""
    manifest = load_manifest(layout, release)
    root = layout.release(release).root().local
    if manifest is None:
        return [f"{root / MANIFEST} not found: {release} is not frozen on this machine"]
    problems = []
    if manifest.get("release") != release:
        problems.append(f"{MANIFEST} is for {manifest.get('release')!r}, not {release!r}")
    for rel, want in manifest["files"].items():
        f = root / rel
        if not f.is_file():
            problems.append(f"{rel}: missing")
        elif f.stat().st_size != want["size"]:
            problems.append(f"{rel}: size {f.stat().st_size}, frozen as {want['size']}")
        elif deep and md5sum(f) != want["md5"]:
            problems.append(f"{rel}: content differs from when it was frozen")
    on_disk = {str(f.relative_to(root)) for f in root.rglob("*") if f.is_file()} - {MANIFEST}
    problems += [f"{x}: not in the manifest (a frozen release does not gain files)" for x in sorted(on_disk - set(manifest["files"]))]
    return problems


def drift(layout, release, cfg=None):
    """-> how the current proteomes and configuration differ from a frozen release (empty if they do not).
    A difference is not an error in the data: it means the release is behind, and a new one is due."""
    manifest = load_manifest(layout, release)
    if manifest is None:
        return []
    cfg = Path(cfg or config_dir())
    out = []
    for prefix, rec in manifest["species"].items():
        faa = layout.proteome(prefix).local
        if not faa.is_file():
            out.append(f"proteome {prefix}: not on this machine (the frozen copy is Releases/{release}/proteomes/{prefix}.faa)")
        elif md5sum(faa) != rec["proteome_md5"]:
            out.append(f"proteome {prefix}: differs from the frozen {release} (a new annotation or assembly: make a new release)")
    for rel, want in manifest["files"].items():
        if rel.startswith("config/") and (cfg / rel[len("config/"):]).is_file() \
                and md5sum(cfg / rel[len("config/"):]) != want["md5"]:
            out.append(f"{rel}: differs from the frozen {release}")
    return out
