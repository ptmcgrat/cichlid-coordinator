"""Who and what produced a transfer or a download; checksums."""
import hashlib
import os
import platform
import subprocess
import time

from . import __version__


def machine_name():
    """Readable machine name (on a Mac, the computer name: platform.node() is the VPN address on VPN)."""
    if platform.system() == "Darwin":
        try:
            r = subprocess.run(["scutil", "--get", "ComputerName"], capture_output=True, text=True, timeout=5)
            if r.returncode == 0 and r.stdout.strip():
                return r.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    return platform.node().split(".")[0] or platform.node()


def git_info(path):
    def git(*a):
        try:
            r = subprocess.run(["git", "-C", str(path)] + list(a), capture_output=True, text=True)
            return r.stdout.strip() if r.returncode == 0 else ""
        except OSError:
            return ""
    return {"branch": git("rev-parse", "--abbrev-ref", "HEAD") or "unknown",
            "commit": git("rev-parse", "--short", "HEAD") or "unknown",
            "uncommitted_changes": bool(git("status", "--porcelain", "--untracked-files=no"))}


def stamp(context="", repo=None):
    """Provenance of an upload: machine, user, time, coordinator version, and the caller's repo state."""
    out = {"context": context, "machine": machine_name(), "user": os.getenv("USER", ""),
           "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "coordinator": __version__}
    if repo:
        out.update(git_info(repo))
    return out


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


# =============================================================================== source records
import json  # noqa: E402
from dataclasses import asdict, dataclass, field  # noqa: E402
from pathlib import Path  # noqa: E402

RECORD = "source.json"
STAGE1_RECORD = "download.json"
# bookkeeping files, never part of an inventory (assembly_report.json is NCBI metadata Stage 1
# fetches separately and writes beside the download)
SKIP = {RECORD, STAGE1_RECORD, "assembly_report.json", ".DS_Store"}


@dataclass
class SourceRecord:
    """What a downloaded resource is, where it came from, and what it contained."""
    name: str
    type: str = ""                  # ncbi_refseq | ncbi_genbank | ensembl | treefam | authors | url | project | unknown
    version: str = ""
    accession: str = ""
    url: str = ""
    retrieved: str = ""
    licence: str = ""
    citation: str = ""
    notes: str = ""
    files: dict = field(default_factory=dict)      # relative path -> {"bytes": n, "md5": hex or None}
    recorded: str = ""
    recorded_by: str = ""

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})

    def write(self, folder):
        Path(folder).mkdir(parents=True, exist_ok=True)
        (Path(folder) / RECORD).write_text(json.dumps(self.to_dict(), indent=2) + "\n")

    @classmethod
    def from_download_json(cls, path):
        """Stage 1's download.json (a species' genome and annotation) as a SourceRecord."""
        d = json.loads(Path(path).read_text())
        row, meta = d.get("manifest_row", {}), d.get("ncbi_metadata", {}) or {}
        acc = row.get("accession", "")
        typ = ("ncbi_refseq" if acc.startswith("GCF_") else "ncbi_genbank" if acc.startswith("GCA_")
               else "authors" if row.get("source") == "local" else row.get("source", ""))
        files = {}
        raw = d.get("files")
        for info in (raw.values() if isinstance(raw, dict) else []):
            if isinstance(info, dict) and info.get("file"):
                files[info["file"]] = {"bytes": info.get("bytes"), "md5": info.get("md5")}
        return cls(name=row.get("species", ""), type=typ, accession=acc,
                   version=meta.get("assembly_name", "") or "", url=(f"https://www.ncbi.nlm.nih.gov/datasets/genome/{acc}/"
                                                                     if acc else ""),
                   retrieved=d.get("downloaded", ""), notes=row.get("notes", "") or "", files=files,
                   recorded=d.get("downloaded", ""), recorded_by="pipeline stage 1")


def is_stage1(path):
    """Stage 1's download.json (a species download), as opposed to other scripts' download.json."""
    try:
        d = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return False
    return isinstance(d, dict) and isinstance(d.get("manifest_row"), dict)


def download_date(path):
    """The download time recorded in any script's download.json, if it has one."""
    try:
        d = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return ""
    if not isinstance(d, dict):
        return ""
    for k in ("downloaded", "retrieved", "date", "time"):
        if isinstance(d.get(k), str) and d[k]:
            return d[k]
    return ""


def read_record(folder):
    """The source record of a folder: its source.json, or Stage 1's download.json for a species
    download. None if neither (other scripts' download.json files are not source records)."""
    folder = Path(folder)
    if (folder / RECORD).is_file():
        return SourceRecord.from_dict(json.loads((folder / RECORD).read_text()))
    if (folder / STAGE1_RECORD).is_file() and is_stage1(folder / STAGE1_RECORD):
        return SourceRecord.from_download_json(folder / STAGE1_RECORD)
    return None


def inventory(folder, md5_limit=1 << 30):
    """{relative path: {bytes, md5}} for every file under folder (md5 only up to md5_limit bytes;
    larger files are recorded by size alone, as checksumming e.g. 200 GB of VCFs takes hours)."""
    folder, out = Path(folder), {}
    if folder.is_file():
        st = folder.stat()
        return {folder.name: {"bytes": st.st_size, "md5": md5sum(folder) if st.st_size <= md5_limit else None}}
    for f in sorted(p for p in folder.rglob("*") if p.is_file() and p.name not in SKIP):
        size = f.stat().st_size
        out[str(f.relative_to(folder))] = {"bytes": size, "md5": md5sum(f) if size <= md5_limit else None}
    return out


def verify(folder, record, deep=False):
    """-> list of problems: missing files, size changes, and (deep) checksum changes. Files present
    but not in the record are reported as notes ('new file')."""
    folder, problems = Path(folder), []
    base = folder.parent if folder.is_file() else folder
    for rel, info in record.files.items():
        f = base / rel
        if not f.is_file():
            problems.append(f"missing: {rel}")
            continue
        if info.get("bytes") is not None and f.stat().st_size != info["bytes"]:
            problems.append(f"size changed: {rel}")
        elif deep and info.get("md5") and md5sum(f) != info["md5"]:
            problems.append(f"content changed: {rel}")
    if folder.is_dir():
        for f in sorted(p for p in folder.rglob("*") if p.is_file() and p.name not in SKIP):
            if str(f.relative_to(folder)) not in record.files:
                problems.append(f"new file (not in record): {f.relative_to(folder)}")
    return problems


def load_sources(path):
    """config/sources.toml -> list of dicts (each with path, name, type, ...)."""
    import tomllib
    return tomllib.loads(Path(path).read_text()).get("source", [])


