"""cichlid-data provenance write | check | show: source records for every input."""
import time
from pathlib import Path

from . import __version__
from .config import config_dir
from .provenance import STAGE1_RECORD, SourceRecord, download_date, inventory, load_sources, read_record, verify

AREAS = ("Inputs/Families", "Inputs/Annotations", "Inputs/Expression", "Inputs/Populations",
         "Inputs/Proteins/local_sources")


def write_records(layout, overwrite=False, md5_limit=1 << 30, say=print, sources_file=None):
    root = Path(layout.local_root)
    entries = load_sources(sources_file or config_dir() / "sources.toml")
    written = 0
    for e in entries:
        folder = root / e["path"]
        if not folder.exists():
            say(f"  {e['path']}: not on this machine (pull it first), skipped")
            continue
        if read_record(folder) and not overwrite:
            say(f"  {e['path']}: already recorded")
            continue
        rec = SourceRecord.from_dict({k: v for k, v in e.items() if k != "path"})
        rec.files = inventory(folder, md5_limit)
        rec.recorded = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        rec.recorded_by = f"cichlid-data {__version__}"
        if not rec.retrieved and (folder / STAGE1_RECORD).is_file():
            rec.retrieved = download_date(folder / STAGE1_RECORD)          # the download script's own date
            rec.notes = (rec.notes + "; " if rec.notes else "") + "download details in download.json"
        if not rec.retrieved:
            rec.retrieved = "before " + rec.recorded[:10] + " (recorded later)"
        rec.write(folder)
        written += 1
        say(f"  {e['path']}: recorded {len(rec.files)} files")
    downloads = root / "Inputs/Proteins/downloads"
    if downloads.is_dir():
        for d in sorted(p for p in downloads.iterdir() if p.is_dir()):
            say(f"  {d.relative_to(root)}: {'recorded by Stage 1 (download.json)' if read_record(d) else 'NO RECORD'}")
    for area, item in uncovered(root, entries):
        say(f"  {area}/{item}: no source description; add it to config/sources.toml")
    return written


def uncovered(root, entries):
    """Entries in the input areas not covered by any source description."""
    paths = [e["path"] for e in entries]
    out = []
    for area in AREAS:
        d = root / area
        if not d.is_dir():
            continue
        for item in sorted(p.name for p in d.iterdir() if not p.name.startswith(".") and p.name != "source.json"):
            rel = f"{area}/{item}"
            if not any(rel == p or rel.startswith(p + "/") or p.startswith(rel + "/") for p in paths):
                out.append((area, item))
    return out


def check_records(layout, deep=False, say=print, sources_file=None):
    root = Path(layout.local_root)
    folders = [root / e["path"] for e in load_sources(sources_file or config_dir() / "sources.toml")]
    downloads = root / "Inputs/Proteins/downloads"
    if downloads.is_dir():
        folders += sorted(p for p in downloads.iterdir() if p.is_dir())
    bad = 0
    for f in folders:
        if not f.exists():
            continue
        rec = read_record(f)
        if rec is None:
            say(f"  {f.relative_to(root)}: NO RECORD")
            bad += 1
            continue
        problems = verify(f, rec, deep=deep)
        real = [p for p in problems if not p.startswith("new file")]
        bad += bool(real)
        say(f"  {f.relative_to(root)}: " + ("ok" if not problems else "; ".join(problems[:6]) +
                                           (f" (+{len(problems) - 6} more)" if len(problems) > 6 else "")))
    return bad
