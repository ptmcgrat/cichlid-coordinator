"""cichlid-data download populations [options]   (moved from the pipeline; default folder Inputs/Populations/Blumer2025)


Download the Blumer et al. 2025 data deposit (Apollo, University of Cambridge) and back it up
to the lab Dropbox. Run on mzebra (the files are large), ideally inside screen.

    python download_blumer.py --list                 # what would be downloaded, with sizes
    python download_blumer.py                        # download VCFs + ancestral sequence, then back up
    python download_blumer.py --only chr7            # just matching files (testing)
    python download_blumer.py --all                  # also the two assemblies the paper used
    python download_blumer.py --no-backup            # download only

Files go to <local master>/Data/Blumer2025/ (see sync.py env), and are copied to
ptm_dropbox:/.../TreeFamData/Data/Blumer2025/. File names, sizes, MD5 checksums and download links
come from the repository's DSpace REST API, so nothing is hard-coded. Downloads resume after an
interruption (rerun the same command) and every file is checked against the repository's MD5.

Data terms: CC BY-NC 4.0, under an Access and Benefit Sharing agreement with the Government of
Malawi; research use only. A notice with the citation is written next to the data.
"""
import argparse
import fnmatch
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import requests

from ..backup import backup, cloud_path, settings  # noqa: E402
from ..layout import get_layout  # noqa: E402

REPO = "https://www.repository.cam.ac.uk"              # web pages
API_BASES = ["https://api.repository.cam.ac.uk/server/api",   # Apollo's DSpace REST API
             "https://www.repository.cam.ac.uk/server/api"]
ITEM = "349319a6-6c36-4e49-8838-595c963a5583"
DEFAULT = ["chr*.vcf.gz", "*.ancestral.fa.gz"]          # what the gene browser needs
NOTICE = """Blumer et al. 2025, "Introgression dynamics of sex-linked chromosomal inversions shape the
Malawi cichlid radiation", Science. https://www.science.org/doi/10.1126/science.adr9961
Data: https://www.repository.cam.ac.uk/items/349319a6-6c36-4e49-8838-595c963a5583

The Malawi cichlid samples were collected under permits, and the data are published under an Access
and Benefit Sharing agreement with the Government of Malawi. Licence: CC BY-NC 4.0, research use only.
Commercial use requires a licensing and benefit sharing arrangement with the Government of Malawi.
Genotypes are against fAstCal1.2 (NCBI RefSeq GCF_900246225.1). Contact: Richard Durbin, rd109@cam.ac.uk.
"""


# --------------------------------------------------------------------------- repository listing
def get_json(session, url, params=None, tries=5):
    for i in range(tries):
        try:
            r = session.get(url, params=params, timeout=60, headers={"Accept": "application/json"})
            if r.status_code == 200:
                return r.json()
            if r.status_code < 500:
                r.raise_for_status()
        except (requests.ConnectionError, requests.Timeout):
            pass
        time.sleep(2 ** i)
    raise IOError(f"repository API did not answer: {url}")


def list_files(session, apis=None, item=ITEM):
    """-> [{name, size, md5, url}] for the item's ORIGINAL bundle (the deposited files).
    Tries each API base in turn (Apollo serves its API from a separate host)."""
    bundles, errors = None, []
    for api in (apis or API_BASES):
        try:
            bundles = get_json(session, f"{api}/core/items/{item}/bundles")
            break
        except (IOError, requests.HTTPError) as e:
            errors.append(f"{api}: {e}")
    if bundles is None:
        raise IOError("could not reach the repository API:\n  " + "\n  ".join(errors))
    files = []
    for b in bundles.get("_embedded", {}).get("bundles", []):
        if b.get("name") != "ORIGINAL":
            continue
        url, page = b["_links"]["bitstreams"]["href"], 0
        while url:
            data = get_json(session, url, params={"size": 100, "page": page})
            for bs in data.get("_embedded", {}).get("bitstreams", []):
                cs = bs.get("checkSum") or {}
                files.append({"name": bs["name"], "size": int(bs.get("sizeBytes") or 0),
                              "md5": (cs.get("value") or "").lower() if cs.get("checkSumAlgorithm", "MD5").upper() == "MD5" else "",
                              "url": bs["_links"]["content"]["href"]})
            pg = data.get("page", {})
            page += 1
            url = url if pg and page < pg.get("totalPages", 1) else None
    if not files:
        raise IOError("the repository listed no files; check the item ID or try again later")
    return files


def select(files, patterns):
    return [f for f in files if any(fnmatch.fnmatch(f["name"], p) or p in f["name"] for p in patterns)]


# --------------------------------------------------------------------------- download
def md5sum(path, block=8 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def verified(path, f):
    if not path.is_file() or (f["size"] and path.stat().st_size != f["size"]):
        return False
    return not f["md5"] or md5sum(path) == f["md5"]


def download(session, f, outdir, tries=8):
    """Resumable download to <name>.part, verified, then renamed. -> (name, status)."""
    dest, part = outdir / f["name"], outdir / (f["name"] + ".part")
    if verified(dest, f):
        return f["name"], "already downloaded"
    for attempt in range(tries):
        have = part.stat().st_size if part.is_file() else 0
        if f["size"] and have >= f["size"]:
            break
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with session.get(f["url"], headers=headers, stream=True, timeout=120) as r:
                if r.status_code == 200 and have:             # server ignored the range: start over
                    have = 0
                    part.unlink(missing_ok=True)
                elif r.status_code not in (200, 206):
                    raise IOError(f"HTTP {r.status_code}")
                with open(part, "ab" if have else "wb") as fh:
                    for chunk in r.iter_content(8 << 20):
                        fh.write(chunk)
            if not f["size"] or part.stat().st_size >= f["size"]:
                break
        except (requests.ConnectionError, requests.Timeout, IOError) as e:
            print(f"  {f['name']}: interrupted ({e}); resuming", file=sys.stderr)
            time.sleep(min(60, 2 ** attempt))
    if f["md5"] and md5sum(part) != f["md5"]:
        part.unlink(missing_ok=True)
        return f["name"], "FAILED: checksum mismatch (removed; rerun to download again)"
    if f["size"] and part.stat().st_size != f["size"]:
        return f["name"], f"FAILED: incomplete ({part.stat().st_size} of {f['size']} bytes; rerun to resume)"
    part.replace(dest)
    return f["name"], "downloaded and verified"


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024


# --------------------------------------------------------------------------- main
def main(argv=None, session=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="show files and sizes, download nothing")
    ap.add_argument("--only", help="comma-separated name patterns, e.g. chr7 or 'chr1.vcf.gz,chr2.vcf.gz'")
    ap.add_argument("--all", action="store_true", help="include the fAstCal1.6 and Diplotaxodon assemblies")
    ap.add_argument("--jobs", type=int, default=4, help="parallel downloads")
    ap.add_argument("--no-backup", action="store_true", help="skip the Dropbox copy")
    ap.add_argument("--outdir", type=Path, help="default: <local master>/Data/Blumer2025")
    ap.add_argument("--api", help="DSpace REST API base, e.g. https://api.repository.cam.ac.uk/server/api")
    args = ap.parse_args(argv)

    outdir = args.outdir or get_layout().population("Blumer2025").local
    session = session or requests.Session()
    files = list_files(session, [args.api] if args.api else None)
    patterns = [p.strip() for p in args.only.split(",")] if args.only else (["*"] if args.all else DEFAULT)
    chosen = select(files, patterns)
    total = sum(f["size"] for f in chosen)
    print(f"{len(chosen)} of {len(files)} files, {human(total)}, to {outdir}", file=sys.stderr)
    for f in chosen:
        print(f"  {f['name']:<32} {human(f['size']):>10}", file=sys.stderr)
    if args.list or not chosen:
        return

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "DATA_NOTICE.txt").write_text(NOTICE)
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(download, session, f, outdir) for f in chosen]
        for fut in as_completed(futs):
            name, status = fut.result()
            print(f"  {name}: {status}", file=sys.stderr)
            if status.startswith("FAILED"):
                failed.append(name)
    with open(outdir / "files.tsv", "w") as fh:
        fh.write("name\tbytes\tmd5\tsource_url\tverified\n")
        for f in chosen:
            fh.write(f"{f['name']}\t{f['size']}\t{f['md5']}\t{f['url']}\t{f['name'] not in failed}\n")
    (outdir / "download.json").write_text(json.dumps(
        {"repository_item": f"{REPO}/items/{ITEM}", "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "files": len(chosen), "failed": failed}, indent=2))
    if failed:
        sys.exit(f"{len(failed)} files failed: {', '.join(failed)}. Rerun the same command to resume; "
                 f"nothing was backed up yet.")

    if args.no_backup:
        return
    cloud = cloud_path(outdir)
    print(f"\nBacking up to {cloud} (rclone shows progress; already-copied files are skipped)", file=sys.stderr)
    r = subprocess.run([settings().rclone, "copy", str(outdir), cloud, "--exclude", "*.part",
                        "--transfers", "4", "--checkers", "8", "--stats", "30s", "--stats-one-line", "-v"])
    if r.returncode != 0:
        sys.exit(f"rclone copy failed (exit {r.returncode}); rerun with the same command to resume the backup")
    backup(outdir)                              # records it in the cloud index (nothing left to copy)
    print("Backup complete.", file=sys.stderr)


if __name__ == "__main__":
    main()


