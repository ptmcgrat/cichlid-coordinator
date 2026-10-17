"""cichlid-data download annotations [options]   (moved from the pipeline; default folder Inputs/Annotations)


Download Gene Ontology and PANTHER annotations for human and zebrafish (for go_enrichment.py).

    python download_annotations.py              # into <local master>/Data/Annotations/{GO,PANTHER}
    python download_annotations.py --inspect    # print the first lines of each downloaded file

GO (GO Consortium current release): go-basic.obo, goa_human.gaf.gz, zfin.gaf.gz.
PANTHER: the human and zebrafish sequence classification files from the current release folder
(the release number is part of the file names, so the folder listing is read to find them).
Files already present are kept unless --force. Backed up to Dropbox like other inputs.
"""
import argparse
import gzip
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

GO_FILES = {"go-basic.obo": ["http://purl.obolibrary.org/obo/go/go-basic.obo",
                             "https://current.geneontology.org/ontology/go-basic.obo"],
            "goa_human.gaf.gz": ["https://current.geneontology.org/annotations/goa_human.gaf.gz"],
            "zfin.gaf.gz": ["https://current.geneontology.org/annotations/zfin.gaf.gz"]}
PANTHER_DIRS = ["https://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/",
                "http://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/"]
PANTHER_SPECIES = {"human": re.compile(r'href="(PTHR[\d.]+_human)"'), "zebrafish": re.compile(r'href="(PTHR[\d.]+_zebrafish)"')}


def fetch(session, urls, dest, tries=3):
    last = None
    for url in urls:
        for _ in range(tries):
            try:
                with session.get(url, stream=True, timeout=300) as r:
                    if r.status_code != 200:
                        last = f"{url}: HTTP {r.status_code}"
                        break
                    part = dest.with_suffix(dest.suffix + ".part")
                    with open(part, "wb") as fh:
                        for chunk in r.iter_content(1 << 20):
                            fh.write(chunk)
                    part.replace(dest)
                    return url
            except Exception as e:                       # network hiccup: retry
                last = f"{url}: {e}"
    raise RuntimeError(f"could not download {dest.name}: {last}")


def panther_files(session, dirs=PANTHER_DIRS):
    """Find the current release's human and zebrafish file URLs from the folder listing."""
    for d in dirs:
        try:
            r = session.get(d, timeout=120)
        except Exception:
            continue
        if r.status_code != 200:
            continue
        out = {}
        for sp, pat in PANTHER_SPECIES.items():
            m = pat.findall(r.text)
            if m:
                out[sp] = d + sorted(m)[-1]
        if len(out) == 2:
            return out
    raise RuntimeError("could not find the PANTHER classification files; pass --panther-human/--panther-zebrafish URLs")


def head(path, n=3):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        lines = []
        for line in fh:
            if line.startswith("!"):
                continue
            lines.append(line.rstrip("\n")[:300])
            if len(lines) >= n:
                break
    return lines


def main(argv=None, session=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, help="default: <local master>/Data/Annotations")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--panther-human", help="URL of the PANTHER human classification file")
    ap.add_argument("--panther-zebrafish", help="URL of the PANTHER zebrafish classification file")
    ap.add_argument("--no-backup", action="store_true")
    args = ap.parse_args(argv)
    from ..backup import backup
    from ..layout import get_layout
    out = args.outdir or get_layout().annotations().local
    go, pan = out / "GO", out / "PANTHER"
    go.mkdir(parents=True, exist_ok=True)
    pan.mkdir(parents=True, exist_ok=True)
    if args.inspect:
        for p in sorted(list(go.iterdir()) + list(pan.iterdir())):
            if p.is_file() and not p.name.endswith(".json"):
                print(f"== {p.relative_to(out)}")
                for l in head(p):
                    print("   " + l)
        return
    if session is None:
        import requests
        session = requests.Session()
    record = {"downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"), "files": {}}
    for name, urls in GO_FILES.items():
        dest = go / name
        if dest.is_file() and not args.force:
            print(f"  GO/{name}: already here", file=sys.stderr)
            continue
        record["files"][f"GO/{name}"] = fetch(session, urls, dest)
        print(f"  GO/{name}: {dest.stat().st_size / 1e6:.1f} MB", file=sys.stderr)
    urls = {"human": args.panther_human, "zebrafish": args.panther_zebrafish}
    if not all(urls.values()) and (args.force or not all((pan / f"{sp}.tsv").is_file() for sp in urls)):
        found = panther_files(session)
        urls = {sp: urls[sp] or found[sp] for sp in urls}
    for sp, url in urls.items():
        dest = pan / f"{sp}.tsv"
        if dest.is_file() and not args.force:
            print(f"  PANTHER/{sp}.tsv: already here", file=sys.stderr)
            continue
        record["files"][f"PANTHER/{sp}.tsv"] = fetch(session, [url], dest)
        print(f"  PANTHER/{sp}.tsv ({url.rsplit('/', 1)[-1]}): {dest.stat().st_size / 1e6:.1f} MB", file=sys.stderr)
    (out / "download.json").write_text(json.dumps(record, indent=2))
    (out / "DATA_NOTICE.txt").write_text(
        "Gene Ontology annotations (CC BY 4.0): The Gene Ontology Consortium, Nucleic Acids Res (2023).\n"
        "PANTHER classification (http://www.pantherdb.org): Thomas PD et al. (2022) Protein Sci 31:8-22.\n")
    if not args.no_backup:
        print(f"Backed up to {backup(out)}", file=sys.stderr)


if __name__ == "__main__":
    main()


