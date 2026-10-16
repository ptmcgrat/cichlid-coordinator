"""cichlid-data download expression --source hpa|bgee|daniocell   (moved from the pipeline; default folders Inputs/Expression/...)


Download expression reference data for the gene browser and back it up to the lab Dropbox.

    python download_expression.py              # Human Protein Atlas single cell type data
    python download_expression.py --no-backup
    python download_expression.py --base-url https://v24.proteinatlas.org/download   # a specific release
    python download_expression.py --source bgee     # Bgee zebrafish expression calls, adult tissues
    python download_expression.py --source daniocell  # Daniocell Seurat object + annotations (raw/daniocell)

Bgee (Bastian et al. 2021, Nucleic Acids Res 49:D831), into <local master>/Data/Expression/Bgee/:
    Danio_rerio_expr_simple_all_conditions.tsv.gz   expression calls per gene, anatomical entity,
                                                    stage, sex and strain (as downloaded)
    bgee_adult.tsv.gz      expressed calls at adult stages, one row per gene and tissue: highest
                           expression score, best call quality, sexes seen in
    download.json

Human Protein Atlas (Karlsson et al. 2021, Sci Adv 7:eabh2169), into <local master>/Data/Expression/HPA/:
    rna_single_cell_type.tsv.zip   expression (nTPM) per gene and cell type
    proteinatlas.tsv.zip           gene summary, incl. single cell type specificity classes
    download.json                  source URLs, date, sizes, MD5s
HPA has moved its download paths between releases, so each file is tried at the known locations
in turn. Stage 9 (--expression) loads these into the site database.
"""
import argparse
import hashlib
import json
import os
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import requests

from ..backup import backup, cloud_path, settings  # noqa: E402
from ..layout import get_layout  # noqa: E402

HPA_BASES = ["https://www.proteinatlas.org/download/tsv", "https://www.proteinatlas.org/download"]
DANIOCELL_BASE = "https://daniocell.nichd.nih.gov/assets/obj-download"
DANIOCELL_FILES = ["Daniocell2023_SeuratV4.rds", "cluster_annotations.csv", "daniocell_load.R"]
BGEE_FILE = "Danio_rerio_expr_simple_all_conditions.tsv.gz"
BGEE_BASES = ["https://www.bgee.org/ftp/current/download/calls/expr_calls",
              "https://bgee.org/ftp/current/download/calls/expr_calls"]
ADULT_STAGE_IDS = {"UBERON:0000113"}          # post-juvenile adult stage
QUALITY_RANK = {"gold": 2, "gold quality": 2, "high quality": 2, "silver": 1, "silver quality": 1, "low quality": 1}
HPA_FILES = ["rna_single_cell_type.tsv.zip", "proteinatlas.tsv.zip"]
NOTICE = """Human Protein Atlas, single cell type data (https://www.proteinatlas.org/about/download).
Karlsson M et al. (2021) A single-cell type transcriptomics map of human tissues. Sci Adv 7:eabh2169.
HPA data are available under CC BY-SA 4.0; cite the HPA and the paper above when showing them.
"""


BGEE_NOTICE = """Bgee gene expression calls for zebrafish (https://www.bgee.org).
Bastian FB et al. (2021) The Bgee suite: integrated curated expression atlas and comparative transcriptomics
in animals. Nucleic Acids Res 49:D831-D847. Bgee data are available under CC0.
"""


def fetch(session, urls, dest, tries=4, check_zip=True):
    """Download the first URL that works to dest (via .part). -> the URL used."""
    errors = []
    for url in urls:
        for attempt in range(tries):
            try:
                with session.get(url, stream=True, timeout=300) as r:
                    if r.status_code == 404:
                        errors.append(f"{url}: not found")
                        break
                    if r.status_code != 200:
                        raise IOError(f"HTTP {r.status_code}")
                    part = dest.with_suffix(dest.suffix + ".part")
                    with open(part, "wb") as fh:
                        for chunk in r.iter_content(4 << 20):
                            fh.write(chunk)
                if check_zip and not zipfile.is_zipfile(part):
                    part.unlink(missing_ok=True)
                    errors.append(f"{url}: not a zip file")
                    break
                part.replace(dest)
                return url
            except (requests.ConnectionError, requests.Timeout, IOError) as e:
                errors.append(f"{url}: {e}")
                time.sleep(2 ** attempt)
    raise IOError(f"could not download {dest.name}:\n  " + "\n  ".join(errors))


def _col(fields, *wanted):
    low = {f.lower(): f for f in fields}
    for w in wanted:
        if w in low:
            return low[w]
    for w in wanted:
        for f in fields:
            if w in f.lower():
                return f
    return None


def is_adult(stage_name, stage_id=""):
    """Bgee names the zebrafish adult stage 'post-juvenile' (UBERON:0000113, post-juvenile adult stage)."""
    name = (stage_name or "").strip().lower()
    return (stage_id or "").strip() in ADULT_STAGE_IDS or name.startswith("post-juvenile") or "adult" in name


def summarize_bgee(src, dest):
    """Bgee all-conditions expression calls -> adult tissues: per gene and anatomical entity, the
    highest expression score among expressed calls at adult stages. -> (genes, rows)"""
    import csv
    import gzip
    best = {}
    with gzip.open(src, "rt", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        f = rd.fieldnames or []
        c = {k: _col(f, *v) for k, v in {
            "gene": ("gene id",), "name": ("gene name",), "aid": ("anatomical entity id", "anat. entity id"),
            "aname": ("anatomical entity name", "anat. entity name"),
            "stage": ("developmental stage name", "stage name"),
            "stage_id": ("developmental stage id", "stage id"), "sex": ("sex",),
            "expr": ("expression",), "quality": ("call quality",), "score": ("expression score",)}.items()}
        missing = [k for k in ("gene", "aname", "stage", "expr", "score") if not c[k]]
        if missing:
            raise ValueError(f"unexpected Bgee columns, missing {missing}: {f}")
        for r in rd:
            if not is_adult(r[c["stage"]], r.get(c["stage_id"]) if c["stage_id"] else ""):
                continue
            if (r[c["expr"]] or "").strip().lower() not in ("present", "expressed"):
                continue
            try:
                score = float(r[c["score"]])
            except (TypeError, ValueError):
                continue
            key = (r[c["gene"]], r[c["aname"]])
            q = QUALITY_RANK.get((r.get(c["quality"]) or "").strip().lower(), 0) if c["quality"] else 0
            sex = (r.get(c["sex"]) or "").strip() if c["sex"] else ""
            cur = best.get(key)
            if cur is None:
                best[key] = [r.get(c["name"], "") if c["name"] else "", r.get(c["aid"], "") if c["aid"] else "",
                             score, q, {sex} if sex and sex != "any" else set()]
            else:
                if score > cur[2]:
                    cur[2] = score
                cur[3] = max(cur[3], q)
                if sex and sex != "any":
                    cur[4].add(sex)
    with gzip.open(dest, "wt") as out:
        out.write("gene_id\tsymbol\tanatomy_id\ttissue\tscore\tquality\tsexes\n")
        for (g, t), (name, aid, score, q, sexes) in sorted(best.items()):
            out.write(f"{g}\t{name}\t{aid}\t{t}\t{score:.2f}\t{['', 'silver', 'gold'][q]}\t{','.join(sorted(sexes))}\n")
    return len({g for g, _ in best}), len(best)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main(argv=None, session=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, help="default: <local master>/Data/Expression/HPA")
    ap.add_argument("--base-url", help="download directory of a specific HPA release")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--source", choices=["hpa", "bgee", "daniocell"], default="hpa")
    ap.add_argument("--force", action="store_true", help="download again even if the file is already here (Bgee)")
    args = ap.parse_args(argv)

    L = get_layout()
    if args.source == "daniocell":
        out = args.outdir or L.expression("raw/daniocell").local
        out.mkdir(parents=True, exist_ok=True)
        session = session or requests.Session()
        for name in DANIOCELL_FILES:
            if (out / name).is_file() and not args.force:
                print(f"  {name}: already here", file=sys.stderr)
                continue
            fetch(session, [f"{(args.base_url or DANIOCELL_BASE).rstrip('/')}/{name}"], out / name, check_zip=False)
            print(f"  {name}: {(out / name).stat().st_size / 1e6:.1f} MB", file=sys.stderr)
        (out / "DATA_NOTICE.txt").write_text(
            "Daniocell (https://daniocell.nichd.nih.gov), Farrell lab, NICHD. Sur A et al. (2023) Single-cell analysis\n"
            "of shared signatures and transcriptional diversity during zebrafish development. Dev Cell 58:3028.\n")
        print(f"Daniocell files in {out}. Next: Rscript export_seurat.R {out / DANIOCELL_FILES[0]} {out / 'export'}",
              file=sys.stderr)
        return
    if args.source == "bgee":
        out = args.outdir or L.expression("Bgee").local
        out.mkdir(parents=True, exist_ok=True)
        session = session or requests.Session()
        bases = [args.base_url.rstrip("/")] if args.base_url else BGEE_BASES
        prev = json.loads((out / "download.json").read_text()) if (out / "download.json").is_file() else {}
        reused = (out / BGEE_FILE).is_file() and not args.force
        if reused:
            url = prev.get("files", {}).get(BGEE_FILE, {}).get("url", "previous download")
            print(f"Using the existing {BGEE_FILE} (--force to download again)", file=sys.stderr)
        else:
            url = fetch(session, [f"{b}/{BGEE_FILE}" for b in bases], out / BGEE_FILE, check_zip=False)
        genes, rows = summarize_bgee(out / BGEE_FILE, out / "bgee_adult.tsv.gz")
        (out / "download.json").write_text(json.dumps({
            "source": "Bgee", "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "files": {BGEE_FILE: {"url": url, "bytes": (out / BGEE_FILE).stat().st_size, "md5": md5(out / BGEE_FILE)}},
            "adult_genes": genes, "adult_rows": rows}, indent=2))
        (out / "DATA_NOTICE.txt").write_text(BGEE_NOTICE)
        print(f"Bgee: {genes} zebrafish genes expressed in adult tissues ({rows} gene x tissue rows)", file=sys.stderr)
        if not args.no_backup:
            print(f"Backed up to {backup(out)}", file=sys.stderr)
        return
    out = args.outdir or L.expression("HPA").local
    out.mkdir(parents=True, exist_ok=True)
    session = session or requests.Session()
    bases = [args.base_url.rstrip("/")] if args.base_url else HPA_BASES
    record = {"source": "Human Protein Atlas", "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "files": {}}
    for name in HPA_FILES:
        url = fetch(session, [f"{b}/{name}" for b in bases], out / name)
        record["files"][name] = {"url": url, "bytes": (out / name).stat().st_size, "md5": md5(out / name)}
        print(f"  {name}: {record['files'][name]['bytes'] / 1e6:.1f} MB from {url}", file=sys.stderr)
    (out / "download.json").write_text(json.dumps(record, indent=2))
    (out / "DATA_NOTICE.txt").write_text(NOTICE)
    if not args.no_backup:
        print(f"Backed up to {backup(out)}", file=sys.stderr)


if __name__ == "__main__":
    main()
