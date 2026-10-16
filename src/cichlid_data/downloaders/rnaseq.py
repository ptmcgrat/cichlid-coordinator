"""cichlid-data summarize rnaseq --remote|--local DIR --dataset NAME [options]   (moved from the pipeline; default folder
Inputs/Expression/<dataset>)


Combine per-tissue nf-core/rnaseq results (e.g. the Brawand et al. 2014 cichlid tissues) into one
gene x tissue expression table for the browser and inversion_tissues.py.

    # straight from Dropbox (only the small merged tables are copied)
    python summarize_tissue_rnaseq.py --remote \
        "ptm_dropbox:/CoS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/2026_Nikesh/Brawand_data" \
        --dataset Brawand --species Mzeb
    # or from a local copy with the same layout (<tissue>/output/star_salmon/salmon.merged.gene_tpm.tsv)
    python summarize_tissue_rnaseq.py --local /path/to/Brawand_data --dataset Brawand --species Mzeb

Each tissue folder's star_salmon/salmon.merged.gene_tpm.tsv (columns gene_id, gene_name, then one
TPM column per sample) is read; tissues without it yet (e.g. still uploading) are skipped and listed.
Folders named like pools (all_tissues_pool, unknown_pool) are skipped unless --keep-pools.
When a tissue has several samples, the tissue value is their mean; the sample values are kept.
Writes <local master>/Data/Expression/<dataset>/ (or --outdir):
    tissue_tpm.tsv.gz   gene_id, gene_name, tissue, tpm, n_samples, sample_tpm
    samples.tsv         tissue, sample, source file
    summary.json        dataset, species, tissues, skipped, genes, date
"""
import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

TPM_FILE = "salmon.merged.gene_tpm.tsv"
REL = Path("output") / "star_salmon" / TPM_FILE


def fetch_remote(remote, rclone, dest):
    """Copy only the merged gene TPM tables from an rclone remote; -> local root."""
    cmd = [rclone, "copy", remote.rstrip("/"), str(dest), "--include", f"*/output/star_salmon/{TPM_FILE}"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"rclone copy failed: {r.stderr.strip()[:500]}")
    return dest


def read_tpm(path):
    """-> (samples, {gene_id: (gene_name, [tpm per sample])})"""
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        lower = [h.lower() for h in head]
        gi = lower.index("gene_id") if "gene_id" in lower else 0
        gn = lower.index("gene_name") if "gene_name" in lower else None
        samples = [(i, h) for i, h in enumerate(head) if i not in (gi, gn)]
        out = {}
        for t in rd:
            if not t or len(t) < len(head):
                continue
            vals = []
            for i, _ in samples:
                try:
                    vals.append(float(t[i]))
                except ValueError:
                    vals.append(0.0)
            out[t[gi]] = (t[gn] if gn is not None else t[gi], vals)
    return [h for _, h in samples], out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--remote", help="rclone path of the folder holding one subfolder per tissue")
    src.add_argument("--local", type=Path, help="local folder holding one subfolder per tissue")
    ap.add_argument("--dataset", required=True, help="name, e.g. Brawand")
    ap.add_argument("--species", default="Mzeb", help="browser prefix of the species whose genes these are")
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--keep-pools", action="store_true")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--rclone", default=os.environ.get("RCLONE", "rclone"))
    args = ap.parse_args(argv)
    from ..backup import backup
    from ..layout import get_layout
    out = args.outdir or get_layout().expression(args.dataset).local
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        root = fetch_remote(args.remote, args.rclone, Path(tmp)) if args.remote else args.local
        tissues, skipped, table, samples_rows = [], [], {}, []
        for d in sorted(p for p in Path(root).iterdir() if p.is_dir() and not p.name.startswith(".")):
            if "pool" in d.name.lower() and not args.keep_pools:
                skipped.append(f"{d.name} (pool)")
                continue
            f = d / REL
            if not f.is_file():
                skipped.append(f"{d.name} (no {TPM_FILE} yet)")
                continue
            samples, data = read_tpm(f)
            tissue = d.name
            tissues.append(tissue)
            for s in samples:
                samples_rows.append({"tissue": tissue, "sample": s, "source": str(Path(d.name) / REL)})
            for g, (name, vals) in data.items():
                table.setdefault(g, [name, {}])[1][tissue] = vals
    if not tissues:
        raise SystemExit("no tissue has a merged gene TPM table yet")
    rows = 0
    with gzip.open(out / "tissue_tpm.tsv.gz", "wt") as fh:
        fh.write("gene_id\tgene_name\ttissue\ttpm\tn_samples\tsample_tpm\n")
        for g, (name, by_t) in sorted(table.items()):
            for t, vals in by_t.items():
                mean = sum(vals) / len(vals) if vals else 0.0
                fh.write(f"{g}\t{name}\t{t}\t{mean:.4g}\t{len(vals)}\t{','.join(f'{v:.4g}' for v in vals)}\n")
                rows += 1
    with open(out / "samples.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, ["tissue", "sample", "source"], delimiter="\t")
        w.writeheader()
        w.writerows(samples_rows)
    info = {"dataset": args.dataset, "species": args.species, "tissues": tissues, "skipped": skipped,
            "genes": len(table), "rows": rows, "source": args.remote or str(args.local),
            "made": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (out / "summary.json").write_text(json.dumps(info, indent=2))
    print(f"{args.dataset}: {len(tissues)} tissues ({', '.join(tissues)}), {len(table)} genes -> {out}", file=sys.stderr)
    if skipped:
        print(f"  skipped: {', '.join(skipped)}", file=sys.stderr)
    if not args.no_backup:
        print(f"Backed up to {backup(out)}", file=sys.stderr)


if __name__ == "__main__":
    main()
