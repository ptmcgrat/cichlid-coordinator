"""cichlid-data populations samples [options]   (moved from the pipeline; default folder Inputs/Populations/Blumer2025)


Snapshot the sample metadata for the Blumer et al. 2025 VCF samples, so the variant analysis uses
a fixed, recorded table.

    python blumer_samples.py                  # Table S1 of the paper (Data/Blumer2025/Table_S1.tsv)
    python blumer_samples.py --sheet          # also fill gaps from the lab Sample Database sheet
    python blumer_samples.py --no-backup

Source of metadata, in order: Table S1 of Blumer et al. 2025 (sequence_id = VCF sample name; clade
becomes Ecogroup), then, with --sheet, the lab's Sample Database for any sample Table S1 lacks.

Steps:
  1. sample names from the header of one Blumer VCF (all 22 share the same samples)
  2. the SampleDatabase and DNAReadsDatabase tabs of the Sample Database sheet, read with the lab's
     service-account key (downloaded with rclone from the CichlidSequencingData Dropbox folder to a
     private file, ~/.config/cichlid-trees/sheets_creds.json, mode 600; never uploaded)
  3. Data/Blumer2025/samples.tsv: one row per VCF sample with SampleID, Species, Ecogroup,
     Subgroup, Sex, ProjectIDs and in_database; plus samples_report.txt listing VCF samples missing
     from the sheet. Both are backed up next to the VCFs in Dropbox.
"""
import argparse
import csv
import gzip
import os
import stat
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from ..backup import backup, cloud_path, settings  # noqa: E402
from ..layout import get_layout  # noqa: E402

SHEET_ID = "1NmgB_TWoO01Qz2ufvECuZFkxXayhUsyu8wQGStVB_8k"            # Sample Database
CREDS_CLOUD = "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/cichlidsrsequencing_api_creds.json"
CREDS_LOCAL = Path.home() / ".config" / "cichlid-trees" / "sheets_creds.json"
COLUMNS = ["SampleID", "Species", "Ecogroup", "Subgroup", "Sex", "Location", "Sublocation", "SpecimenID",
           "BioSample", "BioProject", "Inversion10", "ProjectIDs", "source", "in_database"]


def read_table_s1(path):
    """Blumer et al. 2025 Table S1 -> {sequence_id: metadata in our column names}."""
    text = Path(path).read_text(encoding="utf-8-sig").replace("\r", "")
    out = {}
    for r in csv.DictReader(text.splitlines(), delimiter="\t"):
        sid = (r.get("sequence_id") or "").strip()
        if not sid:
            continue
        clean = lambda v: "" if (v or "").strip() in ("", "NA") else v.strip()
        out[sid] = {"Species": f"{clean(r.get('genus'))} {clean(r.get('species'))}".strip(),
                    "Ecogroup": clean(r.get("clade")), "Subgroup": "", "Sex": clean(r.get("sex")),
                    "Location": clean(r.get("location")), "Sublocation": clean(r.get("sublocation")),
                    "SpecimenID": clean(r.get("specimen_id")), "BioSample": clean(r.get("Biosample")),
                    "BioProject": clean(r.get("Bioproject")), "Inversion10": clean(r.get("Inversion10"))}
    return out


def vcf_samples(path):
    """Sample names from a (b)gzipped VCF header."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    raise ValueError(f"no #CHROM header line in {path}")


def ensure_credentials(rclone="rclone", path=CREDS_LOCAL):
    if path.is_file():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    os.chmod(path.parent, stat.S_IRWXU)
    r = subprocess.run([rclone, "copyto", CREDS_CLOUD, str(path)], capture_output=True, text=True)
    if r.returncode != 0 or not path.is_file():
        raise RuntimeError(f"could not download the Google Sheets key ({r.stderr.strip()[-200:]})")
    os.chmod(path, stat.S_IRUSR | stat.S_IWUSR)
    return path


def read_sheet(creds):
    """-> (sample rows, reads rows) as lists of dicts, from the SampleDatabase and DNAReadsDatabase tabs."""
    import gspread
    from google.oauth2.service_account import Credentials
    scopes = ["https://www.googleapis.com/auth/spreadsheets.readonly", "https://www.googleapis.com/auth/drive.readonly"]
    gc = gspread.authorize(Credentials.from_service_account_file(str(creds), scopes=scopes))
    book = gc.open_by_key(SHEET_ID)
    return (book.worksheet("SampleDatabase").get_all_records(),
            book.worksheet("DNAReadsDatabase").get_all_records())


def build_table(names, samples, reads, table_s1=None):
    by_id = {}
    for r in samples:
        sid = str(r.get("SampleID", "")).strip()
        if sid:
            by_id.setdefault(sid, r)
    projects = {}
    for r in reads:
        sid = str(r.get("SampleID", "")).strip()
        p = str(r.get("ProjectID", "")).strip()
        if sid and p:
            projects.setdefault(sid, set()).add(p)
    table_s1 = table_s1 or {}
    rows, missing = [], []
    for n in names:
        row = {c: "" for c in COLUMNS}
        row["SampleID"] = n
        t = table_s1.get(n)
        r = by_id.get(n)
        if t:
            row.update(t)
            row["source"] = "Table S1"
        elif r:
            row.update({c: str(r.get(c, "")).strip() for c in ("Species", "Ecogroup", "Subgroup", "Sex", "Inversion10")})
            row["source"] = "lab sheet"
        else:
            missing.append(n)
        row["ProjectIDs"] = ",".join(sorted(projects.get(n, ())))
        row["in_database"] = bool(t or r)
        rows.append(row)
    return rows, missing


def main(argv=None, sheet_reader=read_sheet):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, help="default: <local master>/Data/Blumer2025")
    ap.add_argument("--vcf", help="VCF to read sample names from (default: chr1.vcf.gz in --dir)")
    ap.add_argument("--table-s1", type=Path, help="Blumer et al. Table S1 (default: Table_S1.tsv in --dir)")
    ap.add_argument("--sheet", action="store_true", help="also use the lab Sample Database for samples not in Table S1")
    ap.add_argument("--creds", type=Path, help=f"service-account key (default: {CREDS_LOCAL}, downloaded if missing)")
    ap.add_argument("--no-backup", action="store_true")
    args = ap.parse_args(argv)

    d = args.dir or get_layout().population("Blumer2025").local
    vcf = Path(args.vcf) if args.vcf else d / "chr1.vcf.gz"
    if not vcf.is_file():
        sys.exit(f"{vcf} not found; run download_blumer.py first")
    names = vcf_samples(vcf)
    print(f"{len(names)} samples in {vcf.name}", file=sys.stderr)

    t1_path = args.table_s1 or d / "Table_S1.tsv"
    table_s1 = read_table_s1(t1_path) if t1_path.is_file() else {}
    if table_s1:
        print(f"Table S1: {len(table_s1)} samples ({t1_path.name})", file=sys.stderr)
    samples, reads = [], []
    if args.sheet or not table_s1:
        if not table_s1:
            print("No Table S1 found; using the lab Sample Database only", file=sys.stderr)
        creds = args.creds or ensure_credentials(settings().rclone)
        samples, reads = sheet_reader(creds)
    rows, missing = build_table(names, samples, reads, table_s1)

    with open(d / "samples.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    eco = {}
    for r in rows:
        if r["in_database"]:
            eco[r["Ecogroup"] or "(none)"] = eco.get(r["Ecogroup"] or "(none)", 0) + 1
    report = [f"Sample metadata for {vcf.name} samples (Table S1{' + lab sheet' if samples else ''}), "
              f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}",
              f"VCF samples: {len(names)}; found in the sheet: {len(names) - len(missing)}; missing: {len(missing)}",
              f"Species: {len({r['Species'] for r in rows if r['Species']})}",
              "Ecogroups: " + ", ".join(f"{k} {v}" for k, v in sorted(eco.items(), key=lambda kv: -kv[1]))]
    if missing:
        report += ["", "VCF samples not in the sheet:"] + missing
    (d / "samples_report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report[:4]), file=sys.stderr)
    if missing:
        print(f"  ({len(missing)} missing samples listed in samples_report.txt)", file=sys.stderr)

    if not args.no_backup:
        backup(d / "samples.tsv")                  # pushes the population folder (files are copied only if new or changed)
        print("Backed up samples.tsv and samples_report.txt to Dropbox.", file=sys.stderr)


if __name__ == "__main__":
    main()


