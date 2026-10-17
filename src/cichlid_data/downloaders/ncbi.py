"""Genomes and annotations (the pipeline's Stage 1): GFF3/GTF, protein FASTA and sequence report
from NCBI Datasets, or from local files, for every species in config/species.tsv.

    cichlid-data download genomes [--species Bspl,Gacu] [--check-only] [--force] [--keep-zip] [--push]

Output per species, in Inputs/Proteins/downloads/<prefix>/:
    genomic.gff | genomic.gtf   annotation
    protein.faa                 all annotated proteins (isoform selection is in proteomes.py)
    sequence_report.jsonl       per-sequence roles; used to find mitochondrial contigs
    assembly_report.json        NCBI metadata (annotation release, BUSCO, gene counts)
    download.json               provenance: what was fetched, when, md5 of each file
and Inputs/Proteins/downloads/download_summary.tsv (or check_summary.tsv), one row per species.

source=ncbi rows are fetched from NCBI (verified against NCBI's checksums); source=local rows
copy the files named in gff_path and protein_path (and a sequence_report.jsonl beside the
annotation, if present). A species already downloaded and verified is skipped; one missing
locally but present in the cloud is pulled from there instead of being fetched again.

Set NCBI_API_KEY in the environment to raise the NCBI rate limit.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import requests

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]{1,15}$")
ACC_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")
COLUMNS = ["prefix", "species", "accession", "source", "annotation_format",
           "group", "gff_path", "protein_path", "notes"]

# our key -> (Datasets include_annotation_type, dataset_catalog fileType, output filename)
FILE_TYPES = {
    "gff3":      ("GENOME_GFF",      "GFF3",            "genomic.gff"),
    "gtf":       ("GENOME_GTF",      "GTF",             "genomic.gtf"),
    "protein":   ("PROT_FASTA",      "PROTEIN_FASTA",   "protein.faa"),
    "seqreport": ("SEQUENCE_REPORT", "SEQUENCE_REPORT", "sequence_report.jsonl"),
    "cds":       ("CDS_FASTA",       "CDS_NUCLEOTIDE_FASTA", "cds_from_genomic.fna"),   # variant stages
}


class ManifestError(Exception):
    pass


# --------------------------------------------------------------------------- manifest
def read_manifest(path, base=None):
    """Parse and validate the manifest. Collects every problem before failing. Relative gff_path
    and protein_path (source=local) are taken relative to base (the data root), when given."""
    rows, errors, seen = [], [], set()
    with open(path, newline="") as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    missing = {"prefix", "species", "source"} - set(reader.fieldnames or [])
    if missing:
        raise ManifestError(f"manifest missing required columns: {sorted(missing)}")

    for n, raw in enumerate(reader, start=2):
        r = {c: (raw.get(c) or "").strip() for c in COLUMNS}
        r["annotation_format"] = (r["annotation_format"] or "gff3").lower()
        r["source"] = r["source"].lower()
        where = f"row {n} ({r['prefix'] or '?'})"

        if not PREFIX_RE.match(r["prefix"]):
            errors.append(f"{where}: prefix must be 2-16 letters/digits, starting with a letter, no underscores")
        elif r["prefix"].lower() in seen:
            errors.append(f"{where}: duplicate prefix (case-insensitive)")
        seen.add(r["prefix"].lower())

        if not r["species"]:
            errors.append(f"{where}: species is empty")
        if r["annotation_format"] not in ("gff3", "gtf"):
            errors.append(f"{where}: annotation_format must be gff3 or gtf")

        if r["source"] == "ncbi":
            if not ACC_RE.match(r["accession"]):
                errors.append(f"{where}: accession '{r['accession']}' is not like GCF_000000000.1")
        elif r["source"] == "local":
            for col in ("gff_path", "protein_path"):
                if r[col] and base is not None and not Path(r[col]).is_absolute():
                    r[col] = str(Path(base) / r[col])
                if not r[col]:
                    errors.append(f"{where}: source=local needs {col}")
                elif not Path(r[col]).is_file():
                    errors.append(f"{where}: {col} not found: {r[col]}")
        else:
            errors.append(f"{where}: source must be ncbi or local")
        rows.append(r)

    if errors:
        raise ManifestError("manifest problems:\n  " + "\n  ".join(errors))
    if not rows:
        raise ManifestError("manifest has no species rows")
    return rows


# --------------------------------------------------------------------------- HTTP
def make_session():
    s = requests.Session()
    s.headers["accept"] = "application/json"
    key = os.environ.get("NCBI_API_KEY")
    if key:
        s.headers["api-key"] = key
    return s


def _request(session, url, retries=5, timeout=120, **kw):
    """GET with backoff on 429/5xx and connection errors."""
    err = None
    for attempt in range(retries):
        try:
            r = session.get(url, timeout=timeout, **kw)
        except requests.ConnectionError as e:
            err = e
        else:
            if r.status_code < 400:
                return r
            if r.status_code not in (429, 500, 502, 503, 504):
                r.raise_for_status()
            err = requests.HTTPError(f"HTTP {r.status_code} for {url}")
        time.sleep(2 ** attempt)
    raise err


def fetch_report(session, accession):
    r = _request(session, f"{API}/genome/accession/{accession}/dataset_report")
    reports = r.json().get("reports", [])
    if not reports:
        raise LookupError(f"NCBI returned no record for {accession}")
    return reports[0]


def download_package(session, accession, types, dest):
    """Stream the Datasets zip to dest (atomic: .part then rename)."""
    params = [("include_annotation_type", FILE_TYPES[t][0]) for t in types]
    url = f"{API}/genome/accession/{accession}/download"
    part = dest.with_suffix(dest.suffix + ".part")
    with _request(session, url, params=params, stream=True, timeout=600,
                  headers={"accept": "application/zip"}) as r, open(part, "wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    if not zipfile.is_zipfile(part):
        part.unlink()
        raise IOError(f"download for {accession} is not a zip file")
    part.rename(dest)


# --------------------------------------------------------------------------- metadata
def _dig(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def summarize_report(rep):
    return {
        "organism": _dig(rep, "organism", "organism_name", default=""),
        "tax_id": _dig(rep, "organism", "tax_id", default=""),
        "assembly_name": _dig(rep, "assembly_info", "assembly_name", default=""),
        "assembly_level": _dig(rep, "assembly_info", "assembly_level", default=""),
        "assembly_status": _dig(rep, "assembly_info", "assembly_status", default=""),
        "contig_n50": _dig(rep, "assembly_stats", "contig_n50", default=""),
        "annotation_name": _dig(rep, "annotation_info", "name", default=""),
        "annotation_provider": _dig(rep, "annotation_info", "provider", default=""),
        "annotation_release_date": _dig(rep, "annotation_info", "release_date", default=""),
        "busco_complete": _dig(rep, "annotation_info", "busco", "complete", default=""),
        "busco_duplicated": _dig(rep, "annotation_info", "busco", "duplicated", default=""),
        "busco_lineage": _dig(rep, "annotation_info", "busco", "busco_lineage", default=""),
        "protein_coding_genes": _dig(rep, "annotation_info", "stats", "gene_counts",
                                     "protein_coding", default=""),
    }


def check_report(row, s):
    """Return (errors, warnings) for one species' NCBI metadata."""
    errors, warnings = [], []
    if not s["annotation_name"]:
        errors.append("no annotation on this assembly")
    if s["assembly_status"] and s["assembly_status"].lower() != "current":
        warnings.append(f"assembly status is '{s['assembly_status']}', not current")
    if s["organism"] and s["organism"].lower() != row["species"].lower():
        warnings.append(f"NCBI organism is '{s['organism']}', manifest says '{row['species']}'")
    if row["accession"].startswith("GCA_") and s["annotation_name"]:
        warnings.append(f"GenBank annotation by '{s['annotation_provider']}', not RefSeq; QC before use")
    return errors, warnings


# --------------------------------------------------------------------------- files
def md5sum(path, block=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def count_fasta(path):
    with open(path, "rb") as fh:
        return sum(1 for line in fh if line.startswith(b">"))


def extract_package(zip_path, accession, types, outdir):
    """Pull the wanted files out of a Datasets zip, verify md5, rename to standard names."""
    with zipfile.ZipFile(zip_path) as z:
        names = set(z.namelist())
        # locate files via the catalog, fall back to the standard filenames
        located = {}
        catalog_name = "ncbi_dataset/data/dataset_catalog.json"
        if catalog_name in names:
            catalog = json.loads(z.read(catalog_name))
            for asm in catalog.get("assemblies", []):
                if asm.get("accession") != accession:
                    continue
                for f in asm.get("files", []):
                    located[f.get("fileType")] = "ncbi_dataset/data/" + f["filePath"]
        expected_md5 = {}
        if "md5sum.txt" in names:
            for line in z.read("md5sum.txt").decode().splitlines():
                parts = line.split()
                if len(parts) == 2:
                    expected_md5[parts[1]] = parts[0]

        out = {}
        for t in types:
            _, ftype, outname = FILE_TYPES[t]
            member = located.get(ftype) or f"ncbi_dataset/data/{accession}/{outname}"
            if member not in names:
                raise FileNotFoundError(f"{accession}: package has no {ftype} ({member})")
            dest = outdir / outname
            with z.open(member) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            digest = md5sum(dest)
            if member in expected_md5 and expected_md5[member] != digest:
                raise IOError(f"{accession}: md5 mismatch for {member}")
            out[t] = {"file": outname, "md5": digest, "bytes": dest.stat().st_size}

        report_member = "ncbi_dataset/data/assembly_data_report.jsonl"
        if report_member in names:
            first = z.read(report_member).decode().splitlines()[0]
            (outdir / "assembly_report.json").write_text(
                json.dumps(json.loads(first), indent=2))
    return out


def already_done(outdir):
    """True if download.json exists and every listed file is present with the recorded md5."""
    prov = outdir / "download.json"
    if not prov.is_file():
        return False
    try:
        files = json.loads(prov.read_text())["files"]
        return all((outdir / f["file"]).is_file() and md5sum(outdir / f["file"]) == f["md5"]
                   for f in files.values())
    except (KeyError, ValueError):
        return False


# --------------------------------------------------------------------------- per species
def process_ncbi(session, row, outdir, keep_zip):
    types = [row["annotation_format"], "protein", "seqreport"]
    zip_path = outdir / f"{row['accession']}.zip"
    download_package(session, row["accession"], types, zip_path)
    files = extract_package(zip_path, row["accession"], types, outdir)
    if not keep_zip:
        zip_path.unlink()
    return files


def process_local(row, outdir):
    files = {}
    for t, col in ((row["annotation_format"], "gff_path"), ("protein", "protein_path")):
        dest = outdir / FILE_TYPES[t][2]
        shutil.copyfile(row[col], dest)
        files[t] = {"file": dest.name, "md5": md5sum(dest), "bytes": dest.stat().st_size,
                    "source_path": str(Path(row[col]).resolve())}
    # a sequence report next to the annotation (e.g. from prepare_local_species.py) says which
    # sequences are chromosomes; copy it like NCBI's
    rep = Path(row["gff_path"]).parent / "sequence_report.jsonl"
    if rep.is_file():
        dest = outdir / "sequence_report.jsonl"
        shutil.copyfile(rep, dest)
        files["sequence_report"] = {"file": dest.name, "md5": md5sum(dest), "bytes": dest.stat().st_size,
                                    "source_path": str(rep.resolve())}
    return files


# =============================================================================== driver
SUMMARY_COLS = ["prefix", "species", "accession", "source", "status", "organism", "assembly_name",
                "assembly_level", "assembly_status", "contig_n50", "annotation_name",
                "annotation_provider", "annotation_release_date", "busco_lineage", "busco_complete",
                "busco_duplicated", "protein_coding_genes", "protein_records", "messages"]


def merge_summary(path, records):
    """Write the summary table, replacing the given species' rows and keeping the others."""
    old = []
    if path.is_file():
        with open(path, newline="") as fh:
            old = list(csv.DictReader(fh, delimiter="\t"))
    new = {r["prefix"]: r for r in records}
    rows = [r for r in old if r.get("prefix") not in new] + list(new.values())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, SUMMARY_COLS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["prefix"]))


def download(layout, manifest, species=None, check_only=False, force=False, keep_zip=False,
             session=None, store=None, say=None):
    """Download every manifest species (or the given prefixes). -> (summary records, failed prefixes)"""
    say = say or (lambda m: print(m, file=sys.stderr))
    rows = read_manifest(manifest, base=layout.local_root)   # raises ManifestError, listing every problem
    if species:
        wanted = {p.strip().lower() for p in species}
        rows = [r for r in rows if r["prefix"].lower() in wanted]
        if not rows:
            raise ManifestError(f"none of {', '.join(sorted(wanted))} is in {manifest}")
    session = session or make_session()
    root = layout.path("Inputs/Proteins/downloads").local
    root.mkdir(parents=True, exist_ok=True)
    summary, failed = [], []
    for row in rows:
        p = row["prefix"]
        sp_dir = layout.download(p).local
        rec = {"prefix": p, "species": row["species"], "accession": row["accession"],
               "source": row["source"], "status": "", "messages": ""}
        try:
            meta = {}
            if row["source"] == "ncbi":
                rep = fetch_report(session, row["accession"])
                meta = summarize_report(rep)
                errs, warns = check_report(row, meta)
                for w in warns:
                    say(f"[{p}] WARNING: {w}")
                if errs:
                    raise RuntimeError("; ".join(errs))
                rec["messages"] = "; ".join(warns)
            rec.update(meta)

            if not check_only and store is not None and not sp_dir.is_dir() and not force:
                try:                                    # already downloaded on another machine?
                    store.ensure(layout.download(p))
                    say(f"[{p}] pulled from the cloud")
                except Exception:
                    pass
            if check_only:
                rec["status"] = "checked"
                say(f"[{p}] {row['accession'] or 'local files'}: " + (
                    f"{meta.get('assembly_name', '')} ({meta.get('assembly_level', '')}), "
                    f"{meta.get('annotation_name', '')}, {meta.get('protein_coding_genes', '')} protein-coding genes"
                    if meta else "local files listed in species.tsv found"))
            elif already_done(sp_dir) and not force:
                rec["status"] = "up-to-date"
                say(f"[{p}] up to date, skipping")
            else:
                sp_dir.mkdir(parents=True, exist_ok=True)
                say(f"[{p}] downloading {row['accession'] or 'local files'} ...")
                files = (process_ncbi(session, row, sp_dir, keep_zip)
                         if row["source"] == "ncbi" else process_local(row, sp_dir))
                n_prot = count_fasta(sp_dir / FILE_TYPES["protein"][2])
                prov = {"manifest_row": row, "ncbi_metadata": meta, "files": files,
                        "protein_records": n_prot,
                        "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "api": API if row["source"] == "ncbi" else None}
                (sp_dir / "download.json").write_text(json.dumps(prov, indent=2))
                rec["status"] = "downloaded"
                rec["protein_records"] = n_prot
                say(f"[{p}] done: {n_prot} protein records")
        except Exception as e:  # keep going; report at the end
            rec["status"] = "FAILED"
            rec["messages"] = str(e)
            failed.append(p)
            say(f"[{p}] FAILED: {e}")
        summary.append(rec)
        time.sleep(0.2 if session.headers.get("api-key") else 0.4)  # stay under NCBI rate limits

    path = root / ("check_summary.tsv" if check_only else "download_summary.tsv")
    merge_summary(path, summary)
    counts = {}
    for r in summary:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    say(f"{len(summary)} species: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + f"; wrote {path}")
    return summary, failed


