"""Proteomes: one representative protein per protein-coding gene, for every species.

    cichlid-data proteomes [--species Mzeb,Bspl] [--select complete-longest|longest] [--force] [--push]

For each species reads its download (Inputs/Proteins/downloads/<prefix>/: genomic.gff or
genomic.gtf, protein.faa, sequence_report.jsonl) and writes, through the schemas:

    Inputs/Proteins/proteomes/<prefix>.faa        one protein per gene, IDs like Onil_100534530
    Inputs/Proteins/proteomes/<prefix>.genes.tsv  per-gene table (symbol, location, chosen isoform, flags)
    Inputs/Proteins/proteomes/stage2_summary.tsv  per-species QC across the run

then rebuilds the catalog. Species whose proteome is newer than their download are skipped
unless --force.

Selection (default complete-longest), unchanged from the pipeline's Stage 2:
    1. complete models before partial ones (partial only if a gene has no complete model)
    2. longest protein
    3. models NCBI does not flag as LOW QUALITY PROTEIN
    4. protein ID, so reruns are deterministic
--select longest drops rule 1.

CDS features are traced up their Parent chain to the top-level gene, so CDSs attached
directly to a gene (mitochondrial genes in RefSeq) are handled without special cases.
Works on NCBI GFF3, NCBI GTF, and GFF3/GTF without protein_id (e.g. Braker), where the
transcript ID is used to find the protein in the FASTA.

Sequences are written for downstream tools: terminal '*' removed, selenocysteine U and any
other non-standard letter replaced with X (counted in the table).
"""
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from . import schemas
from .catalog import read_manifest

SKIP_TYPES = {"exon", "cDNA_match", "match", "region", "match_part", "five_prime_UTR",
              "three_prime_UTR", "start_codon", "stop_codon", "intron"}
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWYX")
SAFE_ID = re.compile(r"[^A-Za-z0-9.\-]")   # IDs must survive hmmscan, MAFFT, RAxML, Newick


# --------------------------------------------------------------------------- FASTA
def read_fasta(path):
    """-> {id: (description, sequence)}"""
    seqs, name, desc, buf = {}, None, "", []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = (desc, "".join(buf))
                head = line[1:].rstrip("\n")
                name, _, desc = head.partition(" ")
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = (desc, "".join(buf))
    return seqs


def clean_sequence(seq):
    """Strip terminal stop, replace non-standard residues with X. -> (seq, n_replaced, internal_stops)"""
    seq = seq.upper().rstrip("*")
    internal_stops = seq.count("*")
    out = [c if c in STANDARD_AA else "X" for c in seq]
    n_replaced = sum(1 for a, b in zip(seq, out) if a != b)
    return "".join(out), n_replaced, internal_stops


# --------------------------------------------------------------------------- annotation parsing
def parse_gff3_attrs(s):
    attrs = {}
    for part in s.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = unquote(v.strip())
    return attrs


def parse_gtf_attrs(s):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', s):
        k, v = m.group(1), m.group(2)
        if k == "db_xref":                      # repeated key; keep all
            attrs[k] = attrs[k] + "," + v if k in attrs else v
        else:
            attrs.setdefault(k, v)
    return attrs


def geneid_from_dbxref(attrs):
    m = re.search(r"GeneID:(\d+)", attrs.get("Dbxref", "") + "," + attrs.get("db_xref", ""))
    return m.group(1) if m else None


def is_partial(attrs):
    return attrs.get("partial", "").lower() == "true" or "start_range" in attrs or "end_range" in attrs


class Annotation:
    """Genes, transcripts and CDS->protein links from one annotation file."""

    def __init__(self):
        self.genes = {}                       # gene key -> info dict
        self.proteins = defaultdict(lambda: {"partial": False, "transcript": None, "gene": None,
                                             "product": ""})
        self.n_transcripts = defaultdict(set)

    # ---- GFF3
    def load_gff3(self, path):
        feats = {}          # ID -> (type, parent, attrs, seqid, start, end, strand)
        cds_lines = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                t = line.rstrip("\n").split("\t")
                if len(t) < 9 or t[2] in SKIP_TYPES:
                    continue
                attrs = parse_gff3_attrs(t[8])
                if t[2] == "CDS":   # keep only what selection needs
                    cds_lines.append((attrs.get("Parent", "").split(",")[0], attrs.get("protein_id"),
                                      attrs.get("ID"), attrs.get("product", ""), is_partial(attrs)))
                    continue
                fid = attrs.get("ID")
                if fid and fid not in feats:
                    parent = attrs.get("Parent", "").split(",")[0] or None
                    feats[fid] = (t[2], parent, attrs, t[0], int(t[3]), int(t[4]), t[6])

        def top(fid):
            seen = 0
            while fid in feats and feats[fid][1] and seen < 10:
                fid, seen = feats[fid][1], seen + 1
            return fid

        for fid, (ftype, parent, attrs, seqid, start, end, strand) in feats.items():
            if parent is None and ftype in ("gene", "pseudogene"):
                self._add_gene(fid, ftype, attrs, seqid, start, end, strand)

        for parent, protein_id, cds_id, product, cds_partial in cds_lines:
            gene = top(parent)
            if gene not in self.genes:
                continue
            transcript = parent if parent != gene else None
            pkey = protein_id or transcript or cds_id
            if not pkey:
                continue
            p = self.proteins[pkey]
            p["gene"], p["transcript"] = gene, transcript
            p["product"] = p["product"] or product
            tattrs = feats[transcript][2] if transcript in feats else {}
            p["partial"] = p["partial"] or cds_partial or is_partial(tattrs)
            self.n_transcripts[gene].add(transcript or pkey)
        return self

    # ---- GTF
    def load_gtf(self, path):
        gene_rows, cds_rows = {}, []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                t = line.rstrip("\n").split("\t")
                if len(t) < 9:
                    continue
                attrs = parse_gtf_attrs(t[8])
                gid = attrs.get("gene_id")
                if not gid:
                    continue
                if t[2] == "gene":
                    gene_rows[gid] = ("gene", attrs, t[0], int(t[3]), int(t[4]), t[6])
                elif t[2] == "CDS":
                    cds_rows.append(attrs)
                    # genes without a 'gene' line (e.g. Braker GTF): build extent from CDS lines
                    old = gene_rows.get(gid)
                    if old is None:
                        gene_rows[gid] = ("implied", attrs, t[0], int(t[3]), int(t[4]), t[6])
                    elif old[0] == "implied":
                        gene_rows[gid] = ("implied", old[1], old[2], min(old[3], int(t[3])),
                                          max(old[4], int(t[4])), old[5])
        for gid, (kind, attrs, seqid, start, end, strand) in gene_rows.items():
            self._add_gene(gid, "gene", attrs, seqid, start, end, strand)
        for attrs in cds_rows:
            gid, tid = attrs["gene_id"], attrs.get("transcript_id")
            pkey = attrs.get("protein_id") or tid
            if gid not in self.genes or not pkey:
                continue
            p = self.proteins[pkey]
            p["gene"], p["transcript"] = gid, tid
            p["product"] = p["product"] or attrs.get("product", "")
            p["partial"] = p["partial"] or is_partial(attrs)
            self.n_transcripts[gid].add(tid or pkey)
        return self

    def _add_gene(self, key, ftype, attrs, seqid, start, end, strand):
        biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or attrs.get("biotype")
        if ftype == "pseudogene":
            biotype = biotype or "pseudogene"
        self.genes[key] = {
            "key": key, "type": ftype, "biotype": biotype,   # None = unknown (non-NCBI)
            "ncbi_gene_id": geneid_from_dbxref(attrs),
            "symbol": attrs.get("gene") or attrs.get("gene_name") or attrs.get("Name") or "",
            "description": attrs.get("description", ""),
            "seqid": seqid, "start": start, "end": end, "strand": strand,
        }


# --------------------------------------------------------------------------- selection
def select_proteins(ann, fasta, mode="complete-longest", biotypes=("protein_coding",)):
    """-> list of per-gene result dicts, plus counters."""
    by_gene = defaultdict(list)
    stats = defaultdict(int)
    for pkey, p in ann.proteins.items():
        if pkey not in fasta:
            stats["cds_protein_not_in_fasta"] += 1
            continue
        desc, seq = fasta[pkey]
        low_q = "LOW QUALITY PROTEIN" in desc
        partial = p["partial"] or bool(re.search(r"\bpartial\b", desc))
        by_gene[p["gene"]].append({"protein": pkey, "transcript": p["transcript"],
                                   "length": len(seq.rstrip("*")), "partial": partial,
                                   "low_quality": low_q, "seq": seq,
                                   "product": p["product"] or desc})

    linked = {pk for pk, p in ann.proteins.items() if pk in fasta}
    stats["fasta_proteins_not_in_annotation"] = sum(1 for k in fasta if k not in linked)

    def rank(c):
        complete_first = (0 if not c["partial"] else 1) if mode == "complete-longest" else 0
        return (complete_first, -c["length"], c["low_quality"], c["protein"])

    results = []
    for gkey, g in ann.genes.items():
        cands = by_gene.get(gkey, [])
        # NCBI genes carry gene_biotype; genes from other annotators count if they have a CDS
        coding = (g["biotype"] in biotypes) if g["biotype"] else bool(cands)
        if not coding:
            stats[f"skipped_biotype:{g['biotype']}"] += 1
            continue
        stats["protein_coding_genes"] += 1
        if not cands:
            stats["coding_genes_without_protein"] += 1
            results.append({**g, "chosen": None, "n_proteins": 0,
                            "n_transcripts": len(ann.n_transcripts.get(gkey, ()))})
            continue
        best = min(cands, key=rank)
        results.append({**g, "chosen": best, "n_proteins": len(cands),
                        "n_transcripts": len(ann.n_transcripts.get(gkey, ())),
                        "all_partial": all(c["partial"] for c in cands)})
    return results, stats


def output_id(prefix, g):
    raw = g["ncbi_gene_id"] or g["key"]
    raw = re.sub(r"^gene[-:]", "", raw)
    return f"{prefix}_{SAFE_ID.sub('-', raw)}"


def mito_seqids(seqreport_path):
    ids = set()
    if seqreport_path and Path(seqreport_path).is_file():
        with open(seqreport_path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                r = json.loads(line)
                loc = r.get("assignedMoleculeLocationType") or r.get("assigned_molecule_location_type", "")
                if loc.lower() == "mitochondrion":
                    for k in ("refseqAccession", "genbankAccession", "refseq_accession",
                              "genbank_accession", "sequenceName", "sequence_name"):
                        if r.get(k):
                            ids.add(r[k])
    return ids


# --------------------------------------------------------------------------- per species
GENE_COLS = ["gene_id", "ncbi_gene_id", "annotation_gene_key", "symbol", "description", "biotype",
             "seqid", "start", "end", "strand", "mitochondrial", "n_transcripts", "n_proteins",
             "chosen_protein", "chosen_transcript", "protein_length", "partial", "low_quality",
             "residues_replaced", "internal_stops", "product"]



# =============================================================================== driver
SUMMARY = "stage2_summary.tsv"


def build_species(prefix, fmt, download_dir, out_dir, species_name="", mode="complete-longest"):
    """One species: select proteins, write the genes table and proteome through the schemas.
    -> summary dict. Nothing is left behind if the result does not validate."""
    sp = Path(download_dir)
    ann_path = sp / ("genomic.gff" if fmt == "gff3" else "genomic.gtf")
    for f in (ann_path, sp / "protein.faa"):
        if not f.is_file():
            raise FileNotFoundError(f"{f} missing; download the species first")
    ann = Annotation().load_gff3(ann_path) if fmt == "gff3" else Annotation().load_gtf(ann_path)
    fasta = read_fasta(sp / "protein.faa")
    results, stats = select_proteins(ann, fasta, mode)
    mito = mito_seqids(sp / "sequence_report.jsonl")

    ids_seen, rows, proteins = {}, [], []
    lengths, n_partial, n_lowq, n_mito, n_repl = [], 0, 0, 0, 0
    for g in sorted(results, key=lambda g: (g["seqid"], g["start"])):
        gid = output_id(prefix, g)
        if gid in ids_seen:                      # e.g. same GeneID on two loci; keep unique
            ids_seen[gid] += 1
            gid = f"{gid}.{ids_seen[gid]}"
            stats["duplicate_gene_ids"] += 1
        else:
            ids_seen[gid] = 1
        c = g["chosen"]
        is_mito = g["seqid"] in mito
        rec = {"gene_id": gid, "ncbi_gene_id": g["ncbi_gene_id"] or "",
               "annotation_gene_key": g["key"], "symbol": g["symbol"],
               "description": g["description"], "biotype": g["biotype"] or "",
               "seqid": g["seqid"], "start": g["start"], "end": g["end"], "strand": g["strand"],
               "mitochondrial": is_mito, "n_transcripts": g["n_transcripts"],
               "n_proteins": g["n_proteins"]}
        if c:
            seq, n_r, n_stop = clean_sequence(c["seq"])
            rec.update(chosen_protein=c["protein"], chosen_transcript=c["transcript"] or "",
                       protein_length=len(seq), partial=c["partial"], low_quality=c["low_quality"],
                       residues_replaced=n_r, internal_stops=n_stop, product=c["product"])
            proteins.append((gid, seq))
            lengths.append(len(seq))
            n_partial += c["partial"]
            n_lowq += c["low_quality"]
            n_mito += is_mito
            n_repl += n_r > 0
            stats["internal_stop_proteins"] += n_stop > 0
        rows.append(rec)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genes_path, faa_path = out_dir / f"{prefix}.genes.tsv", out_dir / f"{prefix}.faa"
    schemas.GENES.write(genes_path, rows, prefix=prefix)       # refuses invalid rows
    tmp = faa_path.with_suffix(".faa.part")
    with open(tmp, "w") as faa:
        for gid, seq in proteins:
            faa.write(f">{gid}\n")
            for i in range(0, len(seq), 60):
                faa.write(seq[i:i + 60] + "\n")
    tmp.replace(faa_path)
    problems = schemas.validate_proteome(faa_path, prefix, genes_path)
    if problems:
        faa_path.unlink()
        genes_path.unlink()
        raise schemas.SchemaError(f"{prefix}: proteome failed validation: " + "; ".join(problems[:5]))

    lengths.sort()
    n = len(lengths)
    summary = {
        "prefix": prefix, "species": species_name, "mode": mode,
        "protein_coding_genes": stats["protein_coding_genes"],
        "proteins_written": n,
        "coding_genes_without_protein": stats["coding_genes_without_protein"],
        "chosen_partial": n_partial, "chosen_low_quality": n_lowq,
        "mitochondrial_proteins": n_mito, "proteins_with_replaced_residues": n_repl,
        "internal_stop_proteins": stats["internal_stop_proteins"],
        "median_length": lengths[n // 2] if n else 0,
        "short_lt_50aa": sum(1 for x in lengths if x < 50),
        "cds_protein_not_in_fasta": stats["cds_protein_not_in_fasta"],
        "fasta_proteins_not_in_annotation": stats["fasta_proteins_not_in_annotation"],
        "duplicate_gene_ids": stats["duplicate_gene_ids"],
        "skipped_biotypes": "; ".join(f"{k.split(':', 1)[1]}={v}" for k, v in sorted(stats.items())
                                      if k.startswith("skipped_biotype:")),
    }
    try:                                         # NCBI's own count, from the download record
        meta = json.loads((sp / "download.json").read_text()).get("ncbi_metadata", {})
        summary["ncbi_protein_coding_genes"] = meta.get("protein_coding_genes", "")
    except (FileNotFoundError, ValueError):
        summary["ncbi_protein_coding_genes"] = ""
    return summary


def up_to_date(layout, prefix):
    """The proteome and genes table exist and are newer than the species' download."""
    faa, genes, dl = layout.proteome(prefix).local, layout.genes_table(prefix).local, layout.download(prefix).local
    if not (faa.is_file() and genes.is_file() and dl.is_dir()):
        return False
    newest_input = max((f.stat().st_mtime for f in dl.iterdir() if f.is_file()), default=0)
    return min(faa.stat().st_mtime, genes.stat().st_mtime) >= newest_input


def build(layout, manifest, species=None, mode="complete-longest", force=False, store=None, say=print):
    """Build proteomes for the manifest's species (or the given prefixes). With a store, missing
    downloads are fetched first. -> (summaries, failed prefixes)"""
    rows = read_manifest(manifest)
    if species:
        wanted = {s.strip().lower() for s in species}
        rows = {p: r for p, r in rows.items() if p.lower() in wanted}
        missing = wanted - {p.lower() for p in rows}
        for m in sorted(missing):
            say(f"[{m}] not in {manifest}; add it there to build its proteome")
    out_dir = layout.proteomes().local
    summaries, failed = [], []
    for p, row in rows.items():
        try:
            if store is not None and not layout.download(p).local.is_dir():
                store.ensure(layout.download(p))
            if up_to_date(layout, p) and not force:
                say(f"[{p}] up to date, skipping")
                continue
            s = build_species(p, row.get("annotation_format") or "gff3", layout.download(p).local, out_dir,
                              row.get("species", ""), mode)
        except Exception as e:                   # keep going; report at the end
            say(f"[{p}] FAILED: {e}")
            failed.append(p)
            continue
        summaries.append(s)
        warn = []
        if s["ncbi_protein_coding_genes"] and s["protein_coding_genes"] != s["ncbi_protein_coding_genes"]:
            warn.append(f"NCBI reports {s['ncbi_protein_coding_genes']} protein-coding genes")
        if s["protein_coding_genes"] and s["coding_genes_without_protein"] / s["protein_coding_genes"] > 0.005:
            warn.append(f"{s['coding_genes_without_protein']} coding genes have no protein")
        say(f"[{p}] {s['proteins_written']} proteins from {s['protein_coding_genes']} coding genes; "
            f"{s['chosen_partial']} partial, median {s['median_length']} aa" + ("" if not warn else "  WARNING: " + "; ".join(warn)))
    if summaries:
        merge_summary(out_dir / SUMMARY, summaries)
    return summaries, failed


def merge_summary(path, summaries):
    """Update the per-species QC table: new rows replace a species' old row, others are kept."""
    old = []
    if path.is_file():
        with open(path, newline="") as fh:
            old = list(csv.DictReader(fh, delimiter="\t"))
    new = {s["prefix"]: s for s in summaries}
    rows = [r for r in old if r.get("prefix") not in new] + list(new.values())
    cols = list(summaries[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["prefix"]))
