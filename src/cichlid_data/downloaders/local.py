"""
Turn a published coding-sequence set into files the pipeline reads like an NCBI download, for a
species whose annotation is not in RefSeq (first use: Etroplus suratensis, Katneni et al. 2024,
Sci Data 11:1238; figshare 26303968, CC BY 4.0).

    python prepare_local_species.py --prefix Esur --cds Genes.fasta.gz --genome Genome.fa.gz \
        --functional Etroplus_annotation.gff --outdir ~/Temp/TreeFamData/Data/Esur/prepared --threads 16

Inputs: the coding sequences (one per gene, nucleotides), the genome FASTA, optionally a
functional-annotation table (Blast2GO-style GFF: one line per gene with Description= and Gene=), and
optionally the authors' protein FASTA (--proteins), which then replaces our translation. The coding
sequences are still used to place genes on the genome, and each supplied protein is compared with
our translation (printed, and saved in prepare_summary.json).
Outputs in --outdir:
    protein.faa             translations, IDs = cleaned gene IDs. Per coding sequence: frame 1 if it has
                            no internal stop; else frame 1 cut at its first stop if >= --min-fraction
                            remains (the sequence runs on past the true stop); else another frame
                            (forward 2/3 or reverse) that is clean throughout; else the gene is left
                            out (excluded.tsv, with the reason)
                            ('evm.model.Chr_1.1' -> 'Chr1-1', so browser IDs read Esur_Chr1-1)
    genomic.gff             gene / mRNA / CDS / exon features, from spliced alignment of each coding
                            sequence to the genome (minimap2 -x splice:hq); genes that do not align
                            cleanly go on a placeholder sequence 'unplaced_genes', in their original order
    sequence_report.jsonl   one line per genome sequence: Chr_N -> chromosome N (assembled-molecule),
                            anything else unplaced
    prepare_summary.json    counts: genes, aligned, unplaced, internal stops, symbols kept
Then add a species.tsv line with source=local, gff_path=<outdir>/genomic.gff,
protein_path=<outdir>/protein.faa (stage 1 copies the sequence report from the same folder).
"""
import argparse
import gzip
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
CHR_RE = re.compile(r"^(?:chr|Chr|CHR)_?(\d+)$")
UNIPROT_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9\-\.]{1,15}$")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    name, desc, seq = None, "", []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(seq)
                parts = line[1:].split(None, 1)
                name, desc, seq = parts[0], (parts[1] if len(parts) > 1 else ""), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, desc, "".join(seq)


def clean_id(raw):
    """'evm.model.Chr_1.1' -> 'Chr1-1': letters, digits and hyphens only, the characters every
    pipeline stage keeps unchanged (underscores separate the species prefix in browser IDs)."""
    s = re.sub(r"^(?:evm\.model\.|evm\.TU\.)", "", raw)
    s = re.sub(r"_", "", s)
    return re.sub(r"[^A-Za-z0-9]+", "-", s).strip("-")


def chrom_of_id(raw):
    m = re.search(r"(Chr_?\d+)", raw)
    return m.group(1) if m else ""


def translate(cds):
    cds = cds.upper()
    prot = "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - len(cds) % 3, 3))
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def protein_for(cds, min_fraction=0.75):
    """Coding sequence -> (protein or None, how). In order: frame 1 clean; frame 1 cut at its first
    stop if that keeps >= min_fraction of the protein (the sequence runs on past the true stop; a
    wrong frame meets a stop every ~20 codons, so a late first stop means frame 1 is right); another
    frame (forward 2/3, reverse 1-3) clean through the whole sequence; otherwise None."""
    s = cds.upper()
    p1 = translate(s)
    if "*" not in p1:
        return p1, "clean"
    cut = p1.index("*")
    if cut >= min_fraction * len(p1):
        return p1[:cut], "cut at first stop"
    rc = s.translate(COMP)[::-1]
    for how, seq in (("frame 2", s[1:]), ("frame 3", s[2:]), ("reverse 1", rc), ("reverse 2", rc[1:]), ("reverse 3", rc[2:])):
        p = translate(seq)
        if p and "*" not in p and len(p) >= min_fraction * len(p1):
            return p, f"reframed ({how})"
    return None, "first stop too early, no clean frame"


def protein_keys(name, desc=""):
    """IDs a protein header may refer to its gene by: the first word, the part after its last '|'
    ('protein|evm.model.Chr_1.1' -> 'evm.model.Chr_1.1'), and an ID= field in the description."""
    keys = [name]
    if "|" in name:
        keys.append(name.rsplit("|", 1)[-1])
    m = re.search(r"(?:^|[\s|;])ID=([^\s|;]+)", desc or "")
    if m:
        keys.append(m.group(1))
    return [k for k in keys if k]


def clean_description(d):
    d = re.sub(r"^PREDICTED:\s*", "", d or "").strip()
    d = re.sub(r"\s*\[[^\]]+\]\s*$", "", d).strip()
    return "" if d in ("---NA---", "NA") else d


def read_functional(path):
    """Blast2GO-style GFF -> {raw id: (description, symbol)}"""
    out = {}
    if not path or not Path(path).is_file():
        return out
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in t[8].split(";") if "=" in kv)
            sym = (attrs.get("Gene") or "").strip()
            if sym and (UNIPROT_RE.match(sym) or "_" in sym or not SYMBOL_RE.match(sym)):
                sym = ""
            out[t[0]] = (clean_description(attrs.get("Description", "")), sym)
    return out


def exons_from_cg(tstart, cg):
    """PAF target start (0-based) and cg:Z CIGAR -> list of (start, end) exons, 1-based inclusive."""
    exons, pos, cur = [], tstart, None
    for n, op in re.findall(r"(\d+)([MIDN=X])", cg):
        n = int(n)
        if op in "M=XD":
            if cur is None:
                cur = pos
            pos += n
        elif op == "N":
            if cur is not None:
                exons.append((cur + 1, pos))
                cur = None
            pos += n
        # I: insertion in the query, no genome span
    if cur is not None:
        exons.append((cur + 1, pos))
    return exons


def run_minimap2(minimap2, genome, cds_fa, threads, paf):
    cmd = [minimap2, "-x", "splice:hq", "-uf", "--secondary=no", "-c", "-t", str(threads), str(genome), str(cds_fa)]
    with open(paf, "w") as out:
        r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise SystemExit(f"minimap2 failed: {r.stderr[-800:]}")


def best_alignments(paf, min_cov=0.9, min_ident=0.95):
    """PAF -> {query: (target, strand, exons)} for one clean primary alignment per query."""
    best = {}
    with open(paf) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            q, qlen, qs, qe, strand, tname = t[0], int(t[1]), int(t[2]), int(t[3]), t[4], t[5]
            matches, block = int(t[9]), int(t[10])
            tags = dict(x.split(":", 2)[::2] for x in t[12:] if x.count(":") >= 2)
            if tags.get("tp", "P") != "P" or "cg" not in tags:
                continue
            cov, ident = (qe - qs) / qlen if qlen else 0, matches / block if block else 0
            if cov < min_cov or ident < min_ident:
                continue
            score = matches
            if q not in best or score > best[q][0]:
                best[q] = (score, tname, strand, exons_from_cg(int(t[7]), tags["cg"]))
    return {q: v[1:] for q, v in best.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--cds", required=True, type=Path)
    ap.add_argument("--genome", type=Path, help="genome FASTA; without it every gene goes on 'unplaced_genes'")
    ap.add_argument("--functional", type=Path, help="Blast2GO-style table for descriptions and symbols")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--paf", type=Path, help="reuse an existing alignment instead of running minimap2")
    ap.add_argument("--proteins", type=Path,
                    help="protein FASTA from the authors: used instead of translating the coding sequences")
    ap.add_argument("--min-fraction", type=float, default=0.75,
                    help="keep a protein cut at its first stop only if this much of it remains")
    args = ap.parse_args(argv)
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    func = read_functional(args.functional)

    given = {}
    if args.proteins:
        n_read = 0
        for name, desc, p in read_fasta(args.proteins):
            p = re.sub(r"\s", "", p).rstrip("*")
            n_read += 1
            for key in protein_keys(name, desc):
                given.setdefault(key, p)
                given.setdefault(clean_id(key), p)
        print(f"Read {n_read} proteins from {args.proteins.name}", file=sys.stderr)
    compare = Counter()
    genes, stops, hows, excluded = [], 0, Counter(), []
    with open(out / "protein.faa", "w") as pf, open(out / "cds.fa", "w") as cf:
        for raw, desc, seq in read_fasta(args.cds):
            gid = clean_id(raw)
            ours = translate(seq)
            if "*" in ours:
                stops += 1
            if args.proteins:
                prot = given.get(raw) or given.get(gid)
                if prot is None:
                    how = "not in the protein file"
                elif "*" in prot:
                    prot, how = None, "internal stop in the supplied protein"
                else:
                    how = "supplied protein"
                    if "*" not in ours:
                        rel = "ours clean: " + ("identical" if prot == ours else "different")
                    else:
                        cut = ours.split("*")[0]
                        rel = "ours had stops: " + ("matches ours cut at first stop" if prot == cut else
                                                    "same start as ours" if prot[:30] == ours[:30] else "different")
                    compare[rel] += 1
            else:
                prot, how = protein_for(seq, args.min_fraction)
            hows[how] += 1
            d, sym = func.get(raw, (clean_description(desc), ""))
            if prot is None:
                excluded.append((gid, raw, how, len(seq)))
                continue
            genes.append({"raw": raw, "id": gid, "desc": d, "symbol": sym, "len": len(seq), "how": how})
            pf.write(f">{gid} {d}\n" + "\n".join(prot[i:i + 60] for i in range(0, len(prot), 60)) + "\n")
            cf.write(f">{gid}\n{seq}\n")

    aligned, seq_lengths = {}, {}
    if args.genome:
        for name, _, s in read_fasta(args.genome):
            seq_lengths[name] = len(s)
        paf = args.paf or out / "cds_to_genome.paf"
        if not args.paf:
            run_minimap2(args.minimap2, args.genome, out / "cds.fa", args.threads, paf)
        aligned = best_alignments(paf)

    unplaced_pos = 1
    with open(out / "genomic.gff", "w") as gf:
        gf.write("##gff-version 3\n")
        for g in genes:
            a = aligned.get(g["id"])
            attrs_g = f"ID=gene-{g['id']};Name={g['symbol'] or g['id']};gene_biotype=protein_coding"
            if g["symbol"]:
                attrs_g += f";gene={g['symbol']}"
            if g["desc"]:
                attrs_g += f";description={g['desc'].replace(';', ',').replace('=', ' ')}"
            attrs_g += f";Note=source id {g['raw']}"
            if a:
                seqid, strand, exons = a
            else:
                seqid, strand = "unplaced_genes", "+"
                exons = [(unplaced_pos, unplaced_pos + g["len"] - 1)]
                unplaced_pos += g["len"] + 1000
            start, end = min(e[0] for e in exons), max(e[1] for e in exons)
            gf.write(f"{seqid}\tlocal\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs_g}\n")
            gf.write(f"{seqid}\tlocal\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID=rna-{g['id']};Parent=gene-{g['id']}\n")
            for k, (s, e) in enumerate(exons, 1):
                gf.write(f"{seqid}\tlocal\texon\t{s}\t{e}\t.\t{strand}\t.\tID=exon-{g['id']}-{k};Parent=rna-{g['id']}\n")
                gf.write(f"{seqid}\tlocal\tCDS\t{s}\t{e}\t.\t{strand}\t.\tID=cds-{g['id']};Parent=rna-{g['id']};"
                         f"protein_id={g['id']}\n")

    with open(out / "excluded.tsv", "w") as fh:
        fh.write("gene_id\tsource_id\treason\tcds_length\n")
        for row in excluded:
            fh.write("\t".join(map(str, row)) + "\n")
    n_unplaced = sum(1 for g in genes if g["id"] not in aligned)
    with open(out / "sequence_report.jsonl", "w") as fh:
        for name, length in list(seq_lengths.items()) + ([("unplaced_genes", unplaced_pos)] if n_unplaced else []):
            m = CHR_RE.match(name)
            rec = {"refseqAccession": name, "genbankAccession": name, "sequenceName": name, "length": length,
                   "role": "assembled-molecule" if m else "unplaced-scaffold"}
            if m:
                rec["chrName"] = m.group(1)
            fh.write(json.dumps(rec) + "\n")
    off_chr_ids = sum(1 for g in genes if not chrom_of_id(g["raw"]))
    summary = {"prefix": args.prefix, "genes": len(genes), "aligned": sum(1 for g in genes if g["id"] in aligned),
               "unplaced": n_unplaced, "excluded": len(excluded), "translation": dict(hows),
               "cds_with_internal_stops_in_frame_1": stops, "symbols": sum(1 for g in genes if g["symbol"]),
               "ids_off_chromosomes": off_chr_ids, "genome_sequences": len(seq_lengths)}
    (out / "prepare_summary.json").write_text(json.dumps(summary, indent=2))
    (out / "cds.fa").unlink()
    print(f"{args.prefix}: {len(genes) + len(excluded)} coding sequences -> " + "; ".join(f"{k}: {v}" for k, v in hows.most_common()),
          file=sys.stderr)
    if compare:
        print("  supplied proteins against our translation of the coding sequences:", file=sys.stderr)
        for k, v in sorted(compare.items()):
            print(f"    {v:>7}  {k}", file=sys.stderr)
        summary["comparison_with_cds_translation"] = dict(compare)
        (out / "prepare_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  kept {len(genes)} genes ({summary['aligned']} placed on the genome, {n_unplaced} on 'unplaced_genes'), "
          f"excluded {len(excluded)} (see excluded.tsv); {summary['symbols']} gene symbols kept -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()


