"""cichlid-data summarize h5ad FILE... --dataset NAME [options]   (moved from the pipeline; needs the [expression] extra:
numpy and scipy; default folder Inputs/Expression/<dataset>)


Summarize an annotated single-cell atlas (.h5ad) per cell type, for the gene browser's
Expression tab. Works for any atlas (Zebrahub, ZMAP, Daniocell conversions, ...).

    python summarize_h5ad.py zebrahub.h5ad --inspect          # list annotation columns first
    python summarize_h5ad.py zebrahub.h5ad --dataset Zebrahub \
        --celltype-col zebrafish_anatomy_ontology_class --stage-col timepoint
    python summarize_h5ad.py zf_atlas_*_v1_release.h5ad --dataset Zebrahub --celltype-col <col>
        # several files, e.g. one per stage: the stage comes from each file name
        # (zf_atlas_2dpf_v1_release.h5ad -> 2dpf) unless --stage-col is given, and genes are
        # lined up across files by Ensembl ID (or symbol)

For each gene and cell type:
    mean_cp10k       average expression in counts per 10,000 (each cell's counts scaled to 10,000)
    pct_expressing   share of the cell type's cells with any counts for the gene
(the two numbers of a dot plot). Rows where no cell expresses the gene are left out.
Raw counts are used: the --layer if given, else a 'counts' layer if present, else X when it
holds whole numbers; an X that is already log-normalized is turned back with expm1.

Writes to <local master>/Data/Expression/<dataset>/ (or --outdir):
    expression.tsv.gz   gene_id, symbol, cell_type, mean_cp10k, pct_expressing
    cell_types.tsv      cell_type, n_cells, stages (stages with >= 10 cells), zfa_id, system
    summary.json        source file, columns used, counts, date
  with --by-stage (needs a stage, from --stage-col or file names):
    expression_by_stage.tsv.gz   gene_id, symbol, cell_type, stage, mean_cp10k, pct_expressing
    cell_types_by_stage.tsv      cell_type, stage, n_cells
  with --zfa-col (a column of Zebrafish Anatomy Ontology IDs, e.g. zebrafish_anatomy_ontology_id):
    each cell type is assigned to an organ system by walking up the ontology (is_a and part_of),
    used to group rows in the browser's dot plot; the ontology is downloaded once (--zfa-obo to
    give a local copy).
then backs the folder up to Dropbox (--no-backup to skip). Stage 9 loads every dataset folder
under Data/Expression that has these files.
"""
import argparse
import gzip
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

from collections import defaultdict

import numpy as np
from scipy import sparse


ENSEMBL_GENE = re.compile(r"^ENS[A-Z]*G\d{11}")


# --------------------------------------------------------------------------- core (tested)
def aggregate(chunks, celltypes, stages=None, already_log=False, min_stage_cells=10):
    """chunks: iterable of (start, CSR matrix cells x genes) covering all cells in order.
    celltypes: array of labels per cell (None/'' = unlabelled, skipped).
    -> (types, mean_cp10k [types x genes], pct [types x genes], n_cells per type, stages per type)"""
    labels = np.asarray([("" if c is None else str(c)) for c in celltypes], dtype=object)
    types = sorted({c for c in labels if c and c.lower() != "nan"})
    tindex = {t: i for i, t in enumerate(types)}
    code = np.array([tindex.get(c, -1) for c in labels])
    n_types = len(types)
    sums = expr = None
    for start, block in chunks:
        X = sparse.csr_matrix(block, dtype=np.float64)
        if already_log:
            X.data = np.expm1(X.data)
        if sums is None:
            sums = np.zeros((n_types, X.shape[1]))
            expr = np.zeros((n_types, X.shape[1]))
        c = code[start:start + X.shape[0]]
        keep = c >= 0
        if not keep.any():
            continue
        X = X[keep]
        lib = np.asarray(X.sum(axis=1)).ravel()
        lib[lib == 0] = 1
        Xn = sparse.diags(1e4 / lib) @ X
        G = sparse.csr_matrix((np.ones(keep.sum()), (c[keep], np.arange(keep.sum()))), shape=(n_types, keep.sum()))
        sums += np.asarray((G @ Xn).todense())
        expr += np.asarray((G @ (X > 0).astype(np.float64)).todense())
    n_cells = np.bincount(code[code >= 0], minlength=n_types).astype(float)
    denom = np.maximum(n_cells, 1)[:, None]
    stage_lists = [[] for _ in types]
    if stages is not None:
        st = np.asarray([("" if s is None else str(s)) for s in stages], dtype=object)
        for i, t in enumerate(types):
            vals, counts = np.unique(st[code == i], return_counts=True)
            stage_lists[i] = [v for v, k in zip(vals, counts) if v and k >= min_stage_cells]
    if sums is None:
        sums = expr = np.zeros((n_types, 0))
    return types, sums / denom, expr / denom, n_cells.astype(int), stage_lists


def write_outputs(outdir, genes, symbols, types, mean, pct, n_cells, stage_lists, info=None, by_stage=None,
                  systems=None):
    outdir.mkdir(parents=True, exist_ok=True)
    rows = 0
    with gzip.open(outdir / "expression.tsv.gz", "wt") as fh:
        fh.write("gene_id\tsymbol\tcell_type\tmean_cp10k\tpct_expressing\n")
        for j, (g, s) in enumerate(zip(genes, symbols)):
            for i, t in enumerate(types):
                if pct[i, j] > 0:
                    fh.write(f"{g}\t{s}\t{t}\t{mean[i, j]:.4g}\t{100 * pct[i, j]:.3g}\n")
                    rows += 1
    systems = systems or {}
    with open(outdir / "cell_types.tsv", "w") as fh:
        fh.write("cell_type\tn_cells\tstages\tzfa_id\tsystem\n")
        for t, n, st in zip(types, n_cells, stage_lists):
            z, sy = systems.get(t, ("", ""))
            fh.write(f"{t}\t{n}\t{', '.join(st)}\t{z}\t{sy}\n")
    if by_stage:
        pairs, m2, p2, n2 = by_stage
        with gzip.open(outdir / "expression_by_stage.tsv.gz", "wt") as fh:
            fh.write("gene_id\tsymbol\tcell_type\tstage\tmean_cp10k\tpct_expressing\n")
            for j, (g, s) in enumerate(zip(genes, symbols)):
                for i, (ct, st) in enumerate(pairs):
                    if p2[i, j] > 0:
                        fh.write(f"{g}\t{s}\t{ct}\t{st}\t{m2[i, j]:.4g}\t{100 * p2[i, j]:.3g}\n")
        with open(outdir / "cell_types_by_stage.tsv", "w") as fh:
            fh.write("cell_type\tstage\tn_cells\n")
            for (ct, st), n in zip(pairs, n2):
                fh.write(f"{ct}\t{st}\t{n}\n")
    info = {**info, "rows": rows, "cell_types": len(types), "genes": len(genes),
            "made": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (outdir / "summary.json").write_text(json.dumps(info, indent=2))
    return info


def looks_ensembl(values):
    sample = [str(x) for x in list(values)[:200]]
    return bool(sample) and sum(bool(ENSEMBL_GENE.match(v)) for v in sample) > len(sample) / 2


ZFA_URL = "http://purl.obolibrary.org/obo/zfa.obo"
# organ systems for grouping rows, matched against a term's own name and its ancestors' names, in order
SYSTEMS = [("Eye", ("visual system", "eye", "retina", "lens")),
           ("Ear and lateral line", ("auditory system", "lateral line", "ear", "otic", "neuromast")),
           ("Nervous system", ("nervous system", "neural", "neuron", "brain", "spinal cord", "glial")),
           ("Neural crest and pigment", ("neural crest", "pigment", "melanocyte", "iridophore", "xanthophore")),
           ("Blood and vasculature", ("cardiovascular system", "hematopoietic", "blood", "vascular", "heart")),
           ("Immune system", ("immune system", "macrophage", "leukocyte")),
           ("Muscle", ("musculature system", "muscle", "myotome", "somite")),
           ("Skeleton and fins", ("skeletal system", "cartilage", "bone", "fin", "notochord")),
           ("Digestive and endoderm", ("digestive system", "endoderm", "liver", "pancrea", "intestin", "gut")),
           ("Kidney", ("renal system", "pronephr", "kidney")),
           ("Skin", ("integument", "epiderm", "periderm", "skin", "ionocyte", "mucous")),
           ("Germ line", ("reproductive system", "germ", "gonad")),
           ("Mesoderm and progenitors", ("mesoderm", "mesenchyme", "progenitor", "tail bud", "segmental plate")),
           ("Ectoderm", ("ectoderm",))]


def parse_obo(text):
    """-> ({id: name}, {id: [parent ids via is_a and part_of]})"""
    names, parents, cur = {}, defaultdict(list), None
    for line in text.splitlines():
        line = line.strip()
        if line == "[Term]":
            cur = None
        elif line.startswith("id: "):
            cur = line[4:].strip()
        elif cur and line.startswith("name: "):
            names[cur] = line[6:].strip()
        elif cur and line.startswith("is_a: "):
            parents[cur].append(line[6:].split("!")[0].strip())
        elif cur and line.startswith("relationship: part_of "):
            parents[cur].append(line[len("relationship: part_of "):].split("!")[0].strip())
        elif line.startswith("[") and line != "[Term]":
            cur = None
    return names, parents


def system_of(term, names, parents):
    """Organ system of an ontology term: the first SYSTEMS entry matching the term or an ancestor
    (nearest first, so a 'retinal neuron' goes to Eye, not Nervous system)."""
    seen, frontier, order = {term}, [term], [term]
    while frontier:
        nxt = []
        for t in frontier:
            for p in parents.get(t, []):
                if p not in seen:
                    seen.add(p)
                    nxt.append(p)
                    order.append(p)
        frontier = nxt
    labels = [names.get(t, "").lower() for t in order]
    for lab in labels:
        for system, keys in SYSTEMS:
            if any(k in lab for k in keys):
                return system
    return "Other"


STAGE_RE = re.compile(r"\s*(\d+(?:\.\d+)?)(?:\s*[-–]\s*\d+(?:\.\d+)?)?\s*(hpf|dpf|h|d)\b", re.I)


BARE_STAGE = re.compile(r"\s*\d+(?:\.\d+)?(?:\s*[-–]\s*\d+(?:\.\d+)?)?\s*")


def stage_label(stage):
    """Bare numbers or ranges (Daniocell's stage groups, e.g. '  3-4') are hours post fertilization:
    '  3-4' -> '3-4 hpf'. Other labels are only trimmed."""
    s = str(stage).strip()
    return f"{s} hpf" if s and BARE_STAGE.fullmatch(s) else s


def celltype_values(obs, spec):
    """--celltype-col: one column, or 'a+b' to label cells 'a (b)' (b omitted when empty or equal to a),
    e.g. Daniocell's identity.super + identity.sub."""
    cols = spec.split("+")
    if len(cols) == 1:
        return list(obs[spec])
    a, b = obs[cols[0]], obs[cols[1]]
    out = []
    for x, y in zip(a, b):
        x, y = str(x).strip(), str(y).strip()
        bad = lambda v: v in ("", "nan", "NA", "None")
        if bad(x):
            out.append("")
        elif bad(y) or y == x:
            out.append(x)
        else:
            out.append(f"{x} ({y})")
    return out


def stage_hours(stage):
    """'10hpf' -> 10, '2dpf' -> 48, '24-34hpf' -> 24 (a range by its start), '5 somites' -> None."""
    m = STAGE_RE.match(str(stage))
    if not m:
        return None
    v = float(m.group(1))
    return v * 24 if m.group(2).lower().startswith("d") else v


def gene_columns(var_names, var):
    """-> (gene_ids, symbols): Ensembl IDs where the atlas has them, and symbols."""
    names = [str(x) for x in var_names]
    if looks_ensembl(names):
        ids = names
        sym_col = next((c for c in var.columns if c.lower() in ("gene_symbols", "gene_symbol", "symbol",
                                                                 "feature_name", "gene_name", "gene_names")), None)
        syms = [str(x) for x in var[sym_col]] if sym_col else [""] * len(ids)
        return ids, syms
    id_col = next((c for c in var.columns if looks_ensembl(var[c])), None)
    return ([str(x) for x in var[id_col]] if id_col else [""] * len(names)), names


# --------------------------------------------------------------------------- input access (thin)
class SeuratExport:
    """A directory written by export_seurat.R, presented like an AnnData object (cells x genes)."""

    def __init__(self, d, annotations=None, join=None):
        import pandas as pd
        d = Path(d)
        n_genes, n_cells, nnz = (int(x) for x in (d / "dims.txt").read_text().split())
        x = np.fromfile(d / "counts_x.bin", dtype=np.float64, count=nnz)
        i = np.fromfile(d / "counts_i.bin", dtype=np.int32, count=nnz)
        p = np.fromfile(d / "counts_p.bin", dtype=np.int32, count=n_cells + 1)
        # genes x cells in compressed columns == cells x genes in compressed rows
        self.X = sparse.csr_matrix((x, i, p), shape=(n_cells, n_genes))
        self.n_obs, self.n_vars = n_cells, n_genes
        genes = (d / "genes.txt").read_text().splitlines()
        self.var = pd.DataFrame(index=pd.Index(genes))
        self.var_names = self.var.index
        meta = pd.read_csv(d / "meta.tsv", sep="\t", dtype=str, keep_default_na=False)
        cells = (d / "cells.txt").read_text().splitlines()
        meta = meta.set_index("cell").reindex(cells).reset_index()
        if annotations is not None:
            ann = pd.read_csv(annotations, dtype=str, keep_default_na=False)
            left, right = join.split("=", 1) if join and "=" in join else (join, join)
            if not left or left not in meta.columns or right not in ann.columns:
                raise SystemExit(f"--join {join!r}: need a cell metadata column ({', '.join(meta.columns)}) "
                                 f"and an annotation column ({', '.join(ann.columns)})")
            ann = ann.drop_duplicates(right).rename(columns={c: f"ann_{c}" for c in ann.columns if c != right})
            meta = meta.merge(ann, how="left", left_on=left, right_on=right, suffixes=("", "_ann"))
        self.obs = meta
        self.layers = {}
        src = (d / "source.txt").read_text().split() if (d / "source.txt").is_file() else []
        self.log_normalized = len(src) > 1 and src[1] == "data"


def open_h5ad(path, annotations=None, join=None):
    if Path(path).is_dir() and (Path(path) / "dims.txt").is_file():
        return SeuratExport(path, annotations, join)
    import anndata
    return anndata.read_h5ad(path, backed="r")


def inspect(ad):
    print(f"{ad.n_obs} cells x {ad.n_vars} genes; layers: {list(ad.layers.keys()) if ad.layers is not None else []}")
    print("Cell annotation columns (obs): name, number of distinct values, examples")
    for c in ad.obs.columns:
        v = ad.obs[c]
        u = v.unique()
        print(f"  {c:<38} {len(u):>7}   {', '.join(map(str, list(u)[:4]))[:80]}")
    print("Gene columns (var):", ", ".join(map(str, ad.var.columns)) or "(none)", "| index example:", ad.var_names[0])


def matrix_chunks(ad, layer, chunk=20000):
    X = ad.layers[layer] if layer else ad.X
    for start in range(0, ad.n_obs, chunk):
        block = X[start:min(start + chunk, ad.n_obs)]
        yield start, (block if sparse.issparse(block) else sparse.csr_matrix(np.asarray(block)))


STAGE_IN_NAME = re.compile(r"(\d+(?:\.\d+)?\s*(?:hpf|dpf|somites?|ss))", re.I)


def stage_from_name(path):
    m = STAGE_IN_NAME.search(Path(path).name)
    return m.group(1) if m else Path(path).name.split(".")[0]


def raw_chunks(ad, layer, chunk=20000):
    """Chunks of raw counts for one file (expm1 applied if X is log-normalized)."""
    first = next(matrix_chunks(ad, layer, chunk=500))[1]
    sample = first.data[:10000] if first.nnz else np.array([0.0])
    is_log = getattr(ad, "log_normalized", False) or (layer is None and not np.allclose(sample, np.round(sample)))
    def gen():
        for start, block in matrix_chunks(ad, layer, chunk):
            b = sparse.csr_matrix(block, dtype=np.float64)
            if is_log:
                b.data = np.expm1(b.data)
            yield start, b
    return is_log, gen()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad", nargs="+", help="one or more .h5ad files (e.g. one per stage)")
    ap.add_argument("--inspect", action="store_true", help="list annotation columns and stop")
    ap.add_argument("--dataset", help="short name, e.g. Zebrahub (becomes the folder name)")
    ap.add_argument("--celltype-col", help="column with cell type labels, or 'a+b' for labels 'a (b)'")
    ap.add_argument("--stage-col")
    ap.add_argument("--layer", help="layer with raw counts (default: 'counts' if present, else X)")
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--by-stage", action="store_true", help="also summarize each cell type at each stage")
    ap.add_argument("--zfa-col", help="column of Zebrafish Anatomy Ontology IDs, to group cell types by system")
    ap.add_argument("--zfa-obo", type=Path, help="local copy of zfa.obo (default: download)")
    ap.add_argument("--annotations", type=Path, help="cluster annotation CSV to join onto cells (e.g. Daniocell's)")
    ap.add_argument("--join", help="META_COL=ANNOTATION_COL for --annotations (annotation columns get an ann_ prefix)")
    ap.add_argument("--system-col", help="column giving each cell type's tissue or system (used instead of --zfa-col)")
    args = ap.parse_args(argv)

    paths = args.h5ad
    if args.inspect:
        for p in paths:
            print(f"== {Path(p).name} (stage from file name: {stage_from_name(p)})")
            inspect(open_h5ad(p, args.annotations, args.join))
        return
    if not (args.dataset and args.celltype_col):
        sys.exit("give --dataset and --celltype-col (run with --inspect to see the columns)")

    # pass 1: gene lists of every file -> one combined gene list
    ads, keys_per_file, gene_ids, symbols, index = [], [], [], [], {}
    for p in paths:
        ad = open_h5ad(p, args.annotations, args.join)
        if any(c not in ad.obs.columns for c in args.celltype_col.split("+")):
            sys.exit(f"{Path(p).name}: no column {args.celltype_col}; run --inspect")
        ids, syms = gene_columns(ad.var_names, ad.var)
        keys = [i or s_ for i, s_ in zip(ids, syms)]
        for k, i, s_ in zip(keys, ids, syms):
            if k not in index:
                index[k] = len(gene_ids)
                gene_ids.append(i)
                symbols.append(s_)
        ads.append(ad)
        keys_per_file.append(keys)
    n_genes = len(gene_ids)

    # pass 2: counts from every file, re-indexed onto the combined gene list
    celltypes, stages, zfa, tissue, any_log, layers = [], [], [], [], False, {}
    for p, ad in zip(paths, ads):
        celltypes += celltype_values(ad.obs, args.celltype_col)
        if args.system_col and args.system_col in ad.obs.columns:
            tissue += list(ad.obs[args.system_col])
        stages += [stage_label(x) for x in ad.obs[args.stage_col]] if args.stage_col else [stage_from_name(p)] * ad.n_obs
        if args.zfa_col and args.zfa_col in ad.obs.columns:
            zfa += list(ad.obs[args.zfa_col])

    def all_chunks():
        nonlocal any_log
        offset = 0
        for p, ad, keys in zip(paths, ads, keys_per_file):
            layer = args.layer or ("counts" if ad.layers is not None and "counts" in ad.layers.keys() else None)
            layers[Path(p).name] = layer
            cols = np.array([index[k] for k in keys])
            P = sparse.csr_matrix((np.ones(len(cols)), (np.arange(len(cols)), cols)), shape=(len(cols), n_genes))
            is_log, gen = raw_chunks(ad, layer)
            if is_log:
                any_log = True
                print(f"{Path(p).name}: X is not whole-number counts; treated as log1p-normalized", file=sys.stderr)
            for start, block in gen:
                yield offset + start, block @ P
            offset += ad.n_obs

    print(f"Summarizing {len(celltypes)} cells from {len(paths)} file(s), {n_genes} genes, by {args.celltype_col} ...",
          file=sys.stderr)
    types, mean, pct, n_cells, stage_lists = aggregate(all_chunks(), celltypes, stages)
    genes = gene_ids
    by_stage = None
    if args.by_stage:
        SEP = "\x1f"
        combo = [f"{c}{SEP}{s_}" if c not in (None, "") and str(c).lower() != "nan" else "" for c, s_ in zip(celltypes, stages)]
        print("Summarizing each cell type at each stage ...", file=sys.stderr)
        t2, m2, p2, n2, _ = aggregate(all_chunks(), combo)
        by_stage = ([tuple(t.split(SEP, 1)) for t in t2], m2, p2, n2)
    systems = {}
    if tissue:
        for c, t in zip(celltypes, tissue):
            if c not in (None, "") and str(c).lower() != "nan" and t not in (None, ""):
                systems.setdefault(str(c), ("", str(t)))
        print(f"Tissues from {args.system_col} assigned to {len(systems)} cell types", file=sys.stderr)
    elif zfa:
        id_of = {}
        for c, z in zip(celltypes, zfa):
            if c not in (None, "") and z not in (None, "") and str(z).lower() != "nan":
                id_of.setdefault(str(c), str(z))
        try:
            text = args.zfa_obo.read_text() if args.zfa_obo else __import__("requests").get(ZFA_URL, timeout=120).text
            names, parents = parse_obo(text)
            systems = {c: (z, system_of(z, names, parents)) for c, z in id_of.items()}
            print(f"Organ systems assigned to {len(systems)} cell types from the anatomy ontology", file=sys.stderr)
        except Exception as e:                     # rows are just ungrouped then
            print(f"Could not read the anatomy ontology ({e}); cell types will not be grouped", file=sys.stderr)
    from ..backup import backup
    from ..layout import get_layout
    out = args.outdir or get_layout().expression(args.dataset).local
    info = write_outputs(out, genes, symbols, types, mean, pct, n_cells, stage_lists,
                         by_stage=by_stage, systems=systems, info=
                         dict(**{"dataset": args.dataset,
                          "source_file": ", ".join(Path(p).name for p in paths) if len(paths) <= 3
                          else f"{len(paths)} files ({Path(paths[0]).name} ...)",
                          "celltype_col": args.celltype_col,
                          "stage_col": args.stage_col or "from file names", "layers": layers,
                          "log_normalized_input": bool(any_log), "by_stage": bool(by_stage)}))
    print(f"{info['cell_types']} cell types, {info['genes']} genes, {info['rows']} expressed gene x cell type rows "
          f"-> {out}", file=sys.stderr)
    if not args.no_backup:
        print(f"Backed up to {backup(out)}", file=sys.stderr)


if __name__ == "__main__":
    main()


