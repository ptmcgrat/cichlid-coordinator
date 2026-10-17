"""cichlid-data prepare compara [options]   (moved from the pipeline; defaults: the dumps in Inputs/Families/Compara116/dumps,
config/compara_panel.tsv, output Inputs/Families/Compara116/prepared)


Prepare Ensembl Compara gene trees as reference families (replaces TreeFam for --source compara).

    python compara_prepare.py \
        --nhx Data/Compara116/Compara.116.protein_default.nhx.emf.gz \
        --aln Data/Compara116/Compara.116.protein_default.aa.fasta.gz \
        --panel compara_panel.tsv --outdir Data/Compara116/prepared

Streams both dumps once. Tree n in the NHX dump pairs with alignment block n in the FASTA dump
(each block ends with '//'); every pair is checked to contain exactly the same protein IDs.
Families are numbered in dump order: C116_000001, C116_000002, ...

For each family, keeps the panel's 'reference' and 'panel' species and writes
    trees/<family>.aln     alignment of the kept proteins (all-gap columns removed)
    trees/<family>.nhx     the Compara tree pruned to the kept proteins (branch lengths kept)
in the same layout as TreeFam's files, so Stages 4 and 5 run unchanged. Also writes:
    family_audit.tsv       status per family in the Stage 3a format: ok / tiny_tree (de novo) /
                           cichlid_only (only Ensembl cichlid copies: de novo from the manifest
                           cichlids alone) / no_panel_members / mismatch
    protein_to_species.txt kept + validation proteins -> species (for Stages 6 and 6a)
    members.tsv            every kept or validation protein: gene, species, role, family, label
    search_db.fa           unaligned sequences of kept + validation proteins (Stage 3c database)
    prepare.json           release, panel and counts
'validation' species (Ensembl copies of manifest species) are removed from the trees, because a
RefSeq gene would otherwise sit next to Ensembl's copy of the same gene and look like a duplicate,
but they go into the search database so Stage 3c can record the matching Ensembl gene.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter
from operator import itemgetter
from pathlib import Path

from .. import newick  # noqa: E402

PREPARE_VERSION = 2              # bump when output rules change, so prepared trees are rebuilt
KEEP_ROLES = {"reference", "panel"}
ROLES = KEEP_ROLES | {"validation", "exclude"}
MIN_LEAVES = 4


def read_panel(path):
    roles = {}
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        role = r["role"].strip()
        if role not in ROLES:
            raise ValueError(f"panel: unknown role '{role}' for {r['species']}")
        roles[r["species"].strip()] = role
    return roles


def nhx_blocks(path):
    """-> (members [(species, protein, gene, label)], tree text) per block."""
    members, tree, in_data = [], [], False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("//"):
                yield members, "".join(tree)
                members, tree, in_data = [], [], False
            elif line.startswith("SEQ "):
                t = line.split()
                label = t[8] if len(t) > 8 and t[8] != "NULL" else ""
                members.append((t[1], t[2], t[7] if len(t) > 7 else "", label))
            elif line.startswith("DATA"):
                in_data = True
            elif in_data and line.strip():
                tree.append(line.strip())


def aln_blocks(path):
    """-> {protein: aligned sequence} per '//'-separated block."""
    block, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("//"):
                if name:
                    block[name] = "".join(buf)
                yield block
                block, name, buf = {}, None, []
            elif line.startswith(">"):
                if name:
                    block[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if block or name:
        if name:
            block[name] = "".join(buf)
        yield block


def drop_gap_columns(rows):
    """Remove columns that are gaps in every kept row (column scans run in C via zip/count)."""
    if not rows:
        return rows
    names, seqs = list(rows), list(rows.values())
    n = len(seqs)
    keep = [i for i, col in enumerate(zip(*seqs)) if col.count("-") != n]
    if len(keep) == len(seqs[0]):
        return rows
    if not keep:
        return {k: "" for k in names}
    pick = itemgetter(*keep)
    out = {}
    for k, s in zip(names, seqs):
        got = pick(s)
        out[k] = got if isinstance(got, str) else "".join(got)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nhx", required=True, type=Path)
    ap.add_argument("--aln", required=True, type=Path)
    ap.add_argument("--panel", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--release", help="Ensembl release (default: read from the file name)")
    ap.add_argument("--force", action="store_true", help="redo even if already prepared with this panel")
    args = ap.parse_args(argv)

    release = args.release or (re.search(r"Compara\.(\d+)\.", args.nhx.name) or [None, "X"])[1]
    roles = read_panel(args.panel)
    panel_md5 = hashlib.md5(args.panel.read_bytes()).hexdigest()
    inputs = {"nhx_size": args.nhx.stat().st_size, "aln_size": args.aln.stat().st_size,
              "panel_md5": panel_md5, "release": release, "prepare_version": PREPARE_VERSION}
    done = args.outdir / "prepare.json"
    if not args.force and done.is_file():
        old = json.loads(done.read_text())
        if old.get("inputs") == inputs and (args.outdir / "family_audit.tsv").is_file():
            print(f"Compara {release} already prepared with this panel; skipping (--force to redo)",
                  file=sys.stderr)
            return
    trees_dir = args.outdir / "trees"
    trees_dir.mkdir(parents=True, exist_ok=True)

    audit, n_fam, stats = [], 0, Counter()
    unknown_species = Counter()
    with open(args.outdir / "members.tsv", "w") as mem, \
         open(args.outdir / "protein_to_species.txt", "w") as p2s, \
         open(args.outdir / "search_db.fa", "w") as db:
        mem.write("protein\tgene\tspecies\trole\tfamily\tlabel\n")
        for (members, tree_text), aln in zip(nhx_blocks(args.nhx), aln_blocks(args.aln)):
            n_fam += 1
            fam = f"C{release}_{n_fam:06d}"
            rec = {"family": fam, "status": "ok", "aln_seqs": 0, "aln_length": "", "tree_leaves": 0,
                   "members_total": len(members), "members_kept": 0, "validation_members": 0,
                   "species_kept": 0, "detail": ""}
            ids = {p for _, p, _, _ in members}
            if ids != set(aln):
                rec.update(status="mismatch",
                           detail=f"{len(ids - set(aln))} tree-only, {len(set(aln) - ids)} alignment-only proteins")
                audit.append(rec)
                stats["mismatch"] += 1
                continue
            keep, sp_kept = set(), set()
            for sp, prot, gene, label in members:
                role = roles.get(sp)
                if role is None:
                    unknown_species[sp] += 1
                    role = "exclude"
                if role in KEEP_ROLES or role == "validation":
                    mem.write(f"{prot}\t{gene}\t{sp}\t{role}\t{fam}\t{label}\n")
                    p2s.write(f"{prot}\t{sp}\n")
                    db.write(f">{prot}\n{aln[prot].replace('-', '')}\n")
                if role in KEEP_ROLES:
                    keep.add(prot)
                    sp_kept.add(sp)
                elif role == "validation":
                    rec["validation_members"] += 1
            rec.update(members_kept=len(keep), species_kept=len(sp_kept))
            if not keep:
                # only Ensembl cichlid copies (validation): Stage 4 builds a de novo tree from the
                # manifest cichlids' genes alone; an empty reference alignment marks it
                rec["status"] = "cichlid_only" if rec["validation_members"] else "no_panel_members"
                if rec["status"] == "cichlid_only":
                    (trees_dir / f"{fam}.aln").write_text("")
                audit.append(rec)
                stats[rec["status"]] += 1
                continue

            rows = drop_gap_columns({p: aln[p] for p in sorted(keep)})
            with open(trees_dir / f"{fam}.aln", "w") as fh:
                for p, s in rows.items():
                    fh.write(f">{p}\n{s}\n")
            rec.update(aln_seqs=len(rows), aln_length=len(next(iter(rows.values()))))
            try:
                root = newick.prune(newick.parse(tree_text + ("" if tree_text.endswith(";") else ";")), keep)
                leaves = [l.name for l in root.leaves()]
                (trees_dir / f"{fam}.nhx").write_text(newick.write(root) + "\n")
                rec["tree_leaves"] = len(leaves)
                if set(leaves) != keep:
                    rec.update(status="mismatch", detail="pruned tree leaves differ from alignment rows")
                elif len(leaves) < MIN_LEAVES:
                    rec["status"] = "tiny_tree"
            except ValueError as e:
                rec.update(status="tree_parse_error", detail=str(e)[:100])
            stats[rec["status"]] += 1
            audit.append(rec)
            if n_fam % 5000 == 0:
                print(f"  {n_fam} trees ...", file=sys.stderr)

    cols = ["family", "status", "aln_seqs", "aln_length", "tree_leaves", "members_total", "members_kept",
            "validation_members", "species_kept", "detail"]
    with open(args.outdir / "family_audit.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(audit)
    info = {"inputs": inputs, "release": release, "families": n_fam, "status": dict(stats),
            "panel": str(args.panel), "species_not_in_panel": dict(unknown_species)}
    (args.outdir / "prepare.json").write_text(json.dumps(info, indent=2))

    print(f"\nCompara {release}: {n_fam} gene trees", file=sys.stderr)
    for k, v in stats.most_common():
        print(f"  {k:<18}{v:>7}", file=sys.stderr)
    if unknown_species:
        print(f"  {len(unknown_species)} species missing from the panel were excluded: "
              f"{', '.join(sorted(unknown_species)[:5])}", file=sys.stderr)


if __name__ == "__main__":
    main()


