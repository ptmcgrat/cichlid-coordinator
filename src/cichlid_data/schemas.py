"""Formats of the files that cross repository boundaries: readers, writers and validators.

    from cichlid_data import schemas
    rows = schemas.GENES.read(path)                 # list of dicts (strings), header checked
    schemas.GENES.write(path, rows)                 # columns in order; refuses invalid rows
    problems = schemas.validate(path)               # the kind is inferred from the file name

Kinds: genes (<prefix>.genes.tsv), assign (<prefix>.assign.tsv), orthologs (<prefix>.orthologs.tsv),
proteome (<prefix>.faa), run (run.json), residue_map (residue_map.tsv), family (a families/<family>/ folder),
families (a whole families/ folder: every family in it).

Family IDs: letters, digits, '.', '_' and '-', starting with a letter or digit (TF101001, C116_012345,
OG0001234, N0.HOG0001234). The family contract (any method's Results/<run>/families/<family>/) is in
validate_family below and in docs/coordinator_spec.md.

Gene IDs: <prefix>_<id>; the prefix is letters and digits (starting with a letter), the id letters,
digits, '.' and '-' (no underscores: the first underscore separates the prefix). Tables are
tab-separated with a header row; empty strings mean missing; booleans are True or False.
"""
import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import newick
from .layout import RELEASE_RE, RESULT_PARTS, RESULTS

PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
GENE_ID_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*)_([A-Za-z0-9.\-]+)$")
FAMILY_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")      # not '_' first: "_top" is a shard name
ORTHOFINDER_INTERNAL_RE = re.compile(r"^\d+_\d+$")                # OrthoFinder's own sequence names (0_123)
ALIGNMENT_FILE = "combined.aln"          # aligned FASTA; the placement methods write it beside queries.fa, ref.tree, ...
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYX")
MAX_PROBLEMS = 25


class SchemaError(ValueError):
    pass


@dataclass(frozen=True)
class Column:
    name: str
    kind: str = "str"            # str | int | float | bool | enum | id | family
    required: bool = False       # must be non-empty
    values: tuple = ()           # allowed values, for enum

    def check(self, v):
        """-> error message, or None."""
        if v == "" or v is None:
            return f"{self.name} is empty" if self.required else None
        v = str(v)
        if self.kind == "int":
            try:
                int(v)
            except ValueError:
                return f"{self.name} is not an integer: {v!r}"
        elif self.kind == "float":
            try:
                float(v)
            except ValueError:
                return f"{self.name} is not a number: {v!r}"
        elif self.kind == "bool" and v not in ("True", "False"):
            return f"{self.name} is not True/False: {v!r}"
        elif self.kind == "enum" and v not in self.values:
            return f"{self.name} {v!r} is not one of {', '.join(self.values)}"
        elif self.kind == "id" and not GENE_ID_RE.match(v):
            return f"{self.name} {v!r} is not a gene ID (<prefix>_<id>; id without underscores)"
        elif self.kind == "family" and not FAMILY_ID_RE.match(v):
            return f"{self.name} {v!r} is not a family ID (letters, digits, '.', '_' and '-')"
        return None


def valid_family_id(name):
    return bool(FAMILY_ID_RE.match(str(name)))


def gene_prefix(gene_id):
    m = GENE_ID_RE.match(gene_id)
    return m.group(1) if m else None


@dataclass
class TableSchema:
    kind: str
    columns: list
    key: str = "gene_id"
    file_suffix: str = ""
    doc: str = ""

    @property
    def names(self):
        return [c.name for c in self.columns]

    def check_header(self, header):
        if header != self.names:
            missing = [c for c in self.names if c not in header]
            extra = [c for c in header if c not in self.names]
            msg = []
            if missing:
                msg.append(f"missing columns: {', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")
            if extra:
                msg.append(f"unexpected columns: {', '.join(extra[:8])}{' ...' if len(extra) > 8 else ''}")
            if not msg:
                msg.append("columns are in the wrong order")
            return ["header: " + "; ".join(msg)]
        return []

    def check_rows(self, rows, prefix=None):
        problems, seen = [], set()
        cols = self.columns
        for i, r in enumerate(rows, start=2):                  # line numbers, after the header
            for c in cols:
                err = c.check(r.get(c.name, ""))
                if err:
                    problems.append(f"line {i}: {err}")
            k = r.get(self.key, "")
            if k in seen:
                problems.append(f"line {i}: duplicate {self.key} {k!r}")
            seen.add(k)
            if prefix and k and gene_prefix(k) not in (None, prefix):
                problems.append(f"line {i}: {k!r} does not start with the file's prefix {prefix}_")
            if len(problems) >= MAX_PROBLEMS:
                problems.append("(stopping after too many problems)")
                break
        return problems

    def read(self, path):
        with open(path, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            problems = self.check_header(rd.fieldnames or [])
            if problems:
                raise SchemaError(f"{path}: {problems[0]}")
            return list(rd)

    def validate_file(self, path, prefix=None):
        with open(path, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            problems = self.check_header(rd.fieldnames or [])
            if problems:
                return problems
            return self.check_rows(rd, prefix or _prefix_from_name(path))

    def write(self, path, rows, prefix=None):
        rows = [{c: ("" if r.get(c) is None else r.get(c)) for c in self.names} for r in rows]
        for r in rows:                                         # Python booleans -> True/False text
            for c in self.columns:
                if c.kind == "bool" and isinstance(r[c.name], bool):
                    r[c.name] = str(r[c.name])
        problems = self.check_rows(rows, prefix)
        if problems:
            raise SchemaError(f"not writing {path}: " + "; ".join(problems[:5]))
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(path) + ".part")
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, self.names, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        tmp.replace(path)


def _prefix_from_name(path):
    p = Path(path).name.split(".")[0]
    return p if PREFIX_RE.match(p) else None


# =============================================================================== genes table
GENES = TableSchema("genes", [
    Column("gene_id", "id", True), Column("ncbi_gene_id"), Column("annotation_gene_key"), Column("symbol"),
    Column("description"), Column("biotype"), Column("seqid", required=True), Column("start", "int", True),
    Column("end", "int", True), Column("strand", "enum", values=("+", "-", ".", "?")),
    Column("mitochondrial", "bool", True), Column("n_transcripts", "int"), Column("n_proteins", "int"),
    Column("chosen_protein"), Column("chosen_transcript"), Column("protein_length", "int"),
    Column("partial", "bool"), Column("low_quality", "bool"), Column("residues_replaced", "int"),
    Column("internal_stops", "int"), Column("product")],
    file_suffix=".genes.tsv", doc="one row per gene; proteins only for protein-coding genes")

# =============================================================================== family assignments
ASSIGN = TableSchema("assign", [
    Column("gene_id", "id", True), Column("protein_length", "int"), Column("best_family", "family"), Column("score", "float"),
    Column("evalue", "float"), Column("protein_coverage", "float"), Column("family_coverage", "float"),
    Column("protein_region"), Column("n_families", "int"), Column("second_family", "family"), Column("second_score", "float"),
    Column("fusion_partner"), Column("fusion_partner_score", "float"), Column("fusion_partner_region"),
    Column("fusion_partner_family_coverage", "float"), Column("flags")],
    file_suffix=".assign.tsv", doc="each protein's best family (same format for every method; for OrthoFinder the "
    "family is the orthogroup and the score columns may be empty)")

# =============================================================================== orthology calls
ORTHO_BASE = [
    Column("gene_id", "id", True), Column("ncbi_gene_id"), Column("ncbi_symbol"), Column("ncbi_description"),
    Column("seqid"), Column("start", "int"), Column("end", "int"), Column("strand", "enum", values=("+", "-", ".", "?")),
    Column("mitochondrial", "bool"), Column("protein"), Column("protein_length", "int"), Column("treefam_family"),
    Column("tree_mode"), Column("placement_status"), Column("placement_lwr", "float")]
REF_FIELDS = [
    ("relationship", "enum", ("Ortholog", "Paralog", "NoHomolog")), ("ortholog_type", "enum", ("OneToOne", "OneToTwo", "OneToMany")),
    ("symbol", "str", ()), ("official_id", "str", ()), ("ensembl_gene", "str", ()), ("treefam_protein", "str", ()),
    ("description", "str", ()), ("name_status", "str", ()), ("ref_in_subtree", "str", ()),
    ("own_species_in_subtree", "str", ()), ("median_cichlid_copies", "float", ()), ("distance", "float", ()),
    ("call_support", "str", ())]
ORTHO_END = [Column("human_uniprot"), Column("ncbi_symbol_agrees"),
             Column("confidence", "enum", values=("high", "medium", "low")), Column("caveats")]
DEFAULT_REFS = ("human", "mouse", "zebrafish")


def orthologs_schema(refs=DEFAULT_REFS):
    cols = list(ORTHO_BASE)
    for ref in refs:
        cols += [Column(f"{ref}_{name}", kind, values=vals) for name, kind, vals in REF_FIELDS]
    return TableSchema("orthologs", cols + ORTHO_END, file_suffix=".orthologs.tsv",
                       doc="orthology calls of one method, one block of columns per reference species")


ORTHOLOGS = orthologs_schema()


def refs_in_header(header):
    """Reference names, in order, from an orthologs header (blocks of <ref>_<field> columns)."""
    refs, i = [], len(ORTHO_BASE)
    n = len(REF_FIELDS)
    while i + n <= len(header) and header[i].endswith("_relationship"):
        ref = header[i][:-len("_relationship")]
        refs.append(ref)
        i += n
    return refs


# =============================================================================== proteome FASTA
def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0] if line[1:].strip() else "", []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def validate_proteome(path, prefix=None, genes_table=None):
    """Headers are unique gene IDs with the file's prefix; sequences use the 20 amino acids and X
    (no stop codons); with a genes table, every protein-coding gene has exactly one protein."""
    prefix = prefix or _prefix_from_name(path)
    problems, ids = [], set()
    for name, seq in read_fasta(path):
        if not GENE_ID_RE.match(name):
            problems.append(f"{name!r} is not a gene ID")
        elif prefix and gene_prefix(name) != prefix:
            problems.append(f"{name!r} does not start with {prefix}_")
        if name in ids:
            problems.append(f"duplicate protein {name!r}")
        ids.add(name)
        if not seq:
            problems.append(f"{name}: empty sequence")
        bad = set(seq.upper()) - AMINO_ACIDS
        if bad:
            problems.append(f"{name}: characters outside the amino-acid alphabet: {''.join(sorted(bad))}")
        if len(problems) >= MAX_PROBLEMS:
            return problems + ["(stopping after too many problems)"]
    if genes_table and Path(genes_table).is_file():
        coding = {r["gene_id"] for r in GENES.read(genes_table) if r["chosen_protein"]}
        for g in sorted(coding - ids)[:10]:
            problems.append(f"{g}: has a chosen protein in the genes table but none in the proteome")
        for g in sorted(ids - coding)[:10]:
            problems.append(f"{g}: in the proteome but has no chosen protein in the genes table")
    return problems


# =============================================================================== run record and family
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
TREE_MODES = ("place", "denovo")             # how a family's tree was made: placed on a reference tree, or built from scratch
RESIDUE_MAP_FILE = "residue_map.tsv"
RESIDUE_MAP_HEADER = ["gene_id", "protein_length", "aln_length", "residues_not_in_alignment", "columns"]
# run.json: besides RUN_KEYS, an optional "proteomes_md5" (Catalog.proteomes_md5 of the run's species): the checksum of exactly
# the proteomes the run used, which `release validate` compares with a frozen release
RUN_KEYS = ("release", "method", "method_version", "parameters", "started", "finished", "coordinator_version",
            "catalog_md5", "species")


def release_of_location(path, kind):
    """The release a results file or folder sits in, from its path: Results/<release>/<part>/run.json
    (kind "run") or Results/<release>/<part>/families/<family> (kind "family"); None if it sits elsewhere."""
    p = Path(path).resolve()
    if kind == "family":
        if p.parent.name != "families":
            return None
        p = p.parent.parent / "x"                                        # as if a file in the part's folder
    chain = [p.parent, p.parent.parent, p.parent.parent.parent]          # part, release, Results
    if chain[0].name in RESULT_PARTS and RELEASE_RE.match(chain[1].name) and chain[2].name == RESULTS:
        return chain[1].name
    return None


def validate_run(path):
    try:
        d = json.loads(Path(path).read_text())
    except ValueError as e:
        return [f"not valid JSON: {e}"]
    problems = [f"missing key: {k}" for k in RUN_KEYS if k not in d]
    if "release" in d:
        if not RELEASE_RE.match(str(d["release"])):
            problems.append(f"release {d['release']!r} must be like v1")
        else:
            where = release_of_location(path, "run")
            if where and where != d["release"]:
                problems.append(f"release is {d['release']} but the file is in the results of {where}")
    if "proteomes_md5" in d and not MD5_RE.match(str(d["proteomes_md5"])):
        problems.append("proteomes_md5 must be a 32-character MD5 (Catalog.proteomes_md5)")
    if "species" in d and not (isinstance(d["species"], list) and all(PREFIX_RE.match(str(s)) for s in d["species"])):
        problems.append("species must be a list of prefixes")
    if "parameters" in d and not isinstance(d["parameters"], dict):
        problems.append("parameters must be an object")
    return problems


def write_run(path, **record):
    missing = [k for k in RUN_KEYS if k not in record]
    if missing:
        raise SchemaError(f"run record missing: {', '.join(missing)}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(record, indent=2) + "\n")


def read_done(folder):
    """A family's done.json as a dict; None if it is missing or is not a JSON object."""
    path = Path(folder) / "done.json"
    try:
        d = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    return d if isinstance(d, dict) else None


def family_rooting(folder):
    """-> (rooted, how): what done.json says about the root of tree.nwk.

    `rooted` is True only when done.json says "rooted": true; `how` is its "rooting" text (for example
    "orthofinder species tree"). A family that says nothing is not rooted as far as consumers are
    concerned: they midpoint-root it. Consumers use an existing root and root only otherwise.
    """
    d = read_done(folder) or {}
    if d.get("rooted") is True:
        return True, str(d.get("rooting", ""))
    return False, ""


def _alignment_problems(path, leaves):
    rows, lengths, problems = set(), set(), []
    for name, seq in read_fasta(path):
        if name in rows:
            problems.append(f"{ALIGNMENT_FILE}: duplicate row {name!r}")
        rows.add(name)
        lengths.add(len(seq))
    if not rows:
        return [f"{ALIGNMENT_FILE} is empty"]
    if len(lengths) > 1:
        problems.append(f"{ALIGNMENT_FILE}: rows have different lengths ({min(lengths)} to {max(lengths)}): not aligned")
    absent = sorted(set(leaves) - rows)
    if absent:
        more = f" (and {len(absent) - 5} more)" if len(absent) > 5 else ""
        problems.append(f"{ALIGNMENT_FILE}: tree leaves without a row: {', '.join(absent[:5])}{more}")
    return problems


def _leaf_problems(leaves, prefixes):
    problems, seen = [], set()
    for name in leaves:
        if not name:
            problems.append("tree.nwk has a leaf with no name")
        elif name in seen:
            problems.append(f"tree.nwk: duplicate leaf {name!r}")
        elif any(ch.isspace() for ch in name):
            problems.append(f"tree.nwk: leaf {name!r} contains whitespace")
        elif ORTHOFINDER_INTERNAL_RE.match(name):
            problems.append(f"tree.nwk: leaf {name!r} looks like an OrthoFinder internal sequence name; "
                            "map it back to our gene ID")
        elif prefixes and "_" in name and name.split("_", 1)[0] in prefixes and not GENE_ID_RE.match(name):
            problems.append(f"tree.nwk: leaf {name!r} starts with a species prefix and '_' but is not a gene ID "
                            "(<prefix>_<id>; id without underscores)")
        seen.add(name)
        if len(problems) >= MAX_PROBLEMS:
            break
    return problems


def validate_family(folder, prefixes=None, check_alignment=True):
    """A finished family folder, Results/<run>/families/<family>/, for any method:

    folder name   the family ID: letters, digits, '.', '_' and '-' (not starting with '_')
    tree.nwk      a Newick tree; leaf names are unique and are our gene IDs (<prefix>_<id>) or reference
                  protein IDs. With `prefixes` (the catalog's species prefixes) a leaf whose text before
                  its first '_' is one of them must be a well-formed gene ID; any other name (a bare
                  prefix, a longer name that only starts with the same letters) is a reference ID and
                  is left alone. OrthoFinder's internal sequence names are refused.
    done.json     a JSON object (its inputs hash, and so on). If it has "release" (v1, v2, ...) it matches the
                  release of the folder it sits in. If it has "rooted", that is true or false;
                  "rooted": true needs a "rooting" text and a tree whose root has two children.
    combined.aln  required: aligned FASTA (equal row lengths) with a row for every leaf of the tree
                  (rows without a leaf are allowed). Other files in the folder (queries.fa, ref.tree,
                  epa_result.jplace, placements.tsv, residue_map.tsv, ...) belong to the method.
    """
    folder, problems = Path(folder), []
    name = folder.resolve().name
    if not FAMILY_ID_RE.match(name):
        problems.append(f"family ID {name!r} must use letters, digits, '.', '_' and '-' (not starting with '_')")
    tree, root = folder / "tree.nwk", None
    if not tree.is_file():
        problems.append("tree.nwk is missing")
    else:
        t, tree_problems = tree.read_text().strip(), []
        if not t.endswith(";"):
            tree_problems.append("tree.nwk does not end with ';'")
        depth = 0
        for ch in t:
            depth += (ch == "(") - (ch == ")")
            if depth < 0:
                break
        if depth != 0:
            tree_problems.append("tree.nwk has unbalanced parentheses")
        if not tree_problems:                                  # only a well-formed tree is parsed and checked
            try:
                root = newick.parse(t)
            except ValueError as e:
                tree_problems.append(f"tree.nwk does not parse: {e}")
        problems += tree_problems
    leaves = [n.name for n in root.leaves()] if root is not None else []
    if root is not None:
        problems += _leaf_problems(leaves, prefixes)
    done = folder / "done.json"
    if not done.is_file():
        problems.append("done.json is missing (the family is not finished)")
    else:
        try:
            d = json.loads(done.read_text())
        except ValueError:
            d = None
            problems.append("done.json is not valid JSON")
        if d is not None and not isinstance(d, dict):
            problems.append("done.json must be a JSON object")
        elif isinstance(d, dict) and "release" in d:
            if not RELEASE_RE.match(str(d["release"])):
                problems.append(f"done.json: release {d['release']!r} must be like v1")
            else:
                where = release_of_location(folder, "family")
                if where and where != d["release"]:
                    problems.append(f"done.json says release {d['release']} but the family is in the results of {where}")
        if isinstance(d, dict) and "tree_mode" in d and d["tree_mode"] not in TREE_MODES:
            problems.append(f"done.json: tree_mode {d['tree_mode']!r} must be one of {', '.join(TREE_MODES)}")
        if isinstance(d, dict) and "rooted" in d:
            if not isinstance(d["rooted"], bool):
                problems.append(f"done.json: rooted must be true or false, not {d['rooted']!r}")
            elif d["rooted"]:
                if not (isinstance(d.get("rooting"), str) and d["rooting"].strip()):
                    problems.append('done.json: "rooted": true needs a "rooting" text saying how it was rooted')
                if root is not None and len(root.children) != 2:
                    problems.append(f"done.json says rooted but the root of tree.nwk has {len(root.children)} "
                                    "children (a rooted tree has two)")
    aln = folder / ALIGNMENT_FILE
    if check_alignment:
        if not aln.is_file():
            problems.append(f"{ALIGNMENT_FILE} is missing")
        elif root is not None:
            problems += _alignment_problems(aln, leaves)
    rmap = folder / RESIDUE_MAP_FILE
    if rmap.is_file():                                   # optional; checked when present
        problems += [f"{RESIDUE_MAP_FILE}: {pr}" for pr in
                     validate_residue_map(rmap, _alignment_length(aln) if check_alignment else None)]
    return problems


def _alignment_length(path):
    """The common row length of an aligned FASTA; None if it is missing, empty or ragged."""
    try:
        lengths = {len(seq) for _, seq in read_fasta(path)}
    except OSError:
        return None
    return lengths.pop() if len(lengths) == 1 else None


def _residue_row_problem(f, aln_length):
    """The first problem with one residue_map.tsv row (already split into 5 fields), or None."""
    gid, plen, alen, notin, cols = f
    if not GENE_ID_RE.match(gid):
        return f"{gid!r} is not a gene ID (<prefix>_<id>; id without underscores)"
    try:
        plen, alen, notin = int(plen), int(alen), int(notin)
    except ValueError:
        return f"{gid}: protein_length, aln_length and residues_not_in_alignment must be whole numbers"
    if plen < 1 or alen < 1 or notin < 0:
        return f"{gid}: protein_length and aln_length must be at least 1, residues_not_in_alignment 0 or more"
    entries = cols.split(",") if cols else []
    if len(entries) != plen:
        return f"{gid}: columns has {len(entries)} entries for a protein of length {plen} (one per residue)"
    last, dashes = 0, 0
    for e in entries:
        if e == "-":
            dashes += 1
        elif not e.isdigit() or not 1 <= int(e) <= alen:
            return f"{gid}: column {e!r} must be '-' or a whole number from 1 to aln_length ({alen})"
        elif int(e) <= last:
            return f"{gid}: columns must increase along the sequence ({e} after {last})"
        else:
            last = int(e)
    if dashes != notin:
        return f"{gid}: residues_not_in_alignment is {notin} but columns has {dashes} '-' entries"
    return None


def validate_residue_map(path, aln_length=None):
    """A residue_map.tsv: for each protein, the column of the family's combined.aln that each of its residues is in.

    Tab-separated, with the header gene_id, protein_length, aln_length, residues_not_in_alignment, columns; one row
    per protein, sorted by gene_id, no duplicates. `columns` is a comma-separated list with one entry per residue in
    sequence order: the 1-based column, or '-' for a residue left out of the alignment (an insertion MAFFT
    --keeplength dropped). Columns increase along the sequence; aln_length is the same on every row (and is
    `aln_length`, the length of combined.aln, if given). Rows can be very long, so no csv module is used."""
    try:
        fh = open(path, encoding="utf-8")
    except OSError as e:
        return [str(e)]
    problems, previous, lengths = [], None, set()
    with fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        if header != RESIDUE_MAP_HEADER:
            return [f"header must be {', '.join(RESIDUE_MAP_HEADER)}"]
        for n, line in enumerate(fh, start=2):
            line = line.rstrip("\r\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) != len(RESIDUE_MAP_HEADER):
                problems.append(f"line {n}: expected {len(RESIDUE_MAP_HEADER)} columns, found {len(f)}")
            else:
                pr = _residue_row_problem(f, aln_length)
                if pr:
                    problems.append(f"line {n}: {pr}")
                else:
                    lengths.add(int(f[2]))
                    if previous is not None and f[0] <= previous:
                        problems.append(f"line {n}: {f[0]} is not after {previous} (rows are sorted by gene_id, with no duplicates)")
                    previous = f[0]
            if len(problems) >= MAX_PROBLEMS:
                return problems + ["(stopping after too many problems)"]
    if len(lengths) > 1:
        problems.append(f"aln_length differs between rows ({min(lengths)} to {max(lengths)})")
    elif lengths and aln_length is not None and lengths != {aln_length}:
        problems.append(f"aln_length is {next(iter(lengths))} but combined.aln is {aln_length} columns long")
    return problems


def validate_families(folder, prefixes=None, check_alignment=True):
    """Every family folder inside a families/ folder -> problems, each starting with the family's ID."""
    folder, problems, n = Path(folder), [], 0
    if not folder.is_dir():
        return [f"{folder} is not a folder"]
    for fam in sorted(e for e in folder.iterdir() if e.is_dir() and not e.name.startswith(".")):
        n += 1
        problems += [f"{fam.name}: {pr}" for pr in validate_family(fam, prefixes, check_alignment)]
        if len(problems) >= MAX_PROBLEMS:
            return problems[:MAX_PROBLEMS] + [f"(stopping after too many problems; {n} families checked)"]
    return problems or ([] if n else [f"{folder} has no family folders"])


# =============================================================================== dispatch
def kind_of(path):
    p, name = Path(path), Path(path).name
    if p.is_dir():
        return "families" if p.resolve().name == "families" else "family"
    for kind, suffix in (("genes", ".genes.tsv"), ("assign", ".assign.tsv"), ("orthologs", ".orthologs.tsv")):
        if name.endswith(suffix):
            return kind
    if name.endswith((".faa", ".fa", ".fasta")):
        return "proteome"
    if name == "run.json":
        return "run"
    if name == RESIDUE_MAP_FILE:
        return "residue_map"
    return None


def validate(path, kind=None, prefix=None, prefixes=None):
    """-> list of problems (empty if the file is valid). `prefix` is the one species a table or proteome
    must belong to; `prefixes` the set of all species prefixes, for the leaf names of family trees."""
    kind = kind or kind_of(path)
    if kind == "genes":
        return GENES.validate_file(path, prefix)
    if kind == "assign":
        return ASSIGN.validate_file(path, prefix)
    if kind == "orthologs":
        with open(path, newline="") as fh:
            header = next(csv.reader(fh, delimiter="\t"), [])
        return orthologs_schema(refs_in_header(header) or DEFAULT_REFS).validate_file(path, prefix)
    if kind == "proteome":
        genes = Path(str(path)[: -len(Path(path).suffix)] + ".genes.tsv")
        return validate_proteome(path, prefix, genes if genes.is_file() else None)
    if kind == "run":
        return validate_run(path)
    if kind == "residue_map":
        return validate_residue_map(path)
    if kind == "family":
        return validate_family(path, prefixes)
    if kind == "families":
        return validate_families(path, prefixes)
    raise SchemaError(f"can't tell what kind of file {path} is; give the kind")


SCHEMAS = {"genes": GENES, "assign": ASSIGN, "orthologs": ORTHOLOGS}


