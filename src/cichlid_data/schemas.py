"""Formats of the files that cross repository boundaries: readers, writers and validators.

    from cichlid_data import schemas
    rows = schemas.GENES.read(path)                 # list of dicts (strings), header checked
    schemas.GENES.write(path, rows)                 # columns in order; refuses invalid rows
    problems = schemas.validate(path)               # the kind is inferred from the file name

Kinds: genes (<prefix>.genes.tsv), assign (<prefix>.assign.tsv), orthologs (<prefix>.orthologs.tsv),
proteome (<prefix>.faa), run (run.json), family (a families/<family>/ folder).

Gene IDs: <prefix>_<id>; the prefix is letters and digits (starting with a letter), the id letters,
digits, '.' and '-' (no underscores: the first underscore separates the prefix). Tables are
tab-separated with a header row; empty strings mean missing; booleans are True or False.
"""
import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
GENE_ID_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*)_([A-Za-z0-9.\-]+)$")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYX")
MAX_PROBLEMS = 25


class SchemaError(ValueError):
    pass


@dataclass(frozen=True)
class Column:
    name: str
    kind: str = "str"            # str | int | float | bool | enum | id
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
        return None


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
    Column("gene_id", "id", True), Column("protein_length", "int"), Column("best_family"), Column("score", "float"),
    Column("evalue", "float"), Column("protein_coverage", "float"), Column("family_coverage", "float"),
    Column("protein_region"), Column("n_families", "int"), Column("second_family"), Column("second_score", "float"),
    Column("fusion_partner"), Column("fusion_partner_score", "float"), Column("fusion_partner_region"),
    Column("fusion_partner_family_coverage", "float"), Column("flags")],
    file_suffix=".assign.tsv", doc="each protein's best family (same format for every method)")

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
RUN_KEYS = ("method", "method_version", "parameters", "started", "finished", "coordinator_version",
            "catalog_md5", "species")


def validate_run(path):
    try:
        d = json.loads(Path(path).read_text())
    except ValueError as e:
        return [f"not valid JSON: {e}"]
    problems = [f"missing key: {k}" for k in RUN_KEYS if k not in d]
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


def validate_family(folder):
    """A finished family folder: tree.nwk (a Newick tree) and done.json (its inputs hash)."""
    folder, problems = Path(folder), []
    tree = folder / "tree.nwk"
    if not tree.is_file():
        problems.append("tree.nwk is missing")
    else:
        t = tree.read_text().strip()
        if not t.endswith(";"):
            problems.append("tree.nwk does not end with ';'")
        depth = 0
        for ch in t:
            depth += (ch == "(") - (ch == ")")
            if depth < 0:
                break
        if depth != 0:
            problems.append("tree.nwk has unbalanced parentheses")
    done = folder / "done.json"
    if not done.is_file():
        problems.append("done.json is missing (the family is not finished)")
    else:
        try:
            json.loads(done.read_text())
        except ValueError:
            problems.append("done.json is not valid JSON")
    return problems


# =============================================================================== dispatch
def kind_of(path):
    p, name = Path(path), Path(path).name
    if p.is_dir():
        return "family"
    for kind, suffix in (("genes", ".genes.tsv"), ("assign", ".assign.tsv"), ("orthologs", ".orthologs.tsv")):
        if name.endswith(suffix):
            return kind
    if name.endswith((".faa", ".fa", ".fasta")):
        return "proteome"
    if name == "run.json":
        return "run"
    return None


def validate(path, kind=None, prefix=None):
    """-> list of problems (empty if the file is valid)."""
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
    if kind == "family":
        return validate_family(path)
    raise SchemaError(f"can't tell what kind of file {path} is; give the kind")


SCHEMAS = {"genes": GENES, "assign": ASSIGN, "orthologs": ORTHOLOGS}
