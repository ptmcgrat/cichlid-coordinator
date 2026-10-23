"""The names of the reference proteins of TreeFam and Compara: each human, mouse and zebrafish reference protein, its current
gene, symbol and official ID (HGNC, MGI or ZFIN).

    RN = ReferenceNames.load(get_layout(), "TreeFam")        # or "Compara"
    r = RN.get("ENSP00000409362")
    r.symbol, r.official_id, r.status, r.reliable            # "FAM138A", "HGNC:32334", "legacy_id", True

The tables are a frozen snapshot made by the earlier pipeline (its stage 6a), kept under
Inputs/Annotations/ReferenceNames/<part>/reference_names.tsv. TreeFam 9's Ensembl protein IDs date from about 2013 and many are
retired, so each protein is resolved by up to five routes, in this order, and the status says which one succeeded:

    current        the Ensembl REST service gave the protein's current gene
    grch37         (human only) the frozen GRCh37 service did, and the gene is still current
    legacy_id      the 2018 tables gave an HGNC, ZFIN or MGI accession, a UniProt accession or an NCBI gene ID, matched to the
                   nomenclature authority's current record
    replaced       the best-scoring replacement protein from Ensembl's ID history
    legacy_symbol  the 2018 symbol, unverified
    unmapped       none worked

The earlier pipeline tried the 2018 tables before the ID history on purpose: the 2018 files record which gene the protein belonged
to, while a replacement can occasionally point to a different gene. So `reliable` is true for current, grch37 and legacy_id, and
false for replaced, legacy_symbol and unmapped. This is a judgement taken from the pipeline's own comments, not a measurement.
The symbol comes from the nomenclature authority (HGNC, ZFIN or MGI); only when the authority has no entry is Ensembl's display
name used (`symbol_source` says which). Compara 116's IDs are current, so nearly all its rows are `current`.

OrthoFinder has no such table: its reference genes are our own proteomes (`Hsap_<gene ID>`, `Drer_<gene ID>`), whose symbols are in
the genes tables. Their HGNC and ZFIN IDs are a later addition.
"""
import re
from dataclasses import dataclass
from pathlib import Path

from .schemas import MAX_PROBLEMS

FILE = "reference_names.tsv"
HEADER = ["protein_id", "species", "status", "route", "mapped_protein", "transcript_id", "gene_id", "symbol", "symbol_source",
          "official_id", "description"]
STATUSES = ("current", "grch37", "legacy_id", "replaced", "legacy_symbol", "unmapped")
RELIABLE = ("current", "grch37", "legacy_id")
HAS_GENE = ("current", "grch37", "legacy_id", "replaced")           # a mapped protein has a gene
SPECIES = ("homo_sapiens", "mus_musculus", "danio_rerio")
ID_PATTERNS = {"homo_sapiens": re.compile(r"^HGNC:\d+$"), "mus_musculus": re.compile(r"^MGI:\d+$"),
               "danio_rerio": re.compile(r"^ZDB-[A-Z0-9]+-\d{6}-\d+$")}
PARTS = ("TreeFam", "Compara")


def table_path(layout, part):
    return Path(layout.annotations(f"ReferenceNames/{part}").local) / FILE


@dataclass(frozen=True)
class Record:
    protein_id: str
    species: str
    status: str
    route: str
    mapped_protein: str
    transcript_id: str
    gene_id: str
    symbol: str
    symbol_source: str
    official_id: str
    description: str

    @property
    def reliable(self):
        """The protein was resolved by a route the earlier pipeline trusted (current, grch37 or legacy_id)."""
        return self.status in RELIABLE

    @property
    def named(self):
        return bool(self.symbol)


class ReferenceNames:
    def __init__(self, records):
        self._by = records                                    # protein_id -> Record

    @classmethod
    def load(cls, layout, part, only=None):
        """The table of one part. `only`: a set of protein IDs to keep (to bound memory)."""
        if part not in PARTS:
            raise KeyError(f"no reference-names table for {part!r}: only {', '.join(PARTS)} have one (OrthoFinder's reference "
                           "genes are our own proteomes: use the genes tables)")
        path = table_path(layout, part)
        if not path.is_file():
            raise FileNotFoundError(f"{path} not found (run `cichlid-data check --fetch`)")
        keep = set(only) if only is not None else None
        out = {}
        with open(path, encoding="utf-8") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            if header != HEADER:
                raise ValueError(f"{path.name}: the header must be {', '.join(HEADER)}")
            for n, line in enumerate(fh, start=2):
                line = line.rstrip("\r\n")
                if not line:
                    continue
                f = line.split("\t")
                if len(f) != len(HEADER):
                    raise ValueError(f"{path.name} line {n}: expected {len(HEADER)} columns, found {len(f)}")
                if keep is not None and f[0] not in keep:
                    continue
                if f[0] in out:
                    raise ValueError(f"{path.name} line {n}: {f[0]} appears twice")
                out[f[0]] = Record(*f)
        return cls(out)

    def get(self, protein):
        """The record of a reference protein, or None."""
        return self._by.get(protein)

    def __contains__(self, protein):
        return protein in self._by

    def __len__(self):
        return len(self._by)

    def proteins(self, species=None):
        return [p for p, r in self._by.items() if species is None or r.species == species]

    def status_counts(self, species=None):
        out = {}
        for r in self._by.values():
            if species is None or r.species == species:
                out[r.status] = out.get(r.status, 0) + 1
        return out


def validate_reference_names(path):
    """Problems with a reference_names.tsv: the header, repeated proteins, an unknown species or status, an official ID that is not an
    HGNC, MGI or ZFIN ID of the right species, and a mapped protein with no gene."""
    problems = []
    try:
        fh = open(path, encoding="utf-8")
    except OSError as e:
        return [str(e)]
    with fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        if header != HEADER:
            return [f"header must be {', '.join(HEADER)}"]
        seen = set()
        for n, line in enumerate(fh, start=2):
            line = line.rstrip("\r\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) != len(HEADER):
                problems.append(f"line {n}: expected {len(HEADER)} columns, found {len(f)}")
            else:
                r = dict(zip(HEADER, f))
                if not r["protein_id"]:
                    problems.append(f"line {n}: protein_id is empty")
                elif r["protein_id"] in seen:
                    problems.append(f"line {n}: {r['protein_id']} appears twice")
                seen.add(r["protein_id"])
                if r["species"] not in SPECIES:
                    problems.append(f"line {n}: species {r['species']!r} must be one of {', '.join(SPECIES)}")
                elif r["official_id"] and not ID_PATTERNS[r["species"]].match(r["official_id"]):
                    problems.append(f"line {n}: {r['official_id']!r} is not an official ID of {r['species']}")
                if r["status"] not in STATUSES:
                    problems.append(f"line {n}: status {r['status']!r} must be one of {', '.join(STATUSES)}")
                elif r["status"] in HAS_GENE and not r["gene_id"]:
                    problems.append(f"line {n}: {r['protein_id']} has status {r['status']} but no gene_id")
            if len(problems) >= MAX_PROBLEMS:
                return problems + ["(stopping after too many problems)"]
    return problems


def coverage(layout, part):
    """{species: (proteins of the part's mapping file, how many of them have a row)} for human, mouse and zebrafish.
    The earlier pipeline made one row for every reference protein of the mapping file, so each pair should be equal."""
    from .leaf_species import MAPPINGS
    from .master_tree import canonical
    folder, name = MAPPINGS[part]
    mapping = Path(layout.family_resource(folder).local) / name
    have = set(ReferenceNames.load(layout, part)._by)
    total, covered = {s: 0 for s in SPECIES}, {s: 0 for s in SPECIES}
    with open(mapping, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) == 2 and canonical(f[1]) in total:
                sp = canonical(f[1])
                total[sp] += 1
                covered[sp] += f[0] in have
    return {s: (total[s], covered[s]) for s in SPECIES}
