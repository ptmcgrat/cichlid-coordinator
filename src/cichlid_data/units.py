"""Named, transferable pieces of the data tree.

Policy:
    immutable    written once; a push refuses to replace a different cloud copy (unless forced)
    regenerable  pipeline output; replaced by newer runs
    on_demand    large; skipped by bulk pulls, fetched only by name or by Store.ensure
Packing:
    files        copied file by file (large files: HMM libraries, dumps, VCFs)
    tar          one tar, <rel>.tar, holding paths relative to the unit's parent folder
    shards       one tar per shard, <rel>/<shard>.tar, for folders with many entries
                 (families: by ID prefix; downloads: by species); loose files go to _top

Tars hold paths relative to the unit's parent folder (e.g. "families/TF101001/..." for
Results/<run>/families), so a unit can move without being repacked.
"""
import fnmatch
import json
from dataclasses import dataclass, field

from .layout import ANALYSES, FAMILIES, INPUTS, PROTEINS, RESULTS, ResultsLayout

FAMILY_EXCLUDES = ("*/stage5/*.fa", "*.raxml.rba", "*.raxml.ckp", "*.raxml.startTree", "*.raxml.mlTrees",
                   "*.part", "residues.tsv.gz", "residues.done.json")


def family_shard(entry, is_dir):
    """TF101001 -> TF101, C116_012345 -> C116_012 (about 1,000 families per shard); files -> _top."""
    return entry[:-3] if is_dir and len(entry) > 4 else "_top"


def species_shard(entry, is_dir):
    return entry if is_dir else "_top"


@dataclass(frozen=True)
class DataUnit:
    name: str                       # e.g. "results/Compara2026/families"
    rel: str                        # folder or file, relative to the root
    policy: str = "regenerable"     # immutable | regenerable | on_demand
    packing: str = "tar"            # files | tar | shards
    excludes: tuple = ()
    shard_by: object = field(default=None, compare=False)
    description: str = ""

    @property
    def parent(self):
        return self.rel.rsplit("/", 1)[0] if "/" in self.rel else ""

    @property
    def base(self):
        return self.rel.rsplit("/", 1)[-1]

    def cloud_tar(self, shard=None):
        return f"{self.rel}/{shard}.tar" if self.packing == "shards" else f"{self.rel}.tar"

    def key(self, shard=None):
        return f"{self.name}/{shard}" if shard else self.name


class Registry:
    def __init__(self, layout, index=None):
        self.layout = layout
        self.index = index or {}
        self._units = {u.name: u for u in self._build()}

    def _build(self):
        U = DataUnit
        yield U("proteins/downloads", f"{PROTEINS}/downloads", "immutable", "shards", shard_by=species_shard,
                description="genomes and annotations as downloaded, one shard per species")
        yield U("proteins/local_sources", f"{PROTEINS}/local_sources", "immutable", "files",
                description="raw files for species not annotated by NCBI")
        yield U("proteins/proteomes", f"{PROTEINS}/proteomes", "regenerable", "tar",
                description="one protein per gene, with gene tables")
        yield U("proteins/catalog", f"{PROTEINS}/catalog.json", "regenerable", "files", description="species catalog")
        yield U("families/TreeFam9/hmm", f"{FAMILIES}/TreeFam9/hmm", "immutable", "files", description="TreeFam 9 HMMs")
        yield U("families/TreeFam9/trees", f"{FAMILIES}/TreeFam9/trees", "immutable", "tar",
                description="TreeFam 9 alignments and trees")
        yield U("families/TreeFam9/protein_to_species", f"{FAMILIES}/TreeFam9/protein_to_species.txt", "immutable", "files")
        yield U("families/Compara116/dumps", f"{FAMILIES}/Compara116/dumps", "on_demand", "files",
                description="raw Ensembl Compara dumps")
        yield U("families/Compara116/prepared", f"{FAMILIES}/Compara116/prepared", "regenerable", "tar",
                description="Compara trees pruned to the panel")
        yield U("annotations", f"{INPUTS}/Annotations", "regenerable", "files", description="GO, PANTHER, names")
        yield U("expression", f"{INPUTS}/Expression", "regenerable", "files", excludes=("raw/*",),
                description="expression summaries")
        yield U("expression/raw", f"{INPUTS}/Expression/raw", "on_demand", "files", description="raw expression data")
        yield U("populations", f"{INPUTS}/Populations", "on_demand", "files", description="population data (VCFs)")
        yield U("archive", "Archive", "on_demand", "files", description="records of earlier projects (2018 resource)")
        for run in self._runs():
            for part in ResultsLayout.PARTS:
                if part == "families":
                    yield U(f"results/{run}/families", f"{RESULTS}/{run}/families", "regenerable", "shards",
                            excludes=FAMILY_EXCLUDES, shard_by=family_shard)
                else:
                    yield U(f"results/{run}/{part}", f"{RESULTS}/{run}/{part}", "regenerable", "tar")
            yield U(f"results/{run}/run", f"{RESULTS}/{run}/run.json", "regenerable", "files")
        for name in self._analyses():
            yield U(f"analyses/{name}", f"{ANALYSES}/{name}", "regenerable", "tar")

    def _discover(self, top, prefix):
        found = set()
        d = self.layout.local_root / top
        if d.is_dir():
            found |= {e.name for e in d.iterdir() if e.is_dir() and not e.name.startswith(".")}
        found |= {k.split("/")[1] for k in self.index if k.startswith(prefix) and k.count("/") >= 1}
        return sorted(found)

    def _runs(self):
        return self._discover(RESULTS, "results/")

    def _analyses(self):
        return self._discover(ANALYSES, "analyses/")

    def all(self):
        return list(self._units.values())

    def get(self, name):
        if name not in self._units:
            raise KeyError(f"unknown unit {name!r}; see `cichlid-data units`")
        return self._units[name]

    def matching(self, pattern):
        """Units whose name matches a glob (e.g. 'results/Compara2026/*') or equals it."""
        return [u for u in self.all() if u.name == pattern or fnmatch.fnmatch(u.name, pattern)]

    def for_path(self, rel):
        """-> (unit, shard or None) owning a path relative to the root."""
        rel = rel.strip("/")
        best = None
        for u in self.all():
            if rel == u.rel or rel.startswith(u.rel + "/"):
                if best is None or len(u.rel) > len(best.rel):
                    best = u
        if best is None:
            raise KeyError(f"no unit owns {rel!r}")
        if best.packing != "shards":
            return best, None
        rest = rel[len(best.rel):].strip("/")
        if not rest:
            return best, None
        entry = rest.split("/")[0]
        is_dir = "/" in rest or (self.layout.local_root / best.rel / entry).is_dir() or "." not in entry
        return best, best.shard_by(entry, is_dir)


def load_index_text(text):
    return json.loads(text) if text and text.strip() else {}
