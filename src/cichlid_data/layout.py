"""Every data path in the project. Callers ask for things; they never build paths.

    L = get_layout()
    L.proteome("Mzeb").local      # Path on this machine
    L.proteome("Mzeb").cloud      # rclone path in the cloud
    open(L.proteome("Mzeb"))      # DataPath works wherever a path does (it is the local path)

All paths are relative to one root, identical locally and in the cloud.
"""
import os
from dataclasses import dataclass
from pathlib import Path

INPUTS, RESULTS, ANALYSES = "Inputs", "Results", "Analyses"
PROTEINS = f"{INPUTS}/Proteins"
FAMILIES = f"{INPUTS}/Families"


@dataclass(frozen=True)
class DataPath:
    rel: str
    layout: "Layout"

    @property
    def local(self):
        return self.layout.local_root / self.rel

    @property
    def cloud(self):
        return f"{self.layout.cloud_remote}/{self.rel}"

    def __fspath__(self):
        return str(self.local)

    def __str__(self):
        return str(self.local)

    def __truediv__(self, other):
        return DataPath(f"{self.rel.rstrip('/')}/{other}", self.layout)

    def exists(self):
        return self.local.exists()


class ResultsLayout:
    """Identical structure for every method's results: Results/<run>/...

    native/ holds a method's own outputs (OrthoFinder's result folders, say), kept beside the standard
    parts so they are backed up and pulled like everything else."""
    PARTS = ("assign", "families", "names", "native", "orthology", "residues", "site")

    def __init__(self, layout, run):
        self.layout, self.run = layout, run
        self.base = f"{RESULTS}/{run}"

    def _p(self, rel):
        return DataPath(f"{self.base}/{rel}", self.layout)

    def root(self):
        return DataPath(self.base, self.layout)

    def part(self, part):
        if part not in self.PARTS:
            raise KeyError(f"unknown results part {part!r} (one of {', '.join(self.PARTS)})")
        return self._p(part)

    def assign(self, prefix):
        return self._p(f"assign/{prefix}.assign.tsv")

    def family(self, family):
        return self._p(f"families/{family}")

    def names(self):
        return self._p("names")

    def native(self):
        return self._p("native")

    def orthology(self, prefix):
        return self._p(f"orthology/{prefix}.orthologs.tsv")

    def residues(self, prefix):
        return self._p(f"residues/{prefix}.residues.tsv.gz")

    def site(self):
        return self._p("site")

    def run_record(self):
        return self._p("run.json")

    def unit(self, part):
        from .units import Registry
        return Registry(self.layout).get(f"results/{self.run}/{part}")


class Layout:
    def __init__(self, local_root, cloud_remote):
        self.local_root = Path(local_root)
        self.cloud_remote = cloud_remote.rstrip("/")

    def path(self, rel):
        return DataPath(rel.strip("/"), self)

    # ---- proteins
    def proteins(self):
        return self.path(PROTEINS)

    def catalog(self):
        return self.path(f"{PROTEINS}/catalog.json")

    def download(self, prefix):
        return self.path(f"{PROTEINS}/downloads/{prefix}")

    def proteomes(self):
        return self.path(f"{PROTEINS}/proteomes")

    def proteome(self, prefix):
        return self.path(f"{PROTEINS}/proteomes/{prefix}.faa")

    def genes_table(self, prefix):
        return self.path(f"{PROTEINS}/proteomes/{prefix}.genes.tsv")

    def local_source(self, prefix):
        """Raw files for a species whose annotation is not from NCBI (before conversion)."""
        return self.path(f"{PROTEINS}/local_sources/{prefix}")

    # ---- other inputs
    def family_resource(self, name):
        return self.path(f"{FAMILIES}/{name}")

    def annotations(self, name=""):
        return self.path(f"{INPUTS}/Annotations/{name}")

    def expression(self, name=""):
        return self.path(f"{INPUTS}/Expression/{name}")

    def population(self, name=""):
        return self.path(f"{INPUTS}/Populations/{name}")

    def archive(self, name=""):
        """Records of earlier projects, e.g. archive("Project2018")."""
        return self.path(f"Archive/{name}")

    # ---- results and analyses
    def results(self, run):
        return ResultsLayout(self, run)

    def runs(self):
        d = self.local_root / RESULTS
        return sorted(e.name for e in d.iterdir() if e.is_dir()) if d.is_dir() else []

    def analysis(self, name):
        return self.path(f"{ANALYSES}/{name}")

    # ---- bookkeeping
    def state_dir(self):
        return self.local_root / ".cichlid-data"

    def env(self):
        """Shell variables for scripts."""
        return {"CICHLID_DATA_ROOT": str(self.local_root), "CICHLID_DATA_REMOTE": self.cloud_remote,
                "CICHLID_PROTEOMES": str(self.proteomes().local),
                "CICHLID_TREEFAM9": str(self.family_resource("TreeFam9").local),
                "CICHLID_COMPARA": str(self.family_resource("Compara116").local)}


def get_layout(settings=None):
    from .config import Settings
    s = settings or Settings.load()
    return Layout(s.local_root, s.cloud_remote)


