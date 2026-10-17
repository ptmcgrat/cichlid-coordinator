"""Every data path in the project. Callers ask for things; they never build paths.

    L = get_layout()
    L.proteome("Mzeb").local      # Path on this machine
    L.proteome("Mzeb").cloud      # rclone path in the cloud
    open(L.proteome("Mzeb"))      # DataPath works wherever a path does (it is the local path)

All paths are relative to one root, identical locally and in the cloud.

Results belong to a release (one number shared by every project: v1, v2, ...) and a part (the project
that wrote them): Results/<release>/<part>/. Releases/<release>/ is the frozen archive of a release.

    L.results("v1", "OrthoFinder").family("OG0001234").local
    L.release("v1").proteome("Mzeb").local
"""
import os
import re
from dataclasses import dataclass
from pathlib import Path

INPUTS, RESULTS, ANALYSES, RELEASES = "Inputs", "Results", "Analyses", "Releases"
PROTEINS = f"{INPUTS}/Proteins"
FAMILIES = f"{INPUTS}/Families"
RELEASE_RE = re.compile(r"^v[0-9]+$")
RESULT_PARTS = ("TreeFam", "Compara", "OrthoFinder", "Homology", "Analysis")


def check_release(release):
    if not RELEASE_RE.match(str(release)):
        raise ValueError(f"release {release!r} must be 'v' and a number (v1, v2, ...)")
    return release


def check_part(part):
    if part not in RESULT_PARTS:
        raise KeyError(f"unknown part {part!r} (one of {', '.join(RESULT_PARTS)})")
    return part


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
    """Identical structure for every part's results: Results/<release>/<part>/...

    native/ holds a part's own outputs (OrthoFinder's result folders, say), kept beside the standard
    sections so they are backed up and pulled like everything else."""
    SECTIONS = ("assign", "families", "names", "native", "orthology", "residues", "site")

    def __init__(self, layout, release, part):
        self.layout, self.release, self.part = layout, check_release(release), check_part(part)
        self.base = f"{RESULTS}/{release}/{part}"

    def _p(self, rel):
        return DataPath(f"{self.base}/{rel}", self.layout)

    def root(self):
        return DataPath(self.base, self.layout)

    def section(self, section):
        if section not in self.SECTIONS:
            raise KeyError(f"unknown results section {section!r} (one of {', '.join(self.SECTIONS)})")
        return self._p(section)

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

    def unit(self, section):
        from .units import Registry
        return Registry(self.layout).get(f"results/{self.release}/{self.part}/{section}")


class ReleaseLayout:
    """Releases/<release>/: the frozen archive of one release (written by `cichlid-data release freeze`):
    release.json (what the release is, with checksums), catalog.json, proteomes/, genes/ and config/
    (the species trees), exactly as they were when it was frozen."""

    def __init__(self, layout, release):
        self.layout, self.release = layout, check_release(release)
        self.base = f"{RELEASES}/{release}"

    def _p(self, rel):
        return DataPath(f"{self.base}/{rel}", self.layout)

    def root(self):
        return DataPath(self.base, self.layout)

    def manifest(self):
        return self._p("release.json")

    def catalog(self):
        return self._p("catalog.json")

    def proteome(self, prefix):
        return self._p(f"proteomes/{prefix}.faa")

    def genes_table(self, prefix):
        return self._p(f"genes/{prefix}.genes.tsv")

    def species_tree(self, name):
        return self._p(f"config/{name}")

    def unit(self):
        from .units import Registry
        return Registry(self.layout).get(f"releases/{self.release}")


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
    def results(self, release, part):
        return ResultsLayout(self, release, part)

    def release(self, release):
        return ReleaseLayout(self, release)

    def runs(self):
        """(release, part) of every results folder on this machine."""
        d = self.local_root / RESULTS
        out = []
        for r in sorted(e for e in d.iterdir() if e.is_dir()) if d.is_dir() else []:
            out += [(r.name, p.name) for p in sorted(r.iterdir()) if p.is_dir() and p.name in RESULT_PARTS]
        return out

    def releases(self):
        """The releases that have results or a frozen archive on this machine."""
        found = {r for r, _ in self.runs()}
        d = self.local_root / RELEASES
        if d.is_dir():
            found |= {e.name for e in d.iterdir() if e.is_dir() and RELEASE_RE.match(e.name)}
        return sorted(found)

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


