"""The species-to-clade table (config/species_clades.tsv): where each species sits on the ladder.

    CL = Clades.load()                                  # reads config/species_clades.tsv
    CL.level("Mzeb")                                    # -> "cichlid"  (a prefix, or a species name, in any case)
    CL.common_level(["Mzeb", "danio_rerio"])            # -> "teleost"  (see below)
    CL.within("percomorph")                             # species at that level or inside it

The ladder has seven nested levels, from outside to inside, and every species has the narrowest one it belongs to:

    outside_vertebrates   no vertebrate whole-genome duplication (invertebrates, tunicates, amphioxus; plants, fungi)
    cyclostome            lamprey and hagfish. DISPUTED: whether they share the second vertebrate round is disputed
    gnathostome           jawed vertebrates outside the teleosts: both vertebrate rounds (elephant shark, coelacanth,
                          tetrapods, gar, bowfin, reedfish; the `lineage` column tells them apart)
    teleost               Teleostei, but not Percomorpha: both rounds plus the teleost-specific round
    percomorph            Percomorpha, but not Ovalentaria
    ovalentaria           Ovalentaria, but not cichlid
    cichlid               Cichlidae

The levels are nested, so the last common ancestor of a set of species sits at the OUTERMOST of their levels:
two cichlids meet at `cichlid`, a cichlid and a percomorph at `percomorph`, a cichlid and a zebrafish at `teleost`,
a cichlid and a human at `gnathostome`. This is what dates a duplication: the common level of the species below it.
"""
from pathlib import Path

from .config import config_dir

LADDER = ("outside_vertebrates", "cyclostome", "gnathostome", "teleost", "percomorph", "ovalentaria", "cichlid")
COLUMNS = ["species", "prefix", "clade", "lineage", "basis", "notes"]
FILE = "species_clades.tsv"
DISPUTED = ("cyclostome",)           # levels whose position is disputed in the literature


class Clades:
    def __init__(self, rows):
        self.rows = rows                                         # list of dicts with the COLUMNS
        self._by = {}
        for r in rows:
            for k in (r["species"], r["prefix"]):
                if k:
                    self._by.setdefault(k.lower(), r)

    @classmethod
    def load(cls, cfg=None):
        path = Path(cfg or config_dir()) / FILE
        if not path.is_file():
            raise FileNotFoundError(f"{path} not found")
        lines = [ln.rstrip("\r\n") for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        header = lines[0].split("\t") if lines else []
        if header != COLUMNS:
            raise ValueError(f"{path.name}: the header must be {', '.join(COLUMNS)}")
        rows = []
        for n, ln in enumerate(lines[1:], start=1):
            f = ln.split("\t")
            f += [""] * (len(COLUMNS) - len(f))
            if len(f) != len(COLUMNS):
                raise ValueError(f"{path.name}: row {n} has {len(f)} columns, not {len(COLUMNS)}")
            rows.append(dict(zip(COLUMNS, f)))
        return cls(rows)

    def __contains__(self, name):
        return str(name).lower() in self._by

    def _row(self, name):
        try:
            return self._by[str(name).lower()]
        except KeyError:
            raise KeyError(f"{name!r} is not in {FILE}") from None

    def level(self, name):
        """The narrowest ladder level of a species, by prefix or by species name (any case)."""
        return self._row(name)["clade"]

    def lineage(self, name):
        """The finer group of a chordate or non-teleost vertebrate (chondrichthyan, sarcopterygian, ...); '' if none."""
        return self._row(name)["lineage"]

    def rank(self, level):
        """0 for the outermost level (outside_vertebrates), 6 for the innermost (cichlid)."""
        try:
            return LADDER.index(level)
        except ValueError:
            raise KeyError(f"{level!r} is not a level of the ladder ({', '.join(LADDER)})") from None

    def common_level(self, names):
        """The level where the last common ancestor of these species sits: the outermost of their levels."""
        names = list(names)
        if not names:
            raise ValueError("no species given")
        return min((self.level(n) for n in names), key=self.rank)

    def within(self, level):
        """The species whose level is `level` or inside it (for example `teleost` gives every teleost)."""
        r = self.rank(level)
        return [x["species"] for x in self.rows if self.rank(x["clade"]) >= r]

    def validate(self, catalog=None, cfg=None):
        """-> problems: bad levels, repeated names or prefixes, organisms of the catalog with no row, and (given
        `cfg`) species of the Compara panel's tree leaves with no row."""
        problems, seen_sp, seen_pre = [], {}, {}
        for n, r in enumerate(self.rows, start=1):
            if r["clade"] not in LADDER:
                problems.append(f"{r['species'] or 'row ' + str(n)}: clade {r['clade']!r} must be one of {', '.join(LADDER)}")
            if not r["species"] or r["species"] != r["species"].lower() or " " in r["species"]:
                problems.append(f"row {n}: species {r['species']!r} must be lower case with underscores")
            for key, seen, what in ((r["species"].lower(), seen_sp, "species"), (r["prefix"].lower(), seen_pre, "prefix")):
                if key and key in seen:
                    problems.append(f"{what} {key!r} appears twice (rows {seen[key]} and {n}); case does not count")
                elif key:
                    seen[key] = n
        if catalog is not None:
            for p in catalog.organisms:
                if p.lower() not in seen_pre:
                    problems.append(f"{p} ({catalog.organisms[p].species}) is in the catalog but has no row")
        if cfg is not None:
            panel = Path(cfg) / "compara_panel.tsv"
            if panel.is_file():
                for ln in panel.read_text().splitlines():
                    f = ln.split("\t")
                    if ln.startswith("#") or len(f) < 3 or f[0] == "species" or f[1] in ("exclude", "validation"):
                        continue                          # only the species that are leaves of the pruned Compara trees
                    if f[0].lower() not in seen_sp:
                        problems.append(f"{f[0]} (Compara panel, {f[1]}) has no row")
        return problems
