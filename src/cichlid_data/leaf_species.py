"""Which species a leaf of a family tree belongs to.

    LS = LeafSpecies.load(get_layout(), "TreeFam")        # or "Compara", or "OrthoFinder"
    LS.species("ENSGALP00000012291")                      # -> "gallus_gallus"   (a reference leaf)
    LS.species("Mzeb_101467881")                          # -> "maylandia_zebra" (one of our genes)
    LS.species("no such leaf")                            # -> None
    LS.of_leaves(tree_leaf_names)                         # {leaf: species}; KeyError listing the unknown ones
    LS.shadowed()                                         # reference species that are also project species

A tree of TreeFam or Compara has two kinds of leaf. **Reference leaves** are the proteins of the reference
families (TreeFam 9's, or Compara 116's pruned to the panel); their species comes from the part's own mapping file
(`protein_to_species.txt`: TreeFam's, or `prepared/protein_to_species.txt` of Compara), one `protein<TAB>species` line each.
**Our genes** are `<prefix>_<id>` with a catalog prefix. The reference mapping is looked up first. OrthoFinder has only our
genes (every species of its 35 is a catalog organism), so it needs no mapping file.

Species names are the master species tree's leaf names (`master_tree.canonical`: lower case, underscores, and the few
alternative spellings folded together).

**Shadowed species.** A reference species that is also a project species appears twice in a tree: as Ensembl's reference
leaves and as our placed genes (a different annotation of the same genome). Species overlap would call every such pair a
duplication, so a consumer ignores the reference leaves of these species when it calls duplications. Release 1: TreeFam
shadows `gasterosteus_aculeatus` and `oreochromis_niloticus`; Compara shadows `betta_splendens` and `gasterosteus_aculeatus`
(its other Ensembl copies of project species, the `validation` ones, were removed from the trees when they were prepared).
"""
from pathlib import Path

from .catalog import Catalog
from .master_tree import canonical

MAPPINGS = {"TreeFam": ("TreeFam9", "protein_to_species.txt"),
            "Compara": ("Compara116", "prepared/protein_to_species.txt"),
            "OrthoFinder": None}


class LeafSpecies:
    def __init__(self, reference, prefixes, project):
        self.reference = reference            # reference leaf -> canonical species
        self.prefixes = prefixes              # catalog prefix -> canonical species
        self.project = frozenset(project)     # canonical names of the project species

    @classmethod
    def load(cls, layout, part, catalog=None, only=None):
        """The mapping for one part of the release. `only`: a set of leaf names to keep from the (large) mapping file,
        to bound memory when only some families are read."""
        if part not in MAPPINGS:
            raise KeyError(f"unknown part {part!r} (one of {', '.join(MAPPINGS)})")
        catalog = catalog or Catalog.load(layout)
        prefixes = {p: canonical(s.species) for p, s in catalog.organisms.items()}
        project = {canonical(s.species) for s in catalog.species(role="project")}
        reference = {}
        if MAPPINGS[part]:
            folder, name = MAPPINGS[part]
            path = Path(layout.family_resource(folder).local) / name
            if not path.is_file():
                raise FileNotFoundError(f"{path} not found (run `cichlid-data check --fetch`)")
            keep = set(only) if only is not None else None
            with open(path, encoding="utf-8") as fh:
                for n, line in enumerate(fh, start=1):
                    line = line.rstrip("\r\n")
                    if not line:
                        continue
                    f = line.split("\t")
                    if len(f) != 2 or not f[0] or not f[1]:
                        raise ValueError(f"{path.name} line {n}: expected protein<TAB>species")
                    if keep is not None and f[0] not in keep:
                        continue
                    sp = canonical(f[1])
                    if reference.setdefault(f[0], sp) != sp:
                        raise ValueError(f"{path.name} line {n}: {f[0]} is given two species")
        return cls(reference, prefixes, project)

    def species(self, leaf):
        """The canonical species of a leaf, or None if it is neither a reference leaf nor one of our genes."""
        sp = self.reference.get(leaf)
        if sp is not None:
            return sp
        prefix, sep, rest = leaf.partition("_")
        return self.prefixes.get(prefix) if sep and rest else None

    def is_project_gene(self, leaf):
        """True for one of our genes: not a reference leaf, and `<prefix>_<id>` with a catalog prefix."""
        if leaf in self.reference:
            return False
        prefix, sep, rest = leaf.partition("_")
        return bool(sep and rest and prefix in self.prefixes)

    def of_leaves(self, leaves):
        """{leaf: species} for these leaves. KeyError names the leaves that no mapping gives a species."""
        out, unknown = {}, []
        for leaf in leaves:
            sp = self.species(leaf)
            if sp is None:
                unknown.append(leaf)
            else:
                out[leaf] = sp
        if unknown:
            more = f" (and {len(unknown) - 5} more)" if len(unknown) > 5 else ""
            raise KeyError(f"no species for {len(unknown)} leaves: {', '.join(unknown[:5])}{more}")
        return out

    def reference_species(self):
        return set(self.reference.values())

    def shadowed(self):
        """Reference species that are also project species (see the module notes)."""
        return self.reference_species() & self.project


def reference_species_in_mapping(layout, part):
    """The set of canonical species a part's mapping file names (reads the whole file, but keeps only the names)."""
    folder, name = MAPPINGS[part]
    path = Path(layout.family_resource(folder).local) / name
    out = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) == 2 and f[1]:
                out.add(canonical(f[1]))
    return out
