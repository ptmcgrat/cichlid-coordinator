"""The master species tree (config/species_tree_master.nwk): one rooted tree over every species that is a leaf in any
method's trees. The homology decider roots gene trees with it. Topology only: no branch lengths.

    T = load_master_tree()               # a newick.Node
    canonical("Canis_familiaris")        # -> "canis_lupus_familiaris": the leaf name of a species

`build()` makes it from three sources, and the committed file must equal what `build()` gives (a test checks):

  * the tree agreed for OrthoFinder (config/species_tree_orthofinder.nwk, 35 species) is the backbone for the fish and the
    main vertebrate branches, and is kept exactly: the master, pruned to those 35 species, equals it;
  * TreeFam 9's species tree (config/species_tree_treefam9.phy, from the TreeFam FTP archive) gives everything else
    (mammals, birds, invertebrates, ...). It is a taxonomy tree with polytomies and an older fish classification (it has no
    Ovalentaria), so its 8 fish are replaced by the agreed teleosts; its other clades are kept;
  * the fish that only the Compara panel (or TreeFam) has are attached: a TreeFam fish next to its TreeFam relative, a
    Compara fish at its ladder level, as an unresolved polytomy at that level's clade. Polytomies do little harm: the
    tree only chooses the root of a gene tree, and duplications are called by species overlap.

The shark, gar, reedfish and amphioxus (not in TreeFam's tree) go where the agreed tree has them.
`check()` tests the invariants: every species once, the agreed tree preserved, every ladder level a clade.
"""
import sys
import tomllib
from pathlib import Path

from . import newick
from .clades import LADDER, Clades
from .config import config_dir

FILE = "species_tree_master.nwk"
TREEFAM_FILE = "species_tree_treefam9.phy"
AGREED_FILE = "species_tree_orthofinder.nwk"
TELEOST_LEVELS = ("teleost", "percomorph", "ovalentaria", "cichlid")

# alternative spellings of one species -> the leaf name used in the master tree
SYNONYMS = {
    "canis_familiaris": "canis_lupus_familiaris",
    "spermophilus_tridecemlineatus": "ictidomys_tridecemlineatus",
    "caenorhabditis_briggsae_af16": "caenorhabditis_briggsae",
    "drosophila_pseudoobscura_pseudoobscura": "drosophila_pseudoobscura",
    "pediculus_humanus_corporis": "pediculus_humanus",
    "proterospongia": "proterospongia_sp",
    "saccharomyces_cerevisiae_s288c": "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe_972h-": "schizosaccharomyces_pombe",
}


def canonical(name):
    n = str(name).strip().lower().replace(" ", "_")
    return SYNONYMS.get(n, n)


def master_tree_path(cfg=None):
    return Path(cfg or config_dir()) / FILE


def load_master_tree(cfg=None):
    return newick.parse(master_tree_path(cfg).read_text())


# ------------------------------------------------------------------ small tree helpers
def _copy(node, rename=None):
    new = newick.Node(rename(node.name) if (rename and not node.children) else node.name)
    for c in node.children:
        new.add(_copy(c, rename))
    return new


def _leafset(node):
    return {l.name for l in node.leaves()}


def _find(root, name):
    return next(n for n in root.iter() if n.name == name)


def _mrca(root, names):
    names = set(names)
    best = None
    for n in root.iter():
        ls = _leafset(n)
        if names <= ls and (best is None or len(ls) < best[0]):
            best = (len(ls), n)
    return best[1]


def _replace(node, new):
    """Put `new` where `node` is in its parent (the caller then adds `node` below `new`, if it wants it kept)."""
    parent = node.parent
    parent.children[parent.children.index(node)] = new
    new.parent, node.parent = parent, None


def _new(name, *children):
    n = newick.Node(name)
    for c in children:
        n.add(c)
    return n


def clade_sets(root, keep=None):
    """The clades of a rooted tree as a set of frozensets of leaf names, optionally restricted to `keep`."""
    keep = set(keep) if keep is not None else None
    out = set()
    for n in root.iter():
        if n.children:
            ls = _leafset(n)
            ls = ls & keep if keep is not None else ls
            if len(ls) >= 2 and (keep is None or len(ls) < len(keep)):
                out.add(frozenset(ls))
    return out


# ------------------------------------------------------------------ sources and build
def sources(cfg=None):
    """-> (agreed tree, TreeFam tree, the species that must be leaves, clade table, prefix -> species name).
    All names are canonical and lower case."""
    cfg = Path(cfg or config_dir())
    sets = tomllib.loads((cfg / "species_sets.toml").read_text())
    prefix_name = {o["prefix"]: canonical(o["species"]) for o in sets["organism"]}
    agreed = _copy(newick.parse((cfg / AGREED_FILE).read_text()), lambda p: prefix_name[p])
    treefam = _copy(newick.parse((cfg / TREEFAM_FILE).read_text()), canonical)
    wanted = _leafset(agreed) | _leafset(treefam)
    for ln in (cfg / "compara_panel.tsv").read_text().splitlines():            # the leaves of the pruned Compara trees
        f = ln.split("\t")
        if not ln.startswith("#") and len(f) >= 3 and f[0] != "species" and f[1] not in ("exclude", "validation"):
            wanted.add(canonical(f[0]))
    return agreed, treefam, wanted, Clades.load(cfg), prefix_name


def build(cfg=None):
    """The master tree (a newick.Node)."""
    agreed, tf, wanted, clades, prefix_name = sources(cfg)
    in_agreed = _leafset(agreed)
    level = lambda sp: clades.level(sp)

    # 1. the teleosts: the agreed teleost clade, with a node for each ladder level inside it
    teleosts = _copy(_mrca(agreed, [s for s in in_agreed if level(s) in TELEOST_LEVELS]))
    teleosts.name = "Teleostei"
    node_of = {"teleost": teleosts}
    for lv, label in (("percomorph", "Percomorpha"), ("ovalentaria", "Ovalentaria"), ("cichlid", "Cichlidae")):
        members = [s for s in in_agreed if clades.rank(level(s)) >= clades.rank(lv)]
        node_of[lv] = _mrca(teleosts, members)
        node_of[lv].name = label

    # 2. the fish that are not in the agreed tree
    extras = sorted(s for s in wanted - in_agreed if level(s) in TELEOST_LEVELS)
    tf_leaf = {l.name: l for l in tf.leaves()}
    for sp in extras:
        if sp in tf_leaf:                      # a TreeFam fish: next to its relatives in TreeFam's tree
            anc = tf_leaf[sp].parent
            while not (_leafset(anc) & in_agreed):
                anc = anc.parent
            anchor = _mrca(teleosts, _leafset(anc) & in_agreed)
            pair = newick.Node("")
            _replace(anchor, pair)
            pair.add(anchor)
            pair.add(newick.Node(sp))
        else:                                  # a Compara-only fish: at its ladder level
            node_of[level(sp)].add(newick.Node(sp))

    # 3. TreeFam's tree, with its fish replaced, and the taxa that TreeFam lacks put where the agreed tree has them
    master = _copy(tf)
    gar, reedfish, shark, amphioxus = (prefix_name[p] for p in ("Locu", "Ecal", "Cmil", "Bflo"))
    actinopterygii = _new("Actinopterygii", newick.Node(reedfish), _new("", newick.Node(gar), teleosts))
    _replace(_find(master, "Clupeocephala"), actinopterygii)
    euteleostomi = _find(master, "Euteleostomi")
    gnathostomata = newick.Node("Gnathostomata")
    _replace(euteleostomi, gnathostomata)
    gnathostomata.add(newick.Node(shark))
    gnathostomata.add(euteleostomi)
    chordata = _find(master, "Chordata")
    olfactores = _new("Olfactores", *list(chordata.children))
    chordata.children = []
    chordata.add(newick.Node(amphioxus))
    chordata.add(olfactores)
    return master


# ------------------------------------------------------------------ checks
def check(tree, cfg=None):
    """-> problems with a master tree: names, the species it must hold, the agreed tree, the ladder."""
    agreed, tf, wanted, clades, _ = sources(cfg)
    problems, leaves = [], [l.name for l in tree.leaves()]
    dup = sorted({x for x in leaves if leaves.count(x) > 1})
    if dup:
        problems.append(f"leaves that appear twice: {', '.join(dup)}")
    ls = set(leaves)
    if wanted - ls:
        problems.append(f"species missing from the tree: {', '.join(sorted(wanted - ls))}")
    if ls - wanted:
        problems.append(f"leaves that no method has: {', '.join(sorted(ls - wanted))}")
    unknown = sorted(x for x in ls if x not in clades)
    if unknown:
        problems.append(f"leaves with no row in species_clades.tsv: {', '.join(unknown)}")
    if len(tree.children) != 2:
        problems.append(f"the root has {len(tree.children)} children; it must be rooted (two)")
    if any(len(n.children) == 1 for n in tree.iter()):
        problems.append("the tree has a node with a single child")
    agreed_species = _leafset(agreed)
    if not agreed_species <= ls:
        problems.append("the agreed OrthoFinder species are not all in the tree")
    elif clade_sets(tree, agreed_species) != clade_sets(agreed, agreed_species):
        problems.append("pruned to the 35 OrthoFinder species the tree is not the agreed tree")
    if not problems:
        sets = clade_sets(tree)
        everything = frozenset(ls)
        for lv in LADDER:
            members = frozenset(canonical(s) for s in clades.within(lv)) & everything
            if members and members != everything and members not in sets:
                problems.append(f"the species at level {lv} and inside it do not form a clade")
    return problems


def main(argv=None):
    """python -m cichlid_data.master_tree [--write]: build the tree, check it, and (with --write) save it to config/."""
    argv = list(sys.argv[1:] if argv is None else argv)
    tree = build()
    problems = check(tree)
    for p in problems:
        print("problem:", p)
    text = newick.write(tree, internal_names=True) + "\n"
    if "--write" in argv and not problems:
        master_tree_path().write_text(text)
        print(f"wrote {master_tree_path()}")
    print(f"{len(tree.leaves())} species, {sum(1 for n in tree.iter() if n.children)} internal nodes, {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
