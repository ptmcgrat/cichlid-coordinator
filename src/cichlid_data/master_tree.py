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
  * the fish that only the Compara panel (or TreeFam) has are placed by config/species_tree_placements.tsv: one row per
    fish, with its published source and a support level. A placement is resolved only where a published phylogeny
    supports it; where the published work disagrees or is silent, the fish joins an unresolved polytomy on purpose
    (an invented resolution would turn a speciation into a duplication in the homology decider, which reconciles
    gene trees with this tree). A fish with no row falls back to the older rule: next to its TreeFam relative, or as a
    child of its ladder level's clade.

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
PLACEMENTS_FILE = "species_tree_placements.tsv"
PLACEMENT_COLUMNS = ["species", "op", "anchor", "group", "support", "source", "note"]
PLACEMENT_OPS = ("sister_of", "join_group", "group_under", "child_of")
SUPPORT = ("robust", "contested", "by_default")
LABELS = ("Teleostei", "Percomorpha", "Ovalentaria", "Cichlidae")        # nodes that an anchor can name as @Label
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
    # the owner's 2018 cichlid entries in TreeFam's mapping file: older names of project species. They are not leaves in
    # any tree; the names are folded in so that none is ever an orphan (the meaning of the _b is to confirm)
    "astatotilapia_burtoni": "haplochromis_burtoni",
    "metriaclima_zebra": "maylandia_zebra",
    "metriaclima_zebra_b": "maylandia_zebra",
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


# ------------------------------------------------------------------ placements
def parse_placements(text):
    """The text of species_tree_placements.tsv -> a list of dicts (comments and blank lines skipped)."""
    rows, header = [], False
    for n, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if not header:
            if f != PLACEMENT_COLUMNS:
                raise ValueError(f"{PLACEMENTS_FILE} line {n}: the header must be {', '.join(PLACEMENT_COLUMNS)}")
            header = True
            continue
        if len(f) != len(PLACEMENT_COLUMNS):
            raise ValueError(f"{PLACEMENTS_FILE} line {n}: expected {len(PLACEMENT_COLUMNS)} columns, found {len(f)}")
        rows.append(dict(zip(PLACEMENT_COLUMNS, (x.strip() for x in f))))
    if not header:
        raise ValueError(f"{PLACEMENTS_FILE}: no header line")
    return rows


def load_placements(cfg=None):
    path = Path(cfg or config_dir()) / PLACEMENTS_FILE
    return parse_placements(path.read_text()) if path.is_file() else []


def _anchor_tokens(anchor):
    return [a.strip() for a in anchor.split(",") if a.strip()]


def placement_problems(cfg=None):
    """Problems with the placement table: its columns, the vocabulary, that each row is a fish missing from the agreed
    tree, that its anchors exist when it is applied, and that groups are used consistently."""
    cfg = Path(cfg or config_dir())
    try:
        rows = load_placements(cfg)
    except ValueError as e:
        return [str(e)]
    agreed, _tf, wanted, clades, _ = sources(cfg)
    in_agreed = _leafset(agreed)
    extras = {sp for sp in wanted - in_agreed if clades.level(sp) in TELEOST_LEVELS}
    problems, seen, known, group_anchor = [], set(), set(in_agreed), {}
    for i, r in enumerate(rows, start=1):
        tag = f"{PLACEMENTS_FILE} row {i} ({r['species']})"
        sp, op = r["species"], r["op"]
        if sp in seen:
            problems.append(f"{tag}: the species appears twice")
        elif sp not in extras:
            problems.append(f"{tag}: not a fish that is missing from the agreed tree (the agreed tree cannot be changed)")
        if op not in PLACEMENT_OPS:
            problems.append(f"{tag}: op {op!r} must be one of {', '.join(PLACEMENT_OPS)}")
        if r["support"] not in SUPPORT:
            problems.append(f"{tag}: support {r['support']!r} must be one of {', '.join(SUPPORT)}")
        if not r["source"]:
            problems.append(f"{tag}: no source")
        tokens = _anchor_tokens(r["anchor"])
        if not tokens:
            problems.append(f"{tag}: no anchor")
        for t in tokens:
            if t.startswith("@"):
                if t[1:] not in LABELS:
                    problems.append(f"{tag}: anchor {t} must be one of {', '.join('@' + x for x in LABELS)}")
            elif t not in known:
                problems.append(f"{tag}: anchor {t} is not in the agreed tree or placed by an earlier row")
        if op in ("join_group", "group_under") and not r["group"]:
            problems.append(f"{tag}: {op} needs a group name")
        if op == "child_of" and r["group"]:
            problems.append(f"{tag}: child_of makes no new node, so it takes no group name")
        if op in ("join_group", "group_under") and r["group"]:
            if group_anchor.setdefault(r["group"], (op, r["anchor"])) != (op, r["anchor"]) and op == "group_under":
                problems.append(f"{tag}: group {r['group']} is anchored differently by an earlier row")
        seen.add(sp)
        known.add(sp)
    return problems


def _anchor_node(root, labels, anchor, where):
    if anchor.startswith("@"):
        if anchor[1:] not in labels:
            raise ValueError(f"{where}: unknown anchor {anchor}")
        return labels[anchor[1:]]
    names = _anchor_tokens(anchor)
    missing = [x for x in names if x not in _leafset(root)]
    if missing or not names:
        raise ValueError(f"{where}: anchor species not in the tree: {', '.join(missing) or '(none given)'}")
    return _mrca(root, names)


def apply_placements(root, labels, rows):
    """Place each row's fish in `root` (modified in place), in order. `labels` maps a label to a node, for @Label anchors.
    -> the set of species placed."""
    placed, groups, present = set(), {}, _leafset(root)
    for i, r in enumerate(rows, start=1):
        where = f"{PLACEMENTS_FILE} row {i} ({r['species']})"
        sp, op, group = r["species"], r["op"], r["group"]
        if sp in present:
            raise ValueError(f"{where}: the species is already in the tree")
        anchor = _anchor_node(root, labels, r["anchor"], where)
        leaf = newick.Node(sp)
        if op == "sister_of":
            if anchor.parent is None:
                raise ValueError(f"{where}: cannot be the sister of the root")
            pair = newick.Node(group)
            _replace(anchor, pair)
            pair.add(anchor)
            pair.add(leaf)
        elif op == "join_group":
            if group in groups:
                if anchor.parent is not groups[group]:
                    raise ValueError(f"{where}: its anchor is not inside group {group}")
                groups[group].add(leaf)
            else:
                if anchor.parent is None:
                    raise ValueError(f"{where}: cannot group with the root")
                g = newick.Node(group)
                _replace(anchor, g)
                g.add(anchor)
                g.add(leaf)
                groups[group] = g
        elif op == "group_under":
            if group in groups:
                if groups[group].parent is not anchor:
                    raise ValueError(f"{where}: group {group} hangs under a different node")
                groups[group].add(leaf)
            else:
                g = newick.Node(group)
                g.add(leaf)
                anchor.add(g)
                groups[group] = g
        elif op == "child_of":
            anchor.add(leaf)
        else:
            raise ValueError(f"{where}: unknown op {op!r}")
        present.add(sp)
        placed.add(sp)
    return placed


def polytomies(tree):
    """-> [(name, number of children, leaf set)] for every node with more than two children, largest first."""
    out = [(n.name, len(n.children), frozenset(_leafset(n))) for n in tree.iter() if len(n.children) > 2]
    return sorted(out, key=lambda t: (-len(t[2]), t[0]))


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

    # 2. the fish that are not in the agreed tree: placed by the table, else by the older rule
    extras = sorted(s for s in wanted - in_agreed if level(s) in TELEOST_LEVELS)
    problems = placement_problems(cfg)
    if problems:
        raise ValueError("the placement table is not valid: " + "; ".join(problems[:3]))
    labels = {"Teleostei": teleosts, "Percomorpha": node_of["percomorph"], "Ovalentaria": node_of["ovalentaria"],
              "Cichlidae": node_of["cichlid"]}
    placed = apply_placements(teleosts, labels, load_placements(cfg))
    tf_leaf = {l.name: l for l in tf.leaves()}
    for sp in extras:
        if sp in placed:
            continue
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
    problems += placement_problems(cfg)
    fish = next((n for n in tree.iter() if n.name == "Teleostei"), None)         # names are looked up inside the fish clade;
    names = [n.name for n in fish.iter() if n.children and n.name] if fish else []      # TreeFam's own names elsewhere are not ours
    twice = sorted({x for x in names if names.count(x) > 1})
    if twice:
        problems.append(f"internal node names in the fish clade that appear twice: {', '.join(twice)}")
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
