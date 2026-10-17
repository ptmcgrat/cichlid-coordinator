"""
Minimal Newick/NHX tree handling for preparing TreeFam reference trees.
No dependencies; only what Stage 5 needs: parse, write, prune, collapse unary nodes,
unroot, and resolve multifurcations into a strictly binary unrooted tree.
"""
import re

NHX_COMMENT = re.compile(r"\[[^\]]*\]")


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name="", length=None):
        self.name, self.length, self.children, self.parent = name, length, [], None

    def add(self, child):
        child.parent = self
        self.children.append(child)
        return child

    def is_leaf(self):
        return not self.children

    def iter(self):
        """Pre-order, iterative (TreeFam trees can be thousands of levels deep)."""
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(reversed(n.children))

    def leaves(self):
        return [n for n in self.iter() if n.is_leaf()]


def parse(text):
    """Parse Newick (NHX comments stripped). Internal labels kept as names."""
    s = NHX_COMMENT.sub("", text).strip()
    if not s.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    s = s[:-1]
    root = Node()
    cur, i, n = root, 0, len(s)
    expect_label = False

    def read_token(i):
        if i < n and s[i] == "'":
            j = s.index("'", i + 1)
            return s[i + 1:j], j + 1
        j = i
        while j < n and s[j] not in "(),:;":
            j += 1
        return s[i:j].strip(), j

    while i < n:
        c = s[i]
        if c == "(":
            cur = cur.add(Node())
            i += 1
        elif c == ",":
            if cur.parent is None:
                raise ValueError("unbalanced parentheses")
            cur = cur.parent.add(Node())
            i += 1
        elif c == ")":
            if cur.parent is None:
                raise ValueError("unbalanced parentheses")
            cur = cur.parent
            i += 1
        elif c == ":":
            tok, i = read_token(i + 1)
            cur.length = float(tok) if tok else None
        elif c.isspace():
            i += 1
        else:
            tok, i = read_token(i)
            cur.name = tok
    if cur is not root:
        raise ValueError("unbalanced parentheses")
    return root


def write(root, internal_names=False):
    """Iterative writer (deep trees would overflow recursion)."""
    out = []
    stack = [("open", root)]
    while stack:
        kind, node = stack.pop()
        if kind == "comma":
            out.append(",")
        elif kind == "close":
            out.append(")")
            out.append(_label(node, internal_names))
        elif node.children:
            out.append("(")
            stack.append(("close", node))
            for idx, ch in enumerate(reversed(node.children)):
                stack.append(("open", ch))
                if idx < len(node.children) - 1:
                    stack.append(("comma", None))
        else:
            out.append(_label(node, True))
    return "".join(out) + ";"


def _label(node, with_name):
    s = node.name if (with_name and node.name) else ""
    if node.length is not None and node.parent is not None:
        s += f":{node.length:.6g}"
    return s


def prune(root, keep):
    """Remove leaves not in `keep`; then drop empty internal nodes. Returns root."""
    for leaf in root.leaves():
        if leaf.name not in keep:
            _detach(leaf)
    changed = True
    while changed:              # internal nodes left childless
        changed = False
        for node in list(root.iter()):
            if node is not root and not node.children and node.name not in keep:
                _detach(node)
                changed = True
    return collapse_unary(root)


def _sum(a, b):
    if a is None and b is None:
        return None
    return (a or 0.0) + (b or 0.0)


def _detach(node):
    p = node.parent
    if p is not None:
        p.children.remove(node)
        node.parent = None


def collapse_unary(root):
    """Remove internal nodes with one child, summing branch lengths. Returns (new) root."""
    for node in list(root.iter()):
        if node is not root and len(node.children) == 1:
            child, p = node.children[0], node.parent
            child.length = _sum(child.length, node.length)
            idx = p.children.index(node)
            p.children[idx] = child
            child.parent = p
    while len(root.children) == 1:
        root = root.children[0]
        root.parent, root.length = None, None
    return root


def unroot(root):
    """Turn a bifurcating root into a trifurcation (unrooted representation)."""
    if len(root.children) != 2:
        return root
    a, b = root.children
    inner, other = (a, b) if a.children else (b, a)
    if not inner.children:          # two-leaf tree: nothing to do
        return root
    other.length = _sum(other.length, inner.length)
    inner.parent, inner.length = None, None
    inner.add(other)
    return inner


def resolve_polytomies(root, min_length=1e-6):
    """Make every internal node binary (root: exactly 3 children) with near-zero branches."""
    for node in list(root.iter()):
        limit = 3 if node is root else 2
        while len(node.children) > limit:
            a, b = node.children[-2], node.children[-1]
            node.children = node.children[:-2]
            joint = node.add(Node(length=min_length))
            joint.add(a)
            joint.add(b)
    return root


def fill_lengths(root, default=1e-3):
    for node in root.iter():
        if node.parent is not None and (node.length is None or node.length < 0):
            node.length = default
    return root


def prepare_reference(text, keep):
    """TreeFam NHX text -> (binary unrooted Newick string, report dict) restricted to `keep`."""
    root = parse(text)
    names = [l.name for l in root.leaves()]
    dupes = sorted({x for x in names if names.count(x) > 1}) if len(names) != len(set(names)) else []
    if dupes:
        raise ValueError(f"duplicate leaf names in reference tree, e.g. {dupes[:3]}")
    dropped = sorted(set(names) - set(keep))
    root = prune(root, set(keep))
    root = unroot(root)
    root = resolve_polytomies(root)
    root = fill_lengths(root)
    for n in root.iter():           # internal labels (supports) would confuse some tools
        if n.children:
            n.name = ""
    leaves = [l.name for l in root.leaves()]
    if len(leaves) < 4:
        raise ValueError(f"only {len(leaves)} reference leaves after pruning")
    return write(root), {"leaves": len(leaves), "pruned": dropped,
                         "missing_from_tree": sorted(set(keep) - set(leaves))}


# --------------------------------------------------------------------------- rooting
def _adjacency(root):
    """Undirected view: {node: [(neighbour, branch length)]}."""
    adj = {}
    for n in root.iter():
        adj.setdefault(n, [])
        if n.parent is not None:
            L = n.length or 0.0
            adj[n].append((n.parent, L))
            adj.setdefault(n.parent, []).append((n, L))
    return adj


def _farthest(adj, start):
    """-> (farthest leaf, distance, predecessor map) by iterative DFS."""
    dist, prev, stack = {start: 0.0}, {start: None}, [start]
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v not in dist:
                dist[v], prev[v] = dist[u] + L, u
                stack.append(v)
    leaves = [n for n in dist if n.is_leaf()]
    far = max(leaves, key=lambda n: (dist[n], n.name))
    return far, dist[far], prev


def _rebuild(adj, a, b, len_a, len_b):
    """New rooted tree with the root on edge (a, b): a at len_a below the root, b at len_b."""
    new_root = Node()
    mapping = {}
    for start, length, blocked in ((a, len_a, b), (b, len_b, a)):
        top = Node(start.name, length)
        new_root.add(top)
        mapping[start] = top
        stack = [(start, blocked)]
        while stack:
            u, came_from = stack.pop()
            for v, L in adj[u]:
                if v is came_from:
                    continue
                child = Node(v.name, L)
                mapping[u].add(child)
                mapping[v] = child
                stack.append((v, u))
    for n in new_root.iter():          # old internal labels are meaningless after re-rooting
        if n.children:
            n.name = ""
    return collapse_unary(new_root)


def midpoint_root(root):
    """Root at the midpoint of the longest leaf-to-leaf path. Falls back to the most balanced
    edge (by leaf count) when all branch lengths are zero, as the 2018 code did."""
    leaves = root.leaves()
    if len(leaves) < 3:
        return root
    adj = _adjacency(root)
    a, _, _ = _farthest(adj, leaves[0])
    b, total, prev = _farthest(adj, a)
    if total > 0:
        path = [b]
        while prev[path[-1]] is not None:
            path.append(prev[path[-1]])           # b ... a
        half, walked = total / 2.0, 0.0
        for u, v in zip(path, path[1:]):          # walk from b towards a
            L = next(l for w, l in adj[u] if w is v)
            if walked + L >= half:
                x = half - walked                 # distance from u to the root point
                return _rebuild(adj, u, v, x, L - x)
            walked += L
    # zero-length tree: split leaves as evenly as possible
    n_leaves = len(leaves)
    counts = {}
    for n in reversed(list(root.iter())):
        counts[n] = 1 if n.is_leaf() else sum(counts[c] for c in n.children)
    best = min((n for n in root.iter() if n.parent is not None),
               key=lambda n: (abs(counts[n] - n_leaves / 2), n.name or ""))
    L = best.length or 0.0
    return _rebuild(adj, best, best.parent, L / 2, L / 2)


def depths(root):
    """{node: distance from root}."""
    d = {root: 0.0}
    for n in root.iter():
        for c in n.children:
            d[c] = d[n] + (c.length or 0.0)
    return d


