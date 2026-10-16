#!/usr/bin/env python3
"""
Draw "Where cichlids sit among fishes and other vertebrates": a group-level tree from tunicates to
cichlids, with whole-genome duplications, shaded blocks for the major groups, and colors showing
which groups contain project species, which are references already in our trees, and which are
missing.

    python make_species_phylogeny.py                   # writes species_phylogeny.svg (and .png if it can)
    python make_species_phylogeny.py --outdir docs     # somewhere else
    python make_species_phylogeny.py --no-png          # SVG only

The PNG needs a way to render SVG: Playwright (pip install playwright; playwright install chromium)
or CairoSVG (pip install cairosvg). Without either, only the SVG is written.

HOW TO EDIT
-----------
Everything about content is in the TREE section below; the rest is layout.

* A group (tip) is  T(name, examples, style, wgd).
  style: "project"  groups containing project species (blue, bold)
         "ref"      reference groups already in our trees (black)  [default]
         "missing"  groups with no annotated genome in our trees yet (red, dashed branch)
         "excluded" groups deliberately left out (grey)
         "mammal"   like "ref", but highlighted
* A clade is  N(name, [children...], wgd=...). The name only matters if a shaded block uses it
  (BLOCKS below). More than two children draws an unresolved fork.
* wgd puts a whole-genome-duplication diamond on the branch leading to that group, labeled with
  the text given ("3R", "carps", "?" ...). For tips, TIP_OFFSET staggers markers on nearby rows.
* BLOCKS lists the shaded groups: clade name, label, fill color and opacity. Nested blue blocks
  stack, so the shading deepens toward the cichlids.

Example: when a Ptychochrominae genome is added, change its style from "missing" to "project".
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from xml.sax.saxutils import escape

INK, GREY, RED, LAKE, GOLD, BG = "#16303B", "#6C7A80", "#C8435A", "#14627A", "#E2A31B", "#FFFFFF"
TITLE = "Where cichlids sit among fishes and other vertebrates"


class N:
    def __init__(self, name="", kids=(), ex="", style="ref", wgd=""):
        self.name, self.kids, self.ex, self.style, self.wgd = name, list(kids), ex, style, wgd


def T(name, ex, style="ref", wgd=""):
    return N(name, ex=ex, style=style, wgd=wgd)


# ================================================================================ TREE
cichlids = N("Cichlids", [
    T("Etroplinae", "India, Sri Lanka, Madagascar: Etroplus suratensis (Esur)", "project"),
    N("", [T("Ptychochrominae", "Madagascar (Paratilapia, Ptychochromis): second-earliest cichlid branch", "missing"),
           N("", [T("Neotropical cichlids", "Archocentrus centrarchus (Acen); Amphilophus", "project"),
                  T("African cichlids", "8 project species, from tilapia to the Lake Malawi haplochromines", "project")])])])
ovalentaria = N("Ovalentaria", [
    T("Atherinomorpha", "medaka, killifishes, guppy, platyfish"),
    T("Damselfishes", "clownfish, Acanthochromis"),
    T("Other ovalentarians", "blennies, mullets, silversides", "missing"),
    N("", [T("Pholidichthys", "convict blenny: the sister family of cichlids", "missing"), cichlids])])
percomorphs = N("Percomorphs", [
    T("Eupercaria", "stickleback, fugu, sea bass, seabream, wrasse"),
    T("Anabantaria", "betta, climbing perch"),
    T("Carangaria and others", "turbot, amberjack, seahorse"),
    ovalentaria])
teleosts = N("Teleosts", [
    T("Bonytongues", "elephantfish"),
    N("", [N("", [T("Herrings", "Atlantic herring, denticle herring"),
                  T("Carps, catfishes and relatives", "zebrafish (ortholog calls), cavefish, catfish; not carps", wgd="carps")]),
           N("", [T("Pikes and salmon", "northern pike; salmon and trout not used", wgd="salmonids"),
                  N("", [T("Cods", "Atlantic cod"), percomorphs])])])], wgd="3R")
rayfinned = N("Ray-finned fish", [
    T("Bichirs", "reedfish"),
    N("", [T("Sturgeons and paddlefish", "not used: several extra genome duplications", "excluded", wgd="several"),
           N("", [T("Gars and bowfin", "spotted gar (bowfin: to check)"), teleosts])])])
lobe = N("Lobe-finned fish and tetrapods", [
    T("Coelacanths", "coelacanth"),
    N("", [T("Amphibians", "western clawed frog"),
           N("", [T("Reptiles and birds", "anole, chicken, zebra finch"),
                  T("Mammals", "human and mouse (ortholog calls), rat, dog, cow, opossum, platypus", "mammal")])])])
TREE = N("", [
    T("Tunicates", "Ciona"),
    N("Vertebrates", [T("Lampreys", "sea lamprey", wgd="?"),
                      N("Jawed vertebrates", [T("Cartilaginous fishes", "elephant shark"), N("", [lobe, rayfinned])], wgd="2R")],
      wgd="1R")])

# shaded groups: clade name, label, fill, opacity (blue blocks nest and deepen)
BLOCKS = [("Lobe-finned fish and tetrapods", "Lobe-finned fish and tetrapods", "#E9C98A", 0.30),
          ("Ray-finned fish", "Ray-finned fish", LAKE, 0.07),
          ("Teleosts", "Teleosts", LAKE, 0.07),
          ("Percomorphs", "Percomorphs", LAKE, 0.07),
          ("Ovalentaria", "Ovalentaria", LAKE, 0.08),
          ("Cichlids", "Cichlids", LAKE, 0.10)]
# distance of tip markers from the tips (stagger markers on neighboring rows)
TIP_OFFSET = {"carps": 140, "salmonids": 70, "?": 300, "several": 220}
LEGEND = [("text", LAKE, "700", "Blue bold: groups with project species"),
          ("text", INK, "600", "Black: reference groups already in our trees (examples in grey)"),
          ("line", RED, None, "Missing groups: no annotated genome in our trees yet"),
          ("text", GREY, "600", "Grey: groups we deliberately leave out (extra genome duplications complicate orthology)"),
          ("diamond", None, None, "Whole-genome duplication: 1R and 2R (early vertebrates), 3R (teleosts), "
                                  "and later ones in salmonids, carps and sturgeons; ? = debated")]
FOOTNOTE = ("Branch lengths are not to scale. Forks with more than two branches mark relationships that are "
            "still debated (among percomorph groups, and among ovalentarian families).")

# ================================================================================ LAYOUT
TOP, ROW, X0, STEP = 64, 27, 28, 24      # first row, row height, root x, depth step
LABEL_WIDTH = 600                         # room for tip labels; the blocks end here


def build_svg(tree):
    tips = []

    def walk(n, d):
        n.d = d
        if not n.kids:
            tips.append(n)
        for k in n.kids:
            walk(k, d + 1)
    walk(tree, 0)
    xt = X0 + max(t.d for t in tips) * STEP
    for i, t in enumerate(tips):
        t.y = TOP + i * ROW

    def place(n):
        n.x = X0 + n.d * STEP if n.kids else xt
        if n.kids:
            for k in n.kids:
                place(k)
            n.y = (n.kids[0].y + n.kids[-1].y) / 2
    place(tree)

    out = []

    def line(x1, y1, x2, y2, c=INK, dash=False, w=2):
        d = ' stroke-dasharray="5 4"' if dash else ""
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{c}" stroke-width="{w}"{d} '
                   f'stroke-linecap="round"/>')

    def color(style):
        return RED if style == "missing" else LAKE if style == "project" else GREY if style == "excluded" else INK

    def draw(n):
        if not n.kids:
            return
        line(n.x, n.kids[0].y, n.x, n.kids[-1].y)
        for k in n.kids:
            miss = not k.kids and k.style == "missing"
            line(n.x, k.y, k.x, k.y, RED if miss else INK, dash=miss)
            if miss:
                line(n.x, k.y, n.x, n.y, RED, dash=True, w=2.4)
            if k.wgd:
                mx = (n.x + k.x) / 2 if k.kids else xt - TIP_OFFSET.get(k.wgd, 60)
                out.append(f'<path d="M{mx} {k.y - 6} L{mx + 6} {k.y} L{mx} {k.y + 6} L{mx - 6} {k.y} Z" fill="{GOLD}" '
                           f'stroke="{INK}" stroke-width="1"/>')
                out.append(f'<text x="{mx}" y="{k.y + 19}" font-size="10.5" font-weight="700" fill="#9A6B00" '
                           f'text-anchor="middle">{escape(k.wgd)}</text>')
            draw(k)
    draw(tree)
    line(X0 - 12, tree.y, X0, tree.y)

    lx = xt + 10
    for t in tips:
        c = color(t.style)
        w = "700" if t.style in ("project", "mammal") else "600"
        out.append(f'<text x="{lx}" y="{t.y + 4.5}" font-size="13.5" fill="{c}"><tspan font-weight="{w}">{escape(t.name)}</tspan>'
                   f'<tspan dx="8" font-size="12" font-style="italic" fill="{c if t.style == "missing" else GREY}">'
                   f'{escape(t.ex)}</tspan></text>')

    # shaded blocks and the mammal highlight go underneath everything
    right = lx + LABEL_WIDTH

    def find(n, name):
        if n.name == name:
            return n
        for k in n.kids:
            r = find(k, name)
            if r:
                return r

    def leaves(n):
        return [n] if not n.kids else [l for k in n.kids for l in leaves(k)]
    under = []
    for name, label, fill, op in BLOCKS:
        n = find(tree, name)
        if n is None:
            raise SystemExit(f"BLOCKS refers to '{name}', which is not a clade name in the tree")
        lv = leaves(n)
        y1, y2 = lv[0].y - ROW / 2 + 1.5, lv[-1].y + ROW / 2 - 1.5
        x1 = n.x - 9
        under.append(f'<rect x="{x1}" y="{y1}" width="{right - x1}" height="{y2 - y1}" rx="6" fill="{fill}" fill-opacity="{op}"/>')
        under.append(f'<text x="{right - 8}" y="{lv[0].y + 4.5}" font-size="12" font-weight="700" '
                     f'fill="{"#8A5A00" if fill != LAKE else LAKE}" text-anchor="end" letter-spacing="0.3">{escape(label)}</text>')
    for t in tips:
        if t.style == "mammal":
            under.append(f'<rect x="{xt - 6}" y="{t.y - 12}" width="{right - 8 - (xt - 6)}" height="24" rx="4" '
                         f'fill="#E9C98A" fill-opacity="0.45"/>')
    out[:0] = under

    ly = TOP + len(tips) * ROW + 14
    for i, (kind, c, w, txt) in enumerate(LEGEND):
        yy = ly + i * 20
        if kind == "text":
            out.append(f'<text x="30" y="{yy + 4.5}" font-size="12.5" fill="{c}" font-weight="{w}">{escape(txt)}</text>')
        elif kind == "line":
            line(30, yy, 62, yy, RED, dash=True, w=2.4)
            out.append(f'<text x="72" y="{yy + 4.5}" font-size="12.5" fill="{RED}">{escape(txt)}</text>')
        else:
            out.append(f'<path d="M38 {yy - 6} L44 {yy} L38 {yy + 6} L32 {yy} Z" fill="{GOLD}" stroke="{INK}" stroke-width="1"/>'
                       f'<text x="54" y="{yy + 4.5}" font-size="12.5" fill="{INK}">{escape(txt)}</text>')
    out.append(f'<text x="30" y="{ly + len(LEGEND) * 20 + 8}" font-size="11.5" fill="{GREY}">{escape(FOOTNOTE)}</text>')
    width, height = right + 16, ly + len(LEGEND) * 20 + 22
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
            f'font-family="Helvetica, Arial, sans-serif"><rect width="{width}" height="{height}" fill="{BG}"/>'
            f'<text x="28" y="32" font-size="19" font-weight="700" fill="{INK}">{escape(TITLE)}</text>'
            + "".join(out) + "</svg>")


def write_png(svg_path, png_path, scale=2):
    """Render the SVG to PNG with Playwright or CairoSVG, whichever is installed. -> True if written."""
    try:
        from playwright.sync_api import sync_playwright
        svg = svg_path.read_text()
        with sync_playwright() as p:
            b = p.chromium.launch()
            pg = b.new_page(device_scale_factor=scale, viewport={"width": 1400, "height": 1100})
            pg.set_content(f'<html><body style="margin:0">{svg}</body></html>')
            pg.locator("svg").screenshot(path=str(png_path))
            b.close()
        return True
    except Exception:
        pass
    try:
        import cairosvg
        cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), scale=scale)
        return True
    except Exception:
        pass
    if shutil.which("rsvg-convert"):
        r = subprocess.run(["rsvg-convert", "-z", str(scale), "-o", str(png_path), str(svg_path)])
        return r.returncode == 0
    return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=Path("."))
    ap.add_argument("--name", default="species_phylogeny", help="file name without extension")
    ap.add_argument("--no-png", action="store_true")
    args = ap.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    svg_path = args.outdir / f"{args.name}.svg"
    svg_path.write_text(build_svg(TREE))
    print(f"Wrote {svg_path}", file=sys.stderr)
    if not args.no_png:
        png_path = args.outdir / f"{args.name}.png"
        if write_png(svg_path, png_path):
            print(f"Wrote {png_path}", file=sys.stderr)
        else:
            print("No SVG renderer found (Playwright, CairoSVG or rsvg-convert), so no PNG; the SVG is complete.",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
