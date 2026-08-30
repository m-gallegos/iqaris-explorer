"""fig_style.py -- shared publication ("house") style for the IQARIS figures.

A single, journal-quality look: Helvetica (Nimbus Sans), a colour-blind-safe Okabe--Ito
categorical palette, refined CPK element colours, thin dark-grey spines/ticks, sans-serif
math (STIX), embedded fonts.  Imported by make_results_figures.py and make_chemspace_map.py.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

INK   = "#1a1a1a"     # primary text
AXGREY = "#5a5a5a"    # (legacy) light spines / ticks
AXDARK = "#2b2b2b"    # bold spines / ticks
FAINT = "#9aa0a6"     # secondary annotation

def set_style():
    mpl.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Nimbus Sans", "TeX Gyre Heros", "Helvetica",
                            "Liberation Sans", "Arial", "DejaVu Sans"],
        "mathtext.fontset": "stixsans", "mathtext.default": "regular",
        "font.size": 7.4, "font.weight": "bold",
        "axes.labelsize": 8.5, "axes.labelweight": "bold",
        "axes.titlesize": 8.4, "axes.titleweight": "bold",
        "axes.titlepad": 4.0, "axes.labelpad": 2.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.6,
        "axes.linewidth": 1.1, "axes.edgecolor": AXDARK, "axes.labelcolor": INK,
        "text.color": INK, "xtick.color": AXDARK, "ytick.color": AXDARK,
        "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "xtick.major.width": 1.1, "ytick.major.width": 1.1,
        "xtick.minor.width": 0.7, "ytick.minor.width": 0.7,
        "xtick.major.size": 3.4, "ytick.major.size": 3.4,
        "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.axisbelow": True,
        "axes.grid": False, "grid.color": "#e9e9ee", "grid.linewidth": 0.5,
        "legend.frameon": False, "legend.handletextpad": 0.4,
        "legend.borderaxespad": 0.3, "legend.labelspacing": 0.3,
        "figure.dpi": 150, "savefig.dpi": 500, "savefig.bbox": "tight",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "lines.antialiased": True, "patch.linewidth": 0.0,
    })

# refined CPK element colours (tuned for white background; F is pale blue by house convention,
# kept distinct from N's royal blue, so F/Cl no longer collide as two greens)
ECOL = {"H": "#8a929b", "C": "#3b3b3b", "N": "#2f6fd0", "O": "#e23b2e",
        "F": "#56B4E9", "P": "#ef8a2b", "S": "#c8a01e", "Cl": "#199a57"}
EDARK = {"H": "#5f666e", "C": "#1f1f1f", "N": "#1d4a99", "O": "#a8261c",
         "F": "#2e86c1", "P": "#bd6710", "S": "#937205", "Cl": "#0e6e3c"}
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16, "Cl": 17}
ELEM_ORDER = ["H", "C", "N", "O", "F", "P", "S", "Cl"]

# Okabe--Ito colour-blind-safe qualitative palette
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "vermillion": "#D55E00", "sky": "#56B4E9", "yellow": "#E8C100",
      "purple": "#CC79A7", "grey": "#9aa0a6", "black": "#222222"}
OI_CYCLE = [OI["blue"], OI["orange"], OI["green"], OI["vermillion"], OI["sky"],
            OI["purple"], OI["yellow"], OI["black"], OI["grey"]]

_BOND_PRIO = {"C": 0, "N": 1, "O": 2, "S": 3, "P": 4, "F": 5, "Cl": 6, "H": 7}

def bond_label(e1, e2):
    a, b = (e1, e2) if _BOND_PRIO[e1] <= _BOND_PRIO[e2] else (e2, e1)
    return f"{a}–{b}"          # en-dash

def bond_hetero_elem(e1, e2):
    """The bond's defining (heteroatom) element.  C–H -> H, C–C -> C; the heavier/halogen
    heteroatom wins for mixed-heteroatom bonds (e.g. N–O -> O)."""
    hetero = [e for e in (e1, e2) if e not in ("C", "H")]
    if hetero:
        return max(hetero, key=lambda x: _BOND_PRIO[x])
    if "C" in (e1, e2) and "H" in (e1, e2):
        return "H"
    if e1 == "C":
        return "C"
    return "H"

def bond_hetero_color(e1, e2):
    """CPK colour keyed on the heteroatom, shared by Fig 4 panels B & C so a bond reads the same
    colour in both (red = O-bonds, blue = N-bonds, ...)."""
    return ECOL[bond_hetero_elem(e1, e2)]

# --- per-pair bond colour key + canonical bond set (shared by ALL manuscript figures).
#     Colour is a BIJECTION per pair: the eight carbon bonds take their
#     CPK heteroatom colour, the polar X–H (and other) pairs a distinct SHADE of the same heteroatom
#     family, so no colour ever denotes two different pairs, within or across panels/figures.
PAIRCOL = {
    "C–H": ECOL["H"],  "H–H": "#4b535c",                     # H family (grey)
    "C–C": ECOL["C"],                                        # C
    "C–N": ECOL["N"],  "N–H": "#7d8fe0", "N–N": EDARK["N"],  # N family (blue); N–H periwinkle, kept clear of F's sky
    "C–O": ECOL["O"],  "O–H": EDARK["O"], "N–O": "#f2938a",  # O family (red)
    "C–F": ECOL["F"],                                        # F (sky)
    "C–S": ECOL["S"],  "S–H": EDARK["S"],                    # S family (ochre)
    "C–P": ECOL["P"],  "P–H": EDARK["P"],                    # P family (orange)
    "C–Cl": ECOL["Cl"],                                      # Cl (green)
}
def pair_color(t):
    """Per-pair colour for a bond label like 'C–O' (en-dash); falls back to the heteroatom colour."""
    return PAIRCOL.get(t, bond_hetero_color(*t.split("–")))

# Canonical bond set: the ten most abundant bond types in the DB (~97% of all bonded pairs; a >2x count
# gap separates #10 C–Cl from #11 H–H).  Pinned BY NAME, never by a count threshold -- at full-DB scale
# (~5M rows) many more pairs clear any fixed threshold.  LANDSCAPE_EXTENDED adds the next five (15-pair
# SI variant).
LANDSCAPE_BONDS    = ["C–H", "C–C", "C–N", "C–O", "N–H", "O–H", "C–F", "C–S", "C–P", "C–Cl"]
LANDSCAPE_EXTENDED = LANDSCAPE_BONDS + ["H–H", "N–N", "N–O", "P–H", "S–H"]

def panel(ax, label, title=None, x=-0.085, y=1.02, tsize=8.6, tdx=0.115):
    ax.text(x, y, label.upper(), transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", va="bottom", ha="right", color=INK)
    if title:
        ax.text(x + tdx, y, title, transform=ax.transAxes, fontsize=tsize,
                fontweight="bold", va="bottom", ha="left", color=INK)

def ygrid(ax):
    ax.grid(axis="y", which="major", color="#e9e9ee", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

def style_violin(parts, facecolors, alpha=0.85, edge="white", elw=0.5):
    bodies = parts["bodies"]
    if isinstance(facecolors, str):
        facecolors = [facecolors] * len(bodies)
    for b, c in zip(bodies, facecolors):
        b.set_facecolor(c); b.set_alpha(alpha)
        b.set_edgecolor(edge); b.set_linewidth(elw)
    if "cmedians" in parts:
        parts["cmedians"].set_color(INK); parts["cmedians"].set_linewidth(0.9)
    return parts

def save(fig, path):
    fig.savefig(path)
    fig.savefig(path.replace(".pdf", ".png"), dpi=220)
    plt.close(fig)
    print("  wrote", path.split("/")[-1])
