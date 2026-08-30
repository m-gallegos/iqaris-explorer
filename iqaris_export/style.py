"""style.py -- one object holding every appearance / perception knob for a projection.

`Style` is threaded through `build_scene` (color mapping, atom/bond sizing, bond perception,
hydrogen visibility) and through the emitters (background, labels, sphere quality, window).
The SAME object is what the interactive viewer's key bindings mutate, so CLI, library and live
window all drive identical controls.  Every field has a sensible default that reproduces the
current look, so `Style()` == "what you get today".
"""
from __future__ import annotations
from dataclasses import dataclass


# colormaps the interactive viewer cycles through with the "m" key (name -> matplotlib/seaborn)
CMAP_CYCLE = ("mako", "viridis", "plasma", "inferno", "cividis", "RdBu", "coolwarm", "Spectral")
BOND_MODES = ("qtaim", "di", "distance")


@dataclass
class Style:
    """Renderer-agnostic look + bond-perception settings for one projection."""

    # --- color mapping (override the auto choice) -------------------------------------
    cmap: str = None            # matplotlib/seaborn colormap name; None -> auto (mako/RdBu)
    vmin: float = None          # lower color-scale bound; None -> data min
    vmax: float = None          # upper color-scale bound; None -> data max
    vcenter: float = None       # center for a diverging map; None -> auto (0 if signed)

    # --- atoms -------------------------------------------------------------------------
    atom_scale: float = 1.0     # multiplier on sphere size (spacefill %)
    show_h: bool = True         # draw hydrogens (and bonds to them)
    labels: bool = False        # draw element labels at each atom
    show_values: bool = False   # annotate the projected property's numeric value in 3D
                                # (on atoms / bonds / critical points as appropriate)

    # --- bonds -------------------------------------------------------------------------
    bond_scale: float = 1.0     # multiplier on stick thickness
    bonds: str = "di"           # perception: 'qtaim' (any BCP) | 'di' (DI>cutoff) | 'distance'
    bond_cutoff: float = None   # 'di': DI threshold (default 0.30); 'distance': radii tol (1.20)

    # --- QCT molecular graph (only relevant when graph/paths are drawn) ----------------
    cp_scale: float = 1.0       # multiplier on critical-point dot size
    path_scale: float = 1.0     # multiplier on bond-path tube thickness

    # --- display (emitter only) --------------------------------------------------------
    background: str = "white"
    sphere_res: int = None      # facet resolution override; None -> auto by atom count
    point_size: float = 12.0    # points-only mode marker size
    window_size: tuple = (1280, 860)

    def copy(self, **changes):
        """Return a shallow copy with `changes` applied (used by live key bindings)."""
        from dataclasses import replace
        return replace(self, **changes)
