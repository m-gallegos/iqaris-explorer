"""viz3d.py -- interactive PyVista (VTK) emitter for the backend-neutral :class:`Scene`.

The *interactive* sibling of `project.emit_jmol`: turns the SAME renderer-neutral Scene (atoms /
cylinders / points / disks / colorbar, DB Angstrom frame) into a real depth-sorted VTK window you
rotate with the mouse, and/or a static screenshot.  Self-contained native window -- no browser, no
notebook, no Qt/Tk binding (VTK ships its own interactor).  Jmol stays the publication-grade static
path; this is the exploratory viewer.

Two entry points:
  * `emit_pyvista(scene, ...)` -- render one prebuilt Scene (window or `save=` screenshot).
  * `Viewer(con, sid, ...).show()` -- a *live* window whose key bindings re-project on the fly
    (swap property, colormap, bond perception, atom/bond size, toggle graph/paths/IAS/H/labels).

Everything is driven by a single :class:`~iqaris_export.style.Style`, the same object the CLI fills
and the key bindings mutate, so CLI / library / window all expose identical controls.
"""
from __future__ import annotations
import collections
import os

import numpy as np
import matplotlib.colors as mcolors

# van der Waals radii (Angstrom); sphere radius = vdW * spacefill% so ball-and-stick proportions
# match the Jmol `spacefill N%` convention used for the paper figures.
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "Cl": 1.75}
_VDW_FALLBACK = 1.60


def _mid_ellipsis(s, maxlen):
    """Shorten `s` to `maxlen` chars, eliding the MIDDLE (`ab…yz`) rather than the tail.  IQARIS
    structure ids are `<family>_qm9_<n>`: a shared family prefix (often >12 chars on its own)
    followed by the one part that actually distinguishes entries.  A tail-truncated `s[:n]` chops
    off exactly that distinguishing suffix, so every row in a family shows the identical prefix --
    a middle ellipsis keeps a hint of both ends."""
    if s is None or len(s) <= maxlen:
        return s
    head = (maxlen - 1) // 2
    tail = maxlen - 1 - head
    return f"{s[:head]}…{s[-tail:]}" if tail else f"{s[:head]}…"


def _rgb(hexcolor):
    return (np.array(mcolors.to_rgb(hexcolor)) * 255).astype(np.uint8)


# Value/element labels are 2D on-top overlay labels (add_point_labels) on a rounded pill, so they
# stay legible over atoms/bonds/background, survive zoom, and render headless -- unlike a billboard
# text actor anchored at atom-center depth (occluded by its own sphere, absent from screenshots).
_LABEL_PILL_LIGHTBG = "#243447"   # dark slate pill on a light background (text white)
_LABEL_PILL_DARKBG = "#e8ecf0"    # light pill on a dark background (text near-black)


def _add_point_labels(pl, positions, texts, name, dark_bg=False, font_size=15):
    """Anchor legible labels at each point as always-on-top overlay labels, each on a rounded pill
    for contrast.  Empty strings (e.g. NaN values) are skipped with their positions.  Returns the
    list of actor names added (PyVista registers the label actor as ``f"{name}-labels"``), so a
    Viewer rebuild can remove them by name."""
    pairs = [(p, str(t)) for p, t in zip(positions, texts) if t is not None and str(t) != ""]
    if not pairs:
        return []
    pos = np.array([p for p, _ in pairs], float)
    txt = [t for _, t in pairs]
    pill = _LABEL_PILL_DARKBG if dark_bg else _LABEL_PILL_LIGHTBG
    tcol = "#14202b" if dark_bg else "white"
    pl.add_point_labels(pos, txt, name=name, font_size=int(font_size), bold=True, text_color=tcol,
                        font_family="arial", shape="rounded_rect", shape_color=pill,
                        shape_opacity=0.85, fill_shape=True, margin=3, always_visible=True,
                        show_points=False)
    return [f"{name}-labels"]


def _sphere_res(natoms):
    """Fewer facets per sphere for big systems keeps rotation snappy without looking faceted."""
    if natoms <= 60:
        return 24
    if natoms <= 150:
        return 16
    return 12


def _merge(meshes):
    import pyvista as pv
    meshes = [m for m in meshes if m is not None and m.n_points]
    if not meshes:
        return None
    if len(meshes) == 1:
        return meshes[0]
    return pv.MultiBlock(meshes).combine()


def _cylinder_groups(cylinders):
    """Group cylinders by opacity and merge each group into ONE mesh carrying per-cell rgb, so
    hundreds of bonds/paths cost a handful of actors instead of one actor each."""
    import pyvista as pv
    groups = collections.defaultdict(list)
    for cyl in cylinders:
        p1, p2, col, diam = cyl[0], cyl[1], cyl[2], cyl[3]
        alpha = cyl[4] if len(cyl) > 4 else None
        groups[1.0 if alpha is None else float(alpha)].append((p1, p2, col, diam))
    out = {}
    for alpha, items in groups.items():
        meshes = []
        for p1, p2, col, diam in items:
            p1 = np.asarray(p1, float)
            p2 = np.asarray(p2, float)
            L = float(np.linalg.norm(p2 - p1))
            if L < 1e-6:
                continue
            cy = pv.Cylinder(center=(p1 + p2) / 2, direction=(p2 - p1),
                             radius=diam / 2, height=L, resolution=14, capping=True)
            cy.cell_data["rgb"] = np.tile(_rgb(col), (cy.n_cells, 1))
            meshes.append(cy)
        merged = _merge(meshes)
        if merged is not None:
            out[alpha] = merged
    return out


def _reframe_fn(pos, reframe):
    if reframe and len(pos) >= 2:
        c = pos.mean(0)
        _, _, vt = np.linalg.svd(pos - c, full_matrices=True)
        if np.linalg.det(vt) < 0:
            vt[-1] *= -1
        return lambda p: (np.asarray(p, float) - c) @ vt.T
    return lambda p: np.asarray(p, float)


def add_scene(pl, scene, style, points_only=False, reframe=False, scalar_bar="bottom"):
    """Add all of `scene`'s actors to `pl` under STABLE names, so a later rebuild can replace them
    (via `remove_actor(name)`) without disturbing interactor widgets.  `scalar_bar` places the color
    key: "bottom" (static figures) or "right" (leaves the window's bottom edge free for widgets).
    Returns the list of actor names added."""
    import pyvista as pv
    names = []

    pos0 = np.array([a[0] for a in scene.atoms], float)
    xf = _reframe_fn(pos0, reframe)
    pos = np.array([xf(p) for p in pos0]) if len(pos0) else pos0
    els = [a[3] for a in scene.atoms]
    tcol = "black" if style.background not in ("black", "#000000") else "white"

    # ---- atoms: one glyph actor, exact baked hex ------------------------------------
    atom_opacity = scene.meta.get("atom_translucent") or 1.0
    if len(pos):
        radii = np.array([VDW.get(e, _VDW_FALLBACK) * (a[2] / 100.0)
                          for a, e in zip(scene.atoms, els)])
        cloud = pv.PolyData(pos)
        cloud["radius"] = radii
        cloud["rgb"] = np.array([_rgb(a[1]) for a in scene.atoms])
        if points_only:
            # points are screen-space (fixed pixels), so a big molecule's atoms sit close and read
            # fine while a small molecule's spread out and look like specks -> scale the marker up as
            # the atom count drops (points mode is a fast preview; keep big systems at base size).
            eff_ps = style.point_size * float(np.clip((80.0 / max(len(pos), 1)) ** 0.5, 1.0, 3.0))
            pl.add_mesh(cloud, scalars="rgb", rgb=True, point_size=eff_ps,
                        render_points_as_spheres=True, opacity=atom_opacity, name="atoms")
        else:
            res = style.sphere_res or _sphere_res(len(pos))
            glyph = cloud.glyph(orient=False, scale="radius", factor=1.0,
                                geom=pv.Sphere(theta_resolution=res, phi_resolution=res))
            pl.add_mesh(glyph, scalars="rgb", rgb=True, smooth_shading=True,
                        opacity=atom_opacity, specular=0.3, specular_power=12, name="atoms")
        names.append("atoms")
        if style.labels:
            # shrink the font as atoms crowd, so labels overlap less on big/dense systems
            n = len(pos)
            efs = 16 if n <= 20 else 13 if n <= 45 else 11
            names += _add_point_labels(pl, pos, els, "ellabel", dark_bg=(tcol == "white"),
                                       font_size=efs)

    # ---- bonds + bond paths: merged per opacity -------------------------------------
    if not points_only:
        cyls = [(xf(c[0]), xf(c[1]), c[2], c[3]) + ((c[4],) if len(c) > 4 else ())
                for c in scene.cylinders]
        for k, (alpha, mesh) in enumerate(_cylinder_groups(cyls).items()):
            nm = f"cyl{k}"
            pl.add_mesh(mesh, scalars="rgb", rgb=True, smooth_shading=True, opacity=alpha, name=nm)
            names.append(nm)

    # ---- critical-point dots --------------------------------------------------------
    if scene.points and not points_only:
        cp = pv.PolyData(np.array([xf(p[0]) for p in scene.points], float))
        cp["radius"] = np.array([p[2] / 2 for p in scene.points])
        cp["rgb"] = np.array([_rgb(p[1]) for p in scene.points])
        cpg = cp.glyph(orient=False, scale="radius", factor=1.0,
                       geom=pv.Sphere(theta_resolution=18, phi_resolution=18))
        pl.add_mesh(cpg, scalars="rgb", rgb=True, smooth_shading=True, name="cps")
        names.append("cps")

    # ---- interatomic-surface disks --------------------------------------------------
    if scene.disks and not points_only:
        discs = []
        for center, normal, rad, col, _alpha in scene.disks:
            d = pv.Disc(center=xf(center), inner=0.0, outer=float(rad),
                        normal=np.asarray(normal, float), c_res=40, r_res=1)
            d.cell_data["rgb"] = np.tile(_rgb(col), (d.n_cells, 1))
            discs.append(d)
        merged = _merge(discs)
        if merged is not None:
            pl.add_mesh(merged, scalars="rgb", rgb=True, opacity=float(scene.disks[0][4]), name="ias")
            names.append("ias")

    # ---- projected-value labels (always-on-top pills, anchored on their feature) ------
    if scene.value_labels:
        vpos = [xf(p) for (p, _t) in scene.value_labels]
        vtxt = [t for (_p, t) in scene.value_labels]
        # value labels are wider than element symbols, so they crowd sooner -> shrink the font as the
        # count of drawn (non-empty) labels grows, keeping dense projections legible
        nlab = sum(1 for t in vtxt if t is not None and str(t) != "")
        vfs = 15 if nlab <= 18 else 12 if nlab <= 36 else 10
        names += _add_point_labels(pl, vpos, vtxt, "vallabel", dark_bg=(tcol == "white"),
                                   font_size=vfs)

    # ---- scalar bar (rebuilt from the Scene colorbar) -------------------------------
    if scene.colorbar is not None:
        cmap, nrm, label, unit = scene.colorbar
        vmin, vmax = float(nrm.vmin), float(nrm.vmax)
        bar_title = label
        if scalar_bar == "right":
            # sized/placed to clear the 7-row level radio column below it.
            sb = {"vertical": True, "position_x": 0.925, "position_y": 0.30,
                  "width": 0.045, "height": 0.34}
            # the narrow bar clips long centered titles: keep just the unit here and caption the full
            # "prop (unit) [level]" under the molecule (below), where the canvas width is free.
            bar_title = unit
        elif scalar_bar == "top":
            # horizontal bar near the top, title on top -- gives long property labels the full canvas
            # width (a narrow right-edge vertical bar clips them).
            sb = {"vertical": False, "position_x": 0.27, "position_y": 0.885,
                  "width": 0.46, "height": 0.045}
        else:
            sb = {"position_x": 0.30, "position_y": 0.03, "width": 0.40, "height": 0.06}
        sb.update({"title": bar_title, "n_labels": 5, "color": tcol, "font_family": "arial",
                   "bold": True, "title_font_size": 17, "label_font_size": 13})
        proxy = pv.Line((0, 0, 0), (1e-6, 0, 0))
        proxy["s"] = np.array([vmin, vmax], float)
        pl.add_mesh(proxy, scalars="s", cmap=cmap, clim=(vmin, vmax), opacity=0.0,
                    show_scalar_bar=True, reset_camera=False, name="cbar", scalar_bar_args=sb)
        names.append("cbar")
        if scalar_bar == "right":
            w, h = style.window_size
            cap = pl.add_text(label, position=(w // 2, int(h * 0.235)), font_size=13,
                              color=tcol, name="_cbar_caption")
            try:
                cap.GetTextProperty().SetJustificationToCentered()
                cap.GetTextProperty().SetBold(True)
                cap.GetTextProperty().SetFontFamilyToArial()
            except Exception:
                pass
            names.append("_cbar_caption")
    if scene.title:
        w, h = style.window_size
        _t = pl.add_text(scene.title, position=(14, h - 30), font_size=13, color=tcol, name="_title")
        try:
            _t.GetTextProperty().SetBold(True)
            _t.GetTextProperty().SetFontFamilyToArial()
        except Exception:
            pass
    return names


WINDOW_TITLE = "IQARIS-Explorer"


def build_plotter(scene, off_screen, style, points_only=False, reframe=False):
    """New PyVista ``Plotter`` painted with `scene` (camera reset + gentle tilt)."""
    import pyvista as pv
    pl = pv.Plotter(off_screen=off_screen, window_size=list(style.window_size),
                    title=WINDOW_TITLE)
    pl.background_color = style.background
    pl.enable_anti_aliasing("ssaa")
    add_scene(pl, scene, style, points_only=points_only, reframe=reframe)
    pl.reset_camera()
    pl.camera.azimuth = 20
    pl.camera.elevation = 15
    return pl


def emit_pyvista(scene, interactive=True, save=None, style=None, points_only=False, reframe=False):
    """Render one prebuilt :class:`Scene` with PyVista/VTK.

    interactive : open a native rotatable window (blocks until closed; also writes `save` on close).
    save        : image path; required when ``interactive=False`` (headless off-screen render).
    style       : :class:`~iqaris_export.style.Style` (background, sphere res, window size, labels).
    """
    from .style import Style
    style = style or Style()
    if not interactive and not save:
        raise ValueError("emit_pyvista: pass interactive=True to open a window, or save=PATH "
                         "for a headless screenshot")
    off_screen = not interactive
    pl = build_plotter(scene, off_screen=off_screen, style=style,
                       points_only=points_only, reframe=reframe or off_screen)
    if interactive:
        pl.show(screenshot=save) if save else pl.show()
    else:
        pl.screenshot(save)
        pl.close()
    return save


# ================================================================ LIVE VIEWER ================
# The property dropdown is split by level so it only offers columns that EXIST there: QTAIM (DFT,
# level-free) vs IQA-SQM (level-aware); an impossible pairing like `q_QTAIM` under `PM7` can never be
# selected. Only per-atom/per-bond columns are offered (what a projection can color); molecular
# scalars, labels and flags are excluded. Featured (common) properties float to the top; the dropdown
# still offers every colorable column.
FEATURED_QTAIM = ("q_QTAIM", "N_QTAIM", "K_QTAIM", "G_QTAIM", "L_QTAIM", "LI_QTAIM",
                  "DI_ATOM_QTAIM", "VOL_QTAIM", "DI_QTAIM")
FEATURED_SQM = ("q_SQM", "Eintra_SQM", "N_SQM", "Einter_SQM", "Eexchange_SQM", "Edisp_SQM")
# PBE0 ab initio IQA: the per-atom and per-pair IQA channels worth floating to the top of the dropdown.
FEATURED_IQADFT = ("E_IQA", "EINTRA_IQA", "EINTER_IQA", "V_IQA", "VC_IQA", "VX_IQA",
                   "EINT_IQA", "VXC_IJ_IQA", "VCL_IJ_IQA")
# The PBE0 side carries the whole PBE0 record (ab initio IQA + QTAIM/QCT), so float both headline
# sets (IQA first, then QTAIM).
FEATURED_BY_SIDE = {"qtaim": FEATURED_QTAIM, "sqm": FEATURED_SQM,
                    "iqadft": FEATURED_IQADFT + FEATURED_QTAIM}
# Colormaps offered as clickable buttons (a curated subset of style.CMAP_CYCLE); the "m" key cycles
# this same list so the on-screen selection and the keyboard shortcut never drift apart.
VIEW_CMAPS = ("mako", "viridis", "plasma", "RdBu", "coolwarm")
_GEOM = "(geometry)"
# The level picker shares two DFT records: "DFT" = M06-2X QTAIM density topology (shown "M06-2X"),
# "PBE0" = PBE0 ab initio IQA (the *_IQA props). The five SQM methods keep their names.
LEVEL_DISPLAY = {"DFT": "M06-2X"}


def side_for_level(level):
    """Property side a level of theory admits: 'qtaim' (M06-2X QTAIM density topology),
    'iqadft' (PBE0 ab initio IQA energy decomposition) or 'sqm' (a semi-empirical IQA method)."""
    if level == "DFT":
        return "qtaim"
    if level == "PBE0":
        return "iqadft"
    return "sqm"
# radio-button group names (each group is mutually exclusive; one is a distinct picker)
_G_LEVEL, _G_CMAP, _G_BONDS = "iq_level", "iq_cmap", "iq_bonds"
# property dropdown: rows per column x columns per page (the open menu overlays the scene)
_PROP_ROWS, _PROP_COLS = 16, 2
_PROP_PAGE = _PROP_ROWS * _PROP_COLS


def _level_label(opt):
    """Display name for a level radio option (M06-2X for the DFT sentinel, else the method)."""
    return LEVEL_DISPLAY.get(opt, opt)


def _bind1(fn, arg):
    """Zero-arg callable invoking ``fn(arg)`` -- binds the loop value so every radio button in a
    column carries its own option instead of the last one (classic late-binding trap)."""
    return lambda: fn(arg)


class Viewer:
    """A live IQARIS-Explorer (PyVista/VTK) explorer window.  The **level** of theory (the DFT
    functional M06-2X, or one of the five SQM methods), **colormap** and **bond perception** are
    picked from mutually-exclusive **radio-button** groups; the **property** is chosen from a
    collapsible **dropdown** that pages over EVERY colorable (per-atom / per-pair) column -- not
    just a handful.  Level and property are *coupled*: choosing a level re-scopes the dropdown to
    the properties that exist there (QTAIM under DFT, IQA-SQM under a PM6/PM7 method), so a
    nonsensical pairing such as `q_QTAIM` @ `PM7` is never offered.  Checkbox buttons drive the
    display toggles (graph / paths / IAS / fragment / points / hydrogens / labels / **values** =
    annotate the projected number on its atom, bond or critical point); **color min / color max**
    sliders (with an *auto* reset) set the color-scale range, seeded from the data and hidden for
    a plain-geometry view; sliders size atoms & bonds; a slider + Prev/Next browses a filtered set
    of structures.  The default view is plain CPK geometry until a property is picked.  Every key
    binding still works.  All controls mutate one `Style` + the selection index and re-`build_scene`;
    scene actors are replaced by name so the widgets survive each rebuild.  No notebook, no web,
    no Qt."""

    HELP = ("keys: g/p/i graph paths ias  k fragment  o points  h H  l labels  v values   "
            "m cmap  n bonds  +/- atoms  ,/. sticks   [ ] prev/next   r reset  q quit   "
            "(property: click the dropdown; structure: click the blue button by the title; "
            "level/colormap/bonds: radio buttons; color min/max sliders set the range)")

    def __init__(self, con, sid=None, sids=None, prop=None, method="PM7", graph=False, paths=False,
                 ias=False, fragment=False, style=None, points_only=False, widgets=True, save=None,
                 pbe0=False):
        from .style import Style
        self.con = con
        self.sids = list(sids) if sids else ([sid] if sid else [])
        if not self.sids:
            raise ValueError("Viewer needs a sid or a non-empty sids list")
        self.idx = self.sids.index(sid) if (sid and sid in self.sids) else 0
        # `level` is the single knob that gates which properties are valid: "DFT" (QTAIM, level-free)
        # or one of the SQM methods (level-aware).  `_sqm` remembers the last SQM method chosen so
        # toggling DFT<->SQM restores it (and so fragment engagement, which is always SQM, has one).
        from . import config
        # The interactive window always opens on the plain atoms+bonds view; overlays and coloring
        # are turned on from the controls (a launch-time prop/graph is honored only on the headless
        # emit_pyvista path). `method` is remembered for when an SQM property is later picked.
        self.prop = None
        self._sqm = method if method in config.SQM_METHODS else config.SQM_METHODS[0]
        self.level = "PBE0" if pbe0 else "DFT"   # starting level
        self.graph = self.paths = self.ias = self.fragment = False
        self._frag_prop = "Einter_SQM"   # host<->guest interface channel aggregated in fragment mode
        self._frag_active = False        # True when host<->guest coloring is engaged (nfrag>=2); then
                                         # the property menu picks the channel instead of a color prop
        self.points_only = points_only
        self.style = style or Style()
        self.widgets = widgets
        self._save = save
        self._pl = None
        self._names = []          # scene-actor names added last rebuild
        self._first = True
        self._refit_cam = False   # set when the structure changes -> reframe camera to the new one
        self._building = False    # True while building widget groups -> mute their callbacks
        self._radio_specs = {}    # group -> layout/action, so a key press can resync its highlight
        # collapsible paged MENU, shared by two modes: "prop" (every colorable column) and
        # "entry" (pick another structure).  Only one is ever open, so both reuse ONE fixed
        # button pool + ONE white popup panel; `_menu_mode` says which data it is showing.
        self._menu_mode = None    # None (closed) | "prop" | "entry"
        self._menu_page = 0
        self._mol_open = False    # on-canvas molecular-properties panel visible?
        self._mol_names = []      # its overlay actor names (rebuilt each refresh)
        self._prop_btn = None      # the property toggle button widget
        self._prop_anchor = None   # (x, y) pixel anchor of the property toggle
        self._entry_btn = None     # the "pick another entry" toggle button (by the title)
        self._entry_anchor = None
        self._entry_all = None     # lazily-loaded full (sid) list when not browsing a subset
        self._entry_fmt = {}       # {sid: formula} cache for the entry rows currently shown
        self._panel = None         # white popup background behind the open menu
        self._panel_bd = None      # its thin grey border frame
        # FIXED reusable widget pool for the open menu -- created ONCE, only shown/hidden/relabeled,
        # NEVER torn down.  Destroying a VTK widget mid-session (dropping its last ref while its
        # representation is still in the renderer) is a use-after-free -> heap corruption.
        self._pool_btns = []      # _PROP_PAGE row buttons at fixed grid cells
        self._nav_btns = []       # [prev, next] page buttons
        self._slot_opt = []       # per-slot current option string (None = slot unused)
        self._grid = None         # cached menu grid geometry
        self._radio_widgets = {}  # group -> [(widget, option)]; built ONCE, only re-stated never rebuilt
        # color-range sliders (seeded from the scene's data range; hidden for geometry)
        self._wmin = None
        self._wmax = None
        self._reset_btn = None    # "auto" (reset color range) button widget
        self._wstruct = None      # structure-browser slider (kept in step with Prev/Next)
        self._range_lbls = []     # range-control text-actor names
        self._syncing = False     # mute range-slider callbacks during programmatic value sync
        self._last_range = None   # (min, max) of the scalar shown last render, or None
        self._relayouting = False  # guard against re-entrant resize handling

    @property
    def sid(self):
        return self.sids[self.idx]

    @property
    def method(self):
        """The SQM method to project with: the chosen level when it is an SQM method, else the
        remembered `_sqm` (used for the always-SQM fragment engagement when the level is DFT)."""
        from . import config
        return self.level if self.level in config.SQM_METHODS else self._sqm

    def _side_of(self, name):
        """'qtaim' | 'sqm' | 'iqadft' for a property column, or None for geometry / unknown."""
        if not name:
            return None
        from .registry import get_registry
        specs = get_registry(self.con)["by_name"].get(name)
        return specs[0].side if specs else None

    def _projectable(self, name):
        """True if the column can color a projection -- i.e. it has an atomic (per-atom) or pairs
        (per-bond) reading.  Whole-molecule / CP / angle columns can't be mapped onto atoms."""
        if not name:
            return False
        from .registry import get_registry
        specs = get_registry(self.con)["by_name"].get(name, [])
        return any(s.family in ("atomic", "pairs") for s in specs)

    # -- state mutators --------------------------------------------------------------
    def _set(self, attr, val):
        setattr(self, attr, bool(val))
        self._render_scene()

    def _set_style(self, **kw):
        self.style = self.style.copy(**kw)
        self._render_scene()

    def _set_prop(self, label):
        """Select a property: reset the color range to auto (the new property has its own scale)
        and re-render."""
        self.prop = None if label == _GEOM else label
        self.style = self.style.copy(vmin=None, vmax=None)
        self._render_scene()

    def _set_scale(self, attr, v):
        self.style = self.style.copy(**{attr: round(float(v), 3)})
        self._render_scene()

    def _step_scale(self, attr, factor):
        self._set_scale(attr, max(0.1, getattr(self.style, attr) * factor))

    def _goto(self, i):
        i = max(0, min(len(self.sids) - 1, int(round(i))))
        if i != self.idx:
            self.idx = i
            self._refit_cam = True                   # new molecule -> reframe the camera fresh
            if self._wstruct is not None:            # keep the browser slider in step (set BEFORE
                try:                                 # the render so its handle moves in the same pass)
                    self._wstruct.GetSliderRepresentation().SetValue(self.idx)
                except Exception:
                    pass
            self._render_scene()

    def _step(self, d):
        self._goto(self.idx + d)

    # -- level <-> property coupling -------------------------------------------------
    @property
    def _side(self):
        """The property side the current level admits: 'qtaim' (M06-2X), 'iqadft' (PBE0 ab initio
        IQA) or 'sqm' (an SQM method)."""
        return side_for_level(self.level)

    def _prop_options(self):
        """Geometry + EVERY colorable (atomic or pairs) column that exists at the current level
        (QTAIM under DFT/M06-2X, ab initio IQA under PBE0, IQA-SQM under an SQM method).  Common
        properties are floated to the front
        so they land on the dropdown's first page; the rest follow alphabetically.  Non-numeric
        keys (element labels, boolean flags) are filtered out -- they can't drive a color map."""
        from .registry import get_registry
        reg = get_registry(self.con)
        side = self._side
        names = []
        for fam in ("atomic", "pairs"):
            for s in reg["by_group"].get((fam, side), []):
                if s.unit == "label" or s.dtype.upper().startswith(("VARCHAR", "BOOL", "BLOB")):
                    continue
                names.append(s.name)
        featured = [c for c in FEATURED_BY_SIDE.get(side, ()) if c in names]
        rest = sorted(n for n in names if n not in featured)
        from .project import FRAGMENTS_PROP
        # `fragments` is a structural pseudo-property (color atoms by covalent component); it is
        # level-independent, so it sits right after plain geometry at every level -- parity with the
        # Qt viewer and the `--prop fragments` CLI (host<->guest still lives on the `fragment` toggle).
        return [_GEOM, FRAGMENTS_PROP] + featured + rest

    def _frag_prop_options(self):
        """Interface channels fragment mode can aggregate, restricted to the current level of theory
        (mirrors the Qt viewer): the IQA-SQM decomposition at an SQM method, the ab initio IQA
        interaction split under PBE0, the QTAIM delocalization index at DFT/M06-2X -- so DI never
        shows under PM6, Einter never under M06-2X, and EINT_IQA only under PBE0."""
        from .project import FRAGMENT_PAIR_PROPS
        from .registry import get_registry
        have = {s.name for s in get_registry(self.con)["by_group"].get(("pairs", self._side), [])}
        order = ["Einter_SQM", "Eelstat_SQM", "Eexchange_SQM", "Eresonance_SQM", "Edisp_SQM",
                 "Eelstat_ee_SQM", "Eelstat_en_SQM", "Eelstat_nn_SQM",
                 "EINT_IQA", "VCL_IJ_IQA", "VXC_IJ_IQA", "VNE_IJ_IQA", "VEN_IJ_IQA", "VEE_IJ_IQA",
                 "DI_QTAIM"]
        return [c for c in order if c in FRAGMENT_PAIR_PROPS and c in have]

    def _prop_label(self):
        if self._frag_active:                # fragment mode: the picker names the aggregated channel
            return self._frag_prop
        return self.prop if self.prop else _GEOM

    def _set_level(self, level):
        """Switch level.  If this flips the QTAIM<->SQM side, the current property may no longer
        exist -> remap it to the analogous column (q_QTAIM<->q_SQM, ...) or the level's first
        property, and rebuild the property radio group so only valid choices are shown."""
        from . import config
        if level == self.level:
            return
        old_side = self._side
        self.level = level
        if level in config.SQM_METHODS:
            self._sqm = level
        if self._side != old_side:
            self.prop = self._remap_prop()
            self.style = self.style.copy(vmin=None, vmax=None)   # new side -> auto range
            self._menu_page = 0
            self._refresh_prop_button()
            if self._menu_mode == "prop":
                self._show_menu()
        self._render_scene()

    def _remap_prop(self):
        """Map the current property onto the new level's side: keep geometry as geometry; try the
        obvious suffix swap (`*_QTAIM`<->`*_SQM`); else fall back to the level's first property."""
        from .project import FRAGMENTS_PROP
        options = self._prop_options()                  # already reflects the new level
        if self.prop is None or self.prop == FRAGMENTS_PROP:
            return self.prop                            # geometry / structural: level-independent
        if self.prop in options:
            return self.prop                            # same column exists here (e.g. a QTAIM prop
                                                        # shared by M06-2X and PBE0) -> keep it
        cand = None
        if self._side == "sqm" and self.prop.endswith("_QTAIM"):
            cand = self.prop[:-len("_QTAIM")] + "_SQM"
        elif self._side == "qtaim" and self.prop.endswith("_SQM"):
            cand = self.prop[:-len("_SQM")] + "_QTAIM"
        if cand and cand in options:
            return cand
        # No QTAIM<->SQM suffix analog (e.g. entering the PBE0 ab initio IQA side): fall to the first
        # REAL property, skipping the structural `fragments` pseudo-column.
        rest = [p for p in options if p not in (_GEOM, FRAGMENTS_PROP)]
        return rest[0] if rest else None

    # -- render (replace scene actors by name; leave widgets alone) ------------------
    def _render_scene(self):
        from .project import build_scene
        pl = self._pl
        # keep the user's current rotation across property/toggle changes on the SAME molecule, but
        # reframe fresh whenever the STRUCTURE changed (a new molecule has its own size/center) or on
        # the very first paint -- otherwise a small molecule inherits a big one's zoom and is lost.
        cam = None if (self._first or self._refit_cam) else pl.camera_position
        self._refit_cam = False
        for n in self._names:
            try:
                pl.remove_actor(n, reset_camera=False, render=False)
            except Exception:
                pass
        for t in list(getattr(pl, "scalar_bars", {}).keys()):
            try:
                pl.remove_scalar_bar(t, render=False)
            except Exception:
                pass
        # reconcile the aggregated fragment channel to the CURRENT level BEFORE building, so the render
        # never shows a channel invalid for the level (DI only at DFT, IQA-SQM only at an SQM method,
        # ab initio IQA only at PBE0) -- mirrors the Qt viewer.
        if self.fragment:
            fopts = self._frag_prop_options()
            if fopts and self._frag_prop not in fopts:
                self._frag_prop = fopts[0]
        was_frag = self._frag_active
        scene = build_scene(self.con, self.sid, prop=self.prop, method=self.method,
                            graph=self.graph, paths=self.paths, ias=self.ias,
                            fragment=self.fragment, fragment_prop=self._frag_prop,
                            style=self.style, side=self._side)
        self._frag_active = bool(self.fragment and scene.meta.get("nfrag", 1) >= 2)
        self._last_range = scene.meta.get("data_range")     # drives the color-range sliders
        self._names = add_scene(pl, scene, self.style, points_only=self.points_only,
                                scalar_bar="right")
        pos = f"   [{self.idx + 1}/{len(self.sids)}]" if len(self.sids) > 1 else ""
        w, h = self.style.window_size
        tcol = self._text_color()
        # title is shifted right of the entry-picker toggle button (click it to change structure)
        title = pl.add_text((scene.title or self.sid) + pos, position=(42, h - 30),
                            font_size=13, color=tcol, name="_title")
        try:
            title.GetTextProperty().SetBold(True)
            title.GetTextProperty().SetFontFamilyToArial()
        except Exception:
            pass
        if cam is not None:
            pl.camera_position = cam
        else:
            pl.reset_camera()
            pl.camera.azimuth = 20
            pl.camera.elevation = 15
        self._first = False
        self._sync_range()               # show/hide + reseed the color-range sliders
        self._refresh_prop_button()      # keep the dropdown label in step with the prop / channel
        if self._frag_active != was_frag and self._menu_mode == "prop":
            self._show_menu()            # entering/leaving fragment mode flips the menu list (props
            #                              <-> host-guest channels) -> repaint it if it is open
        if self._mol_open:               # keep the molecular panel current for the new structure/level
            self._refresh_molecular()
        self._sync_overlay_labels()      # re-hide value/element pills if a popup is open (they were
                                         # just re-added above, always visible by default)
        pl.render()

    # -- molecular-properties panel (whole-molecule scalars can't be projected, so list them) --
    def _toggle_molecular(self, on):
        """Show / hide the on-canvas panel of every whole-molecule scalar for the current structure
        and level -- the totals (energies, CP census, dipole, HOMO/LUMO, IQA molecular energy, ...)
        that cannot be colored onto atoms."""
        self._mol_open = bool(on)
        self._refresh_molecular()
        self._sync_overlay_labels()
        try:
            self._pl.render()
        except Exception:
            pass

    def _refresh_molecular(self):
        pl = self._pl
        for n in self._mol_names:
            try:
                pl.remove_actor(n, reset_camera=False, render=False)
            except Exception:
                pass
        self._mol_names = []
        if not self._mol_open or pl is None:
            return
        from .project import molecular_values
        side = self._side
        rows = molecular_values(self.con, self.sid, side=side, method=self.method)
        lvl = (_level_label("DFT") if side == "qtaim"
               else _level_label("PBE0") if side == "iqadft" else self.method)
        w, h = self.style.window_size
        if not rows:
            self._label(f"no molecular data at [{lvl}] for {self.sid}", int(w * 0.30), h - 96,
                        size=12, color="#8a2b2b", name="_mol_hdr")
            self._mol_names.append("_mol_hdr")
            return

        def fmt(v):
            if v is None:
                return "-"
            if isinstance(v, bool):
                return "true" if v else "false"
            try:
                f = float(v)
            except (TypeError, ValueError):
                return str(v)
            return str(int(f)) if (f == int(f) and abs(f) < 1e12) else f"{f:.5g}"

        ncol = 2
        per = (len(rows) + ncol - 1) // ncol
        line_h, pw = 17, 600
        ph = 30 + per * line_h + 14
        x0 = max(170, (w - pw) // 2 - 40)    # centred but nudged left of the top-right pickers
        y1 = h - 74
        y0 = max(16, y1 - ph)
        self._rect2d("_mol_bd_b", (x0 - 2, y0 - 2, x0 + pw + 2, y1 + 2), "#9aa0a6", 1.0, visible=True)
        self._rect2d("_mol_bd", (x0, y0, x0 + pw, y1), "#ffffff", 0.95, visible=True)
        self._mol_names += ["_mol_bd_b", "_mol_bd"]
        self._label(f"molecular  [{lvl}]   {self.sid}   ({len(rows)} totals)",
                    x0 + 12, y1 - 22, size=11, color="#14202b", name="_mol_hdr")
        self._mol_names.append("_mol_hdr")
        colw = (pw - 20) // ncol
        for c in range(ncol):
            chunk = rows[c * per:(c + 1) * per]
            if not chunk:
                continue
            lines = [f"{name} = {fmt(val)}"
                     f"{'' if unit in (None, 'label') else ' ' + str(unit)}"
                     for name, val, unit in chunk]
            nm = f"_mol_col{c}"
            a = pl.add_text("\n".join(lines), position=(x0 + 14 + c * colw, y1 - 30 - per * line_h),
                            font_size=7, color="#22303c", name=nm)
            try:
                a.GetTextProperty().SetFontFamilyToArial()
            except Exception:
                pass
            self._mol_names.append(nm)

    # -- widgets ---------------------------------------------------------------------
    def _add_widgets(self):
        from .style import BOND_MODES
        from . import config
        pl = self._pl
        w, h = self.style.window_size
        tcol = self._text_color()
        self._building = True

        # a small toggle by the title: click it to pick another structure (entry) interactively.
        self._add_entry_toggle(x=16, y=h - 32)

        # left column: display | colormap | bonds. DISPLAY toggles ("values" annotates the numbers).
        self._label("display", 14, h - 74, size=12, color=tcol, name="_hdr_display")
        toggles = [
            ("graph", self.graph, lambda s: self._set("graph", s)),
            ("paths", self.paths, lambda s: self._set("paths", s)),
            ("IAS", self.ias, lambda s: self._set("ias", s)),
            ("fragment", self.fragment, lambda s: self._set("fragment", s)),
            ("points", self.points_only, lambda s: self._set("points_only", s)),
            ("H", self.style.show_h, lambda s: self._set_style(show_h=s)),
            ("labels", self.style.labels, lambda s: self._set_style(labels=s)),
            ("values", self.style.show_values, lambda s: self._set_style(show_values=s)),
            ("molecular", self._mol_open, lambda s: self._toggle_molecular(s)),
        ]
        # On a short window the left column would overlap the color-range slider band (at 0.20 of
        # window height), so compress the inter-row gaps by a scale factor (hit-target size stays
        # fixed) to keep it above the band; floored so text stays legible.
        _LEFT_COL_H = 9 * 27 + 7 * 24 + 104        # natural extent at scale 1.0 (see gaps below)
        # the color-range slider band is pinned at 0.20 of window height, plus ~25px for its value
        # label drawn above the track, plus a 15px safety margin
        headroom = 0.80 * h - 98 - 40
        scale = max(0.65, min(1.0, headroom / _LEFT_COL_H))

        bsz, y, gap = 18, h - 98, round(27 * scale)
        for label, val, cb in toggles:
            pl.add_checkbox_button_widget(cb, value=bool(val), position=(14, y), size=bsz,
                                          color_on="#2c7fb8", color_off="#cfd4d9")
            self._label(label, 14 + bsz + 6, y + 2, size=11, name=f"_lbl_{label}")
            y -= gap

        # a horizontal rule separates the toggle block from the colormap/bond pickers below it,
        # then the two picker groups get their own generous gap + header (was cramped together).
        RS, DY = 13, round(24 * scale)
        sep1 = y - round(4 * scale)
        self._rect2d("_sep1", (12, sep1, 150, sep1 + 2), "#c7cbd1", 1.0, visible=True)
        cmap_top = sep1 - round(26 * scale)
        self._radio_specs[_G_CMAP] = dict(
            options=list(VIEW_CMAPS), action=lambda c: self._set_style(cmap=c),
            x=16, y_top=cmap_top - round(22 * scale), dy=DY, size=RS)
        bonds_top = cmap_top - round(22 * scale) - len(VIEW_CMAPS) * DY - round(30 * scale)
        self._radio_specs[_G_BONDS] = dict(
            options=list(BOND_MODES), action=lambda b: self._set_style(bonds=b),
            x=16, y_top=bonds_top - round(22 * scale), dy=DY, size=RS)
        self._add_radio_column(_G_CMAP, "colormap",
                               self.style.cmap if self.style.cmap in VIEW_CMAPS else VIEW_CMAPS[0],
                               header_y=cmap_top, tcol=tcol)
        sep2 = bonds_top + 6
        self._rect2d("_sep2", (12, sep2, 150, sep2 + 2), "#c7cbd1", 1.0, visible=True)
        self._add_radio_column(_G_BONDS, "bonds", self.style.bonds, header_y=bonds_top, tcol=tcol)

        # top-right: level picker (radios; DFT shown as its functional M06-2X) + a property
        # DROPDOWN.  Choosing a level re-scopes the dropdown to the properties that level offers.
        self._radio_specs[_G_LEVEL] = dict(
            options=["DFT", "PBE0"] + list(config.SQM_METHODS), action=self._set_level,
            x=w - 175, y_top=h - 66, dy=DY, size=RS, label_of=_level_label)
        self._add_radio_column(_G_LEVEL, "level", self.level, header_y=h - 44, tcol=tcol)
        self._add_prop_dropdown(x=w - 360, y=h - 66, header_y=h - 44)

        self._building = False

        # bottom band (two rows): color-range sliders on top (y~0.20), size + browser below (y~0.06)
        sl = dict(style="modern", slider_width=0.018, tube_width=0.004)
        pl.add_slider_widget(lambda v: self._set_scale("atom_scale", v), rng=(0.3, 3.0),
                             value=self.style.atom_scale, title="atom size",
                             pointa=(0.05, 0.06), pointb=(0.24, 0.06), fmt="%.2f", **sl)
        pl.add_slider_widget(lambda v: self._set_scale("bond_scale", v), rng=(0.3, 3.0),
                             value=self.style.bond_scale, title="bond size",
                             pointa=(0.76, 0.06), pointb=(0.95, 0.06), fmt="%.2f", **sl)
        self._add_range_sliders()
        if len(self.sids) > 1:
            self._wstruct = pl.add_slider_widget(
                lambda v: self._goto(v), rng=(0, len(self.sids) - 1), value=self.idx,
                title=f"structure  (of {len(self.sids)})",
                pointa=(0.37, 0.06), pointb=(0.63, 0.06), fmt="%.0f", **sl)
            pl.add_checkbox_button_widget(lambda s: self._step(-1), value=False,
                                          position=(int(0.34 * w), int(0.115 * h)), size=18,
                                          color_on="#dd5577", color_off="#dd5577")
            self._label("< Prev", int(0.34 * w) + 22, int(0.115 * h) + 2, size=10, color="black",
                        name="_prev")
            pl.add_checkbox_button_widget(lambda s: self._step(1), value=False,
                                          position=(int(0.63 * w), int(0.115 * h)), size=18,
                                          color_on="#dd5577", color_off="#dd5577")
            self._label("Next >", int(0.63 * w) + 22, int(0.115 * h) + 2, size=10, color="black",
                        name="_next")
        self._sync_range()      # seed / hide the range sliders for the first scene

    # -- radio-button pickers --------------------------------------------------------
    def _radio_cb(self, fn):
        """Wrap a radio action so it is inert while we are (re)building groups -- VTK fires the
        callback when a button is set on, including our programmatic value=True at build time."""
        def cb():
            if not self._building:
                fn()
        return cb

    def _text_color(self):
        return "black" if self.style.background not in ("black", "#000000") else "white"

    def _label(self, text, x, y, size=11, color=None, name=None):
        """A bold, anti-aliased on-screen text label (crisper than the small default) at a pixel
        position; re-adding with the same `name` replaces it in place."""
        a = self._pl.add_text(str(text), position=(int(x), int(y)), font_size=size,
                              color=color or self._text_color(), name=name)
        try:
            a.GetTextProperty().SetBold(True)
            a.GetTextProperty().SetFontFamilyToArial()
        except Exception:
            pass
        return a

    def _rect2d(self, name, px, color, opacity=1.0, visible=False):
        """A filled 2D rectangle at pixel bounds `px=(x0,y0,x1,y1)` (origin bottom-left), added
        under a STABLE name.  VTK draws overlay 2D actors in insertion order (layer numbers are
        not honored here), so create a background rectangle BEFORE the widgets/labels it must sit
        behind.  Used for the white property/entry popup panel, its border, and the thin
        section-separator rules in the control column.  Returns the `vtkActor2D`."""
        import vtk
        w, h = self.style.window_size
        x0, y0, x1, y1 = px[0] / w, px[1] / h, px[2] / w, px[3] / h
        pts = vtk.vtkPoints()
        for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1)):
            pts.InsertNextPoint(x, y, 0)
        poly = vtk.vtkPolyData()
        poly.SetPoints(pts)
        quad = vtk.vtkCellArray()
        quad.InsertNextCell(4)
        for i in range(4):
            quad.InsertCellPoint(i)
        poly.SetPolys(quad)
        mapper = vtk.vtkPolyDataMapper2D()
        mapper.SetInputData(poly)
        coord = vtk.vtkCoordinate()
        coord.SetCoordinateSystemToNormalizedViewport()
        mapper.SetTransformCoordinate(coord)
        act = vtk.vtkActor2D()
        act.SetMapper(mapper)
        act.GetProperty().SetColor(*mcolors.to_rgb(color))
        act.GetProperty().SetOpacity(float(opacity))
        act.SetVisibility(1 if visible else 0)
        self._pl.add_actor(act, name=name, reset_camera=False, render=False)
        return act

    def _sync_overlay_labels(self):
        """Hide the element/value pill labels while the property/entry popup or the molecular-
        properties panel is showing.  Those pills are `add_point_labels(always_visible=True)` --
        VTK paints them on top of EVERYTHING regardless of insertion order, so left alone they
        float over a popup's text/buttons instead of sitting behind it like every other overlay."""
        hide = bool(self._menu_mode) or self._mol_open
        for nm in ("ellabel-labels", "vallabel-labels"):
            actor = self._pl.actors.get(nm) if self._pl is not None else None
            if actor is not None:
                actor.SetVisibility(0 if hide else 1)

    @staticmethod
    def _set_vis(actor, on):
        if actor is not None:
            try:
                actor.SetVisibility(1 if on else 0)
            except Exception:
                pass

    def _add_radio_buttons(self, group, options, current, action, x, y_top, dy, size, label_of=None):
        """Lay out one vertical stack of mutually-exclusive radio buttons ONCE; `current` is checked.
        The button's built-in title is skipped (huge/unstyleable font) in favor of our own crisp
        label; `label_of` maps an internal option to its display text (DFT sentinel -> "M06-2X").
        The widgets are kept in `self._radio_widgets[group]` so a later keyboard change can move the
        highlight by re-stating them (`_resync_radio`) instead of tearing the group down (destroying
        a live widget corrupts the VTK heap)."""
        pl = self._pl
        pairs = []
        y = y_top
        for i, opt in enumerate(options):
            wdg = pl.add_radio_button_widget(self._radio_cb(_bind1(action, opt)), group,
                                             value=(opt == current), title=None, position=(x, y),
                                             size=size, color_on="#2c7fb8", color_off="#cfd4d9")
            self._label(label_of(opt) if label_of else str(opt), x + size + 6, y + 2,
                        size=11, name=f"_rl_{group}_{i}")
            pairs.append((wdg, opt))
            y -= dy
        self._radio_widgets[group] = pairs

    def _add_radio_column(self, group, header, current, header_y, tcol):
        """Header label + a radio stack for a group whose geometry lives in `self._radio_specs`."""
        sp = self._radio_specs[group]
        self._label(header, sp["x"], header_y, size=12, color=tcol, name=f"_hdr_{group}")
        self._add_radio_buttons(group, sp["options"], current, sp["action"],
                                sp["x"], sp["y_top"], sp["dy"], sp["size"], sp.get("label_of"))

    # -- property dropdown (collapsible, paged over EVERY colorable column) -----------
    def _prop_disp(self, name):
        """Display name for the picker: drop the `_QTAIM`/`_SQM` suffix (the level already says
        which side), so long column names stay short and the two menu columns never collide."""
        if name == _GEOM:
            return name
        for suf in ("_QTAIM", "_SQM"):
            if name.endswith(suf):
                return name[:-len(suf)]
        return name

    def _prop_button_text(self):
        lab = self._prop_label()                            # full name here (the menu is compact)
        lab = lab if len(lab) <= 20 else lab[:19] + ".."
        return lab + "   v"                                 # ASCII only (VTK text stays ASCII-safe)

    # -- structure (entry) picker: a second mode of the SAME shared popup menu ----------
    def _add_entry_toggle(self, x, y):
        """A small toggle button next to the title: clicking it opens the shared popup in "entry"
        mode -- a paged list of structures to jump to, so the entry can be changed IN the window
        (not only from the terminal)."""
        self._entry_anchor = (x, y)
        self._entry_btn = self._pl.add_checkbox_button_widget(
            self._toggle_entry_menu, value=False, position=(x, y), size=18,
            color_on="#2c7fb8", color_off="#2c7fb8")

    def _entry_list(self):
        """The structure ids the entry picker pages over: the browse subset when one was supplied,
        otherwise the whole database (listed once, lazily)."""
        if len(self.sids) > 1:
            return self.sids
        if self._entry_all is None:
            self._entry_all = [r[0] for r in self.con.execute(
                "SELECT structure_id FROM structure ORDER BY structure_id").fetchall()]
        return self._entry_all

    def _fetch_formulas(self, sids):
        """{sid: formula} for a handful of visible ids -- shown beside each id so entries are
        identifiable at a glance (fetched per page, so it stays cheap even over the full DB)."""
        if not sids:
            return {}
        ph = ",".join(["?"] * len(sids))
        rows = self.con.execute(
            f"SELECT structure_id, formula FROM structure WHERE structure_id IN ({ph})",
            list(sids)).fetchall()
        return {r[0]: r[1] for r in rows}

    def _menu_panel_px(self, g, pad=0):
        """Pixel bounds (x0,y0,x1,y1) of the popup panel wrapping the header + both columns + the
        nav row, grown by `pad` (used to draw a slightly larger border rectangle behind it).
        Column 2 gets the SAME width budget as column 1 (`colw`, not a smaller fixed margin) --
        it was narrower before, so column-2 text that fit fine between the two columns could still
        run past the panel's own right edge."""
        return (g["cx0"] - 16 - pad, g["navy"] - 12 - pad,
                g["cx0"] + 2 * g["colw"] + pad, g["hdry"] + 24 + pad)

    def _add_prop_dropdown(self, x, y, header_y):
        """Collapsed picker: a header + a toggle button whose label is the current property.
        Clicking it shows the paged menu of ALL colorable columns for the current level.  The menu
        widgets are a FIXED pool built once here (see `_build_prop_pool`) and only shown/hidden --
        never created/destroyed per open, which would free a live widget and corrupt the heap."""
        self._prop_anchor = (x, y)
        self._label("property", x, header_y, size=12, name="_hdr_prop")
        self._prop_btn = self._pl.add_checkbox_button_widget(
            self._toggle_prop_menu, value=False, position=(x, y), size=20,
            color_on="#7b5cd6", color_off="#7b5cd6")
        self._label(self._prop_button_text(), x + 26, y + 3, size=11, name="_prop_cur")
        self._build_prop_pool()

    def _refresh_prop_button(self):
        """Keep the collapsed dropdown's label in step with `self.prop`."""
        if self._prop_btn is None or self._prop_anchor is None:
            return
        x, y = self._prop_anchor
        self._label(self._prop_button_text(), x + 26, y + 3, size=11, name="_prop_cur")

    # --- fixed reusable menu pool (created ONCE; only shown/hidden/relabeled) ----------
    def _menu_grid(self):
        w, h = self.style.window_size
        return dict(cx0=int(w * 0.20), colw=300, y_top=h - 84, dy=26, size=16,
                    navy=(h - 84) - _PROP_ROWS * 26 - 6, hdry=h - 58)

    def _slot_xy(self, k, g):
        return g["cx0"] + (k // _PROP_ROWS) * g["colw"], g["y_top"] - (k % _PROP_ROWS) * g["dy"]

    def _build_prop_pool(self):
        """Create the _PROP_PAGE row buttons + 2 nav buttons ONCE, at fixed grid cells, hidden.
        They are reused for every page/open; positions never change, so no widget is ever torn
        down mid-session (the operation that was corrupting the VTK heap)."""
        pl, g = self._pl, self._menu_grid()
        self._grid = g
        # white popup panel (+ grey border frame) BEHIND the menu -- created here, BEFORE the pool
        # buttons/labels, so VTK draws it under them (2D overlay actors paint in insertion order).
        # Hidden until a menu opens; gives the property/entry names an opaque backdrop instead of
        # floating over the molecule.
        self._panel_bd = self._rect2d("_menu_panel_bd", self._menu_panel_px(g, pad=4), "#9aa0a6", 1.0)
        self._panel = self._rect2d("_menu_panel", self._menu_panel_px(g, pad=1), "white", 0.97)
        self._slot_opt = [None] * _PROP_PAGE
        was, self._building = self._building, True
        for k in range(_PROP_PAGE):
            x, y = self._slot_xy(k, g)
            bw = pl.add_checkbox_button_widget(self._pool_cb(k), value=False, position=(x, y),
                                               size=g["size"], color_on="#7b5cd6",
                                               color_off="#d9dbe0")
            self._hide_widget(bw)
            self._pool_btns.append(bw)
        for dx, d in ((0, -1), (120, 1)):
            bw = pl.add_checkbox_button_widget(self._nav_cb(d), value=False,
                                               position=(g["cx0"] + dx, g["navy"]), size=g["size"],
                                               color_on="#7b5cd6", color_off="#7b5cd6")
            self._hide_widget(bw)
            self._nav_btns.append(bw)
        self._building = was

    def _pool_cb(self, k):
        """Row-slot callback (bound to slot index k, whose meaning changes per page/mode)."""
        def cb(_state):
            if self._building:
                return
            opt = self._slot_opt[k] if k < len(self._slot_opt) else None
            if opt is not None:
                self._menu_select(opt)
        return cb

    def _nav_cb(self, d):
        def cb(_state):
            if not self._building:
                self._page(d)
        return cb

    @staticmethod
    def _hide_widget(bw):
        try:
            bw.Off()
            bw.GetRepresentation().SetVisibility(0)
        except Exception:
            pass

    @staticmethod
    def _show_widget(bw, checked=None):
        try:
            if checked is not None:
                bw.GetRepresentation().SetState(1 if checked else 0)
            bw.On()
            bw.GetRepresentation().SetVisibility(1)
        except Exception:
            pass

    def _pop_toggle(self, btn):
        """Force a toggle button's representation back to the OFF state (no callback fires)."""
        try:
            btn.GetRepresentation().SetState(0)
        except Exception:
            pass

    def _toggle_prop_menu(self, state):
        if self._building:
            return
        self._open_menu("prop") if state else self._close_menu()

    def _toggle_entry_menu(self, state):
        if self._building:
            return
        self._open_menu("entry") if state else self._close_menu()

    def _open_menu(self, mode):
        """Open the shared popup in `mode` ('prop' | 'entry').  Only one is ever open, so pop the
        other toggle out first.  The entry picker opens on the page holding the current structure,
        so prev/next browses its neighbors."""
        self._pop_toggle(self._entry_btn if mode == "prop" else self._prop_btn)
        self._menu_mode = mode
        self._menu_page = 0
        if mode == "entry":
            try:
                self._menu_page = self._entry_list().index(self.sid) // _PROP_PAGE
            except ValueError:
                self._menu_page = 0
        self._show_menu()

    def _close_menu(self):
        """Hide the popup + pop whichever toggle opened it."""
        self._pop_toggle(self._prop_btn if self._menu_mode == "prop" else self._entry_btn)
        self._menu_mode = None
        self._hide_menu()
        self._pl.render()

    # -- mode dispatch: the shared pool shows properties OR structures --------------------
    def _menu_options(self):
        if self._menu_mode == "entry":
            return self._entry_list()
        return self._frag_prop_options() if self._frag_active else self._prop_options()

    def _menu_current(self):
        return self.sid if self._menu_mode == "entry" else self._prop_label()

    def _menu_disp(self, opt):
        if self._menu_mode == "entry":
            f = self._entry_fmt.get(opt)
            sid = _mid_ellipsis(opt, 9)
            f = _mid_ellipsis(f, 7) if f else f
            return f"{sid} {f}" if f else sid
        return self._prop_disp(opt)

    def _menu_header(self, npages):
        if self._menu_mode == "entry":
            return (f"select structure   (page {self._menu_page + 1}/{npages}, "
                    f"{len(self._entry_list())} total)")
        if self._frag_active:
            return f"select host-guest channel   (page {self._menu_page + 1}/{npages})"
        return f"select property   (page {self._menu_page + 1}/{npages})"

    def _menu_select(self, opt):
        self._choose_entry(opt) if self._menu_mode == "entry" else self._choose_prop(opt)

    def _choose_prop(self, opt):
        self._close_menu()
        if self._frag_active:                               # fragment mode: pick the aggregated channel
            self._frag_prop = opt
            self.style = self.style.copy(vmin=None, vmax=None)
            self._render_scene()
        else:
            self._set_prop(opt)                             # resets range + renders + relabels

    def _choose_entry(self, sid):
        self._close_menu()
        if len(self.sids) > 1 and sid in self.sids:         # browse subset: jump + sync the slider
            self._goto(self.sids.index(sid))
        elif sid != self.sid:                               # whole-DB mode: make it the shown one
            self.sids = [sid]
            self.idx = 0
            self._refit_cam = True                          # new molecule -> reframe the camera fresh
            self._render_scene()

    def _page(self, d):
        self._menu_page += d
        self._show_menu()

    def _show_menu(self):
        """Show the current page by REUSING the fixed pool: relabel + re-check + reveal the needed
        slots, hide the rest; recompute pagination.  Reveals the white backdrop panel too.  No
        widget is created or destroyed.  Works for both the property and the structure list."""
        pl, g = self._pl, self._grid
        opts = self._menu_options()
        npages = max(1, (len(opts) + _PROP_PAGE - 1) // _PROP_PAGE)
        self._menu_page = max(0, min(self._menu_page, npages - 1))
        page = opts[self._menu_page * _PROP_PAGE:(self._menu_page + 1) * _PROP_PAGE]
        self._entry_fmt = self._fetch_formulas(page) if self._menu_mode == "entry" else {}
        cur = self._menu_current()
        self._set_vis(self._panel_bd, True)
        self._set_vis(self._panel, True)
        was, self._building = self._building, True
        self._label(self._menu_header(npages), g["cx0"], g["hdry"], size=12, name="_menu_hdr")
        for k in range(_PROP_PAGE):
            bw = self._pool_btns[k]
            if k < len(page):
                opt = page[k]
                self._slot_opt[k] = opt
                self._show_widget(bw, checked=(opt == cur))
                x, y = self._slot_xy(k, g)
                self._label(self._menu_disp(opt), x + 22, y + 1, size=11, name=f"_menu_lbl_{k}")
            else:
                self._slot_opt[k] = None
                self._hide_widget(bw)
                self._remove_label(f"_menu_lbl_{k}")
        if npages > 1:
            for bw in self._nav_btns:
                self._show_widget(bw)
            self._label("prev", g["cx0"] + 22, g["navy"] + 1, size=11, name="_menu_prev")
            self._label("next", g["cx0"] + 142, g["navy"] + 1, size=11, name="_menu_next")
        else:
            for bw in self._nav_btns:
                self._hide_widget(bw)
            self._remove_label("_menu_prev")
            self._remove_label("_menu_next")
        self._building = was
        self._sync_overlay_labels()
        pl.render()

    def _hide_menu(self):
        """Hide (never destroy) every menu widget + the backdrop panel + remove the text labels."""
        for bw in self._pool_btns + self._nav_btns:
            self._hide_widget(bw)
        self._set_vis(self._panel, False)
        self._set_vis(self._panel_bd, False)
        self._sync_overlay_labels()
        for nm in ([f"_menu_lbl_{k}" for k in range(_PROP_PAGE)]
                   + ["_menu_prev", "_menu_next", "_menu_hdr"]):
            self._remove_label(nm)

    def _remove_label(self, nm):
        try:
            self._pl.remove_actor(nm, reset_camera=False, render=False)
        except Exception:
            pass

    # -- color-range sliders (min/max), seeded from the scene's data range ------------
    def _add_range_sliders(self):
        pl = self._pl
        w, h = self.style.window_size
        sl = dict(style="modern", slider_width=0.018, tube_width=0.004)
        # PyVista fires a slider's callback once at creation with the initial value -- mute that
        # (0/1 are placeholders; the real bounds are seeded by the following _sync_range()).
        self._syncing = True
        try:
            self._wmin = pl.add_slider_widget(self._on_vmin, rng=(0.0, 1.0), value=0.0,
                                              title="color min", pointa=(0.05, 0.20),
                                              pointb=(0.29, 0.20), fmt="%.3g", **sl)
            self._wmax = pl.add_slider_widget(self._on_vmax, rng=(0.0, 1.0), value=1.0,
                                              title="color max", pointa=(0.71, 0.20),
                                              pointb=(0.95, 0.20), fmt="%.3g", **sl)
        finally:
            self._syncing = False
        self._reset_btn = pl.add_checkbox_button_widget(
            lambda s: self._reset_range(), value=False, position=(int(0.475 * w), int(0.195 * h)),
            size=16, color_on="#8a8f98", color_off="#8a8f98")
        self._label("auto", int(0.475 * w) + 22, int(0.195 * h) + 1, size=10, name="_reset_lbl")
        self._range_lbls = ["_reset_lbl"]

    def _range_visible(self, show):
        for wdg in (self._wmin, self._wmax, self._reset_btn):
            if wdg is None:
                continue
            try:
                (wdg.On if show else wdg.Off)()
                wdg.GetRepresentation().SetVisibility(1 if show else 0)
            except Exception:
                pass
        for nm in self._range_lbls:
            try:
                a = self._pl.actors.get(nm)
                if a is not None:
                    a.SetVisibility(1 if show else 0)
            except Exception:
                pass

    def _sync_range(self):
        """Show the color-range sliders when a scalar is on screen, hidden for geometry; reseed
        their bounds to the data range and their handles to the current (or auto) vmin/vmax."""
        if self._wmin is None or self._wmax is None:
            return
        rng = self._last_range
        show = bool(rng) and rng[1] > rng[0]
        self._syncing = True
        try:
            self._range_visible(show)
            if show:
                lo, hi = rng
                for wdg in (self._wmin, self._wmax):
                    rep = wdg.GetSliderRepresentation()
                    rep.SetMinimumValue(lo)
                    rep.SetMaximumValue(hi)
                vmin = lo if self.style.vmin is None else self.style.vmin
                vmax = hi if self.style.vmax is None else self.style.vmax
                self._wmin.GetSliderRepresentation().SetValue(min(max(vmin, lo), hi))
                self._wmax.GetSliderRepresentation().SetValue(min(max(vmax, lo), hi))
        finally:
            self._syncing = False

    def _on_vmin(self, v):
        if self._syncing:
            return
        hi = self.style.vmax if self.style.vmax is not None else (self._last_range or (v, v))[1]
        v = min(float(v), hi - abs(hi) * 1e-4 - 1e-12)
        self.style = self.style.copy(vmin=v)
        self._render_scene()

    def _on_vmax(self, v):
        if self._syncing:
            return
        lo = self.style.vmin if self.style.vmin is not None else (self._last_range or (v, v))[0]
        v = max(float(v), lo + abs(lo) * 1e-4 + 1e-12)
        self.style = self.style.copy(vmax=v)
        self._render_scene()

    def _reset_range(self):
        self.style = self.style.copy(vmin=None, vmax=None)
        self._render_scene()

    def _resync_radio(self, group, current):
        """Move a fixed-option group's highlight to `current` by RE-STATING its held widgets (no
        teardown/rebuild) -- used after an `m`/`n` key changes colormap/bonds so the on-screen
        selection stays honest.  Options and positions never change, so nothing is ever destroyed."""
        if not self.widgets:
            return
        was = self._building
        self._building = True                                # SetState must not fire callbacks
        for wdg, opt in self._radio_widgets.get(group, []):
            try:
                wdg.GetRepresentation().SetState(1 if opt == current else 0)
            except Exception:
                pass
        self._building = was

    def _cycle_cmap(self):
        """'m' key: advance colormap through VIEW_CMAPS (the same list the radio offers) and move
        the radio highlight to match, so keyboard and buttons never disagree."""
        cur = self.style.cmap if self.style.cmap in VIEW_CMAPS else VIEW_CMAPS[-1]
        nxt = VIEW_CMAPS[(VIEW_CMAPS.index(cur) + 1) % len(VIEW_CMAPS)]
        self.style = self.style.copy(cmap=nxt)
        self._resync_radio(_G_CMAP, nxt)
        self._render_scene()

    def _cycle_bonds(self):
        """'n' key: advance the bond-perception mode and resync its radio highlight."""
        from .style import BOND_MODES
        nxt = BOND_MODES[(BOND_MODES.index(self.style.bonds) + 1) % len(BOND_MODES)]
        self.style = self.style.copy(bonds=nxt)
        self._resync_radio(_G_BONDS, nxt)
        self._render_scene()

    # -- key bindings (kept as a fast path alongside the widgets) ---------------------
    def _bind(self):
        p = self._pl
        p.add_key_event("g", lambda: self._set("graph", not self.graph))
        p.add_key_event("p", lambda: self._set("paths", not self.paths))
        p.add_key_event("i", lambda: self._set("ias", not self.ias))
        p.add_key_event("k", lambda: self._set("fragment", not self.fragment))
        p.add_key_event("o", lambda: self._set("points_only", not self.points_only))
        p.add_key_event("h", lambda: self._set_style(show_h=not self.style.show_h))
        p.add_key_event("l", lambda: self._set_style(labels=not self.style.labels))
        p.add_key_event("v", lambda: self._set_style(show_values=not self.style.show_values))
        p.add_key_event("m", self._cycle_cmap)
        p.add_key_event("n", self._cycle_bonds)
        for key in ("plus", "equal", "KP_Add"):
            p.add_key_event(key, lambda: self._step_scale("atom_scale", 1.15))
        for key in ("minus", "KP_Subtract"):
            p.add_key_event(key, lambda: self._step_scale("atom_scale", 1 / 1.15))
        p.add_key_event("period", lambda: self._step_scale("bond_scale", 1.2))
        p.add_key_event("comma", lambda: self._step_scale("bond_scale", 1 / 1.2))
        p.add_key_event("bracketright", lambda: self._step(1))
        p.add_key_event("bracketleft", lambda: self._step(-1))

    def open(self, off_screen=False):
        """Create the plotter, paint the first scene, add widgets + key bindings; return the plotter
        WITHOUT entering the blocking event loop (used by `show` and by headless tests)."""
        import pyvista as pv
        self._pl = pv.Plotter(off_screen=off_screen, window_size=list(self.style.window_size),
                              title=WINDOW_TITLE)
        self._pl.background_color = self.style.background
        self._pl.enable_anti_aliasing("ssaa")
        self._first = True
        self._render_scene()
        if self.widgets:
            self._add_widgets()
        self._bind()
        if self.widgets:
            try:
                self._pl.ren_win.AddObserver("WindowResizeEvent", self._on_window_resize)
            except Exception:
                pass
        return self._pl

    # -- keep the on-canvas layout correct when the OS window is resized/maximized ----
    def _on_window_resize(self, obj, _evname):
        """Every widget/label above is placed in ABSOLUTE PIXELS computed once against
        `style.window_size` -- none of it follows the render surface when the user drags an edge
        or hits maximize (VTK resizes the 3D viewport and the normalized-coordinate scalar bar for
        free, but not pixel-anchored overlays), which is what actually misaligns.  Repaint the
        whole overlay fresh against whatever size VTK just reported."""
        try:
            new_size = tuple(int(v) for v in obj.GetSize())
        except Exception:
            return
        if new_size[0] < 2 or new_size[1] < 2 or new_size == tuple(self.style.window_size):
            return
        self._relayout(new_size)

    def _relayout(self, new_size):
        """Rebuild every widget + label at `new_size`, keeping the current structure/property/
        camera.  Uses pyvista's own bulk widget-clearers (they call `.Off()` on every button/radio/
        slider the plotter is tracking -- including the ones added fire-and-forget, like Prev/Next
        and the atom/bond-size sliders, that this class never kept a handle to) so nothing is ever
        garbage-collected while still live, THEN drops every actor and repaints from scratch --
        the same sequence `open()` runs, just re-entered at the new size instead of the first one."""
        if self._relayouting or self._pl is None:
            return
        self._relayouting = True
        try:
            pl = self._pl
            pl.clear_button_widgets()
            pl.clear_radio_button_widgets()
            pl.clear_slider_widgets()
            pl.clear_actors()
            self._pool_btns, self._nav_btns, self._slot_opt = [], [], []
            self._radio_widgets = {}
            self._prop_btn = self._entry_btn = None
            self._prop_anchor = self._entry_anchor = None
            self._panel = self._panel_bd = None
            self._wmin = self._wmax = self._reset_btn = self._wstruct = None
            self._range_lbls = []
            self._names, self._mol_names = [], []
            self._menu_mode = None      # the popup pool is gone; a resize just closes it
            self._menu_page = 0
            self.style = self.style.copy(window_size=new_size)
            self._render_scene()
            if self.widgets:
                self._add_widgets()
            pl.render()
        finally:
            self._relayouting = False

    def show(self):
        self.open(off_screen=False)
        print(WINDOW_TITLE + " -- " + self.HELP)
        self._pl.show(screenshot=self._save) if self._save else self._pl.show()
        return self._save


_QT_PLATFORM_USABLE = None


def _qt_platform_usable():
    """True only if Qt can actually initialize a windowing platform -- importing PySide6/pyvistaqt is
    NOT enough: a missing system lib (e.g. `libxcb-cursor0`, required by the Qt 6.5+ xcb plugin on
    Linux) makes ``QApplication()`` *abort the process*, which no in-process try/except can catch.  So
    probe it in a throwaway SUBPROCESS (inheriting DISPLAY / QT_QPA_PLATFORM) and cache the verdict."""
    global _QT_PLATFORM_USABLE
    if _QT_PLATFORM_USABLE is not None:
        return _QT_PLATFORM_USABLE
    import subprocess
    import sys
    try:
        r = subprocess.run(
            [sys.executable, "-c", "from qtpy import QtWidgets; QtWidgets.QApplication([])"],
            capture_output=True, timeout=25)
        _QT_PLATFORM_USABLE = (r.returncode == 0)
    except Exception:
        _QT_PLATFORM_USABLE = False
    return _QT_PLATFORM_USABLE


def _qt_backend_wanted(backend):
    """Pick the interactive backend.  The default is the self-contained **VTK** on-canvas `Viewer`
    (no Qt system-lib dependency); the native **Qt** panel is opt-in via ``backend='qt'`` or
    ``IQARIS_VIEWER=qt``.  `backend` = 'auto'/'vtk' (VTK) | 'qt' (Qt).  Returns True only for an
    explicit, working Qt request."""
    backend = (os.environ.get("IQARIS_VIEWER") or backend or "auto").lower()
    if backend != "qt":
        return False                     # 'auto' (default) and 'vtk' -> the VTK on-canvas viewer
    # explicit Qt request: use it only if it both imports AND can initialize a windowing platform
    # (a missing xcb system lib makes QApplication() abort, uncatchable in-process, so _qt_platform_usable
    # probes it in a subprocess).
    try:
        from .viz3d_qt import QT_AVAILABLE
    except Exception:
        QT_AVAILABLE = False
    if not QT_AVAILABLE:
        raise RuntimeError("backend='qt' requested but PySide6/pyvistaqt are not installed")
    if not _qt_platform_usable():
        raise RuntimeError(
            "backend='qt' requested but the Qt platform plugin cannot initialize -- on Linux this "
            "usually means the system lib 'libxcb-cursor0' is missing "
            "(`sudo apt-get install libxcb-cursor0`). Use the default VTK viewer instead.")
    return True


def view(sid=None, sids=None, prop=None, method=None, graph=False, paths=False, ias=False,
         fragment=False, fragment_prop=None, cmap=None, interactive=True, save=None,
         points_only=False, reframe=False, widgets=True, style=None, con=None, root=None,
         backend="auto", pbe0=False, **style_kw):
    """One-call viewer.  `interactive=True` opens the live explorer window -- by default the
    self-contained on-canvas **VTK** widget viewer (no extra system libs); the native **Qt** control
    panel (searchable entry picker + styled controls) is opt-in via ``backend='qt'`` / env
    ``IQARIS_VIEWER=qt``.  Either ALWAYS opens on plain atoms+bonds and browses `sids` with Prev/Next.
    `interactive=False` builds one Scene and renders a headless screenshot (this path DOES honor
    `prop`/`graph`/... for scripted images).  Style may be a `style=` object or loose kwargs (cmap=,
    atom_scale=, bonds=, show_h=, background=, ...).  `backend` = 'auto'/'vtk' (default VTK) | 'qt'."""
    from .db import connect
    from .project import build_scene, DEFAULT_METHOD
    from .style import Style
    method = method or DEFAULT_METHOD
    style = style or Style()
    if cmap:
        style_kw.setdefault("cmap", cmap)
    if style_kw:
        style = style.copy(**style_kw)

    own = con is None
    if own:
        con = connect(root)
    try:
        if interactive:
            if _qt_backend_wanted(backend):     # native Qt window (blocks until closed)
                from .viz3d_qt import open_qt_viewer
                return open_qt_viewer(con, sid=sid, sids=sids, method=method, style=style,
                                      points_only=points_only, save=save,
                                      fragment_prop=fragment_prop, pbe0=pbe0)
            return Viewer(con, sid=sid, sids=sids, prop=prop, method=method, graph=graph,
                          paths=paths, ias=ias, fragment=fragment, style=style,
                          points_only=points_only, widgets=widgets, save=save, pbe0=pbe0).show()
        scene = build_scene(con, sid or (sids[0] if sids else None), prop=prop, method=method,
                            graph=graph, paths=paths, ias=ias, fragment=fragment,
                            fragment_prop=fragment_prop, style=style,
                            side=("iqadft" if pbe0 else None))
        return emit_pyvista(scene, interactive=False, save=save, style=style,
                            points_only=points_only, reframe=reframe)
    finally:
        if own:
            con.close()
