"""project.py -- backend-neutral projection of IQARIS local properties onto a 3D structure.

Given a `structure_id` and (optionally) a property, `build_scene` assembles a *Scene*:
a renderer-independent bag of atoms, bonds, bond paths, critical points and interatomic-
surface disks, all in the DB's Angstrom frame.  A Scene is decoupled from any drawing
backend so the SAME scene can be turned into a publication-grade static PNG (`emit_jmol`,
headless Jmol) or an interactive IQARIS-Explorer window (`viz3d.emit_pyvista` /
`viz3d.Viewer`, PyVista/VTK).

Everything is read LIVE from the read-only IQARIS DuckDB views (geometry / qtaim_m062x_* /
qtaim_pbe0_* / iqa_pbe0_* / iqa_sqm_*); nothing is ever written back to the database.  A
prop's home view is taken from the registry, so ab initio IQA `*_IQA` columns color atoms /
bonds exactly like the QTAIM and IQA-SQM ones.  Paths (DB root, Jmol jar), the SQM level,
camera and output location are all parameters / env-vars.

CLI:  python -m iqaris_export project --sid ID --prop q_QTAIM [--method PM7]
                                      [--graph] [--paths] [--ias] [--fragment] [--out DIR]
"""
from __future__ import annotations
import collections
import os
import re
import subprocess
import textwrap
from dataclasses import dataclass, field

import numpy as np
from matplotlib import cm
from matplotlib.colors import (LinearSegmentedColormap, Normalize, TwoSlopeNorm,
                               to_hex)

from . import config
from .db import connect
from .registry import get_registry

# --- de-hardcoded defaults (overridable by env / kwargs) --------------------------------
JMOL_JAR = os.environ.get("IQARIS_JMOL_JAR", "/usr/share/java/JmolData.jar")
DEFAULT_METHOD = config.SQM_METHODS[0]                       # "PM7"
DEFAULT_SIZE = (2000, 1700)                                 # Jmol frame px (single render)
DEFAULT_VIEW = dict(rx=10, ry=8, rz=0, zoom=82)            # gentle tilt, auto-fit-ish
BOHR2 = config.BOHR_TO_ANGSTROM ** 2                        # bohr^2 -> Angstrom^2
K_KCAL = 627.509                                            # Hartree -> kcal/mol
COV_DI = 0.30                                              # delta(A,B) covalent-bond cutoff


# ========================================================= SCENE ============================
@dataclass
class Scene:
    """Renderer-independent description of one projected structure (Angstrom frame).

    atoms       : list of (pos[3], color_hex, spacefill_pct, label)
    cylinders   : list of (p1[3], p2[3], color_hex, diameter_ang)   -- bonds + bond-path halves
    points      : list of (pos[3], color_hex, diameter_ang, kind)   -- BCP / RCP / CCP dots
    disks       : list of (center[3], normal[3], radius_ang, color_hex, alpha) -- IAS plates
    colorbar    : (cmap, norm, label, unit) for the scalar coloring, or None -- `label` is the full
                  "prop (unit) [level]" description, `unit` just the bare axis unit (short, so a
                  narrow vertical color bar can carry it as a title without clipping)
    value_labels: list of (pos[3], text) -- the projected property's value, anchored on the
                  atom / bond midpoint / critical point it belongs to (only when show_values)
    """
    sid: str
    species: list = field(default_factory=list)
    atoms: list = field(default_factory=list)
    cylinders: list = field(default_factory=list)
    points: list = field(default_factory=list)
    disks: list = field(default_factory=list)
    colorbar: tuple = None
    value_labels: list = field(default_factory=list)
    title: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def natoms(self):
        return len(self.atoms)


# ========================================================= COLOR HELPERS =====================
def hexes(values, cmap, norm):
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    return [to_hex(m.to_rgba(v)) for v in values]


def jmol_hex(h):
    return "[x" + h.lstrip("#") + "]"


def sym_norm(vals):
    """Symmetric diverging norm centred at 0 (sign-aware energy / charge panels)."""
    v = float(np.nanmax(np.abs(vals))) or 1.0
    return TwoSlopeNorm(vcenter=0, vmin=-v, vmax=v)


def channel_style(vals):
    """Per-channel norm + cmap fitted to the data's ACTUAL range.  Sign is encoded by hue
    (blue = stabilizing/<0, red = destabilizing/>0) and the full colour range is used: a
    single-sign channel gets a sequential map (darker = stronger); a channel straddling 0
    keeps a diverging coolwarm centred at 0."""
    import matplotlib.pyplot as plt
    lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
    if lo < 0 < hi:
        return TwoSlopeNorm(vcenter=0, vmin=lo, vmax=hi), plt.cm.coolwarm
    if hi <= 0:
        cmap = LinearSegmentedColormap.from_list("bl", plt.cm.Blues_r(np.linspace(0.0, 0.82, 256)))
    else:
        cmap = LinearSegmentedColormap.from_list("rd", plt.cm.Reds(np.linspace(0.18, 1.0, 256)))
    return Normalize(lo, hi), cmap


def _truncated_viridis():
    import matplotlib.pyplot as plt
    return LinearSegmentedColormap.from_list("vtrunc", plt.cm.viridis(np.linspace(0.30, 1.0, 256)))


def plt_get(name):
    """Resolve a colormap name (or pass a Colormap through)."""
    import matplotlib.pyplot as plt
    return plt.get_cmap(name) if isinstance(name, str) else name


def _fmt_val(v):
    """Compact human string for a projected value shown next to an atom / bond / CP."""
    if v is None:
        return ""
    v = float(v)
    if np.isnan(v):
        return ""
    a = abs(v)
    if a != 0 and (a < 1e-2 or a >= 1e4):
        return f"{v:.2e}"
    if a < 1:
        return f"{v:.3f}"
    if a < 100:
        return f"{v:.2f}"
    return f"{v:.1f}"


# ========================================================= DATA FETCHERS =====================
class _Data:
    """SID-parameterized pulls over the read-only DuckDB views (same SQL as the A11 script,
    but against the catalog views instead of raw parquet globs)."""

    def __init__(self, con, sid, pbe0=False):
        self.con = con
        self.sid = sid
        # pbe0=True routes every QTAIM / QCT read to the PBE0 record (density topology + critical
        # points on the PBE0 density) instead of the M06-2X views; *_IQA and SQM reads are unaffected.
        self.pbe0 = pbe0
        geo = con.execute(
            "SELECT atom_index, element, x_ang, y_ang, z_ang FROM geometry "
            "WHERE structure_id = ? ORDER BY atom_index", [sid]).fetchdf()
        if len(geo) == 0:
            raise ValueError(f"structure_id {sid!r} not found in geometry")
        self.geo = geo
        self.N = len(geo)
        self.xyz = geo[["x_ang", "y_ang", "z_ang"]].values.astype(float)
        self.elem = geo["element"].tolist()
        self.aidx = geo["atom_index"].tolist()               # DB atom_index (row order)
        self.idx_of = {a: k for k, a in enumerate(self.aidx)}  # atom_index -> row

    def _qv(self, view):
        """Resolve a view name to the level-explicit schema.  Already-explicit names
        (qtaim_m062x_*, qtaim_pbe0_*, iqa_*) pass through; a legacy base name `qtaim_<sub>` maps to
        its M06-2X or (under self.pbe0) PBE0 twin; `mopac_<sub>` maps to `iqa_sqm_<sub>`."""
        if view.startswith(("qtaim_m062x_", "qtaim_pbe0_", "iqa_")):
            return view
        if view.startswith("qtaim_"):
            sub = view[len("qtaim_"):]
            return f"qtaim_pbe0_{sub}" if self.pbe0 else f"qtaim_m062x_{sub}"
        if view.startswith("mopac_"):
            return "iqa_sqm_" + view[len("mopac_"):]
        return view

    # -- per-atom scalar (auto-routes QTAIM vs level-aware SQM view) --------------------
    def atom_scalar(self, col, method=DEFAULT_METHOD):
        spec = self._spec(col, families=("atomic",))
        where = "structure_id = ?"
        params = [self.sid]
        if spec.level_aware:
            where += " AND method = ?"
            params.append(method)
        rows = self.con.execute(
            f"SELECT atom_index, {col} FROM {self._qv(spec.view)} WHERE {where}", params).fetchall()
        out = np.full(self.N, np.nan, float)
        for ai, v in rows:
            if ai in self.idx_of and v is not None:
                out[self.idx_of[ai]] = v
        return out, spec

    def _spec(self, col, families=None):
        reg = get_registry(self.con)
        specs = reg["by_name"].get(col)
        if not specs:
            raise ValueError(f"unknown property {col!r} (use `list-props` to discover)")
        if families:
            hit = [s for s in specs if s.family in families]
            if not hit:
                fam = specs[0].family
                raise ValueError(
                    f"property {col!r} is a {fam!r} column; `project --prop` maps per-ATOM "
                    f"scalars (family 'atomic'). Use --graph/--paths/--ias for pair/CP views.")
            specs = hit
        # side selection: under the PBE0 level prefer the 'iqadft' spec (its view is the
        # qtaim_pbe0_* / iqa_pbe0_* PBE0 reading); otherwise the M06-2X / SQM spec.  IQA-only columns
        # have a single ('iqadft') spec and resolve to it regardless.
        if self.pbe0:
            pref = [s for s in specs if s.side == "iqadft"]
        else:
            pref = [s for s in specs if s.side != "iqadft"]
        return (pref or specs)[0]

    # -- per-PAIR scalar over the pair views (auto-routes QTAIM vs level-aware SQM) ------
    def pair_scalar(self, col, method=DEFAULT_METHOD):
        """Return {(i, j): value} for pair column `col` (DB's i<j convention), plus the spec.
        Handles QTAIM and level-aware SQM pair columns alike."""
        spec = self._spec(col, families=("pairs",))
        where = "structure_id = ?"
        params = [self.sid]
        if spec.level_aware:
            where += " AND method = ?"
            params.append(method)
        rows = self.con.execute(
            f"SELECT i, j, {col} FROM {self._qv(spec.view)} WHERE {where}", params).fetchall()
        out = {}
        for i, j, val in rows:
            if val is not None:
                out[(int(i), int(j))] = float(val)
        return out, spec

    # -- QTAIM bonded pairs + delocalization index -------------------------------------
    def bonded_pairs(self):
        df = self.con.execute(
            f"SELECT i, j, DI_QTAIM FROM {self._qv('qtaim_pairs')} "
            "WHERE structure_id = ? AND BONDED_QTAIM", [self.sid]).fetchdf()
        return [(int(i), int(j), float(d)) for i, j, d in zip(df.i, df.j, df.DI_QTAIM)]

    # -- critical points of a given type -----------------------------------------------
    def cps(self, cp_type):
        return self.con.execute(
            "SELECT CP_X_ANG, CP_Y_ANG, CP_Z_ANG, CP_RHO_QTAIM, CP_BCP_I, CP_BCP_J "
            f"FROM {self._qv('qtaim_cp')} WHERE structure_id = ? AND CP_TYPE_QTAIM = ?",
            [self.sid, cp_type]).fetchdf()

    # -- IAS per BCP: position, bond-path normal (HESSEVEC_3), area (from the pair) -----
    def ias(self):
        return self.con.execute(
            "SELECT c.CP_X_ANG x, c.CP_Y_ANG y, c.CP_Z_ANG z, "
            "c.CP_HESSEVEC_3_X_QTAIM nx, c.CP_HESSEVEC_3_Y_QTAIM ny, c.CP_HESSEVEC_3_Z_QTAIM nz, "
            "p.AREA_IAS_002_QTAIM area "
            f"FROM (SELECT * FROM {self._qv('qtaim_cp')} WHERE structure_id = ? AND CP_TYPE_QTAIM='BCP') c "
            f"LEFT JOIN (SELECT i, j, AREA_IAS_002_QTAIM FROM {self._qv('qtaim_pairs')} WHERE structure_id = ?) p "
            "ON c.CP_BCP_I = p.i AND c.CP_BCP_J = p.j", [self.sid, self.sid]).fetchdf()

    # -- host<->guest per-atom engagement across the fragment interface ------------------
    def hostguest_engagement(self, host, guest, col="Einter_SQM", method=DEFAULT_METHOD):
        """Sum an *additive* interface pair property `col` over all host<->guest pairs onto each
        atom (each atom accumulates every interface pair it takes part in).  `col` defaults to the
        IQA total interaction energy but any eligible pair channel works (the IQA decomposition
        Eexchange/Eresonance/Eelstat[/ee/en/nn]/Edisp, or the QTAIM delocalization index DI_QTAIM).
        Energies (unit Eh) are returned in kcal/mol; dimensionless channels (DI) are kept as-is.
        Returns (per-atom array, total, spec)."""
        spec = self._spec(col, families=("pairs",))
        where = "structure_id = ?"
        params = [self.sid]
        if spec.level_aware:
            where += " AND method = ?"
            params.append(method)
        rows = self.con.execute(
            f"SELECT i, j, {col} FROM {spec.view} WHERE {where}", params).fetchall()
        hostset = set(host)
        scale = K_KCAL if spec.unit == "Eh" else 1.0        # energies -> kcal/mol; DI stays in e
        eng = np.zeros(self.N)
        total = 0.0
        for i, j, v in rows:
            if v is None:
                continue
            i, j = int(i), int(j)
            if (i in hostset) == (j in hostset):
                continue                                     # not a host<->guest (interface) pair
            val = float(v) * scale
            if i in self.idx_of:
                eng[self.idx_of[i]] += val
            if j in self.idx_of:
                eng[self.idx_of[j]] += val
            total += val
        return eng, total, spec


def covalent_fragments(aidx, cov_bonds):
    """Split atoms into connected components over the covalent (delta>cutoff) graph."""
    adj = collections.defaultdict(set)
    for i, j, _ in cov_bonds:
        adj[i].add(j)
        adj[j].add(i)
    seen, comps = set(), []
    for a in aidx:
        if a in seen:
            continue
        stack, comp = [a], set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.add(x)
            stack.extend(adj[x] - seen)
        comps.append(comp)
    comps.sort(key=len, reverse=True)
    return comps


# ========================================================= SCENE BUILDER =====================
_CPK = {  # minimal CPK fallback used when no scalar is projected
    "H": "#ffffff", "C": "#909090", "N": "#3050f8", "O": "#ff0d0d", "F": "#90e050",
    "P": "#ff8000", "S": "#ffff30", "Cl": "#1ff01f",
}
GREY_STICK = "#9aa0a6"
NODATA_STICK = "#c9ced4"   # thin neutral stick for a skeleton bond with NO value for the
                           # projected pair property (source data gap); keeps the graph intact
NODATA_ATOM = "#c9ced4"    # neutral atom color for a NULL per-atom scalar (e.g. d_pop on a
                           # non-d element) -- avoids a misleading solid-black "extreme" sphere

# "fragments" pseudo-property: color each atom by which covalent fragment (connected component of
# the DI>0.30 graph) it belongs to, so the user can tell the pieces of a complex apart.  Not a DB
# column -- handled specially in build_scene.  Fragments are ordered largest-first.
FRAGMENTS_PROP = "fragments"
FRAGMENT_PALETTE = [        # qualitative, colorblind-friendly (seaborn "deep"-ish), largest-first
    "#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860",
    "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd",
]

# Pair properties that can be summed across the host<->guest interface onto atoms (additive
# interaction channels, fully populated over the interface): the IQA-SQM interaction-energy
# decomposition + the QTAIM delocalization index.  The IAS-integrated pair columns are excluded --
# they only exist at the few real interface BCPs, so aggregating them yields a near-empty map.
FRAGMENT_PAIR_PROPS = frozenset({
    "Einter_SQM", "Eexchange_SQM", "Eresonance_SQM", "Eelstat_SQM",
    "Eelstat_ee_SQM", "Eelstat_en_SQM", "Eelstat_nn_SQM", "Edisp_SQM", "DI_QTAIM",
    # PBE0 ab initio IQA pair channels: the host<->guest interface can also be aggregated over the
    # IQA interaction energy and its classical / exchange-correlation split.
    "EINT_IQA", "VCL_IJ_IQA", "VXC_IJ_IQA", "VNE_IJ_IQA", "VEN_IJ_IQA", "VEE_IJ_IQA"})


def _formula(elems):
    """Hill-ish molecular formula (C, H first) from a list of element symbols."""
    c = collections.Counter(elems)
    order = ["C", "H", "N", "O", "F", "P", "S", "Cl"]
    keys = [e for e in order if e in c] + sorted(e for e in c if e not in order)
    return "".join(f"{e}{c[e]}" if c[e] > 1 else e for e in keys)


def _prop_short(col):
    """Strip the family suffix from a column name for a compact label (Einter_SQM -> Einter)."""
    for suf in ("_SQM", "_QTAIM"):
        if col.endswith(suf):
            return col[:-len(suf)]
    return col
CP_COLORS = {"BCP": "#141414", "RCP": "#e6262b", "CCP": "#d21fce"}
CP_DIAM = {"BCP": 0.16, "RCP": 0.20, "CCP": 0.22}
COV_RADII = {  # covalent radii (Angstrom) for distance-based bond perception
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02,
}


def perceive_bonds(d, mode="di", cutoff=None):
    """Return the list of drawn bonds `[(i, j, di)]` under a perception `mode`:

    - 'qtaim'    : every QTAIM-bonded pair (a BCP / bond path exists), weak contacts included.
    - 'di'       : QTAIM-bonded pairs with delocalization index DI > `cutoff` (default 0.30) --
                   the covalent skeleton; this is the historical default.
    - 'distance' : classic covalent-radii perception, dist < (r_i + r_j) * `cutoff` (tol 1.20),
                   independent of QTAIM (finds bonds even where no BCP was reported); DI attached
                   when the pair is also in the QTAIM table, else NaN.
    """
    bp = d.bonded_pairs()
    if mode == "qtaim":
        return bp
    if mode == "di":
        c = 0.30 if cutoff is None else cutoff
        return [(i, j, di) for (i, j, di) in bp if di > c]
    if mode == "distance":
        tol = 1.20 if cutoff is None else cutoff
        di_of = {(min(i, j), max(i, j)): di for (i, j, di) in bp}
        out = []
        for a in range(d.N):
            for b in range(a + 1, d.N):
                ra = COV_RADII.get(d.elem[a], 0.77)
                rb = COV_RADII.get(d.elem[b], 0.77)
                if np.linalg.norm(d.xyz[a] - d.xyz[b]) < (ra + rb) * tol:
                    ia, ib = d.aidx[a], d.aidx[b]
                    out.append((ia, ib, di_of.get((min(ia, ib), max(ia, ib)), float("nan"))))
        return out
    raise ValueError(f"unknown bond perception {mode!r} (use qtaim | di | distance)")


def _pad_deg(lo, hi):
    """Widen a degenerate (constant / single-valued) range so the scalar lands mid-colormap instead
    of collapsing to the dark end -- `Normalize(v, v)` maps everything to 0, i.e. black atoms/bonds
    and an all-black colorbar (e.g. two Cl atoms with an identical d population)."""
    if hi > lo:
        return lo, hi
    pad = abs(lo) * 0.5 or 0.5
    return lo - pad, hi + pad


def _scalar_style(vals, style):
    """(cmap, norm) for a scalar channel, honoring `style` overrides (cmap / vmin / vmax /
    vcenter) and otherwise falling back to the auto choice: symmetric diverging RdBu when the
    values straddle zero, sequential mako otherwise."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    arr = np.asarray(vals, float)
    if np.isfinite(arr).any():                              # a property may be all-NaN for a
        dlo, dhi = float(np.nanmin(arr)), float(np.nanmax(arr))  # given structure -> avoid the
    else:                                                   # "All-NaN slice" warning + degenerate
        dlo, dhi = 0.0, 1.0                                 # colorbar; show a flat 0..1 scale
    lo = dlo if style.vmin is None else style.vmin
    hi = dhi if style.vmax is None else style.vmax
    overridden = any(x is not None for x in (style.vmin, style.vmax, style.vcenter))
    signed = (lo < 0 < hi) or style.vcenter is not None
    if style.cmap:
        cmap = plt.get_cmap(style.cmap) if isinstance(style.cmap, str) else style.cmap
    elif signed:
        cmap = plt.cm.RdBu
    else:
        cmap = sns.color_palette("mako", as_cmap=True)
    if signed and not overridden:
        norm = sym_norm(arr)                                # symmetric ±max|v| (current default)
    elif signed:
        center = 0.0 if style.vcenter is None else style.vcenter
        norm = (TwoSlopeNorm(vcenter=center, vmin=lo, vmax=hi)
                if lo < center < hi else Normalize(*_pad_deg(lo, hi)))
    else:
        norm = Normalize(*_pad_deg(lo, hi))
    return cmap, norm


def _set_range(scene, vals):
    """Record the true data (min, max) of the active scalar so an interactive viewer can seed its
    color-range sliders; skipped when everything is NaN."""
    arr = np.asarray(vals, float)
    if arr.size and np.isfinite(arr).any():
        scene.meta["data_range"] = (float(np.nanmin(arr)), float(np.nanmax(arr)))


def molecular_values(con, sid, side="qtaim", method=DEFAULT_METHOD):
    """Every whole-molecule scalar property for `sid` on one property side, as a list of
    ``(name, value, unit)`` rows (registry order).  These per-molecule totals (energies, CP census,
    dipole, HOMO/LUMO, IQA molecular energy, ...) cannot be projected onto atoms, so the viewers
    surface them through a "molecular properties" popup instead of the coloring dropdown.

    side: 'qtaim' -> qtaim_m062x_molecular (M06-2X) | 'iqadft' -> qtaim_pbe0_molecular +
          iqa_pbe0_molecular (PBE0: QTAIM topology + ab initio IQA) | 'sqm' -> iqa_sqm_molecular
          for `method`.
    """
    reg = get_registry(con)
    specs = reg["by_group"].get(("molecular", side), [])
    if not specs:
        return []
    # group by each spec's actual view: the PBE0 (iqadft) side spans qtaim_pbe0_molecular
    # (topology) + iqa_pbe0_molecular (ab initio IQA), so query each view for its own columns.
    from collections import OrderedDict
    by_view = OrderedDict()
    for s in specs:
        by_view.setdefault(s.view, []).append(s)
    out = []
    for view, vspecs in by_view.items():
        where, params = "structure_id = ?", [sid]
        if side == "sqm":
            where += " AND method = ?"
            params.append(method)
        quoted = ", ".join(f'"{s.name}"' for s in vspecs)
        row = con.execute(f"SELECT {quoted} FROM {view} WHERE {where} LIMIT 1", params).fetchone()
        if row is not None:
            out.extend((s.name, v, s.unit) for s, v in zip(vspecs, row))
    return out


def build_scene(con, sid, prop=None, method=DEFAULT_METHOD, graph=False, paths=False,
                ias=False, fragment=False, fragment_prop=None, cmap=None, style=None, side=None):
    """Assemble a backend-neutral :class:`Scene` for `sid`.

    prop     : per-atom scalar column to color atoms by (e.g. q_QTAIM, q_SQM); None -> CPK.
    method   : SQM level for level-aware props / fragment engagement.
    side     : property side driving the QTAIM/QCT functional -- 'iqadft' reads QTAIM density
               topology, critical points and IAS from the PBE0 record and tags scalars `[PBE0]`;
               'qtaim'/'sqm'/None keep the M06-2X QTAIM views (default).
    paths    : draw the molecular-graph bond paths (atom -> BCP -> atom), colored by rho_b.
    graph    : draw ALL critical points (BCP dark / RCP red / CCP magenta) + bond paths.
    ias      : draw interatomic (zero-flux) surfaces as area-scaled disks perpendicular to
               the bond path.
    fragment : if the covalent graph splits into >=2 pieces, color atoms by an aggregated
               host<->guest interface pair property instead of `prop` (default off), and tint
               the two fragment skeletons distinctly.
    fragment_prop : which additive interface channel to aggregate in fragment mode -- any of
               :data:`FRAGMENT_PAIR_PROPS` (the IQA-SQM interaction decomposition Einter / Eexchange
               / Eresonance / Eelstat[/ee/en/nn] / Edisp, or the QTAIM delocalization index
               DI_QTAIM).  None -> use `prop` if it is eligible, else default to Einter_SQM.
    style    : :class:`~iqaris_export.style.Style` bundle of appearance / perception knobs
               (colormap, color range, atom & bond sizing, bond perception, hydrogen/label
               visibility). None -> defaults (reproduces the current look). Legacy `cmap`
               kwarg still honored (folded into `style.cmap`).
    """
    from .style import Style
    style = style or Style()
    if cmap and not style.cmap:                             # back-compat: legacy cmap= kwarg
        style = style.copy(cmap=cmap)

    d = _Data(con, sid, pbe0=(side == "iqadft"))
    show_paths = paths or graph
    scene = Scene(sid=sid, species=list(d.elem))
    scene.meta.update(natoms=d.N, method=method)

    # fragmentation always uses the physical covalent graph (DI>0.30); the DRAWN bonds obey the
    # chosen perception mode so the two concerns stay independent.
    cov_frag = [(i, j, di) for (i, j, di) in d.bonded_pairs() if di > COV_DI]
    comps = covalent_fragments(d.aidx, cov_frag)
    nfrag = len(comps)
    host = comps[0] if nfrag >= 2 else set(d.aidx)
    guest = comps[1] if nfrag >= 2 else set()
    draw_bonds = perceive_bonds(d, style.bonds, style.bond_cutoff)
    scene.meta.update(nfrag=nfrag, nbonds=len(draw_bonds), ncov=len(cov_frag),
                      bond_mode=style.bonds)

    hide = {a for a, e in zip(d.aidx, d.elem) if e == "H"} if not style.show_h else set()

    # ---- classify the requested property: atomic -> atom color, pairs -> bond color ---
    use_fragment = fragment and nfrag >= 2
    frag_color = prop == FRAGMENTS_PROP and not use_fragment   # categorical per-fragment coloring
    atom_prop = bond_prop = None
    if prop is not None and not use_fragment and not frag_color:
        # A column may live in several views (e.g. Eexchange_SQM is both a molecular total and a
        # per-pair term); prefer the projectable reading -- atomic colors atoms, pairs colors bonds.
        fams = {s.family for s in get_registry(d.con)["by_name"].get(prop, [])}
        if not fams:
            raise ValueError(f"unknown property {prop!r} (use `list-props` to discover)")
        if "atomic" in fams:
            atom_prop = prop
        elif "pairs" in fams:
            bond_prop = prop                                # DI_QTAIM, IQA-SQM pair energies, ...
        else:
            fam = sorted(fams)[0]
            raise ValueError(
                f"property {prop!r} is a {fam!r} column; `--prop` colors atoms (atomic) or "
                f"bonds (pairs). Use --graph for the molecular graph or --ias for surfaces.")

    # ---- atom coloring ---------------------------------------------------------------
    atom_values = None                                      # per-row scalar shown by show_values
    if use_fragment:
        # aggregate an eligible interface pair channel across the host<->guest pairs.  Prefer the
        # explicit `fragment_prop`, else the picked `prop` if it is itself eligible, else default to
        # the IQA total interaction energy Einter_SQM.
        chosen = (fragment_prop if fragment_prop in FRAGMENT_PAIR_PROPS
                  else prop if prop in FRAGMENT_PAIR_PROPS else None)
        frag_col = chosen or "Einter_SQM"
        frag_default = chosen is None                        # True when Einter is the fallback
        vals, total, fspec = d.hostguest_engagement(host, guest, col=frag_col, method=method)
        cm_, norm = _scalar_style(vals, style)
        colors = hexes(vals, cm_, norm)
        unit = "kcal mol$^{-1}$" if fspec.unit == "Eh" else fspec.unit
        short = _prop_short(frag_col)
        scene.colorbar = (cm_, norm, rf"host$\leftrightarrow$guest $\sum${short} / {unit}", unit)
        scene.title = f"{sid}: host-guest {short} ($\\Sigma$={total:+.1f})"
        host_els = [d.elem[d.idx_of[a]] for a in host]
        guest_els = [d.elem[d.idx_of[a]] for a in guest]
        scene.meta.update(hg_total=total, hg_col=frag_col, hg_default=frag_default,
                          hg_level=fspec.level_aware, hg_host_n=len(host), hg_guest_n=len(guest),
                          hg_host_formula=_formula(host_els), hg_guest_formula=_formula(guest_els))
        atom_values = vals
        _set_range(scene, vals)
    elif frag_color:
        # categorical: each covalent fragment (connected component) gets a distinct color, largest
        # first.  No scalar -> no colorbar / no range slider; the header reports the component count.
        comp_of = {a: ci for ci, comp in enumerate(comps) for a in comp}
        colors = [FRAGMENT_PALETTE[comp_of.get(a, 0) % len(FRAGMENT_PALETTE)] for a in d.aidx]
        sizes = [len(c) for c in comps]
        formulae = [_formula([d.elem[d.idx_of[a]] for a in c]) for c in comps]
        scene.title = f"{sid}: fragments ({nfrag} component{'s' if nfrag != 1 else ''})"
        scene.meta.update(nfrag=nfrag, frag_sizes=sizes, frag_formulae=formulae)
    elif atom_prop is not None:
        vals, spec = d.atom_scalar(atom_prop, method)
        cm_, norm = _scalar_style(vals, style)
        colors = hexes(vals, cm_, norm)
        # atoms with NO value for this column (e.g. d_pop on a non-d element) map to the colormap's
        # "bad" color = solid black, which reads as a real extreme; recolor them neutral "no-data".
        nan_mask = ~np.isfinite(np.asarray(vals, float))
        if nan_mask.any():
            colors = [NODATA_ATOM if bad else c for c, bad in zip(colors, nan_mask)]
        scene.meta["natom_nodata"] = int(nan_mask.sum())
        lvl = f" [{method}]" if spec.level_aware else (" [PBE0]" if (d.pbe0 or spec.side == "iqadft") else "")
        scene.colorbar = (cm_, norm, f"{atom_prop} ({spec.unit}){lvl}", spec.unit)
        scene.title = f"{sid}: {atom_prop}{lvl}"
        atom_values = vals
        _set_range(scene, vals)
    else:
        colors = [_CPK.get(e, GREY_STICK) for e in d.elem]  # CPK atoms as context for a bond map
        scene.title = f"{sid}: {bond_prop} (bond)" if bond_prop else f"{sid}"

    # translucent "ghost" structure under a molecular graph or IAS overlay (unless a per-atom
    # scalar is the message, in which case keep the atoms opaque and readable)
    ghost = (show_paths or ias) and atom_prop is None and not use_fragment and not frag_color
    spacefill = (14.0 if ghost else 26.0) * style.atom_scale
    # Atom size is purely physical (vdW radius x style.atom_scale); it never encodes a property
    # value -- only color carries the data.
    for k, (pos, col, el) in enumerate(zip(d.xyz, colors, d.elem)):
        if d.aidx[k] in hide:
            continue
        scene.atoms.append((pos, col, spacefill, el))
        if style.show_values and atom_values is not None:   # value on the atom it belongs to
            scene.value_labels.append((pos, _fmt_val(atom_values[k])))
    scene.meta["ghost"] = ghost
    scene.meta["atom_translucent"] = 0.62 if ghost else None

    # ---- bonds / bond paths (skip any touching a hidden hydrogen) --------------------
    # value labels annotate the ONE primary scalar so the 3D view stays legible: atoms win over
    # bonds win over the molecular graph (CP rho) win over interatomic surfaces.
    cps_primary = atom_prop is None and not use_fragment
    ias_primary = cps_primary and bond_prop is None and not show_paths
    vis = [(i, j, x) for (i, j, x) in draw_bonds if i not in hide and j not in hide]
    if bond_prop is not None:                               # color sticks by a per-pair scalar
        _add_bond_prop(scene, d, vis, bond_prop, method, style)
    elif show_paths:
        _add_bond_paths(scene, d, cov_frag, style, hide,
                        label_values=style.show_values and cps_primary)
    elif not ias:                                            # plain grey covalent sticks
        for (i, j, _) in vis:
            if i in d.idx_of and j in d.idx_of:
                scene.cylinders.append(
                    (d.xyz[d.idx_of[i]], d.xyz[d.idx_of[j]], GREY_STICK, 0.14 * style.bond_scale))

    # ---- critical points -------------------------------------------------------------
    if graph:
        for t in ("BCP", "RCP", "CCP"):
            df = d.cps(t)
            for r in df.itertuples():
                if not style.show_h and t == "BCP" and (
                        int(r.CP_BCP_I) in hide or int(r.CP_BCP_J) in hide):
                    continue
                scene.points.append(
                    (np.array([r.CP_X_ANG, r.CP_Y_ANG, r.CP_Z_ANG]),
                     CP_COLORS[t], CP_DIAM[t] * style.cp_scale, t))

    # ---- interatomic surfaces --------------------------------------------------------
    if ias:
        _add_ias(scene, d, style, label_values=style.show_values and ias_primary)

    return scene


def _add_bond_prop(scene, d, cov, col, method, style):
    """Color and width-scale bonds by a per-PAIR scalar (e.g. DI_QTAIM bond order, or an
    IQA-SQM pair energy): straight atom->atom sticks, colormap/range from `style` (or auto:
    sequential mako for a single-sign channel, diverging RdBu centred at 0 for a signed one),
    with the stick thickness growing with the value so strong bonds read at a glance."""
    vals, spec = d.pair_scalar(col, method)
    valued, missing = [], []
    for (i, j, _di) in cov:
        if i not in d.idx_of or j not in d.idx_of:
            continue
        (valued if (i, j) in vals else missing).append((i, j))
    # Bonds with no value for this pair column are drawn as thin neutral "no-data" sticks (not
    # dropped) so the graph stays complete -- IAS-integrated pair columns are populated only for a
    # subset of BCPs.
    for (i, j) in missing:
        scene.cylinders.append(
            (d.xyz[d.idx_of[i]], d.xyz[d.idx_of[j]], NODATA_STICK, 0.09 * style.bond_scale))
    scene.meta.update(nbond_valued=len(valued), nbond_nodata=len(missing))
    if not valued:                                          # column empty for every drawn bond
        return
    v = np.array([vals[(i, j)] for (i, j) in valued], float)
    cm_, norm = _scalar_style(v, style)
    cols = hexes(v, cm_, norm)
    nv = np.clip(norm(v), 0.0, 1.0)
    for (i, j), c, w in zip(valued, cols, nv):
        pi, pj = d.xyz[d.idx_of[i]], d.xyz[d.idx_of[j]]
        diam = (0.10 + 0.34 * float(w)) * style.bond_scale
        scene.cylinders.append((pi, pj, c, diam))
        if style.show_values:                              # value on the interatomic line midpoint
            scene.value_labels.append(((pi + pj) / 2, _fmt_val(vals[(i, j)])))
    lvl = f" [{method}]" if spec.level_aware else (" [PBE0]" if (d.pbe0 or spec.side == "iqadft") else "")
    scene.colorbar = (cm_, norm, f"{col} ({spec.unit}){lvl}", spec.unit)
    _set_range(scene, v)


def _add_bond_paths(scene, d, cov, style, hide=frozenset(), label_values=False):
    """Bond-path cylinders (atom -> BCP -> atom), width & color ~ rho at the BCP; covalent
    opaque, weak/non-covalent paths translucent -- the QTAIM molecular graph skeleton.
    `label_values` annotates each BCP with its rho (when the graph is the projected scalar)."""
    df = d.cps("BCP")
    if len(df) == 0:
        return
    rho = df["CP_RHO_QTAIM"].values.astype(float)
    lo = 0.0 if style.vmin is None else style.vmin
    hi = float(style.vmax if style.vmax is not None else np.nanmax(rho)) or 1.0
    norm = Normalize(lo, hi)
    cmap = plt_get(style.cmap) if style.cmap else _truncated_viridis()
    cols = hexes(rho, cmap, norm)
    nv = np.clip(norm(rho), 0, 1)
    covset = {frozenset((i, j)) for (i, j, _) in cov}
    for k, r in enumerate(df.itertuples()):
        bp = np.array([r.CP_X_ANG, r.CP_Y_ANG, r.CP_Z_ANG])
        bi, bj = int(r.CP_BCP_I), int(r.CP_BCP_J)
        if bi in hide or bj in hide:
            continue
        is_cov = frozenset((bi, bj)) in covset
        diam = (0.04 + 0.28 * float(nv[k])) * style.bond_scale * style.path_scale
        for a in (bi, bj):
            if a in d.idx_of:
                p = d.xyz[d.idx_of[a]]
                # translucency for weak paths is carried as a 5th flag on the tuple
                scene.cylinders.append((p, bp, cols[k], diam, None if is_cov else 0.62))
        if label_values:                                   # rho on the critical point
            scene.value_labels.append((bp, _fmt_val(rho[k])))
    # bond-path colors are baked into the cylinders; only own the scene colorbar when no
    # per-atom scalar already drives it (pure molecular-graph case).
    if scene.colorbar is None:
        scene.colorbar = (cmap, norm, r"$\rho_b$ / a.u.", "a.u.")
        _set_range(scene, rho)


def _add_ias(scene, d, style, label_values=False):
    """Interatomic (zero-flux) surfaces as thin disks at each BCP, oriented perpendicular to
    the bond path (normal = HESSEVEC_3) and sized by the real AIMALL surface area.
    `label_values` annotates each surface with its area (when IAS is the projected scalar)."""
    import seaborn as sns
    df = d.ias()
    df = df[df.area.notna()].copy()
    if len(df) == 0:
        return
    area_ang = df.area.values * BOHR2
    lo = style.vmin if style.vmin is not None else np.nanpercentile(area_ang, 3)
    hi = style.vmax if style.vmax is not None else np.nanpercentile(area_ang, 97)
    norm = Normalize(lo, hi)
    cmap = plt_get(style.cmap) if style.cmap else sns.color_palette("mako", as_cmap=True)
    cols = hexes(area_ang, cmap, norm)
    for k, (r, col) in enumerate(zip(df.itertuples(), cols)):
        c = np.array([r.x, r.y, r.z])
        n = np.array([r.nx, r.ny, r.nz], float)
        n = n / (np.linalg.norm(n) + 1e-12)
        rad = float(np.sqrt(r.area * BOHR2 / np.pi)) * 0.5
        # translucent enough that overlapping surfaces stay readable (many stacked disks at higher
        # opacity compound to a muddy near-opaque mass over the molecule)
        scene.disks.append((c, n, rad, col, 0.4))
        if label_values:                                   # area on the interatomic surface
            scene.value_labels.append((c, _fmt_val(area_ang[k])))
    if scene.colorbar is None:
        scene.colorbar = (cmap, norm, r"IAS area / $\mathrm{\AA}^2$", r"$\mathrm{\AA}^2$")
        _set_range(scene, area_ang)


# ========================================================= JMOL EMITTER ======================
def _principal_frame(xyz):
    """Rotation to the heavy-atom principal-axis frame (deterministic; det +1 keeps chirality)."""
    c = xyz.mean(0)
    _, _, vt = np.linalg.svd(xyz - c, full_matrices=True)
    if np.linalg.det(vt) < 0:
        vt[-1] *= -1
    return c, vt


def _header(view, background="white"):
    return textwrap.dedent(f"""
        background {background}
        set frank off
        set antialiasDisplay true
        set antialiasImages true
        set ambientPercent 45
        set diffusePercent 84
        set specularPercent 22
        set specularExponent 6
        set zoomLarge false
        set autoBond off
        rotate x {view['rx']}
        rotate y {view['ry']}
        rotate z {view['rz']}
        zoom {view['zoom']}
    """).strip("\n").split("\n")


def _run_jmol(lines, png, size, jar):
    sptf = png + ".spt"
    with open(sptf, "w") as f:
        f.write("\n".join(lines).replace("__PNG__", png)
                .replace("__W__", str(size[0])).replace("__H__", str(size[1])))
    subprocess.run(["java", "-Xmx6g", "-Djava.awt.headless=true", "-jar", jar,
                    "-ionx", "-s", sptf], check=True, capture_output=True)
    if not os.path.exists(png):
        raise RuntimeError(f"Jmol produced no image: {png}")


def emit_jmol(scene, out, jar=None, size=None, view=None, reframe=True, style=None):
    """Render a :class:`Scene` to a transparent-background PNG with headless Jmol.

    Reframes to principal axes (static-figure convention) unless `reframe=False`.  Returns
    the PNG path.  De-hardcoded: jar / size / camera all overridable; the output filename is
    derived from the structure id, never the A11-specific `a11_*` path.  `style` (optional)
    supplies the background color and whether to draw element labels.
    """
    from .style import Style
    style = style or Style()
    jar = jar or JMOL_JAR
    size = size or DEFAULT_SIZE
    view = dict(DEFAULT_VIEW if view is None else view)
    if not os.path.exists(jar):
        raise RuntimeError(f"Jmol jar not found: {jar} (set IQARIS_JMOL_JAR or pass jar=...)")

    xyz = np.array([a[0] for a in scene.atoms], float)
    if reframe and len(xyz) >= 2:
        c, vt = _principal_frame(xyz)
        pt = lambda p: (np.asarray(p, float) - c) @ vt.T
        vec = lambda v: np.asarray(v, float) @ vt.T
    else:
        pt = lambda p: np.asarray(p, float)
        vec = lambda v: np.asarray(v, float)

    os.makedirs(out, exist_ok=True)
    tag = _slug(scene.sid) + ("_" + _slug(scene.meta.get("proptag", "")) if scene.meta.get("proptag") else "")
    stem = os.path.join(out, tag or _slug(scene.sid))
    xyzf = stem + ".xyz"
    with open(xyzf, "w") as f:
        f.write(f"{len(scene.atoms)}\n{scene.sid}\n")
        for (pos, _c, _sf, el), pfr in zip(scene.atoms, (pt(a[0]) for a in scene.atoms)):
            f.write(f"{el:2s} {pfr[0]:14.6f} {pfr[1]:14.6f} {pfr[2]:14.6f}\n")

    lines = [f'load "{xyzf}"', *_header(view, style.background)]
    # atoms: group by spacefill %, color per atom
    sf_groups = collections.defaultdict(list)
    for k, (_pos, _col, sf, _el) in enumerate(scene.atoms):
        sf_groups[sf].append(k)
    for sf, ks in sf_groups.items():
        sel = " or ".join(f"atomno={k + 1}" for k in ks)
        lines.append(f"select ({sel}); spacefill {sf:.0f}%")
    lines.append("wireframe off")
    transl = scene.meta.get("atom_translucent")
    if transl:
        lines.append(f"select all; color atoms translucent {transl}")
    for k, (_pos, col, _sf, _el) in enumerate(scene.atoms):
        lines.append(f"select atomno={k + 1}; color atom {jmol_hex(col)}")
    if style.labels:
        lines.append('select all; label "%e"; color labels black; '
                     'set fontSize 16; font label 16 SansSerif Bold')

    # cylinders (bonds + bond paths); optional translucency as a 5th tuple element
    for n, cyl in enumerate(scene.cylinders):
        p1, p2, col, diam = cyl[0], cyl[1], cyl[2], cyl[3]
        alpha = cyl[4] if len(cyl) > 4 else None
        a, b = pt(p1), pt(p2)
        colspec = jmol_hex(col) if not alpha else f"TRANSLUCENT {alpha} {jmol_hex(col)}"
        lines.append(
            f'draw ID c{n} CYLINDER {{{a[0]:.4f} {a[1]:.4f} {a[2]:.4f}}} '
            f'{{{b[0]:.4f} {b[1]:.4f} {b[2]:.4f}}} DIAMETER {diam:.3f} COLOR {colspec}')

    # critical-point dots (drawn spheres -- never dummy atoms, which blacken cylinders)
    for n, (pos, col, diam, _kind) in enumerate(scene.points):
        q = pt(pos)
        lines.append(f'draw ID p{n} DIAMETER {diam:.3f} '
                     f'{{{q[0]:.4f} {q[1]:.4f} {q[2]:.4f}}} COLOR {jmol_hex(col)}')

    # IAS disks (thin fat cylinders oriented along the bond-path normal)
    for n, (center, normal, rad, col, alpha) in enumerate(scene.disks):
        c0 = pt(center)
        nn = vec(normal)
        nn = nn / (np.linalg.norm(nn) + 1e-12)
        p1 = c0 - nn * 0.035
        p2 = c0 + nn * 0.035
        lines.append(
            f'draw ID d{n} CYLINDER {{{p1[0]:.4f} {p1[1]:.4f} {p1[2]:.4f}}} '
            f'{{{p2[0]:.4f} {p2[1]:.4f} {p2[2]:.4f}}} DIAMETER {2 * rad:.3f} '
            f'COLOR TRANSLUCENT {alpha} {jmol_hex(col)}')

    lines.append('write IMAGE __W__ __H__ PNGT "__PNG__"')
    png = stem + ".png"
    _run_jmol(lines, png, size, jar)
    return png


def _slug(s):
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(s)).strip("_")


# ========================================================= HIGH-LEVEL API ====================
def project(sid, prop=None, method=DEFAULT_METHOD, graph=False, paths=False, ias=False,
            fragment=False, fragment_prop=None, cmap=None, out=None, root=None, con=None, jar=None,
            size=None, view=None, style=None, pbe0=False):
    """One-call: connect (read-only) if needed, build the scene, render a static Jmol PNG.

    Returns the PNG path.  `out` defaults to `<exports>/projections`.  `style` (a
    :class:`~iqaris_export.style.Style`) controls appearance / bond perception.  `pbe0=True`
    renders the PBE0 reading of shared QTAIM/QCT properties.
    """
    own = con is None
    if own:
        con = connect(root)
    try:
        scene = build_scene(con, sid, prop=prop, method=method, graph=graph, paths=paths,
                            ias=ias, fragment=fragment, fragment_prop=fragment_prop, cmap=cmap,
                            style=style, side=("iqadft" if pbe0 else None))
        tags = [scene.meta.get("hg_col") if scene.meta.get("hg_col")
                else prop or ("hg" if (fragment and scene.meta.get("nfrag", 1) >= 2) else "geom")]
        if graph:
            tags.append("graph")
        elif paths:
            tags.append("paths")
        if ias:
            tags.append("ias")
        scene.meta["proptag"] = "_".join(t for t in tags if t)
        out = out or str(config.DEFAULT_OUT / "projections")
        return emit_jmol(scene, out, jar=jar, size=size, view=view, style=style)
    finally:
        if own:
            con.close()
