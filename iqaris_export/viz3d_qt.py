"""viz3d_qt.py -- native **Qt** (PySide6 / pyvistaqt) control panel for the IQARIS-Explorer viewer.

The 07-03/07-05 viewer drew its controls with VTK's on-canvas widgets (checkbox squares, radio
dots, on-canvas text).  Those have a hard aesthetic ceiling and -- fatally for browsing 268k
structures -- no text input, so picking an entry meant paging.  This module keeps the exact 3D
representation (same `project.build_scene` + `viz3d.add_scene` on an embedded VTK render window)
but moves EVERY control into a real, styled Qt side panel:

  * a **searchable entry picker** -- type an id or a formula fragment, the matching structures
    appear in a list, click to jump (DB-side `ILIKE`, so it scales to the whole database);
  * native **combo boxes** for property / level / colormap / bond perception (the level rescopes
    the property list to the columns that exist there, QTAIM under M06-2X vs IQA-SQM under a PM
    method); native **checkboxes** for the graph/paths/IAS/fragment/points/H/labels/values
    overlays; native **sliders** for atom & bond size and the color-scale min/max.

The window ALWAYS opens on the plain atoms+bonds view (no property, no overlay); coloring and
overlays are turned on from the panel.  Switching structures reframes the camera to the new
molecule; changing a property/overlay on the same molecule keeps the current rotation.
"""
from __future__ import annotations
import os

from . import config
from .style import Style, BOND_MODES
from .registry import get_registry
from .project import build_scene, FRAGMENTS_PROP
from .viz3d import (add_scene, WINDOW_TITLE, FEATURED_QTAIM, FEATURED_SQM, FEATURED_BY_SIDE,
                    side_for_level, VIEW_CMAPS, LEVEL_DISPLAY, _level_label, _GEOM)

try:                                            # Qt is an OPTIONAL dependency (lazy)
    from qtpy import QtCore, QtGui, QtWidgets
    from qtpy.QtCore import Qt
    from pyvistaqt import QtInteractor
    QT_AVAILABLE = True
except Exception:                               # pragma: no cover - absence handled by the caller
    QT_AVAILABLE = False


# professional flat stylesheet: soft light panels, one blue accent, crisp native font ------------
_ACCENT = "#2c7fb8"
_QSS = f"""
* {{ font-family: "Inter","Segoe UI","Helvetica Neue","DejaVu Sans",sans-serif; font-size: 11pt; }}
QWidget#panel {{ background: #eef1f4; }}
QScrollArea, QScrollArea > QWidget > QWidget {{ background: #eef1f4; border: none; }}
QGroupBox {{
    border: 1px solid #dbe0e6; border-radius: 10px; margin-top: 16px; padding: 10px 10px 6px 10px;
    background: #ffffff;
}}
QGroupBox::title {{
    subcontrol-origin: margin; left: 12px; padding: 0 5px;
    color: #55606b; font-size: 9.5pt; font-weight: 700; text-transform: uppercase;
    letter-spacing: 0.6px;
}}
QLabel {{ color: #22303c; background: transparent; }}
QLabel#title {{ font-size: 15pt; font-weight: 700; color: #14202b; }}
QLabel#subtitle {{ color: #6a7683; font-size: 10pt; }}
QLabel#field {{ color: #55606b; font-size: 9.5pt; font-weight: 600; }}
QComboBox, QLineEdit {{
    border: 1px solid #ccd3da; border-radius: 7px; padding: 5px 9px; background: #ffffff;
    selection-background-color: {_ACCENT};
}}
QComboBox:hover, QLineEdit:hover {{ border-color: #a9b4bf; }}
QComboBox:focus, QLineEdit:focus {{ border-color: {_ACCENT}; }}
QComboBox::drop-down {{ subcontrol-origin: padding; subcontrol-position: center right;
    width: 20px; border: none; }}
QComboBox::down-arrow {{
    width: 0; height: 0; margin-right: 8px;
    border-left: 5px solid transparent; border-right: 5px solid transparent;
    border-top: 6px solid #6a7683;
}}
QComboBox QAbstractItemView {{
    border: 1px solid #ccd3da; background: #ffffff; selection-background-color: {_ACCENT};
    selection-color: #ffffff; outline: none;
}}
QListWidget {{ border: 1px solid #ccd3da; border-radius: 7px; background: #ffffff; padding: 2px; }}
QListWidget::item {{ padding: 4px 6px; border-radius: 4px; }}
QListWidget::item:hover {{ background: #eaf2f8; }}
QListWidget::item:selected {{ background: {_ACCENT}; color: #ffffff; }}
QCheckBox {{ spacing: 7px; color: #22303c; padding: 2px; }}
QCheckBox::indicator {{ width: 16px; height: 16px; border: 1px solid #b8c1cb; border-radius: 4px; background: #fff; }}
QCheckBox::indicator:checked {{ background: {_ACCENT}; border-color: {_ACCENT}; }}
QPushButton {{
    background: {_ACCENT}; color: #ffffff; border: none; border-radius: 7px; padding: 6px 14px;
    font-weight: 600;
}}
QPushButton:hover {{ background: #2a72a4; }}
QPushButton:pressed {{ background: #245f88; }}
QPushButton:disabled {{ background: #cfd6dd; color: #eef1f4; }}
QPushButton#ghost {{ background: #ffffff; color: #55606b; border: 1px solid #ccd3da; }}
QPushButton#ghost:hover {{ border-color: {_ACCENT}; color: {_ACCENT}; }}
QSlider::groove:horizontal {{ height: 4px; background: #d7dde3; border-radius: 2px; }}
QSlider::sub-page:horizontal {{ background: {_ACCENT}; border-radius: 2px; }}
QSlider::handle:horizontal {{
    background: #ffffff; border: 2px solid {_ACCENT}; width: 14px; height: 14px;
    margin: -7px 0; border-radius: 8px;
}}
QSlider:disabled {{ }}
QSlider::sub-page:horizontal:disabled {{ background: #cfd6dd; }}
QSlider::handle:horizontal:disabled {{ border-color: #cfd6dd; }}
"""


def _prop_disp(name):
    """Drop the `_QTAIM`/`_SQM` suffix for the combo label (the level already says the side)."""
    if name == _GEOM:
        return name
    for suf in ("_QTAIM", "_SQM"):
        if name.endswith(suf):
            return name[:-len(suf)]
    return name


if QT_AVAILABLE:

    def _esc(s):
        return (s or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    def _mathtext_safe(latex):
        """Nudge a LaTeX string into the subset matplotlib mathtext understands (so the popup
        can typeset it); the unicode `eq_text` is the fallback when even this fails to parse."""
        return (latex.replace(r"\boldsymbol", r"\mathbf")
                     .replace(r"\tfrac", r"\frac")
                     .replace(r"\!", ""))

    def _equation_pixmap(latex, dpi=170):
        """Typeset a single-line equation via matplotlib mathtext into a QPixmap, or None if it
        cannot be parsed (the caller then shows the plain-text/unicode form instead)."""
        try:
            import io
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            fig = Figure(figsize=(0.01, 0.01))
            FigureCanvasAgg(fig)
            fig.text(0.0, 0.0, f"${_mathtext_safe(latex)}$", fontsize=16, color="#14202b")
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                        pad_inches=0.06, transparent=True)
            buf.seek(0)
            pm = QtGui.QPixmap()
            pm.loadFromData(buf.getvalue(), "PNG")
            return pm if not pm.isNull() else None
        except Exception:
            return None

    class _InfoDialog(QtWidgets.QDialog):
        """Popup describing ONE property: name, symbol, unit, type, level of theory, prose, and a
        typeset equation (with a unicode fallback).  Fed by `explore.property_info`."""

        def __init__(self, info, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Property info")
            self.setStyleSheet(_QSS)
            # Pin the width so every popup is the same size and the word-wrapped labels have a
            # definite height-for-width; `adjustSize()` at the end then computes the height needed at
            # that width, so the description is never clipped (as it was before).
            self._fixed_w = 500
            self.setFixedWidth(self._fixed_w)
            v = QtWidgets.QVBoxLayout(self)
            v.setContentsMargins(20, 18, 20, 16)
            v.setSpacing(8)

            title = QtWidgets.QLabel(info.get("name") or info["column"])
            title.setObjectName("title")
            title.setWordWrap(True)
            v.addWidget(title)

            sym = info.get("symbol") or ""
            sub = QtWidgets.QLabel(
                f"<code>{_esc(info['column'])}</code>"
                + (f" &nbsp;&nbsp;{_esc(sym)}" if sym else "")
                + ("  &middot;  bonus" if info.get("bonus") else ""))
            sub.setObjectName("subtitle")
            sub.setTextFormat(Qt.RichText)
            v.addWidget(sub)

            meta = "  ·  ".join(x for x in (
                f"Unit: {info['unit']}" if info.get("unit") else "",
                f"Type: {info.get('type', '')}" if info.get("type") else "") if x)
            if meta:
                m = QtWidgets.QLabel(meta)
                m.setObjectName("field")
                v.addWidget(m)
            if info.get("method"):
                mm = QtWidgets.QLabel(f"Level of theory: {info['method']}")
                mm.setObjectName("subtitle")
                mm.setWordWrap(True)
                v.addWidget(mm)

            if info.get("description"):
                d = QtWidgets.QLabel(info["description"])
                d.setWordWrap(True)
                d.setStyleSheet("color:#22303c; font-size:10.5pt; padding-top:3px;")
                v.addWidget(d)

            if info.get("eq_latex") or info.get("eq_text"):
                v.addWidget(self._eq_widget(info.get("eq_latex"), info.get("eq_text")))

            if info.get("see_also"):
                sa = QtWidgets.QLabel("See also:  " + ",  ".join(info["see_also"]))
                sa.setObjectName("subtitle")
                sa.setWordWrap(True)
                v.addWidget(sa)

            row = QtWidgets.QHBoxLayout()
            row.addStretch(1)
            ok = QtWidgets.QPushButton("Close")
            ok.clicked.connect(self.accept)
            row.addWidget(ok)
            v.addLayout(row)
            # Width is already pinned (setFixedWidth), so adjustSize computes the height needed at
            # that width; lock that size so the window opens fitting its content, same width every time.
            self.adjustSize()
            self.setFixedSize(self._fixed_w, self.height())

        def _eq_widget(self, eq_latex, eq_text):
            box = QtWidgets.QFrame()
            box.setStyleSheet("background:#f5f7f9; border:1px solid #e3e7ec; border-radius:8px;")
            lay = QtWidgets.QVBoxLayout(box)
            lay.setContentsMargins(12, 10, 12, 10)
            pm = _equation_pixmap(eq_latex) if eq_latex else None
            if pm is not None:
                max_w = 448                       # keep wide integrals inside the 520px popup width
                if pm.width() > max_w:
                    pm = pm.scaledToWidth(max_w, Qt.SmoothTransformation)
                lbl = QtWidgets.QLabel()
                lbl.setPixmap(pm)
                lbl.setAlignment(Qt.AlignCenter)
            else:
                lbl = QtWidgets.QLabel(eq_text or "")
                lbl.setAlignment(Qt.AlignCenter)
                lbl.setWordWrap(True)
                lbl.setStyleSheet("font-family:'DejaVu Sans Mono',monospace; "
                                  "font-size:11pt; color:#14202b; background:transparent; border:none;")
            lay.addWidget(lbl)
            return box

    def _fmt_molval(v):
        """Human value for the molecular table: ints stay ints, floats to ~6 sig figs, None -> dash."""
        if v is None:
            return "—"
        if isinstance(v, bool):
            return "true" if v else "false"
        try:
            f = float(v)
        except (TypeError, ValueError):
            return str(v)
        if f == int(f) and abs(f) < 1e15:
            return str(int(f))
        return f"{f:.6g}"

    class _MolecularDialog(QtWidgets.QDialog):
        """Popup table of every whole-molecule scalar property + value for the current structure and
        level (fed by project.molecular_values).  These per-molecule totals cannot be projected onto
        atoms, so they live here rather than in the coloring dropdown."""

        def __init__(self, title, rows, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Molecular properties")
            self.setStyleSheet(_QSS)
            self.setMinimumWidth(460)
            v = QtWidgets.QVBoxLayout(self)
            v.setContentsMargins(18, 16, 18, 14)
            v.setSpacing(8)
            t = QtWidgets.QLabel(title)
            t.setObjectName("title")
            t.setWordWrap(True)
            v.addWidget(t)
            sub = QtWidgets.QLabel(f"{len(rows)} whole-molecule scalar{'s' if len(rows) != 1 else ''}")
            sub.setObjectName("subtitle")
            v.addWidget(sub)

            table = QtWidgets.QTableWidget(len(rows), 3)
            table.setHorizontalHeaderLabels(["Property", "Value", "Unit"])
            table.verticalHeader().setVisible(False)
            table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            table.setAlternatingRowColors(True)
            for r, (name, val, unit) in enumerate(rows):
                table.setItem(r, 0, QtWidgets.QTableWidgetItem(str(name)))
                vi = QtWidgets.QTableWidgetItem(_fmt_molval(val))
                vi.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                table.setItem(r, 1, vi)
                table.setItem(r, 2, QtWidgets.QTableWidgetItem(
                    "" if unit in (None, "label") else str(unit)))
            table.resizeColumnsToContents()
            table.horizontalHeader().setStretchLastSection(True)
            v.addWidget(table)

            row = QtWidgets.QHBoxLayout()
            row.addStretch(1)
            ok = QtWidgets.QPushButton("Close")
            ok.clicked.connect(self.accept)
            row.addWidget(ok)
            v.addLayout(row)
            self.resize(500, 560)

    class QtViewer(QtWidgets.QMainWindow):
        """The IQARIS-Explorer main window: a Qt control panel docked left of an embedded PyVista/VTK
        render view.  Reuses `build_scene`/`add_scene` verbatim, so the molecule looks identical to
        the static path; only the chrome is native Qt."""

        def __init__(self, con, sid=None, sids=None, method="PM7", style=None,
                     points_only=False, save=None, fragment_prop=None, parent=None, pbe0=False):
            super().__init__(parent)
            self.con = con
            self.sids = list(sids) if sids else ([sid] if sid else [])
            if not self.sids:
                raise ValueError("QtViewer needs a sid or a non-empty sids list")
            self.idx = self.sids.index(sid) if (sid and sid in self.sids) else 0
            self._sqm = method if method in config.SQM_METHODS else config.SQM_METHODS[0]
            # ALWAYS open clean: plain atoms + bonds, no property, no overlay.
            self.prop = None
            from .project import FRAGMENT_PAIR_PROPS
            self._frag_prop = (fragment_prop if fragment_prop in FRAGMENT_PAIR_PROPS
                               else "Einter_SQM")   # interface channel aggregated in fragment mode
            self._frag_active = False         # host<->guest coloring actually engaged (nfrag>=2)
            self._combo_frag_mode = False     # which option set the property combo currently holds
            self.level = "PBE0" if pbe0 else "DFT"   # starting level
            self.graph = self.paths = self.ias = self.fragment = False
            self.points_only = points_only
            self.style = style or Style()
            self._save = save
            self._first = True
            self._refit = True           # reframe camera on the next render (new molecule / first)
            self._syncing = False        # mute control signals during programmatic updates
            self._rng = None             # (lo, hi) data range of the active scalar, or None
            self._last_range = None
            self.setWindowTitle(WINDOW_TITLE)
            self.setStyleSheet(_QSS)
            self.resize(1320, 820)
            self._build_ui()
            self._render()
            self._refresh_prop_combo()
            self._update_nav()
            self._on_search("")               # seed the entry list with the current browse set

        # -- properties ----------------------------------------------------------------
        @property
        def sid(self):
            return self.sids[self.idx]

        @property
        def method(self):
            return self.level if self.level in config.SQM_METHODS else self._sqm

        @property
        def _side(self):
            return side_for_level(self.level)

        def _prop_options(self):
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
            # `fragments` is a structural pseudo-property (color atoms by covalent component); it is
            # level-independent, so it sits right after plain geometry at every level.
            return [_GEOM, FRAGMENTS_PROP] + featured + rest

        def _frag_prop_options(self):
            """Interface channels that fragment mode can aggregate (see project.FRAGMENT_PAIR_PROPS),
            **restricted to the current level of theory** (like the normal property picker): the
            IQA-SQM decomposition exists only at an SQM method, the ab initio IQA interaction split
            only under PBE0, the delocalization index only at the DFT/QTAIM level -- so DI never
            shows up under PM6, Einter never under M06-2X, and EINT_IQA only under PBE0."""
            from .project import FRAGMENT_PAIR_PROPS
            reg = get_registry(self.con)
            have = {s.name for s in reg["by_group"].get(("pairs", self._side), [])}
            order = ["Einter_SQM", "Eelstat_SQM", "Eexchange_SQM", "Eresonance_SQM", "Edisp_SQM",
                     "Eelstat_ee_SQM", "Eelstat_en_SQM", "Eelstat_nn_SQM",
                     "EINT_IQA", "VCL_IJ_IQA", "VXC_IJ_IQA", "VNE_IJ_IQA", "VEN_IJ_IQA", "VEE_IJ_IQA",
                     "DI_QTAIM"]
            return [c for c in order if c in FRAGMENT_PAIR_PROPS and c in have]

        # -- UI ------------------------------------------------------------------------
        def _build_ui(self):
            central = QtWidgets.QWidget()
            root = QtWidgets.QHBoxLayout(central)
            root.setContentsMargins(0, 0, 0, 0)
            root.setSpacing(0)

            # left: scrollable control panel
            panel = QtWidgets.QWidget()
            panel.setObjectName("panel")
            panel.setFixedWidth(330)
            col = QtWidgets.QVBoxLayout(panel)
            col.setContentsMargins(14, 14, 14, 14)
            col.setSpacing(9)
            self._build_header(col)
            self._build_structure_group(col)
            self._build_representation_group(col)
            self._build_display_group(col)
            self._build_range_group(col)
            col.addStretch(1)
            scroll = QtWidgets.QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(panel)
            scroll.setFixedWidth(346)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

            # right: the embedded VTK render view.  (A bare QtInteractor starts with NO lights, so
            # spheres render flat/matte; PyVista's 5-light "light kit" is (re)applied every render in
            # `_render` -- `clear()` drops the lights -- so shading/specular match `emit_pyvista`/Jmol.)
            self.plotter = QtInteractor(central)
            self.plotter.set_background(self.style.background)
            try:
                self.plotter.enable_anti_aliasing("ssaa")
            except Exception:
                pass

            root.addWidget(scroll)
            root.addWidget(self.plotter.interactor, 1)
            self.setCentralWidget(central)

        def _build_header(self, col):
            self.title_lbl = QtWidgets.QLabel(self.sid)
            self.title_lbl.setObjectName("title")
            self.sub_lbl = QtWidgets.QLabel("")
            self.sub_lbl.setObjectName("subtitle")
            col.addWidget(self.title_lbl)
            col.addWidget(self.sub_lbl)

        def _field(self, text):
            lab = QtWidgets.QLabel(text)
            lab.setObjectName("field")
            return lab

        def _build_structure_group(self, col):
            g = QtWidgets.QGroupBox("Structure")
            v = QtWidgets.QVBoxLayout(g)
            v.setSpacing(7)
            # prev / next + position
            nav = QtWidgets.QHBoxLayout()
            self.prev_btn = QtWidgets.QPushButton("‹  Prev")
            self.prev_btn.setObjectName("ghost")
            self.next_btn = QtWidgets.QPushButton("Next  ›")
            self.next_btn.setObjectName("ghost")
            self.pos_lbl = QtWidgets.QLabel("")
            self.pos_lbl.setObjectName("subtitle")
            self.pos_lbl.setAlignment(Qt.AlignCenter)
            self.prev_btn.clicked.connect(lambda: self._step(-1))
            self.next_btn.clicked.connect(lambda: self._step(1))
            nav.addWidget(self.prev_btn)
            nav.addWidget(self.pos_lbl, 1)
            nav.addWidget(self.next_btn)
            v.addLayout(nav)
            # searchable entry picker
            v.addWidget(self._field("jump to entry  (type id or formula)"))
            self.search = QtWidgets.QLineEdit()
            self.search.setPlaceholderText("e.g.  C6H6   or   qm9_1087 ...")
            self.search.setClearButtonEnabled(True)
            self.search.textChanged.connect(self._on_search)
            v.addWidget(self.search)
            self.results = QtWidgets.QListWidget()
            self.results.setFixedHeight(124)
            self.results.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
            self.results.setTextElideMode(Qt.ElideRight)
            self.results.itemActivated.connect(self._on_pick)
            self.results.itemClicked.connect(self._on_pick)
            v.addWidget(self.results)
            col.addWidget(g)

        def _combo(self, items_display_data):
            c = QtWidgets.QComboBox()
            for disp, data in items_display_data:
                c.addItem(disp, data)
            return c

        def _build_representation_group(self, col):
            g = QtWidgets.QGroupBox("Representation")
            form = QtWidgets.QFormLayout(g)
            form.setLabelAlignment(Qt.AlignLeft)
            form.setSpacing(8)
            # property
            self.prop_combo = QtWidgets.QComboBox()
            self.prop_combo.currentIndexChanged.connect(self._on_prop)
            # level
            self.level_combo = self._combo(
                [(_level_label(o), o) for o in (["DFT", "PBE0"] + list(config.SQM_METHODS))])
            self.level_combo.currentIndexChanged.connect(self._on_level)
            # colormap
            self.cmap_combo = self._combo([(c, c) for c in VIEW_CMAPS])
            self.cmap_combo.currentIndexChanged.connect(self._on_cmap)
            # bonds
            self.bond_combo = self._combo([(b, b) for b in BOND_MODES])
            self.bond_combo.setCurrentIndex(list(BOND_MODES).index(self.style.bonds))
            self.bond_combo.currentIndexChanged.connect(self._on_bonds)
            # property row = combo + an "ⓘ" button that pops up the property's definition
            self.info_btn = QtWidgets.QPushButton("ⓘ")
            self.info_btn.setObjectName("ghost")
            self.info_btn.setFixedWidth(34)
            self.info_btn.setCursor(Qt.PointingHandCursor)
            self.info_btn.setToolTip("Definition, units and equation of the selected property")
            self.info_btn.clicked.connect(self._on_prop_info)
            prop_row = QtWidgets.QWidget()
            prl = QtWidgets.QHBoxLayout(prop_row)
            prl.setContentsMargins(0, 0, 0, 0)
            prl.setSpacing(6)
            prl.addWidget(self.prop_combo, 1)
            prl.addWidget(self.info_btn)
            form.addRow(self._field("property"), prop_row)
            form.addRow(self._field("level"), self.level_combo)
            form.addRow(self._field("colormap"), self.cmap_combo)
            form.addRow(self._field("bonds"), self.bond_combo)
            # full-width button -> popup of ALL whole-molecule scalars (the totals that can't color atoms)
            self.mol_btn = QtWidgets.QPushButton("Molecular properties…")
            self.mol_btn.setObjectName("ghost")
            self.mol_btn.setCursor(Qt.PointingHandCursor)
            self.mol_btn.setToolTip("List every whole-molecule scalar (energies, CP census, dipole,\n"
                                    "HOMO/LUMO, IQA molecular energy, …) for this structure at the\n"
                                    "current level -- the totals that can't be colored onto atoms.")
            self.mol_btn.clicked.connect(self._on_molecular)
            form.addRow(self.mol_btn)
            col.addWidget(g)

        def _build_display_group(self, col):
            g = QtWidgets.QGroupBox("Display")
            grid = QtWidgets.QGridLayout(g)
            grid.setHorizontalSpacing(20)
            grid.setVerticalSpacing(9)
            grid.setColumnStretch(0, 1)
            grid.setColumnStretch(1, 1)
            self._checks = {}
            specs = [
                ("graph", "graph"), ("paths", "paths"), ("IAS", "ias"),
                ("fragment", "fragment"), ("points", "points_only"),
                ("H", "show_h"), ("labels", "labels"), ("values", "show_values"),
            ]
            tips = {
                "fragment": "Split into host + guest and color atoms by their summed host↔guest\n"
                            "interaction across the interface; the two fragment skeletons are tinted\n"
                            "distinctly (host slate, guest ochre). Needs ≥2 covalent fragments. While\n"
                            "on, the property menu picks the aggregated channel (Einter / Eelstat /\n"
                            "Eexchange / Eresonance / Edisp, or QTAIM DI).",
            }
            for i, (label, key) in enumerate(specs):
                cb = QtWidgets.QCheckBox(label)
                init = getattr(self.style, key) if key in ("show_h", "labels", "show_values") \
                    else getattr(self, key)
                cb.setChecked(bool(init))
                cb.toggled.connect(lambda st, k=key: self._on_toggle(k, st))
                if key in tips:
                    cb.setToolTip(tips[key])
                grid.addWidget(cb, i // 2, i % 2)
                self._checks[key] = cb
            col.addWidget(g)

        def _slider(self, lo, hi, val):
            s = QtWidgets.QSlider(Qt.Horizontal)
            s.setMinimum(lo)
            s.setMaximum(hi)
            s.setValue(val)
            # Commit (re-render) only when the handle is RELEASED, not on every intermediate value:
            # with tracking OFF, `valueChanged` fires on release / keyboard / click-step, while
            # `sliderMoved` still fires live during the drag (used only to update the number label).
            s.setTracking(False)
            return s

        def _pct(self, x):
            return f"{x / 100.0:.2f}"

        def _range_live(self, x):
            """Live label for a color-range slider: map its 0-1000 position onto the data range."""
            if not self._rng:
                return "--"
            lo, hi = self._rng
            return f"{lo + (hi - lo) * x / 1000.0:.3g}"

        def _divider(self, text):
            """A thin rule + a small muted caption, to visually split a group into subsections."""
            w = QtWidgets.QWidget()
            lay = QtWidgets.QVBoxLayout(w)
            lay.setContentsMargins(0, 5, 0, 1)
            lay.setSpacing(5)
            line = QtWidgets.QFrame()
            line.setFixedHeight(1)
            line.setStyleSheet("background:#e3e7ec; border:none;")
            lay.addWidget(line)
            lay.addWidget(self._field(text))
            return w

        def _build_range_group(self, col):
            g = QtWidgets.QGroupBox("Sizing and color")
            v = QtWidgets.QVBoxLayout(g)
            v.setSpacing(8)
            # atom size / bond size (0.30 - 3.00, stored x100)
            self.atom_slider, self.atom_val = self._labeled_slider(
                v, "atom size", 30, 300, int(self.style.atom_scale * 100),
                lambda x: self._set_scale("atom_scale", x / 100.0), self._pct)
            self.bond_slider, self.bond_val = self._labeled_slider(
                v, "bond size", 30, 300, int(self.style.bond_scale * 100),
                lambda x: self._set_scale("bond_scale", x / 100.0), self._pct)
            # QCT molecular-graph sizing -- ENABLED ONLY while the graph / paths overlay is shown
            self.cp_slider, self.cp_val = self._labeled_slider(
                v, "CP size", 30, 300, int(self.style.cp_scale * 100),
                lambda x: self._set_scale("cp_scale", x / 100.0), self._pct)
            self.path_slider, self.path_val = self._labeled_slider(
                v, "path size", 30, 300, int(self.style.path_scale * 100),
                lambda x: self._set_scale("path_scale", x / 100.0), self._pct)
            # color min / max (0-1000 mapped onto the data range; enabled only for a scalar)
            v.addWidget(self._divider("color range"))
            self.vmin_slider, self.vmin_val = self._labeled_slider(
                v, "color min", 0, 1000, 0, self._on_vmin, self._range_live)
            self.vmax_slider, self.vmax_val = self._labeled_slider(
                v, "color max", 0, 1000, 1000, self._on_vmax, self._range_live)
            self.auto_btn = QtWidgets.QPushButton("auto range")
            self.auto_btn.setObjectName("ghost")
            self.auto_btn.clicked.connect(self._reset_range)
            v.addSpacing(2)
            v.addWidget(self.auto_btn)
            col.addWidget(g)

        def _labeled_slider(self, parent_layout, label, lo, hi, val, cb, live):
            row = QtWidgets.QHBoxLayout()
            name = self._field(label)
            name.setMinimumWidth(74)
            vallbl = QtWidgets.QLabel("")
            vallbl.setObjectName("subtitle")
            vallbl.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            vallbl.setMinimumWidth(52)
            slider = self._slider(lo, hi, val)
            slider.valueChanged.connect(cb)                                # commit (render) on release
            slider.sliderMoved.connect(lambda x: vallbl.setText(live(x)))  # live number, NO render
            row.addWidget(name)
            row.addWidget(slider, 1)
            row.addWidget(vallbl)
            parent_layout.addLayout(row)
            return slider, vallbl

        # -- combos populate / sync ----------------------------------------------------
        def _refresh_prop_combo(self):
            """Repopulate the property combo and reselect the active property.  In fragment mode the
            combo lists the aggregatable interface channels (and drives which one is summed); otherwise
            it lists the projectable atomic/pair columns for the current level."""
            self._syncing = True
            self.prop_combo.blockSignals(True)
            self.prop_combo.clear()
            frag = self._frag_active
            if frag:
                opts = self._frag_prop_options()
                for opt in opts:
                    self.prop_combo.addItem(_prop_disp(opt), opt)
                cur = self._frag_prop if self._frag_prop in opts else (opts[0] if opts else None)
                self._frag_prop = cur          # reconcile: the level may have dropped the old channel
                                               # (e.g. DI_QTAIM at DFT -> switch to PM6), so keep the
                                               # aggregated channel in step with what the combo shows
            else:
                opts = self._prop_options()
                for opt in opts:
                    self.prop_combo.addItem(_prop_disp(opt), opt)
                cur = self.prop if self.prop else _GEOM
            i = self.prop_combo.findData(cur)
            self.prop_combo.setCurrentIndex(max(0, i))
            self.prop_combo.blockSignals(False)
            self._combo_frag_mode = frag
            self._syncing = False

        # -- control callbacks ---------------------------------------------------------
        def _on_prop(self, _i):
            if self._syncing:
                return
            data = self.prop_combo.currentData()
            if self._frag_active:                 # combo drives the aggregated interface channel
                self._frag_prop = data
            else:
                self.prop = None if data == _GEOM else data
            self.style = self.style.copy(vmin=None, vmax=None)
            self._render()

        def _on_prop_info(self):
            """Pop up the definition / units / equation of the currently selected property."""
            data = self.prop_combo.currentData()
            if not data or data == _GEOM:
                QtWidgets.QToolTip.showText(
                    self.info_btn.mapToGlobal(QtCore.QPoint(self.info_btn.width(), 0)),
                    "Pick a property (not plain geometry) to see its definition.")
                return
            if data == FRAGMENTS_PROP:                     # synthetic (not a DB column)
                _InfoDialog({
                    "column": FRAGMENTS_PROP, "name": "Fragments",
                    "type": "structural", "unit": "categorical",
                    "description": "Colors each atom by which covalent fragment (connected component "
                                   "of the delocalization-index > 0.30 bond graph) it belongs to, so "
                                   "the separate pieces of a complex -- e.g. host and guest -- can be "
                                   "told apart. Level-independent; no numeric scale."},
                    parent=self).exec()
                return
            from .explore import property_info
            try:
                info = property_info(self.con, data)
            except Exception:
                info = {"column": data}
            _InfoDialog(info, parent=self).exec()

        def _on_molecular(self):
            """Pop up every whole-molecule scalar for the current structure at the current level."""
            from .project import molecular_values
            side = self._side
            rows = molecular_values(self.con, self.sid, side=side, method=self.method)
            lvl = (_level_label("DFT") if side == "qtaim"
                   else _level_label("PBE0") if side == "iqadft" else self.method)
            if not rows:
                QtWidgets.QToolTip.showText(
                    self.mol_btn.mapToGlobal(QtCore.QPoint(self.mol_btn.width(), 0)),
                    "No molecular data at this level for this structure.")
                return
            _MolecularDialog(f"{self.sid}  —  molecular  [{lvl}]", rows, parent=self).exec()

        def _on_level(self, _i):
            if self._syncing:
                return
            level = self.level_combo.currentData()
            if level == self.level:
                return
            old_side = self._side
            self.level = level
            if level in config.SQM_METHODS:
                self._sqm = level
            if self._side != old_side:
                self.prop = self._remap_prop()
                self.style = self.style.copy(vmin=None, vmax=None)
                self._refresh_prop_combo()
            self._render()

        def _remap_prop(self):
            options = self._prop_options()
            if self.prop is None or self.prop == FRAGMENTS_PROP:
                return self.prop                          # geometry / structural: level-independent
            if self.prop in options:
                return self.prop                          # same column exists here (a QTAIM prop
                                                          # shared by M06-2X and PBE0) -> keep it
            cand = None
            if self._side == "sqm" and self.prop.endswith("_QTAIM"):
                cand = self.prop[:-len("_QTAIM")] + "_SQM"
            elif self._side == "qtaim" and self.prop.endswith("_SQM"):
                cand = self.prop[:-len("_SQM")] + "_QTAIM"
            if cand and cand in options:
                return cand
            # No QTAIM<->SQM suffix analog (e.g. entering the PBE0 ab initio IQA side): fall to the
            # first REAL property, skipping the structural `fragments` pseudo-column.
            rest = [p for p in options if p not in (_GEOM, FRAGMENTS_PROP)]
            return rest[0] if rest else None

        def _on_cmap(self, _i):
            if self._syncing:
                return
            self.style = self.style.copy(cmap=self.cmap_combo.currentData())
            self._render()

        def _on_bonds(self, _i):
            if self._syncing:
                return
            self.style = self.style.copy(bonds=self.bond_combo.currentData())
            self._render()

        def _on_toggle(self, key, state):
            if self._syncing:
                return
            if key in ("show_h", "labels", "show_values"):
                self.style = self.style.copy(**{key: bool(state)})
            else:
                setattr(self, key, bool(state))
            self._render()

        def _set_scale(self, attr, v):
            if self._syncing:
                return
            self.style = self.style.copy(**{attr: round(float(v), 3)})
            self._render()

        def _on_vmin(self, s):
            if self._syncing or not self._rng:
                return
            lo, hi = self._rng
            val = lo + (hi - lo) * s / 1000.0
            top = self.style.vmax if self.style.vmax is not None else hi
            val = min(val, top - abs(top) * 1e-4 - 1e-12)
            self.style = self.style.copy(vmin=val)
            self._render()

        def _on_vmax(self, s):
            if self._syncing or not self._rng:
                return
            lo, hi = self._rng
            val = lo + (hi - lo) * s / 1000.0
            bot = self.style.vmin if self.style.vmin is not None else lo
            val = max(val, bot + abs(bot) * 1e-4 + 1e-12)
            self.style = self.style.copy(vmax=val)
            self._render()

        def _reset_range(self):
            self.style = self.style.copy(vmin=None, vmax=None)
            self._render()

        # -- entry search / navigation -------------------------------------------------
        def _add_result(self, sid, formula, natoms=None):
            txt = f"{sid}    {formula or ''}"
            if natoms is not None:
                txt += f"   ({natoms})"
            it = QtWidgets.QListWidgetItem(txt)
            it.setData(Qt.UserRole, sid)
            self.results.addItem(it)
            return it

        def _hint_item(self, text):
            """A muted, non-selectable placeholder row (empty browse set / no match)."""
            it = QtWidgets.QListWidgetItem(text)
            it.setFlags(Qt.NoItemFlags)
            it.setForeground(QtGui.QColor("#9aa6b2"))
            it.setTextAlignment(Qt.AlignCenter)
            self.results.addItem(it)

        def _on_search(self, text):
            text = text.strip()
            self.results.clear()
            if not text:                                    # no query -> show the current browse set
                if len(self.sids) > 1:
                    sub = self.sids[:500]
                    ph = ",".join(["?"] * len(sub))
                    m = {r[0]: (r[1], r[2]) for r in self.con.execute(
                        f"SELECT structure_id, formula, natoms FROM structure "
                        f"WHERE structure_id IN ({ph})", sub).fetchall()}
                    for sid in sub:
                        f, n = m.get(sid, (None, None))
                        it = self._add_result(sid, f, n)
                        if sid == self.sid:
                            self.results.setCurrentItem(it)
                else:                                       # single structure -> nothing to browse
                    self._hint_item("type above to search entries")
                return
            like = f"%{text}%"
            rows = self.con.execute(
                "SELECT structure_id, formula, natoms FROM structure "
                "WHERE structure_id ILIKE ? OR formula ILIKE ? ORDER BY structure_id LIMIT 400",
                [like, like]).fetchall()
            for sid, formula, natoms in rows:
                self._add_result(sid, formula, natoms)
            if not rows:
                self._hint_item("no match")

        def _on_pick(self, item):
            sid = item.data(Qt.UserRole)
            if sid:
                self._goto_sid(sid)

        def _goto_sid(self, sid):
            if sid in self.sids:
                self._goto(self.sids.index(sid))
            elif sid != self.sid:
                self.sids = [sid]
                self.idx = 0
                self._refit = True
                self._render()
                self._update_nav()

        def _step(self, d):
            self._goto(self.idx + d)

        def _goto(self, i):
            i = max(0, min(len(self.sids) - 1, int(i)))
            if i != self.idx:
                self.idx = i
                self._refit = True
                self._render()
                self._update_nav()

        def _update_nav(self):
            n = len(self.sids)
            self.pos_lbl.setText(f"{self.idx + 1} / {n}" if n > 1 else "single")
            self.prev_btn.setEnabled(self.idx > 0)
            self.next_btn.setEnabled(self.idx < n - 1)
            if not self.search.text().strip():          # keep the browse-list highlight in step
                for i in range(self.results.count()):
                    if self.results.item(i).data(Qt.UserRole) == self.sid:
                        self.results.setCurrentRow(i)
                        break

        # -- render --------------------------------------------------------------------
        def _render(self):
            pl = self.plotter
            cam = None if (self._first or self._refit) else pl.camera_position
            self._refit = False
            pl.clear()
            try:                               # clear() drops the lights -> restore the 5-light kit
                pl.enable_lightkit()           # (glossy shading/specular, same as emit_pyvista/Jmol)
            except Exception:
                pass
            # Reconcile the aggregated fragment channel to the CURRENT level BEFORE building, so the
            # render can never show a channel invalid for the level (DI only at DFT, the IQA energies
            # only at an SQM method).  Done here (not just when repopulating the combo) so the very
            # first frame after toggling fragment / changing level is already consistent -- no lag
            # where the combo says DI but the canvas still shows Einter.
            if self.fragment:
                opts = self._frag_prop_options()
                if opts and self._frag_prop not in opts:
                    self._frag_prop = opts[0]
            scene = build_scene(self.con, self.sid, prop=self.prop, method=self.method,
                                graph=self.graph, paths=self.paths, ias=self.ias,
                                fragment=self.fragment, fragment_prop=self._frag_prop,
                                style=self.style, side=self._side)
            self._last_range = scene.meta.get("data_range")
            # fragment coloring only actually engages when the covalent graph splits (>=2 pieces);
            # when it does, the aggregated interface channel drives coloring instead of `prop`.
            self._frag_active = bool(self.fragment and scene.meta.get("nfrag", 1) >= 2)
            hdr = scene.title or self.sid
            scene.title = ""                       # the panel shows the title; keep the canvas clean
            add_scene(pl, scene, self.style, points_only=self.points_only, scalar_bar="top")
            if cam is not None:
                pl.camera_position = cam
            else:
                pl.reset_camera()
                pl.camera.azimuth = 20
                pl.camera.elevation = 15
            self._first = False
            pl.render()
            self._update_header(hdr, scene)
            self._sync_range_controls()
            self._sync_frag_ui()

        def _sync_frag_ui(self):
            """Keep the property combo in step with fragment mode.  When host<->guest coloring engages,
            the combo lists the aggregatable interface channels and *drives* which one is summed (so it
            stays enabled -- no more misleading 'shows Einter but names q_QTAIM'); when it disengages,
            the combo goes back to the normal projectable columns and the prior selection is restored.
            Repopulate only on an actual mode flip to avoid clobbering the user's in-mode choice."""
            active = getattr(self, "_frag_active", False)
            if active != self._combo_frag_mode:
                self._refresh_prop_combo()
            self.prop_combo.setEnabled(True)
            self.info_btn.setEnabled(True)

        def _update_header(self, hdr, scene):
            self.title_lbl.setText(self.sid)
            meta = scene.meta if isinstance(scene.meta, dict) else {}
            bits = []
            if getattr(self, "_frag_active", False):       # host<->guest coloring OWNS the channel
                col = meta.get("hg_col", "Einter_SQM")
                short = col[:-4] if col.endswith("_SQM") else col[:-6] if col.endswith("_QTAIM") else col
                unit = "kcal/mol" if meta.get("hg_level") else "e"
                lvl = f"  [{self.method}]" if meta.get("hg_level") else ""
                lab = f"host↔guest Σ{short}{lvl}"
                tot = meta.get("hg_total")
                if tot is not None:
                    lab += f"   Σ={tot:+.1f} {unit}"
                bits.append(lab)
            elif self.prop == FRAGMENTS_PROP:              # structural, level-independent
                bits.append("fragments")
            elif self.prop:
                lvl = (f"  [{self.method}]" if self._side == "sqm"
                       else "  [PBE0]" if self._side == "iqadft" else "")
                bits.append(f"{self.prop}{lvl}")
            else:
                bits.append("geometry (atoms + bonds)")
            bits.append(f"{scene.natoms} atoms")
            if getattr(self, "_frag_active", False):       # name the two interacting fragments
                hf, gf = meta.get("hg_host_formula"), meta.get("hg_guest_formula")
                hn, gn = meta.get("hg_host_n"), meta.get("hg_guest_n")
                if hf and gf:
                    bits.append(f"host {hf} ({hn})  ·  guest {gf} ({gn})")
            elif self.prop == FRAGMENTS_PROP:              # list each colored component
                forms, sizes = meta.get("frag_formulae"), meta.get("frag_sizes")
                if forms and sizes:
                    bits.append("  ·  ".join(f"{f} ({n})" for f, n in zip(forms, sizes)))
            elif self.fragment:
                bits.append("fragment: single component")  # toggle on, but molecule does not split
            nodata = meta.get("nbond_nodata", 0)
            if nodata:                                     # grey sticks = no value for this column
                bits.append(f"{nodata} bonds: no data")
            a_nodata = meta.get("natom_nodata", 0)
            if a_nodata:                                   # grey atoms = NULL scalar (e.g. d_pop on H)
                bits.append(f"{a_nodata} atoms: no data")
            self.sub_lbl.setText("   ·   ".join(bits))

        def _sync_range_controls(self):
            rng = self._last_range
            active = bool(rng) and rng[1] > rng[0]
            self._syncing = True
            for w in (self.vmin_slider, self.vmax_slider, self.auto_btn):
                w.setEnabled(active)
            # CP-size / path-size sliders are meaningful only while the QCT graph / paths are drawn:
            # CP dots exist with the full graph; bond paths exist with graph OR paths.
            self.cp_slider.setEnabled(self.graph)
            self.path_slider.setEnabled(self.graph or self.paths)
            self.atom_val.setText(f"{self.style.atom_scale:.2f}")
            self.bond_val.setText(f"{self.style.bond_scale:.2f}")
            self.cp_val.setText(f"{self.style.cp_scale:.2f}")
            self.path_val.setText(f"{self.style.path_scale:.2f}")
            if active:
                lo, hi = rng
                self._rng = (lo, hi)
                vmin = lo if self.style.vmin is None else self.style.vmin
                vmax = hi if self.style.vmax is None else self.style.vmax
                span = (hi - lo) or 1.0
                self.vmin_slider.setValue(int(1000 * (min(max(vmin, lo), hi) - lo) / span))
                self.vmax_slider.setValue(int(1000 * (min(max(vmax, lo), hi) - lo) / span))
                self.vmin_val.setText(f"{vmin:.3g}")
                self.vmax_val.setText(f"{vmax:.3g}")
            else:
                self._rng = None
                self.vmin_val.setText("--")
                self.vmax_val.setText("--")
            self._syncing = False

        # -- lifecycle -----------------------------------------------------------------
        def closeEvent(self, ev):
            try:
                self.plotter.close()
            except Exception:
                pass
            super().closeEvent(ev)


def open_qt_viewer(con, sid=None, sids=None, method="PM7", style=None, points_only=False,
                   save=None, fragment_prop=None, block=True, pbe0=False):
    """Create the QApplication (if needed), show the QtViewer, and (optionally) run the event loop.

    Returns the window (for headless/testing when ``block=False``) or `save` after the loop exits.
    """
    if not QT_AVAILABLE:
        raise RuntimeError("Qt viewer requested but PySide6/pyvistaqt are not importable")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    win = QtViewer(con, sid=sid, sids=sids, method=method, style=style,
                   points_only=points_only, save=save, fragment_prop=fragment_prop, pbe0=pbe0)
    win.show()
    print(WINDOW_TITLE + " (Qt) -- pick a property/level from the panel; type an id or formula to "
          "jump to another entry.")
    if not block:
        return win
    app.exec()
    return save
