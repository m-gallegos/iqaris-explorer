"""iqaris_export -- explore, visualize and export the IQARIS database for ML training.

The IQARIS source calculations are single-point, so there are NO atomic forces: exports carry
per-structure total energies (the PBE0 SCF energy, the PBE0 IQA-recovered energy, plus five SQM
totals) as the training target plus the rich per-atom / per-pair / per-CP quantum labels
(QTAIM density topology at M06-2X and PBE0 + ab initio IQA at PBE0 + IQA-SQM) as ML labels.

Quick start (Python)::

    import iqaris_export as ie
    con = ie.connect()
    sel = ie.Selection(elements=("S",), natoms=(1, 40), one_per_molecule=True)
    print(ie.describe(con, sel))
    ie.export(con, sel, props="atomic,charges", level="PM7",
              fmt="extxyz", out="exports/demo")
"""
from __future__ import annotations

from . import config
from .db import connect, filter_ids, molecule_groups
from .selection import Selection, Predicate
from .style import Style
from .properties import resolve as resolve_props, Plan
from .registry import get_registry, ColumnSpec
from .engine import Exporter
from .explore import describe, coverage, column_stats, value_range, list_props, property_info
from .writers import get_writer, TABULAR_FORMATS, ALL_FORMATS, export_tabular
from . import property_docs

__version__ = "0.1.0"

__all__ = [
    "connect", "Selection", "Predicate", "export", "describe", "coverage", "column_stats",
    "value_range", "list_props", "resolve_props", "Plan", "get_registry",
    "config", "ALL_FORMATS", "project", "build_scene", "emit_jmol", "Scene",
    "plot_property_map", "view", "emit_pyvista", "Style",
    "info", "property_info", "glossary", "property_docs",
]


def info(column, con=None, printit=True):
    """Print (and return) a human-readable description of a database property/column.

    Pulls the exact per-column unit from the live registry when ``con`` is given, otherwise
    answers from the shipped property dictionary alone (no database needed)::

        import iqaris_export as ie
        ie.info("q_QTAIM")                 # dictionary only
        ie.info("q_QTAIM", con=ie.connect())   # with the exact DB unit
    """
    from .registry import _resolve_unit
    # The exact per-column unit is a function of the column name (no DB needed), so resolve it even
    # on the --no-db path rather than showing the composite ";" family blob from doc.unit.
    unit = _resolve_unit(column)
    level_aware = None
    extra = ""
    if con is not None:
        specs = get_registry(con)["by_name"].get(column)
        if specs:
            level_aware = any(s.level_aware for s in specs)
            # list EVERY view/level the column lives in, not just specs[0]: a *_QTAIM
            # column appears in both qtaim_m062x_* and qtaim_pbe0_*.
            views = sorted({s.view for s in specs})
            if views:
                extra = "\nViews  : " + ", ".join(views)
    text = property_docs.format_text(column, unit=unit, level_aware=level_aware) + extra
    if printit:
        print(text)
    return text


def glossary(fmt="markdown"):
    """Return the whole property dictionary rendered as a standalone document (markdown|html)."""
    return property_docs.glossary(fmt)


def export(con, selection, props, level="dft", fmt="extxyz", out=None,
           energy_target=config.DEFAULT_ENERGY_TARGET, with_bonus=False,
           batch_size=200, force_stub=False, progress=True, include_problematic=False):
    """Resolve a property plan and stream the selection to `out` in format `fmt`.

    `include_problematic` switches the PBE0 record to its opt-in `_all` views (clean UNION the
    quality-flagged records) and emits a `pbe0_qc_status` column; default is clean data only.
    """
    out = str(out or (config.DEFAULT_OUT / "export"))
    plan = resolve_props(con, props, level=level, with_bonus=with_bonus)
    reg = get_registry(con)
    if fmt in TABULAR_FORMATS:
        return export_tabular(con, selection, plan, reg, out, fmt=fmt,
                              energy_target=energy_target, level=level,
                              include_problematic=include_problematic)
    writer_cls = get_writer(fmt)
    writer = writer_cls(out, plan, reg, selection, energy_target, level=level,
                        args={"force_stub": force_stub})
    exporter = Exporter(con, selection, plan, energy_target=energy_target,
                        batch_size=batch_size, progress=progress,
                        include_problematic=include_problematic)
    return exporter.run(writer)


def plot_distribution(*a, **k):
    from .viz import plot_distribution as _p
    return _p(*a, **k)


def plot_coverage(*a, **k):
    from .viz import plot_coverage as _p
    return _p(*a, **k)


def project(*a, **k):
    from .project import project as _p
    return _p(*a, **k)


def build_scene(*a, **k):
    from .project import build_scene as _p
    return _p(*a, **k)


def emit_jmol(*a, **k):
    from .project import emit_jmol as _p
    return _p(*a, **k)


def plot_property_map(*a, **k):
    from .viz import plot_property_map as _p
    return _p(*a, **k)


def view(*a, **k):
    """Open an interactive PyVista/VTK 3D window for a structure (or save a headless PNG)."""
    from .viz3d import view as _v
    return _v(*a, **k)


def emit_pyvista(*a, **k):
    from .viz3d import emit_pyvista as _e
    return _e(*a, **k)


def __getattr__(name):
    if name == "Scene":
        from .project import Scene
        return Scene
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
