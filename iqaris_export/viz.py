"""Lightweight visualization: property distributions and coverage, in the IQARIS house style
(`fig_style.py`).  Heavy columns are drawn from a seeded sample so plotting stays interactive.
The per-structure 3D property-map renderer lives in `project.py` (backend-neutral Scene +
Jmol emitter); `plot_property_map` builds a Scene for the requested structure and renders it.
"""
from __future__ import annotations
import os
import sys

from . import config
from .explore import _filter, coverage, _pick_spec
from .registry import get_registry

def _fs():
    from . import fig_style
    fig_style.set_style()
    return fig_style


def _use_interactive_backend():
    """Switch matplotlib to a real GUI backend for `--show`, undoing fig_style's forced Agg
    (Agg is headless-only; `plt.show()` is a silent no-op under it)."""
    import matplotlib
    matplotlib.use("QtAgg", force=True)


def plot_distribution(con, column, selection=None, kind="hist", by=None,
                      method=None, out=None, sample=50000, seed=config.DEFAULT_SEED, bins=60,
                      show=False, side=None):
    """Histogram or per-element violin of any registry column across a selection.

    By default renders headlessly and saves a PDF under `out`. With ``show=True``, opens an
    interactive window instead; pass `out` as well to both display AND save.  `side` disambiguates
    a name registered on both the M06-2X ('qtaim') and PBE0 ('iqadft') sides.
    """
    import matplotlib.pyplot as plt
    fs = _fs()
    if show:
        _use_interactive_backend()
    reg = get_registry(con)
    specs = reg["by_name"].get(column)
    if not specs:
        raise ValueError(f"unknown column {column!r}")
    spec = _pick_spec(specs, column, side)
    where, p = _filter(con, selection)
    extra = ""
    if spec.level_aware:
        method = method or config.SQM_METHODS[0]
        extra = f" AND method = '{method}'"
    has_elem = by == "element"
    cols = f"element, {column}" if has_elem else column
    inner = f"SELECT {cols} FROM {spec.view} WHERE {where}{extra} AND {column} IS NOT NULL"
    q = f"SELECT {cols} FROM ({inner}) USING SAMPLE {int(sample)} ROWS (reservoir, {int(seed)})"
    df = con.execute(q, p).fetchdf()

    explicit_out = out is not None
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    title = f"{column}" + (f" [{method}]" if spec.level_aware else "")
    if has_elem and kind == "violin":
        elems = [e for e in fs.ELEM_ORDER if e in set(df["element"])]
        data = [df.loc[df["element"] == e, column].to_numpy() for e in elems]
        parts = ax.violinplot(data, showmedians=True, showextrema=False)
        fs.style_violin(parts, [fs.ECOL[e] for e in elems])
        ax.set_xticks(range(1, len(elems) + 1))
        ax.set_xticklabels(elems)
        ax.set_ylabel(f"{column} ({spec.unit})")
    else:
        ax.hist(df[column].to_numpy(), bins=bins, color=fs.OI["blue"], alpha=0.85)
        ax.set_xlabel(f"{column} ({spec.unit})")
        ax.set_ylabel("count")
    ax.set_title(title)
    if show:
        plt.show()          # blocks until closed; fs.save() below would otherwise close the
                             # figure first via plt.close(), leaving nothing to show
    path = None
    if explicit_out or not show:
        out = out or str(config.DEFAULT_OUT / "plots")
        os.makedirs(out, exist_ok=True)
        path = os.path.join(out, f"{column}_{'violin' if (has_elem and kind=='violin') else 'hist'}.pdf")
        fs.save(fig, path)
    return path


def plot_coverage(con, selection=None, out=None, show=False):
    """Elemental-coverage bar chart. With ``show=True``, opens an interactive window instead of
    (or, with `out` also given, in addition to) saving a PDF."""
    import matplotlib.pyplot as plt
    fs = _fs()
    if show:
        _use_interactive_backend()
    cov = coverage(con, selection)
    elems = [d["element"] for d in cov["elements"]]
    natoms = [d["n_atoms"] for d in cov["elements"]]
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    ax.bar(range(len(elems)), natoms,
           color=[fs.ECOL.get(e, fs.OI["grey"]) for e in elems])
    ax.set_xticks(range(len(elems)))
    ax.set_xticklabels(elems)
    ax.set_ylabel("atoms in selection")
    ax.set_title("Elemental coverage")
    if show:
        plt.show()
    path = None
    if out is not None or not show:
        out = out or str(config.DEFAULT_OUT / "plots")
        os.makedirs(out, exist_ok=True)
        path = os.path.join(out, "coverage_elements.pdf")
        fs.save(fig, path)
    return path


def plot_property_map(structure_id, prop=None, method="PM7", graph=False, paths=False,
                      ias=False, fragment=False, out=None, con=None, root=None):
    """Project a per-atom property onto a structure's 3D geometry (+ optional QTAIM molecular
    graph / interatomic surfaces) and render a static Jmol PNG.

    This builds a backend-neutral Scene for THE REQUESTED structure via `project.build_scene`
    and renders it -- it no longer shells out to the A11-specific script, which fixed a bug
    where the structure id was silently ignored and every call re-rendered A11.  Returns the
    PNG path.
    """
    from .project import project
    return project(structure_id, prop=prop, method=method, graph=graph, paths=paths,
                   ias=ias, fragment=fragment, out=out, con=con, root=root)
