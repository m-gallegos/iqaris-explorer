"""Exploration helpers: summaries, coverage and per-column statistics computed IN SQL (streamed,
never materialized).  Returned as plain dicts / DataFrames; the CLI pretty-prints them.
"""
from __future__ import annotations

from . import config
from .registry import get_registry


def _filter(con, selection):
    """Return (where_clause, params) restricting to a Selection (or 'TRUE' for the whole DB)."""
    if selection is None:
        return "TRUE", []
    idsql, params = selection.ids_query(con)
    return f"structure_id IN ({idsql})", params


def describe(con, selection=None):
    where, p = _filter(con, selection)
    row = con.execute(
        f"SELECT count(*) n, count(DISTINCT inchikey_block1) nmol, "
        f"min(natoms) amin, median(natoms) amed, max(natoms) amax, "
        f"sum(CASE WHEN sane_qtaim THEN 1 ELSE 0 END) nq, "
        f"sum(CASE WHEN sane_mopac_all5 THEN 1 ELSE 0 END) nm, "
        f"sum(CASE WHEN recommended THEN 1 ELSE 0 END) nr "
        f"FROM manifest WHERE {where}", p).fetchone()
    fam = con.execute(
        f"SELECT source_family, count(*) n FROM manifest WHERE {where} "
        f"GROUP BY 1 ORDER BY n DESC", p).fetchall()
    return {
        "n_structures": row[0], "n_molecules": row[1],
        "natoms_min": row[2], "natoms_median": float(row[3]) if row[3] is not None else None,
        "natoms_max": row[4],
        "sane_qtaim": row[5], "sane_mopac_all5": row[6], "recommended": row[7],
        "by_source_family": dict(fam),
    }


def coverage(con, selection=None):
    where, p = _filter(con, selection)
    elems = con.execute(
        f"SELECT element, count(*) n_atoms, count(DISTINCT structure_id) n_struct "
        f"FROM geometry WHERE {where} GROUP BY 1 ORDER BY n_atoms DESC", p).fetchall()
    conf = con.execute(
        f"SELECT n, count(*) n_mol FROM ("
        f"  SELECT inchikey_block1, count(*) n FROM manifest WHERE {where} GROUP BY 1"
        f") GROUP BY 1 ORDER BY n", p).fetchall()
    return {
        "elements": [{"element": e, "n_atoms": a, "n_struct": s} for e, a, s in elems],
        "conformers_per_molecule": dict(conf),
    }


def analysis_of(spec):
    """The kind of analysis a column reports: QTAIM density topology, ab initio IQA, or IQA-SQM.
    Derived from the column's view, so the PBE0 record's topology columns read 'QTAIM' and its
    energy-decomposition columns read 'IQA'."""
    if spec.side == "sqm":
        return "IQA-SQM"
    if spec.view.startswith("iqa_pbe0"):
        return "IQA"                       # ab initio IQA energy decomposition (PBE0 density)
    return "QTAIM"                          # density topology / QCT (M06-2X or PBE0 density)


def level_of(spec, method=None):
    """Human level-of-theory label for a column spec: the SQM method for level-aware columns,
    'PBE0' for the PBE0 record (side 'iqadft'), else 'M06-2X'."""
    if spec.level_aware:
        return method or config.SQM_METHODS[0]
    if spec.side == "iqadft":
        return "PBE0"
    return "M06-2X"


def _pick_spec(specs, column, side):
    """Choose the ColumnSpec for the requested side.  Many QTAIM/QCT names are
    registered on BOTH the 'qtaim' (M06-2X) and 'iqadft' (PBE0) sides; without a side the first
    (M06-2X) is used, preserving prior behavior.  A requested side that lacks the column errors."""
    if side:
        for s in specs:
            if s.side == side:
                return s
        avail = sorted({s.side for s in specs})
        raise ValueError(f"column {column!r} is not on side {side!r} (available: {avail})")
    return specs[0]


def column_stats(con, column, selection=None, method=None, side=None):
    reg = get_registry(con)
    specs = reg["by_name"].get(column)
    if not specs:
        raise ValueError(f"unknown column {column!r}")
    spec = _pick_spec(specs, column, side)
    where, p = _filter(con, selection)
    extra = ""
    if spec.level_aware:
        if method is None:
            method = config.SQM_METHODS[0]
        extra = f" AND method = '{method}'"
    q = (f"SELECT count({column}) n, min({column}) lo, max({column}) hi, "
         f"avg({column}) mean, stddev({column}) std, "
         f"approx_quantile({column}, 0.5) med "
         f"FROM {spec.view} WHERE {where}{extra}")
    r = con.execute(q, p).fetchone()
    return {"column": column, "view": spec.view, "unit": spec.unit,
            "analysis": analysis_of(spec), "level": level_of(spec, method),
            "level_aware": spec.level_aware, "method": method if spec.level_aware else None,
            "n": r[0], "min": r[1], "max": r[2], "mean": r[3], "std": r[4], "median": r[5]}


def value_range(column, level=None, con=None, selection=None, method=None):
    """(min, max) of a column at a given LEVEL of theory.

    `level` = m062x | pbe0 | one of the SQM methods; opens its own connection if `con` is None,
    so viewer colour-scale ranges match the level being shown (not always M06-2X)."""
    own = con is None
    if own:
        from .db import connect
        con = connect()
    try:
        side = None
        if level:
            lv = level.lower()
            side = _LEVEL_TO_SIDE.get(lv)
            if side is None and any(lv == m.lower() for m in config.SQM_METHODS):
                side, method = "sqm", next(m for m in config.SQM_METHODS if m.lower() == lv)
        s = column_stats(con, column, selection, method, side=side)
        return (s["min"], s["max"])
    finally:
        if own:
            con.close()


# Two user-facing catalogue axes: LEVEL (m062x | pbe0 | the five SQM Hamiltonians) and ANALYSIS
# (qtaim | iqa | iqa-sqm). A *_QTAIM name is listed at both M06-2X and PBE0, a *_SQM name once per
# Hamiltonian -- one row per property per real level.
_LEVEL_TO_SIDE = {"m062x": "qtaim", "dft": "qtaim", "pbe0": "iqadft"}
_ANALYSIS_CANON = {"qtaim": "QTAIM", "iqa": "IQA", "iqa-sqm": "IQA-SQM", "iqasqm": "IQA-SQM"}
_LEVEL_ORDER = {lv: i for i, lv in enumerate(["M06-2X", "PBE0"] + list(config.SQM_METHODS))}


def list_props(con, family=None, level=None, analysis=None):
    """List catalogue columns, filtered by the two public axes (never by an internal token).

    level    : 'm062x' | 'pbe0' | one of the five SQM methods (pm7 | pm6 | pm6-d3h4 | pm6-d3h4x |
               pm6-org), or 'all'/None for every level.  ('sqm' is NOT a level -- use analysis.)
    analysis : 'qtaim' | 'iqa' | 'iqa-sqm' (or None for every analysis type).
    """
    reg = get_registry(con)
    want_side = None
    want_method = None                     # a single SQM level to restrict the expansion to
    if level and level.lower() != "all":
        lv = level.lower()
        if lv in _LEVEL_TO_SIDE:
            want_side = _LEVEL_TO_SIDE[lv]
        else:
            m = next((mm for mm in config.SQM_METHODS if mm.lower() == lv), None)
            if m is None:
                raise ValueError(f"unknown level {level!r} (use m062x | pbe0 | "
                                 f"{' | '.join(mm.lower() for mm in config.SQM_METHODS)} | all)")
            want_side, want_method = "sqm", m
    want_analysis = None
    if analysis:
        want_analysis = _ANALYSIS_CANON.get(analysis.lower())
        if want_analysis is None:
            raise ValueError(f"unknown analysis {analysis!r} (use qtaim | iqa | iqa-sqm)")

    def _row(s, a, lvl):
        return {"name": s.name, "family": s.family, "analysis": a, "level": lvl,
                "unit": s.unit, "dtype": s.dtype, "view": s.view, "description": s.description}

    out = []
    for (fam, sd), specs in reg["by_group"].items():
        if family and fam != family:
            continue
        if want_side and sd != want_side:
            continue
        for s in specs:
            a = analysis_of(s)
            if want_analysis and a != want_analysis:
                continue
            if s.level_aware:              # an SQM property -> one row per Hamiltonian (the levels)
                for m in ([want_method] if want_method else list(config.SQM_METHODS)):
                    out.append(_row(s, a, m))
            else:                          # a QTAIM/IQA property -> its single concrete level
                out.append(_row(s, a, level_of(s)))
    out.sort(key=lambda d: (d["family"], _LEVEL_ORDER.get(d["level"], 99), d["analysis"], d["name"]))
    return out


def property_info(con, column, method=None):
    """Merge the LIVE registry spec (exact unit / view / level-awareness) with the
    :mod:`property_docs` dictionary entry (name, symbol, prose, equation) for a column.

    Returns a plain dict; ``con`` may be None to answer from the dictionary alone.
    """
    from . import property_docs
    from .registry import _resolve_unit
    specs = []
    if con is not None:
        specs = get_registry(con)["by_name"].get(column) or []
    doc = property_docs.lookup(column)
    # unit is a pure function of the column name (DB-free) -> never the ";" blob
    unit = _resolve_unit(column) if (specs or doc) else (doc.unit if doc else None)
    # report ALL views/families/levels the column lives in, not just specs[0]
    info = {
        "column": column,
        "known": bool(specs) or bool(doc),
        "unit": unit,
        "view": specs[0].view if specs else None,           # kept for back-compat
        "views": sorted({s.view for s in specs}),
        "family": specs[0].family if specs else None,
        "families": sorted({s.family for s in specs}),
        "levels": sorted({level_of(s) for s in specs}, key=lambda x: (len(x), x)),
        "level_aware": any(s.level_aware for s in specs) if specs else None,
    }
    if doc is not None:
        info.update(name=doc.name, symbol=doc.symbol, type=doc.type, method=doc.method,
                    description=doc.text, eq_text=doc.eq_text, eq_latex=doc.eq_latex,
                    bonus=doc.bonus, see_also=list(doc.see_also), key=doc.key)
    return info
