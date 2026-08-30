"""Tabular fast path: one file per family, written by streaming DuckDB `COPY` (no per-structure
assembly, no OOM even for the 229M-row pair table).  The analyst / DuckDB-reloadable export.
"""
from __future__ import annotations
from pathlib import Path

from .. import config
from ..registry import _is_iqadft_col
from .base import BaseWriter


def export_tabular(con, selection, plan, reg, outdir, fmt="parquet",
                   energy_target=config.DEFAULT_ENERGY_TARGET, level="dft",
                   include_problematic=False):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    idsql, params = selection.ids_query(con)
    con.execute("CREATE OR REPLACE TEMP TABLE _sel AS " + idsql, params)
    total = con.execute("SELECT count(*) FROM _sel").fetchone()[0]
    opt = "(FORMAT parquet)" if fmt == "parquet" else "(FORMAT csv, HEADER)"
    ext = "parquet" if fmt == "parquet" else "csv"
    sel_filter = "structure_id IN (SELECT structure_id FROM _sel)"

    def _write(name, q):
        path = outdir / f"{name}.{ext}"
        con.execute(f"COPY ({q}) TO '{path}' {opt}")
        n = con.execute(f"SELECT count(*) FROM read_{ext if ext=='parquet' else 'csv'}('{path}')").fetchone()[0]
        print(f"  [tabular] {name:22s} {n:>12,} rows -> {path.name}")

    def dump(name, view, cols, methods=None):
        keycols = config.family_keys(view)
        sel = keycols + [c for c in cols if c not in keycols]
        where = sel_filter
        if methods:
            mlist = ",".join(f"'{m}'" for m in methods)
            where += f" AND method IN ({mlist})"
        _write(name, f"SELECT {', '.join(sel)} FROM {view} WHERE {where} "
                     f"ORDER BY {', '.join(keycols)}")

    def dump_pbe0(fam, icols):
        """One file per PBE0 family, JOINing qtaim_pbe0_* (topology) with iqa_pbe0_* (ab initio
        IQA) on the shared keys (cp/angles are topology-only, no join)."""
        pair = config.FAMILY_PBE0_VIEWS.get(fam)
        if not pair:
            return
        qv, iv = pair
        if include_problematic:
            qv = qv + config.PBE0_ALL_SUFFIX
            iv = (iv + config.PBE0_ALL_SUFFIX) if iv else None
        keys = config.family_keys(qv)
        top = [c for c in icols if not _is_iqadft_col(c)]
        iqa = [c for c in icols if _is_iqadft_col(c)]
        # pbe0_qc_status lives in BOTH _all views (each joins the manifest), so qualify it to q
        qc = ["q." + config.PBE0_QC_STATUS_COL] if include_problematic else []
        if iv and iqa:
            sel = keys + top + iqa + qc       # topology/IQA names disjoint; keys deduped by USING
            _write(f"{fam}_pbe0", f"SELECT {', '.join(sel)} FROM {qv} q "
                                  f"JOIN {iv} i USING ({', '.join(keys)}) "
                                  f"WHERE {sel_filter} ORDER BY {', '.join(keys)}")
        else:                                  # cp/angles or IQA-only: single view
            dump(f"{fam}_pbe0", qv, top + ([config.PBE0_QC_STATUS_COL] if include_problematic else []) + iqa)

    # always: a manifest slice + geometry + the energy target
    con.execute(f"COPY (SELECT * FROM structure WHERE {sel_filter} ORDER BY structure_id) "
                f"TO '{outdir / ('manifest.' + ext)}' {opt}")
    con.execute(f"COPY (SELECT g.* FROM geometry g WHERE {sel_filter} "
                f"ORDER BY structure_id, atom_index) TO '{outdir / ('geometry.' + ext)}' {opt}")
    if energy_target:
        con.execute(f"COPY (SELECT structure_id, {energy_target} FROM energies WHERE {sel_filter} "
                    f"ORDER BY structure_id) TO '{outdir / ('energy.' + ext)}' {opt}")

    for fam in plan.families:
        qview, mview = config.FAMILY_VIEWS[fam]
        if plan.qtaim_cols(fam):
            dump(f"{fam}_qtaim", qview, plan.qtaim_cols(fam))
        if plan.sqm_cols(fam) and plan.methods:
            dump(f"{fam}_sqm", mview, plan.sqm_cols(fam), methods=plan.methods)
        if config.FAMILY_PBE0_VIEWS.get(fam) and plan.iqadft_cols(fam):
            dump_pbe0(fam, plan.iqadft_cols(fam))

    con.execute("DROP TABLE IF EXISTS _sel")
    # sidecars
    w = BaseWriter(outdir, plan, reg, selection, energy_target, level)
    w.open(total=total)
    w.n = total
    w.close()
    return {"structures": total, "format": fmt}
