"""Database connection + DuckDB-native selection helpers.

The released database ships no Python next to the data, so the access layer is implemented here
directly on the DuckDB catalog. Read-only throughout; nothing here ever writes into the database.
"""
from __future__ import annotations
import os
import re

from . import config

_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")   # guard column names spliced into SQL


def connect(root: str = None, read_only: bool = True):
    """Open the pre-baked DuckDB catalog (all views ready); fall back to executing `iqaris_views.sql`
    if the catalog file is absent (e.g. a freshly untarred release before `make_catalog.py`)."""
    root = root or config.DB_ROOT
    import duckdb
    cat = os.path.join(root, config.CATALOG)
    if os.path.exists(cat):
        return duckdb.connect(cat, read_only=read_only)
    con = duckdb.connect()
    with open(os.path.join(root, config.VIEWS_SQL)) as fh:
        con.execute(fh.read())
    return con


def filter_ids(root: str = None, *, method: str = "QTAIM", sane: bool = True,
               natoms=None, con=None, **resid_max):
    """Structure ids passing a granular sanity filter (DuckDB-native), with an optional `natoms`
    range.

    method  : the `sanity.method` scope -- QTAIM | PBE0 | PM7 | PM6 | PM6-D3H4 | PM6-D3H4X | PM6-ORG.
    sane    : require `sanity.sane = TRUE`.
    natoms  : optional (lo, hi) inclusive atom-count range (from `manifest`).
    resid_max : per-column upper bounds, e.g. ``iqa_closure_resid_kj=1.0`` -> keep rows with the
                column < 1.0 (NULL is treated as +inf, i.e. excluded), matching the old behaviour.

    Example: ``filter_ids(method="PBE0", iqa_closure_resid_kj=1.0)`` / ``filter_ids(natoms=(3, 5))``.
    """
    own = con is None
    if own:
        con = connect(root)
    try:
        where, params = ["s.method = ?"], [method]
        if sane:
            where.append("s.sane")
        for col, hi in resid_max.items():
            if not _IDENT.match(col):
                raise ValueError(f"bad residual column name {col!r}")
            where.append(f"s.{col} < ?")            # NULL -> NULL -> excluded (== fillna(inf) < hi)
            params.append(hi)
        sql = f"SELECT s.structure_id FROM sanity s WHERE {' AND '.join(where)}"
        if natoms is not None:
            sql += (" AND s.structure_id IN (SELECT structure_id FROM manifest "
                    "WHERE natoms BETWEEN ? AND ?)")
            params += [int(natoms[0]), int(natoms[1])]
        return {r[0] for r in con.execute(sql, params).fetchall()}
    finally:
        if own:
            con.close()


def molecule_groups(root: str = None, key: str = "inchikey_block1", limit: int = None, con=None):
    """Map a molecule key -> list of its structure_ids (DuckDB-native; replaces the retired toolkit
    function). `key` in {inchikey_block1, inchikey, formula}; `limit` caps the number of groups."""
    if key not in ("inchikey_block1", "inchikey", "formula"):
        raise ValueError("key must be inchikey_block1 | inchikey | formula")
    own = con is None
    if own:
        con = connect(root)
    try:
        sql = f"SELECT {key} AS k, list(structure_id) AS ids FROM manifest GROUP BY {key}"
        if limit:
            sql += f" LIMIT {int(limit)}"
        return {k: ids for k, ids in con.execute(sql).fetchall()}
    finally:
        if own:
            con.close()


def recommended_ids(root: str = None, con=None):
    """Structure ids in the `recommended` tier (DuckDB-native)."""
    own = con is None
    if own:
        con = connect(root)
    try:
        return {r[0] for r in
                con.execute("SELECT structure_id FROM manifest WHERE recommended").fetchall()}
    finally:
        if own:
            con.close()


def run_df(con, sql: str, params=None):
    """Execute and return a pandas DataFrame."""
    return con.execute(sql, params or []).fetchdf()
