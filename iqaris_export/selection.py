"""Selection: a declarative filter over the IQARIS structures that compiles to one DuckDB
query returning a sorted `structure_id` list.  Filters: elements present / composition subset,
source family, atom-count range, formula, sanity tier, explicit ids, molecule-level dedup,
arbitrary value predicates (sanity residuals / any property column), and a reproducible random
sample.  Order is always by `structure_id` so the streaming engine hits contiguous on-disk row
groups.
"""
from __future__ import annotations
from dataclasses import dataclass, field

from .config import DEFAULT_SEED, SQM_METHODS

ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl"}
SOURCE_FAMILIES = {"aminoacid", "atom", "formula_qm9", "gdb", "mobley", "ncia",
                   "opt", "other", "qm9"}
TIERS = {"sane_qtaim", "sane_mopac_all5", "sane_all", "recommended"}
MOL_KEYS = {"inchikey_block1", "inchikey", "formula"}
MOL_PICKS = {"min_id", "max_natoms", "recommended"}

# --- value-predicate vocabulary ---------------------------------------------------------
# `manifest` columns are per-structure -> filtered directly on the `m` alias.
MANIFEST_COLUMNS = {
    "structure_id", "source_family", "natoms", "charge", "multiplicity", "formula",
    "canonical_smiles", "inchikey", "inchikey_block1", "inchi",
    "sane_qtaim", "sane_mopac_all5", "sane_all", "recommended",
    "identity_tool", "identity_status",
    # PBE0 record flags + quality-flag columns
    "has_pbe0", "sane_pbe0", "has_pbe0_flagged", "pbe0_qc_status",
}
# `sanity` columns are per (structure_id, method) -> filtered via a correlated EXISTS so several
# methods can AND together in one query.
SANITY_COLUMNS = {
    "status_ok", "sane", "charge_resid", "charge_ok", "pair_count_ok", "poincare_hopf",
    "has_nna", "L_atom_max", "integ_ok", "q_recon_max", "N_LI_DI_max",
    "scf_converged", "etot_recon_resid", "etot_ok", "disp_supported", "disp_selfconsistent",
    "halogen_quarantine",
    # PBE0 QC columns (populated on the PBE0 sanity rows; NULL for QTAIM/SQM methods)
    "iqa_closure_resid_kj", "iqa_ok", "charge_tol_e", "iqa_closure_tol_kj", "L_tol", "qc_verdict",
    "qc_fails",                                   # the semicolon-joined failing-criteria string
    # 4th PBE0 QC criterion: pair-sum reconstruction residual |DE_IQA_RECON| (Eh) and its applied
    # threshold (0.1 Eh); populated on the PBE0 sanity rows only
    "pair_closure_resid_eh", "pair_closure_tol_eh",
}
# the 7 methods present in the `sanity` table (QTAIM + PBE0 + five SQM)
SANITY_METHODS = ("QTAIM", "PBE0") + SQM_METHODS
ALLOWED_OPS = {"<", "<=", ">", ">=", "==", "!=", "between", "in", "is_true", "is_false",
               "contains", "not_contains"}   # contains: substring match on a text column (e.g. qc_fails)
_SQL_CMP = {"<": "<", "<=": "<=", ">": ">", ">=": ">=", "==": "=", "!=": "<>"}

# `@LEVEL` scope tokens accepted on a predicate (case-insensitive).  The M06-2X sanity row is stored
# under method 'QTAIM', so the user-facing spellings M06-2X / M062X / DFT normalise to it.
_METHOD_ALIASES = {"M06-2X": "QTAIM", "M062X": "QTAIM", "DFT": "QTAIM"}


def _normalize_method(m):
    """Canonical sanity/SQM method for a `@LEVEL` token, or None if unknown."""
    u = str(m).strip().upper()
    if u in _METHOD_ALIASES:
        return _METHOD_ALIASES[u]
    for s in SANITY_METHODS:
        if s.upper() == u:
            return s
    return None


# The Selection `level` picks which reading an unscoped predicate on a shared QTAIM/QCT column (or a
# sanity row) binds to.  Accept the paper's level spellings and map them to the two internal sides.
_INVALID_LEVEL = object()
_LEVEL_SIDE = {
    "pbe0": "iqadft",
    "m062x": "qtaim", "m06-2x": "qtaim", "dft": "qtaim",
}


def _normalize_level(level):
    """Map a level token ('pbe0' | 'm062x'/'dft', case-insensitive) to its internal side, pass
    None through, and return the sentinel for anything else."""
    if level is None:
        return None
    key = str(level).strip().lower()
    if key in _LEVEL_SIDE:
        return _LEVEL_SIDE[key]
    if key in ("qtaim", "iqadft"):          # internal side names, accepted verbatim
        return key
    return _INVALID_LEVEL


def _check(name, values, allowed):
    bad = set(values) - allowed
    if bad:
        raise ValueError(f"unknown {name}: {sorted(bad)} (allowed: {sorted(allowed)})")


@dataclass(frozen=True)
class Predicate:
    """One value filter over a structure/property/sanity column.

    column : a `manifest` column (per-structure), a `sanity` residual/flag (per method, matched
             via a correlated EXISTS), or any registry *property* column (matched via EXISTS over
             its data view; for the multi-row families -- atomic/pairs/cp/angles -- this means
             "at least one row satisfies", so it is a containment test, not an all-rows test).
    op     : one of < <= > >= == != between in is_true is_false.
    value  : scalar for the comparisons; a (lo, hi) pair for `between`; a tuple for `in`;
             None for is_true/is_false.
    method : level scope (`@LEVEL` in the CLI): QTAIM / M06-2X / M062X (the M06-2X sanity row, or
             the M06-2X reading of a shared *_QTAIM column), PBE0 (the PBE0 sanity row / PBE0
             reading), or an SQM method (PM7, PM6, ...) for level-aware *_SQM columns.
             Case-insensitive; M06-2X/M062X/DFT normalise to QTAIM.  When None, the Selection's
             `level` / `sqm_method` (the subcommand's --level / --method) applies, then the
             defaults M06-2X and PM7.
    """
    column: str
    op: str
    value: object = None
    method: str = None

    def __post_init__(self):
        if self.op not in ALLOWED_OPS:
            raise ValueError(f"unknown predicate op {self.op!r} "
                             f"(allowed: {sorted(ALLOWED_OPS)})")
        if self.op in ("is_true", "is_false"):
            if self.value is not None:
                raise ValueError(f"{self.op} takes no value")
        elif self.op == "between":
            if not (isinstance(self.value, (tuple, list)) and len(self.value) == 2):
                raise ValueError("between needs a (lo, hi) pair")
        elif self.op == "in":
            if not (isinstance(self.value, (tuple, list)) and len(self.value) >= 1):
                raise ValueError("in needs a non-empty list of values")
        else:
            if self.value is None or isinstance(self.value, (tuple, list)):
                raise ValueError(f"{self.op} needs a single scalar value")
        if self.method is not None:
            norm = _normalize_method(self.method)
            if norm is None:
                raise ValueError(f"unknown @LEVEL scope {self.method!r} (allowed: "
                                 f"{sorted(SANITY_METHODS)}; M06-2X/M062X/DFT mean QTAIM)")
            object.__setattr__(self, "method", norm)      # frozen dataclass: store the canonical token


def _cond_sql(lhs, pred):
    """Compile the comparison itself against an already-qualified column expression `lhs`."""
    op = pred.op
    if op == "is_true":
        return f"CAST({lhs} AS BOOLEAN) = TRUE", []      # CAST covers BOOLEAN and 0/1 DOUBLE flags
    if op == "is_false":
        return f"CAST({lhs} AS BOOLEAN) = FALSE", []
    if op == "between":
        return f"{lhs} BETWEEN ? AND ?", [pred.value[0], pred.value[1]]
    if op == "in":
        ph = ",".join("?" for _ in pred.value)
        return f"{lhs} IN ({ph})", list(pred.value)
    if op == "contains":
        return f"CAST({lhs} AS VARCHAR) LIKE ?", ["%" + str(pred.value) + "%"]
    if op == "not_contains":
        return f"CAST({lhs} AS VARCHAR) NOT LIKE ?", ["%" + str(pred.value) + "%"]
    return f"{lhs} {_SQL_CMP[op]} ?", [pred.value]


def _pick_spec_for_method(specs, method, level=None):
    """Choose the registry spec for a predicate's level scope.  An explicit `@PBE0` -> the PBE0
    (iqadft) reading of a shared *_QTAIM column, `@QTAIM` (= M06-2X) -> the M06-2X (qtaim) reading;
    with no scope the Selection's `level` (the subcommand's --level) decides which reading a shared
    *_QTAIM column binds to, so `--level pbe0` filters on the PBE0 reading; otherwise the first spec
    (M06-2X for shared QTAIM columns; the sole spec for IQA / SQM columns)."""
    want = None
    if method:
        m = method.upper()
        want = "iqadft" if m == "PBE0" else "qtaim" if m == "QTAIM" else None
    elif level in ("qtaim", "iqadft"):
        want = level
    if want:
        for s in specs:
            if s.side == want:
                return s
    return specs[0]


def _predicate_sql(pred, con, level=None, sqm_method=None):
    """Compile one Predicate to (sql_fragment, params) against the `manifest m` core query.

    `level` ('qtaim' | 'iqadft' | None) and `sqm_method` are the Selection-wide defaults (from the
    subcommand's --level / --method) used when the predicate carries no explicit `@LEVEL` scope.
    """
    col = pred.column
    if col in MANIFEST_COLUMNS:
        return _cond_sql(f"m.{col}", pred)
    if col in SANITY_COLUMNS:
        method = pred.method or ("PBE0" if level == "iqadft" else "QTAIM")
        cond, cp = _cond_sql(f"s.{col}", pred)
        return (f"EXISTS (SELECT 1 FROM sanity s WHERE s.structure_id = m.structure_id "
                f"AND s.method = ? AND {cond})"), [method] + cp
    # any other name must be a live property column -> resolve its view via the registry
    if con is None:
        raise ValueError(
            f"predicate on {col!r} needs a database connection to resolve its view; "
            f"call selection.resolve(con)/count(con) (the CLI does this for you)")
    from .registry import get_registry
    specs = get_registry(con)["by_name"].get(col)
    if not specs:
        raise ValueError(
            f"unknown predicate column {col!r} (not a manifest, sanity, or property column)")
    # bind the requested LEVEL, not always M06-2X: a shared *_QTAIM column is registered on both
    # sides, so `@PBE0` selects the PBE0 reading and `@M06-2X`/`@QTAIM` (or no scope) the M06-2X one.
    spec = _pick_spec_for_method(specs, pred.method, level)
    cond, cp = _cond_sql(f"v.{col}", pred)
    mclause, mp = "", []
    if spec.level_aware:
        method = pred.method or sqm_method or SQM_METHODS[0]
        if method not in SQM_METHODS:
            raise ValueError(f"predicate on the SQM column {col!r}: scope {method!r} is not an SQM "
                             f"method (use one of {list(SQM_METHODS)})")
        mclause, mp = " AND v.method = ?", [method]
    elif pred.method in SQM_METHODS:
        raise ValueError(f"predicate on {col!r}: scope @{pred.method} does not apply to a QTAIM/IQA "
                         f"column (use @QTAIM for M06-2X or @PBE0)")
    return (f"EXISTS (SELECT 1 FROM {spec.view} v "
            f"WHERE v.structure_id = m.structure_id{mclause} AND {cond})"), mp + cp


@dataclass(frozen=True)
class Selection:
    elements: tuple = None          # structures containing ALL of these elements
    elements_only: tuple = None     # composition is a SUBSET of these elements
    source_family: tuple = None
    natoms: tuple = None            # (lo, hi) inclusive
    formula: str = None             # exact, or a LIKE pattern if it contains '%'
    tier: str = None                # sane_qtaim | sane_mopac_all5 | sane_all | recommended
    ids: tuple = None               # explicit structure_id allow-list
    predicates: tuple = None        # tuple of Predicate: sanity residuals / arbitrary value ranges
    one_per_molecule: bool = False
    molecule_key: str = "inchikey_block1"
    molecule_pick: str = "min_id"   # which conformer to keep per molecule
    sample: int = None
    seed: int = DEFAULT_SEED
    limit: int = None
    level: str = None               # 'pbe0' (PBE0) | 'm062x'/'dft' (M06-2X, default): the level that
                                    # UNSCOPED predicates on shared *_QTAIM columns / sanity rows bind
                                    # to (the subcommand's --level); None = M06-2X
    sqm_method: str = None          # default SQM method for unscoped *_SQM predicates (--method)

    def __post_init__(self):
        norm = _normalize_level(self.level)
        if norm is _INVALID_LEVEL:
            raise ValueError("level must be 'pbe0' (PBE0), 'm062x'/'dft' (M06-2X) or None")
        if norm != self.level:
            object.__setattr__(self, "level", norm)
        if self.sqm_method is not None and self.sqm_method not in SQM_METHODS:
            raise ValueError(f"sqm_method must be one of {list(SQM_METHODS)}")
        if self.elements:
            _check("elements", self.elements, ELEMENTS)
        if self.elements_only:
            _check("elements_only", self.elements_only, ELEMENTS)
        if self.source_family:
            _check("source_family", self.source_family, SOURCE_FAMILIES)
        if self.tier and self.tier not in TIERS:
            raise ValueError(f"tier must be one of {sorted(TIERS)}")
        if self.molecule_key not in MOL_KEYS:
            raise ValueError(f"molecule_key must be one of {sorted(MOL_KEYS)}")
        if self.molecule_pick not in MOL_PICKS:
            raise ValueError(f"molecule_pick must be one of {sorted(MOL_PICKS)}")
        if self.predicates and not all(isinstance(p, Predicate) for p in self.predicates):
            raise ValueError("predicates must be a tuple of Predicate instances")

    # -- SQL compilation -----------------------------------------------------------------
    def _predicates(self, con=None):
        preds, params = [], []
        if self.source_family:
            fams = ",".join(f"'{f}'" for f in self.source_family)
            preds.append(f"m.source_family IN ({fams})")
        if self.natoms:
            preds.append("m.natoms BETWEEN ? AND ?")
            params += [int(self.natoms[0]), int(self.natoms[1])]
        if self.formula:
            if "%" in self.formula:
                preds.append("m.formula LIKE ?")
            else:
                preds.append("m.formula = ?")
            params.append(self.formula)
        if self.tier == "recommended":
            preds.append("m.recommended")
        elif self.tier:
            preds.append(f"m.{self.tier}")
        if self.ids:
            ph = ",".join("?" for _ in self.ids)
            preds.append(f"m.structure_id IN ({ph})")
            params += list(self.ids)
        if self.elements:
            conds = " AND ".join(
                f"count(*) FILTER (WHERE element = '{e}') > 0" for e in self.elements)
            preds.append(f"m.structure_id IN (SELECT structure_id FROM geometry "
                         f"GROUP BY structure_id HAVING {conds})")
        if self.elements_only:
            allowed = ",".join(f"'{e}'" for e in self.elements_only)
            preds.append(f"m.structure_id IN (SELECT structure_id FROM geometry "
                         f"GROUP BY structure_id HAVING bool_and(element IN ({allowed})))")
        if self.predicates:
            for pred in self.predicates:
                frag, fparams = _predicate_sql(pred, con, self.level, self.sqm_method)
                preds.append(frag)
                params += fparams
        return preds, params

    def ids_query(self, con=None):
        """Return (sql, params) selecting the sorted list of structure_ids.

        `con` is only needed when a value predicate targets a live property column (to resolve
        its data view); manifest/sanity predicates and all other filters compile without it.
        """
        preds, params = self._predicates(con)
        where = (" WHERE " + " AND ".join(preds)) if preds else ""
        core = (f"SELECT m.structure_id AS structure_id, m.natoms AS natoms, "
                f"m.recommended AS recommended, m.{self.molecule_key} AS molkey "
                f"FROM manifest m{where}")
        if self.one_per_molecule:
            order = {"min_id": "structure_id ASC",
                     "max_natoms": "natoms DESC, structure_id ASC",
                     "recommended": "recommended DESC, structure_id ASC"}[self.molecule_pick]
            core = (f"SELECT structure_id FROM ({core}) "
                    f"QUALIFY row_number() OVER (PARTITION BY molkey ORDER BY {order}) = 1")
        else:
            core = f"SELECT structure_id FROM ({core})"
        if self.sample:
            core = (f"SELECT structure_id FROM ({core}) "
                    f"USING SAMPLE {int(self.sample)} ROWS (reservoir, {int(self.seed)})")
        sql = f"SELECT structure_id FROM ({core}) ORDER BY structure_id"
        if self.limit:
            sql += f" LIMIT {int(self.limit)}"
        return sql, params

    def resolve(self, con):
        """Materialize the sorted structure_id list."""
        sql, params = self.ids_query(con)
        return [r[0] for r in con.execute(sql, params).fetchall()]

    def count(self, con):
        sql, params = self.ids_query(con)
        return con.execute(f"SELECT count(*) FROM ({sql})", params).fetchone()[0]
