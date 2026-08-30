"""PropertySet: resolve friendly `--props` tokens + a `--level` into a concrete pull plan
(which columns from which view, and for which SQM method(s)).

Tokens accepted in `props`:
  * a family shorthand: molecular | atomic | pairs | cp | angles  (-> curated defaults)
  * `all`                                                          (-> all families' defaults)
  * a named bundle: charges | iqa | esp | topology | multipoles | energies
  * an explicit DB column name: e.g. DI_QTAIM, CP_RHO_QTAIM, q_SQM
`with_bonus=True` makes a family shorthand expand to ALL of that family's columns, not just the
curated defaults.
"""
from __future__ import annotations
from dataclasses import dataclass, field

from . import config
from .registry import get_registry, DEFAULTS, BUNDLES, _is_iqadft_col


def normalize_level(level: str):
    """(-> (list of SQM methods, set of allowed property sides)).

    The allowed sides ⊆ {'qtaim','sqm','iqadft'} say which property surface a level exposes;
    'iqadft' is the FULL PBE0 record (QTAIM density topology + QCT + ab initio IQA), 'qtaim' is
    M06-2X, 'sqm' is a semi-empirical IQA method.

    '' / None / 'dft'    -> ([],        {'qtaim'})               QTAIM @ M06-2X
    'pbe0'               -> ([],        {'iqadft'})              QTAIM + QCT + IQA @ PBE0
    a single SQM method  -> ([that],    {'qtaim','sqm'})         IQA-SQM + the M06-2X QTAIM reference
    'all-sqm' / 'allsqm' -> (all five,  {'qtaim','sqm'})
    'all'                -> (all five,  {'qtaim','sqm','iqadft'})  everything the DB carries
    """
    if not level or level.lower() in ("dft", "m062x", "m06-2x", "qtaim"):
        return [], {"qtaim"}
    low = level.lower()
    if low in ("pbe0", "iqadft", "iqa-dft"):
        return [], {"iqadft"}
    if low == "all":
        return list(config.SQM_METHODS), {"qtaim", "sqm", "iqadft"}
    if low in ("all-sqm", "allsqm"):
        return list(config.SQM_METHODS), {"qtaim", "sqm"}
    for m in config.SQM_METHODS:
        if low == m.lower():
            return [m], {"qtaim", "sqm"}
    raise ValueError(
        f"unknown level {level!r}; use m062x (alias dft) | pbe0 | "
        f"{' | '.join(config.SQM_METHODS)} | all-sqm | all")


@dataclass
class Plan:
    columns: dict = field(default_factory=dict)   # (family, side) -> [column names]
    methods: list = field(default_factory=list)   # SQM methods to pull ([] = DFT only)
    iqadft: bool = False                           # pull the PBE0 ab initio IQA side too

    @property
    def families(self):
        fams = []
        for (fam, _side) in self.columns:
            if fam not in fams:
                fams.append(fam)
        return [f for f in config.FAMILIES if f in fams]

    @property
    def suffix_methods(self):
        return len(self.methods) > 1

    def qtaim_cols(self, family):
        return self.columns.get((family, "qtaim"), [])

    def sqm_cols(self, family):
        return self.columns.get((family, "sqm"), [])

    def iqadft_cols(self, family):
        return self.columns.get((family, "iqadft"), []) if self.iqadft else []


def resolve(con, props, level="dft", with_bonus=False):
    """Build a Plan from tokens.  `props` is a list of tokens or a comma string."""
    if isinstance(props, str):
        props = [p.strip() for p in props.split(",") if p.strip()]
    reg = get_registry(con)
    by_group, by_name = reg["by_group"], reg["by_name"]
    methods, sides = normalize_level(level)
    sides = set(sides)

    # ordered accumulation of (family, side, column), de-duplicated
    picked = {}         # (family, side) -> list (order-preserving)
    explicit = {}       # col -> set of native sides (only for explicitly-named columns)

    def add(family, side, col):
        picked.setdefault((family, side), [])
        if col not in picked[(family, side)]:
            picked[(family, side)].append(col)

    def add_family(family):
        for side in ("qtaim", "sqm", "iqadft"):
            specs = by_group.get((family, side))
            if not specs:
                continue
            if with_bonus:
                for s in specs:
                    add(family, side, s.name)
            else:
                for col in DEFAULTS.get((family, side), []):
                    add(family, side, col)

    def add_named(col, is_explicit):
        specs = by_name.get(col)
        if not specs:
            raise ValueError(f"unknown property {col!r} (run `list-props` / pass --props to discover)")
        for s in specs:
            add(s.family, s.side, s.name)
        if is_explicit:
            explicit.setdefault(col, set()).update(s.side for s in specs)

    for tok in props:
        if tok in ("geometry", "geom"):
            continue                          # geometry is always exported
        if tok == "all":
            for fam in config.FAMILIES:
                add_family(fam)
        elif tok in config.FAMILIES:
            add_family(tok)
        elif tok in BUNDLES:
            for col in BUNDLES[tok]:
                add_named(col, is_explicit=False)
        else:
            add_named(tok, is_explicit=True)

    # Level consistency for explicitly named columns whose native side this level doesn't expose:
    # an IQA-only column turns the PBE0 side on; an SQM-only column with no method is a user error;
    # QTAIM columns are always providable.
    for col, native in explicit.items():
        if native & sides:
            continue
        if native == {"iqadft"}:
            sides.add("iqadft")
        elif native == {"sqm"}:
            raise ValueError(
                f"SQM property {col!r} requested but level {level!r} has no SQM method. "
                "Pass --level PM7 (or another method / all-sqm).")

    # Keep only the property sides this level exposes (columns pulled implicitly by a family
    # shorthand / bundle for an inactive side are simply dropped).
    picked = {k: v for k, v in picked.items() if k[1] in sides}

    # When both sides are active for a family (e.g. --level all), the PBE0 side would re-emit the
    # shared *_QTAIM columns under the same names and collide, so restrict it to its unique *_IQA keys
    # (--level pbe0 alone keeps the full PBE0 QTAIM + QCT + IQA surface).
    for fam in config.FAMILIES:
        if (fam, "qtaim") in picked and (fam, "iqadft") in picked:
            picked[(fam, "iqadft")] = [c for c in picked[(fam, "iqadft")] if _is_iqadft_col(c)]
            if not picked[(fam, "iqadft")]:
                del picked[(fam, "iqadft")]

    return Plan(columns=picked, methods=methods, iqadft=("iqadft" in sides))
