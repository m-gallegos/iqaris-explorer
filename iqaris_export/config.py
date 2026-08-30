"""Static configuration for the iqaris_export toolkit.

No logic here -- just paths, level names and the family -> view mapping that the rest of
the package builds on.  The IQARIS database is READ-ONLY; nothing in this package ever
writes into DB_ROOT.  All exports land under DEFAULT_OUT (repo `exports/`) or a user path.
"""
from __future__ import annotations
import os
from pathlib import Path

# --- database (read-only) ---------------------------------------------------------------
# Level-explicit schema: qtaim_m062x_* / qtaim_pbe0_* (density topology at M06-2X / PBE0),
# iqa_pbe0_* (ab initio IQA), iqa_sqm_* (five SQM methods). Clean data by default; opt-in *_all
# views add the quality-flagged PBE0 records. Point IQARIS_DB_ROOT at the release directory.
DB_ROOT = os.environ.get("IQARIS_DB_ROOT", "iqaris_db")
CATALOG = "iqaris.duckdb"          # pre-baked DuckDB catalog inside DB_ROOT
VIEWS_SQL = "iqaris_views.sql"     # fallback if the catalog is absent

# --- repo / output ----------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]     # .../IQARIS
DEFAULT_OUT = REPO_ROOT / "exports"

# --- levels of theory -------------------------------------------------------------------
# QTAIM density topology at M06-2X and PBE0 (ab initio IQA rides on the PBE0 record) + five SQM methods.
DFT_LABEL = "M062X"                                  # M06-2X/def2-TZVP (gas phase, QTAIM only)
PBE0_LABEL = "PBE0"                                  # PBE0/def2-TZVP (gas phase, QTAIM + ab initio IQA)
SQM_METHODS = ("PM7", "PM6", "PM6-D3H4", "PM6-D3H4X", "PM6-ORG")
ALL_LEVELS = ("dft", "pbe0") + SQM_METHODS

# Per-structure energy targets (columns of the `energies` view). The DFT total energy is the PBE0
# SCF energy E_PBE0_Eh (also the IQA reference); E_IQA_MOL is the IQA-recovered molecular energy.
ENERGY_TARGETS = (
    "E_PBE0_Eh",                                     # PBE0 total electronic (SCF) energy (default)
    "E_IQA_MOL",                                     # PBE0 IQA-recovered molecular energy
    "Etot_PM7", "Etot_PM6", "Etot_PM6_D3H4", "Etot_PM6_D3H4X", "Etot_PM6_ORG",
    "HOF_PM7", "HOF_PM6", "HOF_PM6_D3H4", "HOF_PM6_D3H4X", "HOF_PM6_ORG",
)
DEFAULT_ENERGY_TARGET = "E_PBE0_Eh"

# --- property families ------------------------------------------------------------------
# Each family maps to an M06-2X QTAIM view and (for mol/atomic/pairs) a level-aware IQA-SQM view.
FAMILIES = ("molecular", "atomic", "pairs", "cp", "angles")
FAMILY_VIEWS = {
    "molecular": ("qtaim_m062x_molecular", "iqa_sqm_molecular"),
    "atomic":    ("qtaim_m062x_atomic",    "iqa_sqm_atomic"),
    "pairs":     ("qtaim_m062x_pairs",     "iqa_sqm_pairs"),
    "cp":        ("qtaim_m062x_cp",         None),
    "angles":    ("qtaim_m062x_angles",     None),
}
# The PBE0 record is split across qtaim_pbe0_* (density topology / QCT) and iqa_pbe0_* (ab initio
# IQA); the engine joins them on the shared keys, so the PBE0 level exposes both. cp/angles are
# pure QTAIM (no iqa_pbe0 view).
FAMILY_PBE0_VIEWS = {
    "molecular": ("qtaim_pbe0_molecular", "iqa_pbe0_molecular"),
    "atomic":    ("qtaim_pbe0_atomic",    "iqa_pbe0_atomic"),
    "pairs":     ("qtaim_pbe0_pairs",     "iqa_pbe0_pairs"),
    "cp":        ("qtaim_pbe0_cp",         None),
    "angles":    ("qtaim_pbe0_angles",     None),
}
# Opt-in `_all` (clean UNION flagged) view suffix for --include-problematic; the flagged rows carry
# a pbe0_qc_status column.  Only the PBE0 families have flagged siblings (M06-2X/SQM cover all).
PBE0_ALL_SUFFIX = "_all"
PBE0_QC_STATUS_COL = "pbe0_qc_status"
# Per-structure key columns of each view (NOT exported as properties; used for grouping).
FAMILY_KEYS = {
    "qtaim_m062x_molecular": ["structure_id"],
    "iqa_sqm_molecular":     ["structure_id", "method"],
    "qtaim_m062x_atomic":    ["structure_id", "atom_index"],
    "iqa_sqm_atomic":        ["structure_id", "method", "atom_index"],
    "qtaim_m062x_pairs":     ["structure_id", "i", "j"],
    "iqa_sqm_pairs":         ["structure_id", "method", "i", "j"],
    "qtaim_m062x_cp":        ["structure_id", "cp_index"],
    "qtaim_m062x_angles":    ["structure_id", "ANGLE_I", "ANGLE_J", "ANGLE_K"],
    # PBE0 record (split): same key columns as the matching M06-2X views.
    "qtaim_pbe0_molecular": ["structure_id"],
    "iqa_pbe0_molecular":   ["structure_id"],
    "qtaim_pbe0_atomic":    ["structure_id", "atom_index"],
    "iqa_pbe0_atomic":      ["structure_id", "atom_index"],
    "qtaim_pbe0_pairs":     ["structure_id", "i", "j"],
    "iqa_pbe0_pairs":       ["structure_id", "i", "j"],
    "qtaim_pbe0_cp":        ["structure_id", "cp_index"],
    "qtaim_pbe0_angles":    ["structure_id", "ANGLE_I", "ANGLE_J", "ANGLE_K"],
}


def family_keys(view):
    """Key columns for a view name, tolerant of the `_all`/`_flagged` opt-in suffixes."""
    base = view
    for suf in ("_all", "_flagged"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return FAMILY_KEYS[base]

BOHR_TO_ANGSTROM = 0.529177210903
DEFAULT_SEED = 20260702
