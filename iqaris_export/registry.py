"""Property registry: every exportable column, its family, level-awareness and (best-effort)
physical unit.  Built LIVE by DESCRIBE-ing the DuckDB views, so it can never drift from the
database.  Units follow the DB convention (atomic units throughout except lengths in Angstrom);
a compact rule set annotates the common columns, defaulting to `a.u.`.
"""
from __future__ import annotations
import re
from dataclasses import dataclass

from . import config
from . import property_docs

# Columns duplicated from the `geometry` table (provided there, not re-exported as properties).
_GEOM_DUP = {"element", "Z", "X_ANG", "Y_ANG", "Z_ANG", "x_ang", "y_ang", "z_ang"}


@dataclass(frozen=True)
class ColumnSpec:
    name: str          # DB column name (the token users type)
    view: str          # DuckDB view it lives in
    family: str        # molecular | atomic | pairs | cp | angles
    side: str          # 'qtaim' (level-free) or 'sqm' (level-aware)
    level_aware: bool
    unit: str
    dtype: str
    description: str = ""


# unit resolution (best-effort, DB-native units)
_LENGTH = {"CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG", "BPL_QTAIM",
           "GBL_I_QTAIM", "GBL_II_QTAIM", "GBL_III_QTAIM", "GBL_IV_QTAIM"}
_ANGLE = {"BPA_QTAIM", "GBA_I_QTAIM", "GBA_II_QTAIM", "GBA_III_QTAIM", "GBA_IV_QTAIM"}
_BOOL = {"PH_QTAIM", "BONDED_QTAIM"}
_LABEL = {"method", "point_group_SQM", "CP_TYPE_QTAIM", "CP_ATOMS_QTAIM",
          "identity_tool", "identity_status", "source_family", "formula",
          "canonical_smiles", "inchikey", "inchikey_block1", "inchi",
          # structural element-symbol labels: VARCHAR, not numeric
          "element_i", "element_j", "ELEMENT_I", "ELEMENT_J", "ELEMENT_K",
          "pbe0_qc_status", "pbe0_qc_note", "qc_verdict", "qc_fails"}
_COUNT = {"natoms", "n_filled_levels", "n_hbonds", "N_NACP_QTAIM", "N_NNACP_QTAIM",
          "N_BCP_QTAIM", "N_RCP_QTAIM", "N_CCP_QTAIM", "charge", "multiplicity",
          "CP_RANK_QTAIM", "CP_SIG_QTAIM", "CP_BCP_I", "CP_BCP_J"}
# AIMAll radial moments R{-2,-1,0,+1,+2}(A) / GR{-2,-1,0,+1,+2}(A) = <rho*r^k>(A) / <Dot*r^k>(A)
# (Dot = r.grad(Rho), same e*bohr^k ladder).  R_0(A) = <rho>(A) = N(A) is an electron count -> 'e'.
_RADIAL_MOMENTS = {
    "R_M2_QTAIM": "e/bohr²", "R_M1_QTAIM": "e/bohr", "R_0_QTAIM": "e",
    "R_P1_QTAIM": "e·bohr", "R_P2_QTAIM": "e·bohr²",
    "GR_M2_QTAIM": "e/bohr²", "GR_M1_QTAIM": "e/bohr", "GR_0_QTAIM": "e",
    "GR_P1_QTAIM": "e·bohr", "GR_P2_QTAIM": "e·bohr²",
}


# Ab initio IQA (PBE0) energies (Hartree): the *_IQA columns plus the recovered molecular energy
# E_IQA_MOL and the closure residuals DE_IQA_WFN / DE_IQA_RECON (which lack the _IQA suffix).
_IQA_EXTRA = {"E_IQA_MOL", "DE_IQA_WFN", "DE_IQA_RECON", "E_SCF"}


def _is_iqadft_col(col: str) -> bool:
    """True for the ab initio IQA (PBE0) property columns exposed on the `iqadft` side."""
    return col.endswith("_IQA") or col in _IQA_EXTRA


def _unit_for(col: str) -> str:
    c = col
    if _is_iqadft_col(c):                 # all ab initio IQA energies + closure residuals are Hartree
        return "Eh"
    if c.endswith("_ANG") or c in _LENGTH:
        return "Angstrom"
    # COSMO cavity descriptors (MOPAC, semiempirical): one PropDoc covers area + volume
    if c == "cosmo_area_SQM":
        return "Angstrom²"
    if c == "cosmo_vol_SQM":
        return "Angstrom³"
    if c in _ANGLE:
        return "degree"
    if c in _BOOL:
        return "bool"
    if c in _LABEL:
        return "label"
    if c in _COUNT:
        return "count"
    if c == "Z":
        return "atomic_number"
    if c == "LOC_PCT_QTAIM" or c.startswith("pct") or "PCT" in c:
        return "percent"
    # radial moments: must precede the generic energy regex below, which otherwise shadows
    # GR_* via its bare "G" alternative (e.g. would mislabel GR_0_QTAIM as Eh).
    if c in _RADIAL_MOMENTS:
        return _RADIAL_MOMENTS[c]
    # critical-point field densities (evaluated at a point, not integrated)
    if c.startswith("CP_"):
        if "RHO" in c and "DELSQ" not in c and "GRAD" not in c:
            return "e/bohr³"
        if "GRADRHO" in c:
            return "e/bohr⁴"
        if "HESSEIG" in c or c == "CP_DELSQRHO_QTAIM":
            return "e/bohr⁵"
        if "ELLIP" in c:
            return "dimensionless"
        if "EVEC" in c:
            return "dimensionless"
        if c.startswith("CP_ESP"):
            return "Eh/e"
        if c in ("CP_V_QTAIM", "CP_G_QTAIM", "CP_K_QTAIM", "CP_L_QTAIM",
                 "CP_VEN_QTAIM", "CP_VREP_QTAIM"):
            return "Eh/bohr³"
        if c == "CP_VNUC_QTAIM":
            return "Eh/e"
        if c in ("CP_DELSQV_QTAIM", "CP_DELSQVEN_QTAIM", "CP_DELSQVREP_QTAIM",
                 "CP_DELSQG_QTAIM", "CP_DELSQK_QTAIM"):
            return "Eh/bohr⁵"          # Laplacian of an Eh/bohr^3 energy density
        if c.startswith("CP_STRESSEIG_"):
            return "Eh/bohr³"          # eigenvalues of the stress tensor (trace = V)
        if c.startswith("CP_EHRENFEST_"):
            return "Eh/bohr⁴"          # -div(stress tensor) = Ehrenfest force density
        return "a.u."
    # dipole / multipole moments
    if c.startswith("MU4"):
        return "e·bohr²"
    if c.startswith("MU_") or c.startswith("dipole"):
        return "e·bohr"
    # interatomic-surface INTEGRALS  X_IAS(A|B) = ∫_IAS X(r) dS : one power of length
    # LESS than the underlying field X.  Must precede the N_/E.. prefix rules below,
    # which otherwise short-circuit N_IAS -> 'e' and G/K/L/V_IAS -> 'Eh', silently
    # dropping the /bohr.  (AREA_IAS is a genuine bohr^2 area, handled further down.)
    if "IAS" in c and not c.startswith("AREA"):
        if c.startswith("N_"):
            return "e/bohr"             # ∫ ρ dS
        if "DELSQRHO" in c:
            return "e/bohr³"           # ∫ ∇²ρ dS
        return "Eh/bohr"                # ∫ G/K/L/V dS
    # charges, populations, delocalization (electrons)
    if re.match(r"^(q_|N_|LI_|DI_|D2_|Q_CONTRIB|Q_BOND|s_pop|p_pop|d_pop)", c) or c in ("N_QTAIM",):
        return "e"
    if c == "VR_QTAIM":
        return "dimensionless"
    if "IAS" in c and c.startswith("AREA"):
        return "bohr²"
    if c == "NVOL_QTAIM":
        return "e"                     # electron count enclosed by the atomic volume
    if c == "NVOL_DENS_QTAIM":
        return "e/bohr³"              # NVOL_QTAIM / VOL_QTAIM
    if c.startswith("AREA"):
        return "bohr²"                # incl. AREA_ESP_*_IDS (positive/negative ESP surface areas)
    if c.startswith("VOL"):
        return "bohr³"
    # atomic ESP-on-isodensity-surface descriptors (ESP*_IDS, ABSESP_INT_IDS): one PropDoc covers
    # the whole family, whose units are HETEROGENEOUS -- resolve the specific column here so the
    # per-column unit is exact (must precede the generic energy regex, which would mislabel ESP*→Eh).
    if "ESP" in c and "IDS" in c:
        if "VAR" in c:
            return "(Eh/e)²"                          # variances of the surface ESP
        if "2_INT" in c:
            return "Eh²·bohr²/e²"                   # 2nd-moment surface integrals (ESP2/NEG2/POS2)
        if "INT" in c:
            return "Eh·bohr²/e"                       # 1st-moment surface integrals (incl. ABSESP)
        return "Eh/e"                                  # value stats: max / min / avg / mad
    # energies (Hartree): E*/K/G/L/V/T/HOF/HOMO...; the `E` alternative is guarded `E(?![lL])` so it
    # does not swallow the ELEMENT_I/J/K labels (handled by _LABEL above).
    if re.match(r"^(E(?![lL])|K|G|L|V|T|HOF|HOMO|LUMO|gap|IP|EE|VEN|VNN|Exr|Eelstat|Eexch|Eres)", c):
        return "Eh"
    if "IAS" in c:
        return "Eh/bohr"
    return "a.u."


# registry construction (live)
_CACHE = {}

# Label/structural columns duplicated as join keys in the iqa_pbe0_* views; skip them there.
_IQA_JOINDUP = _GEOM_DUP | {"element_i", "element_j"}


def _resolve_unit(col: str) -> str:
    """Exact unit for a column. Prefer the property_docs entry when it carries a single unit; a
    PropDoc documenting a whole family has a ';'/','-joined composite, so fall back to the
    per-column `_unit_for` regex."""
    doc = property_docs.lookup(col)
    unit = getattr(doc, "unit", None) if doc is not None else None
    if unit and ";" not in unit and "," not in unit:
        return unit
    return _unit_for(col)


def get_registry(con):
    """Return {(family, side): [ColumnSpec]} plus a flat name index, memoized per connection.

    Sides: 'qtaim' (M06-2X), 'sqm' (level-aware SQM), 'iqadft' (the PBE0 record: topology from
    qtaim_pbe0_* plus ab initio IQA from iqa_pbe0_*; each spec keeps its own `view`).
    """
    key = id(con)
    if key in _CACHE:
        return _CACHE[key]
    by_group: dict = {}
    by_name: dict = {}
    existing = {r[0] for r in
                con.execute("SELECT table_name FROM information_schema.tables").fetchall()}

    def _register(view, family, side, skip):
        if view is None or view not in existing:
            return
        keys = set(config.FAMILY_KEYS.get(view, []))
        for row in con.execute(f"DESCRIBE {view}").fetchall():
            col, dtype = row[0], row[1]
            if col in keys or col in skip:
                continue
            doc = property_docs.lookup(col)
            spec = ColumnSpec(
                name=col, view=view, family=family, side=side,
                level_aware=view.startswith("iqa_sqm"),
                unit=_resolve_unit(col), dtype=str(dtype),
                description=(doc.name if doc else ""),
            )
            by_group.setdefault((family, side), []).append(spec)
            by_name.setdefault(col, []).append(spec)

    # M06-2X (qtaim side) + the five SQM methods (level-aware sqm side)
    for family, (qview, mview) in config.FAMILY_VIEWS.items():
        _register(qview, family, "qtaim", _GEOM_DUP)
        _register(mview, family, "sqm", _GEOM_DUP)

    # PBE0 record under the single 'iqadft' side; a shared *_QTAIM name then appears twice in
    # by_name (M06-2X spec first, PBE0 second), and consumers pick the reading by side.
    for family, (qpbe0, ipbe0) in getattr(config, "FAMILY_PBE0_VIEWS", {}).items():
        _register(qpbe0, family, "iqadft", _GEOM_DUP)
        _register(ipbe0, family, "iqadft", _IQA_JOINDUP)

    reg = {"by_group": by_group, "by_name": by_name}
    _CACHE[key] = reg
    return reg


# curated defaults + bundles (friendly `--props` tokens)
# Central columns per (family, side) -- what a bare family shorthand expands to.
DEFAULTS = {
    # M06-2X molecular defaults = the CP census + Poincare-Hopf (no total electronic energy; that
    # lives on the PBE0/IQA side, since the additive total is IQA's job).
    ("molecular", "qtaim"): ["N_BCP_QTAIM", "N_RCP_QTAIM", "N_CCP_QTAIM", "PH_QTAIM"],
    ("molecular", "sqm"):   ["Etot_SQM", "HOF_SQM", "Eelec_SQM", "Ecore_SQM",
                              "Etot_inter_SQM", "Eresonance_mol_SQM", "Eexchange_mol_SQM",
                              "Eelstat_mol_SQM", "Edisp_mol_SQM", "dipole_SQM",
                              "HOMO_SQM", "LUMO_SQM", "gap_SQM", "IP_SQM"],
    ("atomic", "qtaim"):    ["q_QTAIM", "N_QTAIM", "K_QTAIM", "G_QTAIM", "L_QTAIM",
                              "LI_QTAIM", "DI_ATOM_QTAIM", "VOL_QTAIM", "MU_INTRA_MAG_QTAIM"],
    ("atomic", "sqm"):      ["q_SQM", "Eintra_SQM", "N_SQM", "s_pop_SQM", "p_pop_SQM",
                              "d_pop_SQM"],
    ("pairs", "qtaim"):     ["DI_QTAIM", "BONDED_QTAIM", "AREA_IAS_001_QTAIM", "V_IAS_QTAIM"],
    ("pairs", "sqm"):       ["Einter_SQM", "Eresonance_SQM", "Eexchange_SQM",
                              "Eelstat_SQM", "Edisp_SQM"],
    # PBE0 record: the QTAIM/QCT headline keys (at the PBE0 density) then the ab initio IQA headline
    # keys. (Rollups + AIMAll halves are bonus, via `with_bonus`.)
    ("molecular", "iqadft"): ["E_SCF", "N_BCP_QTAIM", "N_RCP_QTAIM",
                              "E_IQA_MOL", "EINTRA_SUM_IQA", "EINTER_SUM_IQA",
                              "VCL_SUM_IQA", "VXC_SUM_IQA", "DE_IQA_WFN"],
    ("atomic", "iqadft"):    ["q_QTAIM", "K_QTAIM", "L_QTAIM", "DI_ATOM_QTAIM", "VOL_QTAIM",
                              "E_IQA", "EINTRA_IQA", "EINTER_IQA", "V_IQA", "VC_IQA", "VX_IQA"],
    ("pairs", "iqadft"):     ["DI_QTAIM", "BONDED_QTAIM",
                              "EINT_IQA", "VCL_IJ_IQA", "VXC_IJ_IQA",
                              "VNE_IJ_IQA", "VEN_IJ_IQA", "VEE_IJ_IQA"],
    # QCT critical points + bond-path angles at PBE0 (same headline set as the M06-2X side).
    ("cp", "iqadft"):        ["CP_TYPE_QTAIM", "CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG",
                              "CP_RHO_QTAIM", "CP_DELSQRHO_QTAIM", "CP_ELLIP_QTAIM",
                              "CP_V_QTAIM", "CP_G_QTAIM", "CP_K_QTAIM",
                              "CP_BCP_I", "CP_BCP_J"],
    ("angles", "iqadft"):    ["BPA_QTAIM", "GBA_I_QTAIM", "GBA_II_QTAIM"],
    ("cp", "qtaim"):        ["CP_TYPE_QTAIM", "CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG",
                              "CP_RHO_QTAIM", "CP_DELSQRHO_QTAIM", "CP_ELLIP_QTAIM",
                              "CP_V_QTAIM", "CP_G_QTAIM", "CP_K_QTAIM",
                              "CP_BCP_I", "CP_BCP_J"],
    ("angles", "qtaim"):    ["BPA_QTAIM", "GBA_I_QTAIM", "GBA_II_QTAIM"],
}

# Cross-family named bundles -> explicit column lists.
BUNDLES = {
    "charges":    ["q_QTAIM", "q_SQM"],
    "iqa":        ["Eintra_SQM", "Einter_SQM", "Eresonance_SQM", "Eexchange_SQM",
                   "Eelstat_SQM", "Edisp_SQM"],
    "iqa_dft":    ["E_IQA", "EINTRA_IQA", "EINTER_IQA", "VC_IQA", "VX_IQA",   # PBE0 ab initio IQA
                   "EINT_IQA", "VCL_IJ_IQA", "VXC_IJ_IQA", "E_IQA_MOL"],
    "esp":        ["ESP_INT_IDS_QTAIM", "ESP_MAX_IDS_QTAIM", "ESP_MIN_IDS_QTAIM",
                   "ESP_AVG_IDS_QTAIM"],
    "topology":   ["N_BCP_QTAIM", "N_RCP_QTAIM", "N_CCP_QTAIM", "PH_QTAIM",
                   "CP_TYPE_QTAIM", "CP_RHO_QTAIM", "CP_DELSQRHO_QTAIM", "CP_ELLIP_QTAIM"],
    "multipoles": ["MU_INTRA_MAG_QTAIM", "MU_X_QTAIM", "MU_Y_QTAIM", "MU_Z_QTAIM"],
    "energies":   ["E_SCF", "K_QTAIM", "G_QTAIM", "L_QTAIM", "Etot_SQM", "HOF_SQM",
                   "Einter_SQM", "Eintra_SQM"],
}
