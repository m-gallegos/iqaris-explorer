"""matrices.py -- atom-condensed matrix representations of a structure, built on demand.

IQARIS stores every atom pair (i < j, bonded or not), so two matrix representations follow directly
from the released columns, with no new calculation and nothing stored twice:

* ``ldm`` -- the QTAIM localization-delocalization matrix (LDM)::

      L[A, A] = lambda(A)          LI_QTAIM
      L[A, B] = delta(A, B) / 2    DI_QTAIM / 2

  Each row sums to the atomic electron population N(A) (``N_QTAIM``) and the whole matrix to the
  number of electrons.  Levels: M06-2X (``m062x``) and PBE0 (``pbe0``).

* ``iqa`` -- the IQA energy matrix::

      E[A, A] = E_intra(A)         EINTRA_IQA | Eintra_SQM
      E[A, B] = E_int(A, B) / 2    EINT_IQA   | Einter_SQM

  Each row sums to the additive atomic energy (``E_IQA`` | ``Eadd_SQM``) and the whole matrix to the
  molecular energy.  Levels: ab initio IQA at PBE0 (``pbe0``) and IQA-SQM at the five semiempirical
  Hamiltonians (``PM7``, ``PM6``, ``PM6-D3H4``, ``PM6-D3H4X``, ``PM6-ORG``).  The pairwise empirical
  dispersion of the SQM levels is kept out, as in ``Eadd_SQM``, unless ``dispersion=True``.

:func:`condense` sums a matrix over groups of atoms (a ring, a functional group, a fragment), giving
the group-level matrix whose diagonal is the group's localization (or intra-group energy) and whose
off-diagonal entries are the inter-group delocalization (or interaction energy) halves.

Read-only: everything comes from the DuckDB views; nothing is written to the database.

Python::

    import iqaris_export as ie
    con = ie.connect()
    L = ie.ldm(con, "dsgdb9nsd_000013", level="pbe0")
    L.to_dataframe()                     # labeled n x n table (e)
    L.max_rowsum_deviation()             # |sum_B L[A,B] - N(A)|, a built-in consistency check
    E = ie.iqa_matrix(con, "dsgdb9nsd_000013", level="PM7")
    ie.condense(E, [[1, 2, 3], [4, 5]])   # group-level energy matrix

CLI:  iqaris matrix --sid ID [--kind ldm|iqa] [--level m062x|pbe0|PM7|...] [--fragments "1-6;7,8"]
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import config

KINDS = ("ldm", "iqa")
_KIND_ALIAS = {"ldm": "ldm", "iqa": "iqa", "iqam": "iqa", "energy": "iqa"}
_CHUNK = 500                                         # structure ids per SQL round-trip


@dataclass
class AtomMatrix:
    """A symmetric atom- (or group-) condensed matrix of one structure.

    ``reference`` holds the stored quantity each row must sum to (N_QTAIM, E_IQA or Eadd_SQM), so
    :meth:`max_rowsum_deviation` is a direct consistency check against the database; it is ``None``
    when no stored counterpart exists (e.g. SQM matrices with ``dispersion=True``).
    """
    structure_id: str
    kind: str
    level: str
    unit: str
    labels: list
    values: np.ndarray
    atoms: list = field(default_factory=list)        # 1-based atom_index per row (or per group: list)
    reference: np.ndarray | None = None
    reference_name: str | None = None
    dispersion: bool = False                         # SQM energy matrix with Edisp/2 added

    @property
    def n(self):
        return len(self.labels)

    @property
    def row_sums(self):
        return self.values.sum(axis=1)

    @property
    def total(self):
        return float(self.values.sum())

    def max_rowsum_deviation(self):
        """max_A |sum_B M[A,B] - reference(A)| (NaN when there is no stored reference)."""
        if self.reference is None:
            return float("nan")
        return float(np.nanmax(np.abs(self.row_sums - self.reference)))

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame(self.values, index=self.labels, columns=self.labels)

    def to_text(self, max_n=24):
        """Human-readable block: a two-line header (structure, matrix, level, unit, row-sum check)
        and the labeled matrix; matrices larger than ``max_n`` are shown truncated."""
        grouped = any(isinstance(a, (list, tuple)) for a in self.atoms)
        name = "LDM" if self.kind == "ldm" else "IQA energy matrix"
        digits = 3 if self.unit == "e" else 4
        what = "groups" if grouped else "atoms"
        lines = [f"{self.structure_id} | {name} | {self.level} | {self.n} {what} | unit {self.unit}"]
        if self.reference is not None:
            lines.append(f"row sums vs stored {self.reference_name}: max deviation "
                         f"{self.max_rowsum_deviation():.1e} {self.unit} | total "
                         f"{self.total:.{digits}f} {self.unit}")
        else:
            lines.append(f"total {self.total:.{digits}f} {self.unit} (no stored row-sum reference: "
                         f"{'pairwise dispersion included' if self.dispersion else 'partial groups'})")
        if grouped:
            lines.append("groups: " + "; ".join(
                f"{lab}={','.join(str(a) for a in g)}" for lab, g in zip(self.labels, self.atoms)))
        lines.append("")
        lines.append(self.to_dataframe().to_string(
            float_format=lambda v: f"{v:.{digits}f}", max_rows=max_n, max_cols=max_n))
        if self.n > max_n:
            lines.append(f"({self.n} x {self.n} matrix shown truncated; --format csv, json or npz "
                         f"writes it in full)")
        return "\n".join(lines)

    def as_dict(self):
        """JSON-ready dictionary (used by the CLI)."""
        def _f(x):
            return None if x is None or not np.isfinite(x) else float(x)
        return {
            "structure_id": self.structure_id, "kind": self.kind, "level": self.level,
            "unit": self.unit, "n": self.n, "labels": list(self.labels),
            "atoms": [list(a) if isinstance(a, (list, tuple)) else int(a) for a in self.atoms],
            "matrix": [[_f(v) for v in row] for row in self.values],
            "row_sums": [_f(v) for v in self.row_sums],
            "reference_name": self.reference_name, "dispersion": self.dispersion,
            "reference": None if self.reference is None else [_f(v) for v in self.reference],
            "max_rowsum_deviation": _f(self.max_rowsum_deviation()),
            "total": _f(self.total),
        }


# ------------------------------------------------------------------ level / kind handling
def _norm_kind(kind):
    k = _KIND_ALIAS.get(str(kind).strip().lower())
    if k is None:
        raise ValueError(f"unknown matrix kind {kind!r}: use one of {', '.join(KINDS)}")
    return k


def _norm_level(level):
    """User token -> canonical level label: 'M06-2X', 'PBE0' or an SQM method name."""
    low = str(level).strip().lower()
    if low in ("m062x", "m06-2x", "dft", "qtaim"):
        return "M06-2X"
    if low in ("pbe0", "iqadft"):
        return "PBE0"
    for m in config.SQM_METHODS:
        if low == m.lower():
            return m
    raise ValueError(f"unknown level {level!r}: use m062x, pbe0 or one of "
                     f"{', '.join(config.SQM_METHODS)}")


def _plan(kind, level, include_problematic, dispersion):
    """Views and columns for (kind, level). Returns a dict, or raises ValueError for a combination
    the database does not carry."""
    kind, lev = _norm_kind(kind), _norm_level(level)
    sfx = config.PBE0_ALL_SUFFIX if (include_problematic and lev == "PBE0") else ""
    if kind == "ldm":
        if lev in config.SQM_METHODS:
            raise ValueError("the LDM is a QTAIM (real-space) quantity: use level m062x or pbe0; "
                             "the semiempirical partition carries no localization indices")
        stem = "qtaim_m062x" if lev == "M06-2X" else "qtaim_pbe0"
        p = dict(atomic=f"{stem}_atomic{sfx}", pairs=f"{stem}_pairs{sfx}", method=None,
                 diag="LI_QTAIM", off="DI_QTAIM", ref="N_QTAIM", disp=None, unit="e")
    else:
        if lev == "M06-2X":
            raise ValueError("there is no IQA energy decomposition at M06-2X: use level pbe0 "
                             "(ab initio IQA) or a semiempirical level (PM7, PM6, ...)")
        if lev == "PBE0":
            p = dict(atomic=f"iqa_pbe0_atomic{sfx}", pairs=f"iqa_pbe0_pairs{sfx}", method=None,
                     diag="EINTRA_IQA", off="EINT_IQA", ref="E_IQA", disp=None, unit="Eh")
        else:
            p = dict(atomic="iqa_sqm_atomic", pairs="iqa_sqm_pairs", method=lev,
                     diag="Eintra_SQM", off="Einter_SQM", ref="Eadd_SQM",
                     disp="Edisp_SQM" if dispersion else None, unit="Eh")
    if dispersion and p["disp"] is None:
        raise ValueError("dispersion=True applies only to the semiempirical IQA matrix (the SQM "
                         "levels carry a separate pairwise empirical dispersion term)")
    p.update(kind=kind, level=lev)
    return p


# ------------------------------------------------------------------ builders
def _fetch(con, view, cols, sids, method):
    ph = ",".join("?" * len(sids))
    where = f"structure_id IN ({ph})"
    params = list(sids)
    if method is not None:
        where += " AND method = ?"
        params.append(method)
    return con.execute(f"SELECT structure_id, {', '.join(cols)} FROM {view} WHERE {where}",
                       params).fetchdf()


def atom_matrices(con, ids, kind="ldm", level="m062x", include_problematic=False,
                  dispersion=False):
    """Build the matrix for every structure id in ``ids`` (batched SQL).

    Returns ``(matrices, missing)``: ``matrices`` maps structure_id -> :class:`AtomMatrix`, in the
    order of ``ids``; ``missing`` lists the ids with no record at that level (e.g. structures without
    a clean PBE0 record, unless ``include_problematic=True``).
    """
    p = _plan(kind, level, include_problematic, dispersion)
    ids = [str(s) for s in ids]
    out, missing = {}, []
    for k in range(0, len(ids), _CHUNK):
        chunk = ids[k:k + _CHUNK]
        at = _fetch(con, p["atomic"], ["atom_index", "element", p["diag"], p["ref"]], chunk,
                    p["method"])
        pcols = ["i", "j", p["off"]] + ([p["disp"]] if p["disp"] else [])
        pr = _fetch(con, p["pairs"], pcols, chunk, p["method"])
        at_g = {s: g for s, g in at.groupby("structure_id", sort=False)}
        pr_g = {s: g for s, g in pr.groupby("structure_id", sort=False)}
        for sid in chunk:
            a = at_g.get(sid)
            if a is None or len(a) == 0:
                missing.append(sid)
                continue
            a = a.sort_values("atom_index")
            aidx = a["atom_index"].astype(int).tolist()
            row = {ai: r for r, ai in enumerate(aidx)}
            n = len(aidx)
            M = np.full((n, n), np.nan)
            np.fill_diagonal(M, a[p["diag"]].astype(float).to_numpy())
            g = pr_g.get(sid)
            if g is not None:
                ii = g["i"].astype(int).map(row).to_numpy()
                jj = g["j"].astype(int).map(row).to_numpy()
                v = g[p["off"]].astype(float).to_numpy() / 2.0
                if p["disp"]:
                    v = v + np.nan_to_num(g[p["disp"]].astype(float).to_numpy()) / 2.0
                M[ii, jj] = v
                M[jj, ii] = v
            ref = None if p["disp"] else a[p["ref"]].astype(float).to_numpy()
            labels = [f"{e}{ai}" for e, ai in zip(a["element"], aidx)]
            out[sid] = AtomMatrix(structure_id=sid, kind=p["kind"], level=p["level"],
                                  unit=p["unit"], labels=labels, values=M, atoms=aidx,
                                  reference=ref, reference_name=None if p["disp"] else p["ref"],
                                  dispersion=bool(p["disp"]))
    return {s: out[s] for s in ids if s in out}, missing


def _one(con, sid, kind, level, include_problematic, dispersion):
    mats, missing = atom_matrices(con, [sid], kind=kind, level=level,
                                  include_problematic=include_problematic, dispersion=dispersion)
    if missing:
        lev = _norm_level(level)
        if not con.execute("SELECT count(*) FROM manifest WHERE structure_id = ?",
                           [str(sid)]).fetchone()[0]:
            raise ValueError(f"structure_id {sid!r} not found")
        hint = ""
        if lev == "PBE0" and not include_problematic:
            hint = (" (no clean PBE0 record: the structure is quality-flagged or has no PBE0 data; "
                    "include_problematic=True / --include-problematic reads the flagged records)")
        raise ValueError(f"structure_id {sid!r} has no {lev} record{hint}")
    return mats[str(sid)]


def ldm(con, sid, level="m062x", include_problematic=False):
    """QTAIM localization-delocalization matrix of one structure (see module docstring)."""
    return _one(con, sid, "ldm", level, include_problematic, False)


def iqa_matrix(con, sid, level="pbe0", include_problematic=False, dispersion=False):
    """IQA energy matrix of one structure, ab initio (pbe0) or semiempirical (PM7, ...)."""
    return _one(con, sid, "iqa", level, include_problematic, dispersion)


# ------------------------------------------------------------------ group condensation
def parse_groups(spec):
    """Parse ``"1-6;7,8;name=9-12"`` into ``(groups, names)``: groups are ';'-separated, each a
    comma list of 1-based atom indices or ranges, optionally prefixed by ``name=``."""
    groups, names = [], []
    for k, part in enumerate(p for p in str(spec).split(";") if p.strip()):
        name = None
        if "=" in part:
            name, part = part.split("=", 1)
            name = name.strip()
        idx = []
        for tok in part.split(","):
            tok = tok.strip()
            if not tok:
                continue
            if "-" in tok:
                lo, hi = tok.split("-", 1)
                idx.extend(range(int(lo), int(hi) + 1))
            else:
                idx.append(int(tok))
        if not idx:
            raise ValueError(f"empty atom group in {spec!r}")
        groups.append(idx)
        names.append(name or f"G{k + 1}")
    if not groups:
        raise ValueError(f"no atom groups in {spec!r}")
    return groups, names


def condense(m, groups, names=None, rest=True):
    """Sum an :class:`AtomMatrix` over groups of atoms (1-based atom_index lists).

    F[G, H] = sum_{A in G} sum_{B in H} M[A, B].  For the LDM, F[G, G] is the group localization
    lambda(G) and 2 F[G, H] the inter-group delocalization delta(G, H); for the IQA matrix, F[G, G] is
    the intra-group energy and 2 F[G, H] the inter-group interaction energy.  Groups must be disjoint.
    With ``rest=True`` the atoms in no group form a final "rest" group, so the row sums still
    reproduce the group populations (or additive group energies) and the total is unchanged.
    """
    names = list(names) if names else [f"G{k + 1}" for k in range(len(groups))]
    if len(names) != len(groups):
        raise ValueError("names and groups differ in length")
    pos = {ai: r for r, ai in enumerate(m.atoms)}
    seen, rows = set(), []
    for g in groups:
        bad = [a for a in g if a not in pos]
        if bad:
            raise ValueError(f"atom index {bad[0]} is not in {m.structure_id!r} "
                             f"(atoms 1..{len(m.atoms)})")
        dup = seen.intersection(g)
        if dup:
            raise ValueError(f"atom {sorted(dup)[0]} is in more than one group")
        seen.update(g)
        rows.append([pos[a] for a in g])
    groups = [list(g) for g in groups]
    leftover = [a for a in m.atoms if a not in seen]
    if rest and leftover:
        groups.append(leftover)
        rows.append([pos[a] for a in leftover])
        names.append("rest")
    k = len(rows)
    F = np.zeros((k, k))
    for x in range(k):
        for y in range(k):
            F[x, y] = m.values[np.ix_(rows[x], rows[y])].sum()
    ref = None
    if m.reference is not None and (rest or not leftover):     # row sums only close on a partition
        ref = np.array([m.reference[r].sum() for r in rows])
    return AtomMatrix(structure_id=m.structure_id, kind=m.kind, level=m.level, unit=m.unit,
                      labels=names, values=F, atoms=groups, reference=ref,
                      reference_name=m.reference_name if ref is not None else None,
                      dispersion=m.dispersion)


# ------------------------------------------------------------------ heatmap
K_KCAL = 627.509                                     # Hartree -> kcal/mol


def _seq_cmap():
    """Sequential map, light (weak) -> dark (strong): seaborn 'mako' reversed, as in the 3D views."""
    try:
        import seaborn as sns
        return sns.color_palette("mako_r", as_cmap=True)
    except Exception:                                # seaborn is optional: plain matplotlib fallback
        import matplotlib.pyplot as plt
        return plt.get_cmap("Blues")


def matrix_norm(m, scale="auto", energy_unit="kcal/mol"):
    """(plotted off-diagonal array, norm, cmap, colorbar label) for :func:`plot_matrix` and the viewer.

    LDM: delta/2 on a sequential map, logarithmic when the entries span more than ~2 decades.
    IQA: E_int/2 on a diverging map centered at 0 (blue = stabilizing, red = destabilizing),
    linear for small matrices and symmetric-logarithmic (linear within +-1 kcal/mol) for large ones.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize, SymLogNorm, TwoSlopeNorm
    grouped = any(isinstance(a, (list, tuple)) for a in m.atoms)
    pair = "(G,H)" if grouped else "(A,B)"
    off = np.array(m.values, float)
    np.fill_diagonal(off, np.nan)
    if m.kind == "ldm":
        pos = off[np.isfinite(off) & (off > 0)]
        hi = float(pos.max()) if pos.size else 1.0
        lo = max(float(pos.min()) if pos.size else 1e-4, 1e-4, hi * 1e-4)
        use_log = scale == "log" or (scale == "auto" and hi / lo > 100)
        norm = LogNorm(lo, hi) if use_log else Normalize(0.0, hi)
        cmap = _seq_cmap().copy()
        cmap.set_under("white")                      # below the floor (and <= 0) read as "none"
        cmap.set_bad("white")
        label = f"δ{pair}/2  (e{', log' if use_log else ''})"
        return off, norm, cmap, label
    f = K_KCAL if energy_unit == "kcal/mol" else 1.0
    unit = "kcal/mol" if energy_unit == "kcal/mol" else m.unit
    data = off * f
    lim = float(np.nanmax(np.abs(data))) if np.isfinite(data).any() else 1.0
    lim = lim or 1.0
    use_log = scale == "log" or (scale == "auto" and m.n > 40)
    if use_log:
        norm = SymLogNorm(linthresh=1.0 * f / K_KCAL if unit != "kcal/mol" else 1.0,
                          vmin=-lim, vmax=lim, base=10)
    else:
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("white")
    disp = " incl. dispersion" if m.dispersion else ""
    label = f"E_int{pair}/2{disp}  ({unit}{', symlog' if use_log else ''})"
    return data, norm, cmap, label


def plot_matrix(m, out=None, ax=None, scale="auto", energy_unit="kcal/mol", annotate=None,
                title=None, dpi=200, highlight=None):
    """Heatmap of an :class:`AtomMatrix` (returns the matplotlib Figure; saves it when ``out``).

    Only the off-diagonal entries are colored (the diagonal, lambda(A) or E_intra(A), is one to two
    orders of magnitude larger and would wash the map out); diagonal cells are grey and carry their
    value as text for matrices up to 30 x 30.  ``annotate`` prints every off-diagonal value
    (default: matrices up to 8 x 8, e.g. fragment-condensed ones).  ``ax`` draws into an existing axes.
    ``highlight`` (a 1-based atom index of an atom-level matrix) outlines that atom's row and column,
    e.g. the atom picked in the 3D viewer.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    n = m.n
    data, norm, cmap, label = matrix_norm(m, scale=scale, energy_unit=energy_unit)
    if ax is None:
        side = min(4.6 + 0.16 * n, 10.0)
        fig, ax = plt.subplots(figsize=(side + 1.4, side))
    else:
        fig = ax.figure
    im = ax.imshow(data, cmap=cmap, norm=norm, interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label)
    small = n <= 30
    if annotate is None:
        annotate = n <= 8
    fs = 9 if n <= 8 else 7 if n <= 16 else 5.5
    for k in range(n):
        ax.add_patch(Rectangle((k - 0.5, k - 0.5), 1, 1, facecolor="0.86", edgecolor="none",
                               zorder=2))
        if small:
            d = m.values[k, k]
            ax.text(k, k, f"{d:.2f}" if m.kind == "ldm" else f"{d:.1f}", ha="center",
                    va="center", fontsize=fs, zorder=3)
    if annotate:
        for i in range(n):
            for j in range(n):
                if i == j or not np.isfinite(data[i, j]):
                    continue
                r, g, b, _ = cmap(norm(data[i, j]))
                dark = 0.299 * r + 0.587 * g + 0.114 * b < 0.5
                v = data[i, j]
                txt = f"{v:.3f}" if m.kind == "ldm" else (f"{v:.0f}" if abs(v) >= 10 else f"{v:.1f}")
                ax.text(j, i, txt, ha="center", va="center", fontsize=fs,
                        color="white" if dark else "black", zorder=3)
    if highlight is not None and highlight in m.atoms:
        r = m.atoms.index(highlight)
        for xy, wd, ht in (((-0.5, r - 0.5), n, 1), ((r - 0.5, -0.5), 1, n)):
            ax.add_patch(Rectangle(xy, wd, ht, fill=False, edgecolor="#e8288a", linewidth=1.6,
                                   zorder=4))
    step = 1 if n <= 60 else int(np.ceil(n / 30))
    ticks = list(range(0, n, step))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    tfs = 8 if n <= 16 else 6.5 if n <= 40 else 6
    ax.set_xticklabels([m.labels[t] for t in ticks], rotation=90, fontsize=tfs)
    ax.set_yticklabels([m.labels[t] for t in ticks], fontsize=tfs)
    ax.tick_params(length=2)
    name = "LDM" if m.kind == "ldm" else "IQA energy matrix"
    diag = ("diagonal: λ (e)" if m.kind == "ldm"
            else f"diagonal: E_intra ({m.unit})")
    if title is None:
        title = f"{m.structure_id}  |  {name}  |  {m.level}"
    sub = diag + ("" if small else " (not printed)")
    if m.reference is not None:
        sub += f"  |  max |row sum − {m.reference_name}| = {m.max_rowsum_deviation():.1e} {m.unit}"
    ax.set_title(f"{title}\n{sub}", fontsize=9)
    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
    return fig


# ------------------------------------------------------------------ one row, for 3D projection
def row_values(con, sid, atom, kind="ldm", level="m062x", include_problematic=False):
    """One matrix row as per-atom values, for projecting onto the 3D structure.

    Returns ``(atoms, values, matrix)``: 1-based atom indices, the pair quantity between ``atom`` and
    every atom (``ldm``: the full delocalization index delta(A,B) = 2 L[A,B], in e; ``iqa``: the full
    interaction energy E_int(A,B) = 2 E[A,B], in kcal/mol), NaN at ``atom`` itself, and the
    :class:`AtomMatrix` it came from.
    """
    m = _one(con, sid, kind, level, include_problematic, False)
    if atom not in m.atoms:
        raise ValueError(f"atom {atom} is not in {sid!r} (atoms 1..{len(m.atoms)})")
    r = m.atoms.index(atom)
    vals = 2.0 * np.array(m.values[r], float)
    if m.kind == "iqa":
        vals = vals * K_KCAL
    vals[r] = np.nan
    return list(m.atoms), vals, m
