"""Streaming export engine.

Resolves a Selection to a sorted id list, then walks it in batches.  Per batch it issues one
range-scan per requested view (`WHERE structure_id IN (batch)` -> parquet zonemap pruning on the
clustered, sorted ids), groups by structure_id, assembles a per-structure record (dict of arrays)
and hands it to a writer.  Batch frames are dropped each loop, so memory stays ~one batch.
"""
from __future__ import annotations
import time
import numpy as np

from . import config
from .registry import get_registry, _is_iqadft_col

# CP structural columns always pulled with the `cp` family so the graph writers can build the
# auxiliary-node block (position, type, atom back-index) whatever CP features the user selected.
CP_STRUCTURAL = ["CP_TYPE_QTAIM", "CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG",
                 "CP_BCP_I", "CP_BCP_J", "CP_ATOMS_QTAIM"]


def _cp_cols(plan, side):
    """CP columns to pull for one side ('qtaim' = M06-2X, 'iqadft' = PBE0): the requested features
    plus the structural columns the graph writers need.  Empty when that side asks for no CP data,
    so its view is not scanned at all."""
    base = plan.iqadft_cols("cp") if side == "iqadft" else plan.qtaim_cols("cp")
    if not base:
        return []
    cols = list(base)
    for c in CP_STRUCTURAL:
        if c not in cols:
            cols.append(c)
    return cols


def _pbe0_views(fam, include_problematic):
    """(topology_view, iqa_view) for the PBE0 record of a family, with the `_all` opt-in suffix
    applied when include_problematic.  iqa_view is None for cp/angles (pure QTAIM)."""
    pair = config.FAMILY_PBE0_VIEWS.get(fam)
    if not pair:
        return None, None
    q, i = pair
    if include_problematic:
        suf = config.PBE0_ALL_SUFFIX
        q = (q + suf) if q else None
        i = (i + suf) if i else None
    return q, i


def _scan_pbe0(con, fam, icols, batch, include_problematic):
    """Scan the split PBE0 record for one family and MERGE the two views on their shared keys:
    topology (*_QTAIM) columns from qtaim_pbe0_*, ab initio IQA (*_IQA) columns from iqa_pbe0_*.
    Returns {sid: sub} with all requested PBE0 columns in one frame."""
    qv, iv = _pbe0_views(fam, include_problematic)
    if not qv or not icols:
        return {}
    keys = config.family_keys(qv)
    top = [c for c in icols if not _is_iqadft_col(c)]
    iqa = [c for c in icols if _is_iqadft_col(c)]
    dfq = _scan(con, qv, top, batch, keys) if top else None
    dfi = _scan(con, iv, iqa, batch, config.family_keys(iv)) if (iv and iqa) else None
    if dfq is None and dfi is None:
        return {}
    if dfq is None:
        merged = dfi
    elif dfi is None:
        merged = dfq
    else:
        merged = dfq.merge(dfi, on=keys, how="outer")
    return _groups(merged)


class StructureRecord:
    """Everything requested for one structure, as plain arrays/scalars."""
    __slots__ = ("sid", "natoms", "species", "Z", "pos", "atom_index", "atomic", "pairs", "cp",
                 "angles", "molecular", "energy", "energy_name", "meta")

    def __init__(self, sid):
        self.sid = sid
        self.natoms = 0
        self.species = []
        self.Z = None
        self.pos = None
        self.atom_index = None  # DB atom_index per emitted atom (usually 1..N) -> builds the
                                # atom_index/label -> 0-based position map for pair & CP back-refs
        self.atomic = {}     # colname -> np.array (len N)
        self.pairs = {}      # label ('qtaim'|'sqm:PM7') -> {'i','j', col->arr}
        self.cp = {}         # colname -> np.array (len C)
        self.angles = {}     # colname -> np.array (len A)
        self.molecular = {}  # name -> scalar
        self.energy = None
        self.energy_name = None
        self.meta = {}


def _scan(con, view, cols, ids, keycols, methods=None):
    sel = list(keycols) + [c for c in cols if c not in keycols]
    ph = ",".join("?" for _ in ids)
    sql = f"SELECT {', '.join(sel)} FROM {view} WHERE structure_id IN ({ph})"
    params = list(ids)
    if methods:
        mph = ",".join("?" for _ in methods)
        sql += f" AND method IN ({mph})"
        params += list(methods)
    sql += f" ORDER BY {', '.join(keycols)}"
    return con.execute(sql, params).fetchdf()


def _groups(df):
    if df is None or len(df) == 0:
        return {}
    return {sid: sub for sid, sub in df.groupby("structure_id", sort=False)}


class Exporter:
    def __init__(self, con, selection, plan, energy_target=config.DEFAULT_ENERGY_TARGET,
                 batch_size=200, progress=True, include_problematic=False):
        self.con = con
        self.sel = selection
        self.plan = plan
        self.energy_target = energy_target
        self.batch_size = batch_size
        self.progress = progress
        self.include_problematic = include_problematic
        get_registry(con)   # warm the registry

    def run(self, writer):
        ids = self.sel.resolve(self.con)
        total = len(ids)
        writer.open(total=total)
        t0 = time.time()
        done = 0
        for start in range(0, total, self.batch_size):
            batch = ids[start:start + self.batch_size]
            for rec in self._assemble_batch(batch):
                writer.write_record(rec)
            done += len(batch)
            if self.progress:
                dt = time.time() - t0
                rate = done / dt if dt else 0
                print(f"  [export] {done:,}/{total:,} structures "
                      f"({rate:5.0f}/s, {dt:5.1f}s)", flush=True)
        writer.close()
        return {"structures": total, "seconds": round(time.time() - t0, 2)}

    # -- per-batch assembly --------------------------------------------------------------
    def _assemble_batch(self, batch):
        con, plan = self.con, self.plan
        geo = _groups(_scan(con, "geometry",
                            ["atom_index", "element", "Z", "x_ang", "y_ang", "z_ang"],
                            batch, ["structure_id"]))
        en = {}
        if self.energy_target:
            edf = _scan(con, "energies", [self.energy_target], batch, ["structure_id"])
            en = dict(zip(edf["structure_id"], edf[self.energy_target]))

        # per-structure pbe0_qc_status (opt-in flagged mode): one scalar from the manifest
        qc = {}
        if self.include_problematic:
            qdf = _scan(con, "manifest", [config.PBE0_QC_STATUS_COL], batch, ["structure_id"])
            qc = dict(zip(qdf["structure_id"], qdf[config.PBE0_QC_STATUS_COL]))

        fam_groups = {}   # (view) -> {sid: sub}
        for fam in plan.families:
            qview, mview = config.FAMILY_VIEWS[fam]
            qcols = _cp_cols(plan, "qtaim") if fam == "cp" else plan.qtaim_cols(fam)
            scols = plan.sqm_cols(fam)
            if qcols:
                fam_groups[qview] = _groups(
                    _scan(con, qview, qcols, batch, config.FAMILY_KEYS[qview]))
            if scols and plan.methods:
                fam_groups[mview] = _groups(
                    _scan(con, mview, scols, batch, config.FAMILY_KEYS[mview],
                          methods=plan.methods))
            icols = _cp_cols(plan, "iqadft") if fam == "cp" else plan.iqadft_cols(fam)
            if config.FAMILY_PBE0_VIEWS.get(fam) and icols:
                fam_groups[f"pbe0:{fam}"] = _scan_pbe0(
                    con, fam, icols, batch, self.include_problematic)

        for sid in batch:
            if sid not in geo:
                continue
            yield self._assemble_one(sid, geo[sid], en.get(sid), fam_groups, qc.get(sid))

    def _assemble_one(self, sid, gsub, energy, fam_groups, pbe0_qc_status=None):
        rec = StructureRecord(sid)
        gsub = gsub.sort_values("atom_index")
        rec.natoms = len(gsub)
        rec.species = list(gsub["element"])
        rec.Z = gsub["Z"].to_numpy()
        rec.pos = gsub[["x_ang", "y_ang", "z_ang"]].to_numpy(dtype=float)
        rec.energy = None if energy is None else float(energy)
        rec.energy_name = self.energy_target
        order = gsub["atom_index"].to_numpy()
        rec.atom_index = order
        plan = self.plan
        if pbe0_qc_status is not None:
            rec.meta["pbe0_qc_status"] = pbe0_qc_status

        for fam in plan.families:
            qview, mview = config.FAMILY_VIEWS[fam]
            qcols = _cp_cols(plan, "qtaim") if fam == "cp" else plan.qtaim_cols(fam)
            scols = plan.sqm_cols(fam)
            icols = _cp_cols(plan, "iqadft") if fam == "cp" else plan.iqadft_cols(fam)
            iview = config.FAMILY_PBE0_VIEWS.get(fam)
            qsub = fam_groups.get(qview, {}).get(sid)
            msub = fam_groups.get(mview, {}).get(sid)
            isub = fam_groups.get(f"pbe0:{fam}", {}).get(sid) if iview else None

            if fam == "molecular":
                if qsub is not None and len(qsub):
                    row = qsub.iloc[0]
                    for c in qcols:
                        rec.molecular[c] = _scalar(row[c])
                if msub is not None:
                    for meth, mrow in msub.groupby("method"):
                        r = mrow.iloc[0]
                        for c in scols:
                            rec.molecular[_mname(c, meth, plan)] = _scalar(r[c])
                if isub is not None and len(isub):  # PBE0 ab initio IQA (single level, unique names)
                    row = isub.iloc[0]
                    for c in icols:
                        rec.molecular[c] = _scalar(row[c])

            elif fam == "atomic":
                if qsub is not None:
                    self._merge_atomic(rec, qsub, qcols, order, suffix=None)
                if msub is not None:
                    for meth, mrow in msub.groupby("method"):
                        self._merge_atomic(rec, mrow, scols, order,
                                           suffix=(meth if plan.suffix_methods else None))
                if isub is not None:
                    self._merge_atomic(rec, isub, icols, order, suffix=None)

            elif fam == "pairs":
                if qsub is not None and len(qsub):
                    rec.pairs["qtaim"] = _pair_block(qsub, qcols)
                if msub is not None:
                    for meth, mrow in msub.groupby("method"):
                        rec.pairs[f"sqm:{meth}"] = _pair_block(mrow, scols)
                if isub is not None and len(isub):
                    rec.pairs["iqadft"] = _pair_block(isub, icols)

            elif fam == "cp":
                # PBE0 (iqadft) CP data when the PBE0 side is active, else M06-2X (qtaim); resolve
                # keeps only one side per family, so at most one of qsub/isub is populated.
                csub, ccols = (isub, icols) if (isub is not None and len(isub)) else (qsub, qcols)
                if csub is not None and len(csub):
                    csub = csub.sort_values("cp_index")
                    rec.cp["cp_index"] = csub["cp_index"].to_numpy()
                    for c in ccols:
                        rec.cp[c] = csub[c].to_numpy()

            elif fam == "angles":
                asub, acols = (isub, icols) if (isub is not None and len(isub)) else (qsub, qcols)
                if asub is not None and len(asub):
                    for c in ["ANGLE_I", "ANGLE_J", "ANGLE_K"] + acols:
                        if c in asub.columns:
                            rec.angles[c] = asub[c].to_numpy()
        return rec

    def _merge_atomic(self, rec, sub, cols, order, suffix):
        sub = sub.set_index("atom_index")
        aligned = sub.reindex(order)
        if aligned[cols].isna().all(axis=None) if cols else False:
            rec.meta["atomic_missing"] = True
        for c in cols:
            name = f"{c}__{suffix}" if suffix else c
            rec.atomic[name] = aligned[c].to_numpy()


def _mname(col, meth, plan):
    return f"{col}__{meth}" if plan.suffix_methods else col


def _pair_block(sub, cols):
    block = {"i": sub["i"].to_numpy(), "j": sub["j"].to_numpy()}
    for c in cols:
        block[c] = sub[c].to_numpy()
    return block


def _scalar(v):
    try:
        if hasattr(v, "item"):
            return v.item()
    except Exception:
        pass
    return v
