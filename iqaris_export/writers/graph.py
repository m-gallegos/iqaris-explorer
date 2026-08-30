"""Complete-graph writers -- carry ALL four families (molecular, atomic, pairs, CP/graph) with
their native ragged shapes, keyed by structure_id.  Two backends:
  * npz  : one `graphs/<structure_id>.npz` per structure + an `index.csv` (numpy only, no deps).
  * hdf5 : one file, `/<structure_id>/...` groups (lazy `import h5py`; scalable to large exports).
Arrays avoid object dtype so they reload with allow_pickle=False.

Graph conventions (documented in the sidecars):
  * atom order   = geometry sorted by atom_index; all back-references are 0-based positions here.
  * pair block   = upper-triangular (i<j) `pair_<tag>_ij [E,2]` (0-based) + `pair_<tag> [E,F]`,
                   i.e. a directed-graph edge_index/edge_attr (SchNet4AIM M(M-1)/2).  Any *_IJ/*_JI
                   or *_AB/*_BA directional feature columns are kept in the DB's native i->j / j->i
                   orientation, NOT reordered.
  * CP block     = auxiliary nodes: `cp_pos [C,3]` (Angstrom), `cp_type [C]` (NACP/BCP/RCP/CCP),
                   `cp_feat [C,F]` (+`cp_feat_names`), and `cp_atoms [C,k]` = 0-based atom
                   back-index (-1 padded; BCP via CP_BCP_I/J, ring/cage via CP_ATOMS_QTAIM labels).
"""
from __future__ import annotations
import csv
import numpy as np

from .base import BaseWriter

# CP columns consumed structurally (become cp_pos/cp_type/cp_atoms) rather than plain features.
_CP_STRUCT = {"cp_index", "CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG", "CP_TYPE_QTAIM",
              "CP_ATOMS_QTAIM", "CP_BCP_I", "CP_BCP_J"}


def _split_numeric(names, getter):
    """Return (num_names, num_matrix, str_cols) from a set of columns."""
    num_names, num_cols, str_cols = [], [], {}
    for c in names:
        arr = np.asarray(getter(c))
        if arr.dtype.kind in ("U", "S", "O"):
            str_cols[c] = arr.astype("U")
        else:
            num_names.append(c)
            num_cols.append(arr.astype(float))
    mat = np.stack(num_cols, axis=1) if num_cols else np.zeros((0, 0))
    return num_names, mat, str_cols


def _isnull(v):
    if v is None:
        return True
    try:
        return bool(np.isnan(v))
    except (TypeError, ValueError):
        return False


def _index_maps(rec):
    """Build atom_index -> 0-based position and label ('C1') -> position maps for back-references."""
    n = rec.natoms
    aidx = (np.asarray(rec.atom_index) if rec.atom_index is not None
            else np.arange(1, n + 1))          # fall back to 1..N if not provided
    pos_of = {int(a): p for p, a in enumerate(aidx)}
    label_of = {f"{rec.species[p]}{int(aidx[p])}": p for p in range(n)}
    return pos_of, label_of


def _cp_atoms(cp, ncp, pos_of, label_of):
    """0-based atom back-index per CP, -1 padded to the widest CP.  Prefers CP_BCP_I/J (integer
    atom ids) for BCPs; falls back to parsing the CP_ATOMS_QTAIM label string (ring/cage/NACP)."""
    bi, bj = cp.get("CP_BCP_I"), cp.get("CP_BCP_J")
    atoms_str = cp.get("CP_ATOMS_QTAIM")
    rows = []
    for r in range(ncp):
        idxs = []
        if bi is not None and not _isnull(bi[r]) and not _isnull(bj[r]):
            idxs = [pos_of.get(int(bi[r])), pos_of.get(int(bj[r]))]
        elif atoms_str is not None and isinstance(atoms_str[r], str):
            idxs = [label_of.get(lab) for lab in atoms_str[r].split()]
        rows.append([p for p in idxs if p is not None])
    k = max((len(r) for r in rows), default=0)
    arr = np.full((ncp, k), -1, dtype=np.int32)
    for r, idxs in enumerate(rows):
        arr[r, :len(idxs)] = idxs
    return arr


def record_arrays(rec):
    """Flatten one StructureRecord into a dict of numpy arrays (npz/hdf5 payload)."""
    out = {
        "structure_id": np.array(rec.sid),
        "natoms": np.array(rec.natoms),
        "Z": np.asarray(rec.Z, dtype=int),
        "species": np.asarray(rec.species, dtype="U2"),
        "pos": np.asarray(rec.pos, dtype=float),
    }
    pos_of, label_of = _index_maps(rec)
    if rec.energy is not None:
        out["energy"] = np.array(float(rec.energy))
        out["energy_name"] = np.array(rec.energy_name or "")

    # atomic (all numeric, aligned to geometry order)
    if rec.atomic:
        names = list(rec.atomic.keys())
        out["atomic_names"] = np.asarray(names, dtype="U40")
        out["atomic"] = np.stack([np.asarray(rec.atomic[c], dtype=float) for c in names], axis=1)

    # pairs, per label -> upper-tri 0-based edge_index [E,2] + edge_attr [E,F]
    for label, block in rec.pairs.items():
        tag = label.replace(":", "_")
        cols = [c for c in block if c not in ("i", "j")]
        ij = np.array([[pos_of.get(int(i), -1), pos_of.get(int(j), -1)]
                       for i, j in zip(block["i"], block["j"])], dtype=np.int32)
        out[f"pair_{tag}_ij"] = ij.reshape(-1, 2)
        out[f"pair_{tag}_names"] = np.asarray(cols, dtype="U40")
        out[f"pair_{tag}"] = (np.stack([np.asarray(block[c], dtype=float) for c in cols], axis=1)
                              if cols else np.zeros((len(block["i"]), 0)))

    # critical points -> auxiliary-node block (pos / type / feat / atom back-index)
    if rec.cp:
        cp = rec.cp
        ncp = len(cp["cp_index"]) if "cp_index" in cp else len(next(iter(cp.values())))
        if "cp_index" in cp:
            out["cp_index"] = np.asarray(cp["cp_index"], dtype=int)
        if all(c in cp for c in ("CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG")):
            out["cp_pos"] = np.stack([np.asarray(cp["CP_X_ANG"], dtype=float),
                                      np.asarray(cp["CP_Y_ANG"], dtype=float),
                                      np.asarray(cp["CP_Z_ANG"], dtype=float)], axis=1)
        if "CP_TYPE_QTAIM" in cp:
            out["cp_type"] = np.asarray(cp["CP_TYPE_QTAIM"], dtype="U5")
        out["cp_atoms"] = _cp_atoms(cp, ncp, pos_of, label_of)
        feat_names = [c for c in cp if c not in _CP_STRUCT]
        num_names, mat, str_cols = _split_numeric(feat_names, lambda c: cp[c])
        out["cp_feat_names"] = np.asarray(num_names, dtype="U40")
        out["cp_feat"] = mat
        for c, arr in str_cols.items():          # any non-numeric leftover CP feature
            out[f"cp_str_{c}"] = arr

    # angles
    if rec.angles:
        ang_names = list(rec.angles.keys())
        num_names, mat, _ = _split_numeric(ang_names, lambda c: rec.angles[c])
        out["angle_names"] = np.asarray(num_names, dtype="U40")
        out["angles"] = mat

    # molecular scalars (split numeric / string)
    if rec.molecular:
        mnames = list(rec.molecular.keys())
        num_names, num_vals, str_vals = [], [], {}
        for c in mnames:
            v = rec.molecular[c]
            if isinstance(v, str):
                str_vals[c] = v
            else:
                num_names.append(c)
                num_vals.append(float(v) if v is not None else np.nan)
        out["mol_names"] = np.asarray(num_names, dtype="U40")
        out["mol"] = np.asarray(num_vals, dtype=float)
        if str_vals:
            out["mol_str_names"] = np.asarray(list(str_vals), dtype="U40")
            out["mol_str_values"] = np.asarray(list(str_vals.values()), dtype="U40")
    return out


class _GraphWriter(BaseWriter):
    """Shared sidecar semantics for the ragged graph backends (npz + hdf5)."""

    def _conventions(self):
        c = super()._conventions()
        c.update({
            "keyed_by": "structure_id (one npz file / one HDF5 group per structure)",
            "atom_order": "geometry sorted by atom_index; all back-references are 0-based here.",
            "pair_block": ("upper-triangular (i<j): pair_<tag>_ij [E,2] is a 0-based edge_index "
                           "(directed i->j), pair_<tag> [E,F] the edge_attr; tag = 'qtaim' or "
                           "'sqm_<METHOD>'. E = N(N-1)/2 when all pairs are present."),
            "pair_directional": ("directional feature columns (*_IJ/*_JI, *_AB/*_BA) are kept in "
                                 "the DB's native i->j / j->i orientation and are NOT reordered."),
            "cp_block": ("auxiliary nodes: cp_pos [C,3] Angstrom, cp_type [C] (NACP/BCP/RCP/CCP), "
                         "cp_feat [C,F] (+cp_feat_names), cp_atoms [C,k] = 0-based atom back-index "
                         "(-1 padded; BCP=2 atoms via CP_BCP_I/J, ring/cage via CP_ATOMS_QTAIM)."),
        })
        return c

    def _structural_docs(self):
        return [
            {"column": "pair_<tag>_ij", "family": "pairs", "unit": "0-based atom index",
             "dtype": "INT32[E,2]", "note": "upper-triangular i<j edge_index"},
            {"column": "pair_<tag>", "family": "pairs", "unit": "see pair_<tag>_names",
             "dtype": "FLOAT[E,F]", "note": "edge_attr; directional cols un-swapped"},
            {"column": "cp_pos", "family": "cp", "unit": "Angstrom", "dtype": "FLOAT[C,3]"},
            {"column": "cp_type", "family": "cp", "unit": "label", "dtype": "STR[C]",
             "note": "NACP/BCP/RCP/CCP"},
            {"column": "cp_feat", "family": "cp", "unit": "see cp_feat_names", "dtype": "FLOAT[C,F]"},
            {"column": "cp_atoms", "family": "cp", "unit": "0-based atom index", "dtype": "INT32[C,k]",
             "note": "atoms defining each CP; -1 padded"},
        ]


class GraphNpzWriter(_GraphWriter):
    ext = ".npz"

    def _open(self):
        self._gdir = self.outdir / "graphs"
        self._gdir.mkdir(parents=True, exist_ok=True)
        self._index = []

    def write_record(self, rec):
        arrs = record_arrays(rec)
        np.savez_compressed(self._gdir / f"{rec.sid}.npz", **arrs)
        npairs = sum(len(b["i"]) for b in rec.pairs.values())
        ncp = len(rec.cp.get("cp_index", [])) if rec.cp else 0
        self._index.append((rec.sid, rec.natoms, npairs, ncp,
                            "" if rec.energy is None else rec.energy))
        self.n += 1

    def _close(self):
        with open(self.outdir / "index.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["structure_id", "natoms", "npairs", "ncp", "energy"])
            w.writerows(self._index)


class GraphHdf5Writer(_GraphWriter):
    ext = ".h5"

    def _open(self):
        try:
            import h5py
        except ImportError as e:
            raise RuntimeError(
                "graph-hdf5 needs h5py, which is not installed in this environment. "
                "Use --format graph-npz, or `pip install h5py`.") from e
        self._h5py = h5py
        self._h5 = h5py.File(self.outdir / "iqaris_graph.h5", "w")

    def write_record(self, rec):
        g = self._h5.create_group(rec.sid)
        vlen = self._h5py.string_dtype(encoding="utf-8")
        for k, v in record_arrays(rec).items():
            arr = np.asarray(v)
            if arr.dtype.kind in ("U", "S", "O"):        # h5py needs vlen UTF-8, not numpy 'U'
                arr = arr.astype(vlen)
                g.create_dataset(k, data=arr)
            else:                                        # gzip only chunkable numeric arrays
                kw = {"compression": "gzip"} if (arr.ndim > 0 and arr.size > 0) else {}
                g.create_dataset(k, data=arr, **kw)
        self.n += 1

    def _close(self):
        self._h5.close()
