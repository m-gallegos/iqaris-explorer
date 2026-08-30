"""Extended-XYZ writer -- the MLFF lingua franca (read by ASE, MACE, NequIP, SchNetPack...).

One multi-frame `.extxyz` file.  Carries GEOMETRY + ATOMIC labels (as extended `Properties`
columns) + MOLECULAR scalars + the selected energy target (`energy=`).  Pair/CP/angle data do
NOT fit the per-atom model -- use the graph writer for those.  Written by hand (no ASE dep),
generalizing the gallery `write_xyz` idiom.  NO forces unless `--force-stub` (all-zero placeholder).
"""
from __future__ import annotations
import math

from .base import BaseWriter


def _num(v):
    if v is None:
        return "nan"
    if isinstance(v, bool):
        return "1" if v else "0"
    try:
        f = float(v)
        return "nan" if math.isnan(f) else f"{f:.10g}"
    except (TypeError, ValueError):
        return str(v)


def _kv(key, val):
    if isinstance(val, str):
        return f'{key}="{val}"' if (" " in val or val == "") else f"{key}={val}"
    if isinstance(val, bool):
        return f"{key}={'T' if val else 'F'}"
    if isinstance(val, float):
        return f"{key}={'nan' if math.isnan(val) else f'{val:.10g}'}"
    return f"{key}={val}"


class ExtxyzWriter(BaseWriter):
    ext = ".extxyz"

    def _conventions(self):
        c = super()._conventions()
        c["dropped"] = ("pair / CP / angle data cannot ride in the per-atom extXYZ model and are "
                        "dropped here; use --format graph-hdf5 (or graph-npz) for those families.")
        c["atom_properties"] = ("per-atom labels ride as extended Properties columns; molecular "
                                "scalars + the energy target ride in the comment line (energy=).")
        return c

    def _open(self):
        self._fh = open(self.outdir / "iqaris.extxyz", "w")
        self._warned_drop = False

    def write_record(self, rec):
        atom_cols = list(rec.atomic.keys())
        force_stub = bool(self.args.get("force_stub"))
        if (rec.pairs or rec.cp or rec.angles) and not self._warned_drop:
            print("  [extxyz] note: pair/CP/angle data cannot ride in extXYZ; "
                  "use --format graph-npz for those. Dropping them here.")
            self._warned_drop = True

        props = ["species:S:1", "pos:R:3"]
        for c in atom_cols:
            props.append(f"{c}:R:1")
        if force_stub:
            props.append("forces:R:3")

        comment = [f"Properties={':'.join(props)}"]
        if rec.energy is not None:
            comment.append(_kv("energy", rec.energy))
            # Energy is in Hartree (as everywhere in IQARIS); state it so ASE readers, which
            # otherwise assume eV, do not misinterpret it.
            comment.append("energy_unit=Eh")
        for name, val in rec.molecular.items():
            comment.append(_kv(name, val))
        comment.append(_kv("structure_id", rec.sid))
        comment.append('pbc="F F F"')

        fh = self._fh
        fh.write(f"{rec.natoms}\n")
        fh.write(" ".join(comment) + "\n")
        for a in range(rec.natoms):
            row = [rec.species[a],
                   f"{rec.pos[a, 0]:.8f}", f"{rec.pos[a, 1]:.8f}", f"{rec.pos[a, 2]:.8f}"]
            for c in atom_cols:
                row.append(_num(_at(rec.atomic[c], a)))
            if force_stub:
                row += ["0.0", "0.0", "0.0"]
            fh.write(" ".join(row) + "\n")
        self.n += 1

    def _close(self):
        self._fh.close()


def _at(arr, i):
    v = arr[i]
    try:
        return v.item()
    except AttributeError:
        return v
