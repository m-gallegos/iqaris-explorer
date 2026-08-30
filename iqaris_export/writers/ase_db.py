"""Optional ASE writer: build ase.Atoms (info = molecular scalars + energy, arrays = atomic
labels) into an ASE `.db`.  Lazily imports ase; degrades with a clear message if absent (the
repo's DB venv does not ship ase).  Carries GEOMETRY + ATOMIC + MOLECULAR (not pairs/CP).
"""
from __future__ import annotations
import numpy as np

from .base import BaseWriter


def _looks_numeric(s):
    """True if ASE would reject the string as a key_value_pair because it parses as a number.

    ASE refuses a string value in `key_value_pairs` when ``int(s)`` or ``float(s)`` succeeds (to
    keep query types unambiguous), so a purely numeric structure_id such as ``"008"`` cannot be
    stored there.  Such values are routed to the row's `data` blob instead, where they round-trip
    verbatim.
    """
    s = s.strip()
    if not s:
        return False
    try:
        int(s)
        return True
    except ValueError:
        try:
            float(s)
            return True
        except ValueError:
            return False


class AseDbWriter(BaseWriter):
    ext = ".db"

    def _open(self):
        try:
            from ase.db import connect  # noqa
        except ImportError as e:
            raise RuntimeError(
                "ase-db needs ASE, which is not installed in this environment. "
                "Use --format extxyz (also ASE-readable), or `pip install ase`.") from e
        from ase.db import connect
        path = self.outdir / "iqaris.db"
        if path.exists():
            path.unlink()
        self._db = connect(str(path))

    def write_record(self, rec):
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator
        from ase.db.core import reserved_keys
        atoms = Atoms(symbols=rec.species, positions=rec.pos)
        for c, arr in rec.atomic.items():
            try:
                atoms.set_array(c, np.asarray(arr, dtype=float))
            except Exception:
                pass
        # The per-structure energy is attached via a single-point calculator so it is queryable as
        # `row.energy`; ASE reserves "energy" as a key_value_pairs key and would reject it there.
        if rec.energy is not None:
            atoms.calc = SinglePointCalculator(atoms, energy=float(rec.energy))
        info = dict(rec.molecular)
        info["structure_id"] = rec.sid
        # ASE key_value_pairs cannot use a reserved column name (energy, charge, formula, natoms,
        # an element symbol, ...) or a numeric-looking string value (e.g. the structure_id "008");
        # keep the queryable scalars as key_value_pairs and store the rest verbatim in `data`,
        # which round-trips unchanged.
        kvp, data = {}, {}
        for k, v in info.items():
            if k in reserved_keys:
                data[k] = v
            elif isinstance(v, bool) or isinstance(v, (int, float)):
                kvp[k] = v
            elif isinstance(v, str):
                (data if _looks_numeric(v) else kvp)[k] = v
            else:
                data[k] = v
        self._db.write(atoms, key_value_pairs=kvp, data=(data or None))
        self.n += 1

    def _close(self):
        pass
