"""BaseWriter + self-documenting sidecars.

Every export writes two sidecar files into the output dir:
  * manifest.json -- the exact selection, property plan, level, energy target, counts, timestamp
    (so any export is reproducible), and
  * columns.json  -- per emitted column: family, side, source view, best-effort unit, dtype.
"""
from __future__ import annotations
import json
import time
from dataclasses import asdict
from pathlib import Path


class BaseWriter:
    ext = ""

    def __init__(self, outdir, plan, reg, selection, energy_target, level="dft", args=None):
        self.outdir = Path(outdir)
        self.plan = plan
        self.reg = reg
        self.sel = selection
        self.energy_target = energy_target
        self.level = level
        self.args = dict(args or {})
        self.n = 0
        self.total = 0

    # lifecycle ------------------------------------------------------------------
    def open(self, total=0):
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.total = total
        self._open()

    def _open(self):
        pass

    def write_record(self, rec):
        raise NotImplementedError

    def close(self):
        self._close()
        self._write_sidecars()

    def _close(self):
        pass

    # sidecars -------------------------------------------------------------------
    def _column_docs(self):
        by_name = self.reg["by_name"]
        docs = [{"column": "structure_id", "family": "key", "unit": "id", "dtype": "VARCHAR"},
                {"column": "species", "family": "geometry", "unit": "symbol", "dtype": "VARCHAR"},
                {"column": "Z", "family": "geometry", "unit": "atomic_number", "dtype": "INTEGER"},
                {"column": "pos", "family": "geometry", "unit": "Angstrom", "dtype": "DOUBLE[3]"}]
        for (family, side), cols in self.plan.columns.items():
            for col in cols:
                spec = next((s for s in by_name.get(col, [])
                             if s.family == family and s.side == side), None)
                docs.append({
                    "column": col, "family": family, "side": side,
                    "view": spec.view if spec else None,
                    "unit": spec.unit if spec else "a.u.",
                    "dtype": spec.dtype if spec else None,
                })
        docs.extend(self._structural_docs())
        return docs

    # -- overridable hooks: format-specific array layout + conventions -----------
    def _structural_docs(self):
        """Extra column docs for arrays a writer synthesizes beyond the plan (e.g. graph blocks)."""
        return []

    def _conventions(self):
        """Semantics a consumer must know to read this export correctly."""
        return {
            "no_forces": ("IQARIS source calculations are single-point: NO atomic forces are "
                          "available. Exports carry total energy + local quantum labels only."),
            "units": "atomic units throughout except positions in Angstrom.",
            "energy_target": self.energy_target,
        }

    def _write_sidecars(self):
        conventions = self._conventions()
        manifest = {
            "toolkit": "iqaris_export",
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "format": self.__class__.__name__,
            "selection": asdict(self.sel),
            "level": self.level,
            "methods": self.plan.methods,
            "energy_target": self.energy_target,
            "families": self.plan.families,
            "columns_per_group": {f"{f}:{s}": cols for (f, s), cols in self.plan.columns.items()},
            "structures_written": self.n,
            "conventions": conventions,
            "note": conventions.get("no_forces"),
        }
        (self.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
        (self.outdir / "columns.json").write_text(json.dumps(self._column_docs(), indent=2))
