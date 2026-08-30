"""Writer registry.  Record-based writers implement BaseWriter; tabular formats are handled by
`export_tabular` (a streaming COPY fast path)."""
from __future__ import annotations

from .base import BaseWriter
from .extxyz import ExtxyzWriter
from .graph import GraphNpzWriter, GraphHdf5Writer
from .ase_db import AseDbWriter
from .tabular import export_tabular

RECORD_WRITERS = {
    "extxyz": ExtxyzWriter,
    "graph-npz": GraphNpzWriter,
    "graph-hdf5": GraphHdf5Writer,
    "ase-db": AseDbWriter,
}
TABULAR_FORMATS = {"csv", "parquet"}
ALL_FORMATS = tuple(RECORD_WRITERS) + tuple(sorted(TABULAR_FORMATS))


def get_writer(fmt):
    if fmt not in RECORD_WRITERS:
        raise ValueError(f"unknown record format {fmt!r}; choose from {list(RECORD_WRITERS)}")
    return RECORD_WRITERS[fmt]
