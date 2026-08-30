"""Optional convenience loaders over the exported graph format.

These are NOT the storage format (the ragged npz/hdf5 + sidecars are); they are thin adapters for
popular training stacks, kept behind lazy imports so the core toolkit has no heavy ML dependency.
"""
from __future__ import annotations


def load_npz_as_pyg(path, **kw):
    from .pyg import load_npz_as_pyg as _f
    return _f(path, **kw)


def IqarisGraphDataset(*a, **kw):
    from .pyg import IqarisGraphDataset as _c
    return _c(*a, **kw)
