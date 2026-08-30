"""PyTorch-Geometric adapter over one exported `graphs/<structure_id>.npz` (the graph-npz format).

Convenience only -- the ragged npz/hdf5 + sidecars remain the portable storage format; this maps a
single structure into a `torch_geometric.data.Data` so it drops into a PyG training loop.  `torch`
and `torch_geometric` are imported lazily; a clear error is raised if they are absent (they are NOT
core deps of iqaris_export).

Mapping (see the export's `columns.json` conventions):
  pos        <- pos [N,3] (Angstrom)          z         <- Z [N]
  x          <- atomic [N,F_atom]             y         <- energy (scalar target; may be absent)
  edge_index <- pair_<tag>_ij.T [2,E] (0-based, upper-tri i<j)   edge_attr <- pair_<tag> [E,F_pair]
  CP auxiliary nodes ride as extra tensors (cp_pos / cp_type / cp_feat / cp_atoms) -- a consumer can
  fold them into a heterograph; they are not standard PyG `Data` fields so are attached as attrs.

Note on `undirected=True`: edges are duplicated to (j,i); directional edge features (*_IJ/*_JI,
*_AB/*_BA) are NOT swapped, so prefer directed use (default) when those columns are present.
"""
from __future__ import annotations
import numpy as np


def _require_pyg():
    try:
        import torch                       # noqa: F401
        from torch_geometric.data import Data
    except ImportError as e:               # pragma: no cover - optional dep
        raise RuntimeError(
            "load_npz_as_pyg needs torch + torch_geometric, which are not installed. "
            "`pip install torch torch_geometric`, or read the npz directly with numpy.") from e
    return torch, Data


def _pick_pair_tag(keys, prefer=("qtaim",)):
    tags = sorted({k[len("pair_"):-len("_ij")] for k in keys
                   if k.startswith("pair_") and k.endswith("_ij")})
    if not tags:
        return None
    for p in prefer:
        if p in tags:
            return p
    return tags[0]


def load_npz_as_pyg(path, pair_tag=None, undirected=False):
    """Load one exported .npz (path or an already-loaded dict) into a torch_geometric Data."""
    torch, Data = _require_pyg()
    d = path if isinstance(path, dict) else dict(np.load(path, allow_pickle=False))
    t = lambda a, dt: torch.as_tensor(np.asarray(a), dtype=dt)

    data = Data()
    data.pos = t(d["pos"], torch.float)
    data.z = t(d["Z"], torch.long)
    data.num_nodes = int(d["natoms"])
    data.structure_id = str(d["structure_id"])
    if "atomic" in d and d["atomic"].size:
        data.x = t(d["atomic"], torch.float)
    if "energy" in d:
        data.y = t(float(d["energy"]), torch.float).reshape(1)

    tag = pair_tag or _pick_pair_tag(d.keys())
    if tag is not None:
        ij = np.asarray(d[f"pair_{tag}_ij"])              # [E,2], 0-based, i<j
        attr = np.asarray(d[f"pair_{tag}"])               # [E,F]
        if undirected and len(ij):
            ij = np.concatenate([ij, ij[:, ::-1]], axis=0)
            attr = np.concatenate([attr, attr], axis=0)
        data.edge_index = t(ij.T, torch.long) if len(ij) else torch.zeros((2, 0), dtype=torch.long)
        data.edge_attr = t(attr, torch.float)
        if f"pair_{tag}_names" in d:
            data.edge_attr_names = [str(x) for x in d[f"pair_{tag}_names"]]

    # CP auxiliary nodes (attached as attributes; not standard PyG node fields)
    if "cp_pos" in d:
        data.cp_pos = t(d["cp_pos"], torch.float)
        data.cp_atoms = t(d["cp_atoms"], torch.long) if "cp_atoms" in d else None
        if "cp_feat" in d and d["cp_feat"].size:
            data.cp_feat = t(d["cp_feat"], torch.float)
        if "cp_type" in d:
            data.cp_type = [str(x) for x in d["cp_type"]]
    return data


class IqarisGraphDataset:
    """Minimal in-memory PyG dataset over a `graphs/` dir of exported .npz files.

    Thin wrapper (not a full `InMemoryDataset` with on-disk processing) -- reads the `index.csv`
    order and lazily converts each .npz on access.  For large exports, stream the npz/hdf5 yourself.
    """
    def __init__(self, graphs_dir, pair_tag=None, undirected=False):
        from pathlib import Path
        self.dir = Path(graphs_dir)
        gdir = self.dir / "graphs" if (self.dir / "graphs").is_dir() else self.dir
        self.files = sorted(gdir.glob("*.npz"))
        self.pair_tag, self.undirected = pair_tag, undirected

    def __len__(self):
        return len(self.files)

    def __getitem__(self, i):
        return load_npz_as_pyg(self.files[i], pair_tag=self.pair_tag, undirected=self.undirected)
