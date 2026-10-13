"""Verify admitted native arrays through Microgen's public BoxMesh interface."""

import argparse
import importlib.metadata
import json
import os
import sys
from pathlib import Path

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[name] = "1"

import numpy as np
import pyvista as pv
from microgen import BoxMesh, is_periodic

from reference_sensitivity import read

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report


def verify(paths):
    report(0, "Verify admitted meshes through Microgen BoxMesh", eta_s=15)
    rows = []
    for index, path in enumerate(paths):
        artifact, _, provenance = read(path)
        if artifact["workflow"] != "meshers-pass1":
            raise ValueError("This bridge check requires an admitted direct native artifact")
        saved = np.load(path.parent / "mesh.npz")
        points, tets = saved["points"], saved["tetrahedra"]
        grid = pv.UnstructuredGrid(np.c_[np.full(len(tets), 4), tets].ravel(),
                                   np.full(len(tets), pv.CellType.TETRA), points)
        box = BoxMesh.from_pyvista(grid)
        restored = box.to_pyvista()
        if not np.array_equal(restored.points, points) or not np.array_equal(
            restored.cells_dict[pv.CellType.TETRA], tets
        ):
            raise ValueError("Microgen BoxMesh conversion changed native coordinates or tetrahedra")
        if not is_periodic(box.nodes_coords, tol=1e-10):
            raise ValueError("Microgen's public periodic-node check rejected the native mesh")
        rows.append({**provenance, "case": artifact["case"], "cells": artifact["cells"],
                     "tetrahedra": len(tets), "nodes": len(points), "exact_roundtrip": True,
                     "microgen_is_periodic": True, "periodic_tolerance": 1e-10,
                     "face_interior_node_counts": {k: len(v) for k, v in box.faces.items()}})
        report((index + 1) / len(paths), f"Microgen bridge verified: {artifact['case']}",
               eta_s=0 if index + 1 == len(paths) else 10, tetrahedra=len(tets))
    return {"state": "complete", "microgen_version": importlib.metadata.version("microgen"),
            "scope": "Public BoxMesh exact array roundtrip and periodic-node check; no remeshing or new FEA",
            "rows": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", type=Path, nargs="+")
    args = parser.parse_args()
    print(json.dumps(verify(args.artifacts), indent=2))
