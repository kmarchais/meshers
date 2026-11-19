"""Write a Meshers mesh using the optional VTK HDF writer."""

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkFieldData, vtkUnstructuredGrid
from vtkmodules.vtkIOHDF import vtkHDFWriter

from . import _meshers

if TYPE_CHECKING:
    from .mesh import Mesh


_TET4_NODES = 4


def add_arrays(attributes: vtkFieldData, arrays: dict[str, np.ndarray]) -> None:
    """Attach named arrays, retaining component counts for empty arrays."""
    for name, values in arrays.items():
        array = numpy_to_vtk(values, deep=True)
        array.SetName(name)
        attributes.AddArray(array)


def periodic_masters(mesh: "Mesh") -> np.ndarray:
    """Return the lowest point ID in each periodic equivalence class."""
    parent = np.arange(len(mesh.points), dtype=np.int64)

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for pairs in mesh.periodic_pairs:
        for a, b in pairs:
            ra, rb = root(a), root(b)
            parent[max(ra, rb)] = min(ra, rb)
    return np.array([root(i) for i in range(len(mesh.points))], dtype=np.int64)


def volume_data(mesh: "Mesh") -> dict[str, np.ndarray]:
    """Return volume and shape measures on the final tetrahedra."""
    points, cells = mesh.points, mesh.tetrahedra
    p = points[cells]
    determinant = np.einsum(
        "ij,ij->i",
        p[:, 1] - p[:, 0],
        np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]),
    )
    edge_sum = sum(
        np.sum((p[:, i] - p[:, j]) ** 2, axis=1) for i in range(4) for j in range(i + 1, 4)
    )
    quality = np.sqrt(np.maximum(432 * determinant**2 / edge_sum**3, 0))
    if cells.shape[1] == _TET4_NODES:
        return {"Volume": determinant / 6, "MMGQuality": quality}
    return {
        "Volume": _meshers.quadratic_volumes(points, cells),
        "CornerVolume": determinant / 6,
        "CornerMMGQuality": quality,
    }


def write_grid(mesh: "Mesh", path: str | Path, *, surface: bool) -> None:
    """Write either the boundary or the volume with its solver metadata.

    Raises:
        OSError: VTK could not write the file.
    """
    cells = mesh.surface if surface else mesh.tetrahedra
    points = mesh.points
    # Keep all points in surface exports so periodic master IDs still index the same nodes.
    grid = vtkUnstructuredGrid()
    coordinates = vtkPoints()
    coordinates.SetData(numpy_to_vtk(points, deep=True))
    grid.SetPoints(coordinates)
    offsets = np.arange(len(cells) + 1, dtype=np.int64) * cells.shape[1]
    topology = vtkCellArray()
    topology.SetData(
        numpy_to_vtkIdTypeArray(offsets, deep=True),
        numpy_to_vtkIdTypeArray(cells.ravel(), deep=True),
    )
    vtk_type = {3: 5, 4: 10, 6: 22, 10: 24}[cells.shape[1]]
    grid.SetCells(vtk_type, topology)

    add_arrays(
        grid.GetCellData(), {"BoundaryTag": mesh.boundary_tags} if surface else volume_data(mesh)
    )
    masters = periodic_masters(mesh)
    point_data = {"PeriodicMasterId": masters, "PeriodicShift": points - points[masters]}
    if mesh.parameters is not None:
        point_data["ParameterCoordinates"] = mesh.parameters
    if mesh.constraint_masks is not None:
        point_data["SurfaceConstraints"] = mesh.constraint_masks
    add_arrays(grid.GetPointData(), point_data)
    field_data: dict[str, np.ndarray] = {
        f"PeriodicPairs{axis}": pairs for axis, pairs in enumerate(mesh.periodic_pairs)
    }
    if mesh.periodic_transforms is not None:
        field_data["PeriodicTransforms"] = mesh.periodic_transforms.reshape(3, 16)
    if mesh.constraint_names:
        field_data["ConstraintNamesUTF8"] = np.frombuffer(
            json.dumps(mesh.constraint_names, ensure_ascii=False).encode("utf-8"), dtype=np.uint8
        )
        field_data["FeatureEdges"] = mesh.feature_edges
    add_arrays(grid.GetFieldData(), field_data)

    writer = vtkHDFWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    writer.SetChunkSize(min(writer.GetChunkSize(), max(1, len(points))))
    if not writer.Write() or writer.GetErrorCode():
        raise OSError(f"VTK could not write {path}")
