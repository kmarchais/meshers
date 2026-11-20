# Copyright (c) 2026 Kevin Marchais

"""Direct meshing of bounded implicit solids. CPU float64; Python 3.11+."""

from ._meshers import (
    MAX_CONSTRAINTS,
    CancellationToken,
    CancelledError,
    CompiledField,
    MeshingError,
)
from .expression import TraceError as TraceError
from .expression import compile_field
from .intersection import ADAPTIVE_MAPPED_INTERSECTIONS as ADAPTIVE_MAPPED_INTERSECTIONS
from .intersection import CLOSED_MAPPED_SEAMS as CLOSED_MAPPED_SEAMS
from .intersection import adapt_intersection, generate_intersection
from .mesh import Mesh
from .surface import SURFACE_MESHING_AVAILABLE, SurfaceMesh, generate_surface
from .volume import Field, generate

__all__ = [
    "MAX_CONSTRAINTS",
    "SURFACE_MESHING_AVAILABLE",
    "CancellationToken",
    "CancelledError",
    "CompiledField",
    "Field",
    "Mesh",
    "MeshingError",
    "SurfaceMesh",
    "adapt_intersection",
    "compile_field",
    "generate",
    "generate_intersection",
    "generate_surface",
]
