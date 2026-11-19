"""Direct meshing of bounded implicit solids. CPU float64; Python 3.11+."""

from . import _meshers
from ._meshers import CancellationToken, CancelledError, CompiledField, MeshingError
from .expression import TraceError as TraceError
from .expression import compile_field
from .intersection import ADAPTIVE_MAPPED_INTERSECTIONS as ADAPTIVE_MAPPED_INTERSECTIONS
from .intersection import CLOSED_MAPPED_SEAMS as CLOSED_MAPPED_SEAMS
from .intersection import adapt_intersection, generate_intersection
from .mesh import Mesh
from .volume import Field, generate

__all__ = [
    "CancellationToken",
    "CancelledError",
    "CompiledField",
    "Field",
    "Mesh",
    "MeshingError",
    "adapt_intersection",
    "compile_field",
    "generate",
    "generate_intersection",
]

if hasattr(_meshers, "generate_surface"):  # ruff: ignore[non-empty-init-module] - optional native exports
    from .surface import SurfaceMesh, generate_surface

    __all__ += ["SurfaceMesh", "generate_surface"]
