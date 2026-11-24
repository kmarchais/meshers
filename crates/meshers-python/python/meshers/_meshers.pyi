# Copyright (c) 2026 Kevin Marchais

from collections.abc import Sequence
from typing import NotRequired, TypedDict

import numpy as np
from numpy.typing import NDArray

from .volume import Field

class _GeometryArrays(TypedDict):
    points: NDArray[np.float64]
    tetrahedra: NDArray[np.int64]
    surface: NDArray[np.int64]

class _PeriodicArrays(_GeometryArrays):
    periodic_pairs: list[NDArray[np.int64]]

class _QuadraticResult(_PeriodicArrays):
    linear_volume: float
    curved_edges: int
    straightened_edges: int
    partially_curved_edges: int
    minimum_jacobian_ratio: float
    fit: NotRequired[dict[str, int | float]]

class _VolumeResult(_PeriodicArrays):
    boundary_tags: NDArray[np.uint8]
    diagnostics: dict[str, object]

class _GeneratedVolume(_VolumeResult):
    quadratic: NotRequired[_QuadraticResult]

class _AdaptedIntersection(_VolumeResult):
    constraint_masks: NDArray[np.uint8] | NDArray[np.uint32]

class _IntersectionResult(_GeometryArrays):
    parameters: NDArray[np.float64]
    constraint_masks: NDArray[np.uint8] | NDArray[np.uint32]
    diagnostics: dict[str, object]

class CancellationToken:
    """Thread-safe cooperative cancellation. Tokens remain cancelled after use."""
    def __init__(self) -> None: ...
    @property
    def cancelled(self) -> bool:
        """Whether cancellation has been requested."""
    def cancel(self) -> None:
        """Request interruption at the next mesher cancellation checkpoint."""

class MeshingError(RuntimeError):
    """Meshing failed a topology, budget, field or geometry check."""

class CancelledError(RuntimeError):
    """The operation was interrupted by its cancellation token."""

def generate(
    field: str | Field | CompiledField,
    gradient: Field | None,
    bounds: Sequence[Sequence[float]],
    cells: Sequence[int],
    periodic: Sequence[bool],
    band: Sequence[float] | None,
    tolerance: float,
    minimum_quality: float,
    budget: int,
    passes: int,
    interior_refine_nodes: int,
    snap: float,
    threads: int,
    batch_size: int,
    cancel: CancellationToken | None,
    adapt: tuple[float, float, float, float, int, float, float] | None = None,
    order: int = 1,
    size_field: str | Field | CompiledField | None = None,
    halfspaces: Sequence[Sequence[float]] | None = None,
    regions: Sequence[int] | None = None,
    constraints: Sequence[str | Field | CompiledField] | None = None,
    curved_fit: bool = False,
    *,
    lattice_fallback: bool = True,
    defeaturing: bool = True,
) -> _GeneratedVolume: ...

class CompiledField:
    """Immutable float64 evaluator and optional gradient, owning compiled code.

    Create with meshers.compile_field, or let meshers.generate prepare it automatically.
    """
    @property
    def has_gradient(self) -> bool:
        """Whether this evaluator contains a compiled gradient."""
    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate one physical point."""
    @property
    def node_count(self) -> int:
        """Number of interned nodes in the compiled expression graph."""
    def gradient(self, x: float, y: float, z: float) -> list[float]:
        """Return the three first derivatives, or raise ValueError if unavailable."""
    def values(self, points: NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluate an (n, 3) float64 point array without a Python call per point."""
    @property
    def halfspaces(self) -> list[list[float]] | None:
        """Affine planes whose nonpositive sides define this convex field, if recognized."""
    @property
    def polyhedral_regions(self) -> tuple[list[list[float]], list[int]] | None:
        """Planes and region masks of a recognized Boolean polyhedral field."""
    @property
    def sdf_intersection(self) -> bool:
        """Whether the retained branches require a signed-distance consistency check."""
    def intersection_fields(self) -> list[CompiledField] | None:
        """Return cached fields for an intersection, or None for an unsupported expression."""

def compile_graph(
    nodes: Sequence[tuple[str, Sequence[int], float]],
    outputs: Sequence[int],
    constraint_outputs: Sequence[Sequence[int]] | None = None,
    *,
    retain_sdf: bool = False,
) -> CompiledField: ...

MAX_CONSTRAINTS: int
ADAPTIVE_POLYHEDRAL: bool
PERIODIC_POLYHEDRAL: bool
POLYHEDRAL_REGIONS: bool
RETAINED_INTERSECTIONS: bool
SDF_INTERSECTIONS: bool
ADAPTIVE_RETAINED_INTERSECTIONS: bool

def available_threads() -> int: ...
def face_components(tetrahedra: NDArray[np.int64]) -> int: ...
def quadratic_volumes(
    points: NDArray[np.float64],
    tetrahedra: NDArray[np.int64],
) -> NDArray[np.float64]: ...
def generate_intersection(
    fields: Sequence[str | Field | CompiledField],
    mapping: Sequence[str | Field | CompiledField] | None,
    bounds: Sequence[Sequence[float]],
    cells: Sequence[int],
    periodic: Sequence[bool],
    closed: Sequence[bool],
    minimum_quality: float,
    budget: int,
    passes: int,
    snap: float,
    batch_size: int,
    cancel: CancellationToken | None,
) -> _IntersectionResult: ...
def _periodic_affine_traces(
    planes: Sequence[Sequence[float]],
    bounds: Sequence[Sequence[float]],
    periodic: Sequence[bool],
    regions: Sequence[int] | None = None,
) -> bool: ...
def adapt_intersection(
    fields: Sequence[str | Field | CompiledField],
    points: NDArray[np.float64],
    tetrahedra: NDArray[np.int64],
    masks: NDArray[np.uint32],
    bounds: Sequence[Sequence[float]],
    tolerance: float,
    minimum_quality: float,
    budget: int,
    passes: int,
    adapt: tuple[float, float, float, float, int, float, float],
    seed_tolerance: float,
    batch_size: int,
    cancel: CancellationToken | None,
) -> _AdaptedIntersection: ...
