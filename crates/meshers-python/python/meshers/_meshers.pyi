from typing import Any

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
    field: Any,
    gradient: Any,
    bounds: Any,
    cells: Any,
    periodic: Any,
    band: Any,
    tolerance: float,
    minimum_quality: float,
    budget: int,
    passes: int,
    snap: float,
    threads: int,
    batch_size: int,
    cancel: CancellationToken | None,
) -> dict[str, Any]: ...

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

def _compile_expression(
    nodes: list[tuple[str, list[int], float]], outputs: list[int]
) -> CompiledField: ...
def available_threads() -> int: ...
def generate_intersection(
    fields: Any,
    mapping: Any,
    bounds: Any,
    cells: Any,
    periodic: Any,
    minimum_quality: float,
    budget: int,
    passes: int,
    snap: float,
    batch_size: int,
    cancel: CancellationToken | None,
) -> dict[str, Any]: ...
