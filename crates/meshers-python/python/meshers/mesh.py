"""Array-backed tetrahedral meshes and their boundary metadata."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class Mesh:
    """Owned arrays and diagnostics for a tetrahedral mesh.

    Attributes:
        points: Float64 coordinates, shape ``(N, 3)``.
        tetrahedra: Int64 connectivity, shape ``(M, 4)`` or ``(M, 10)``.
        surface: Boundary triangles, shape ``(K, 3)`` or ``(K, 6)``.
        boundary_tags: Zero for implicit faces, 1 and 2 for x-/x+, then y and z.
        periodic_pairs: Three arrays of paired point IDs, one for each box axis.
            For scalar fields, every material-cap node is matched. Existing
            edge/corner images form complete constraint classes. A wall-only
            contact can lack an opposite image when that image has no
            tetrahedral material support.
        diagnostics: Integrated volume, corner-based MMG quality, timing,
            sampled error, evaluator and worker count. For Tet10, quadratic
            diagnostics also retain the linear parent's volume.
        parameters: Parameter coordinates for mapped intersections, otherwise None.
        constraint_masks: Bit j marks membership in constraint j, otherwise None.
        constraint_names: Names corresponding to the constraint-mask bits.
        periodic_transforms: Three low-to-high rigid 4x4 transforms for intersections.
            For mapped meshes, use these to rotate vector degrees of freedom; the
            exported PeriodicShift alone only describes point displacement.

    Arrays remain mutable even though the container is frozen. Editing coordinates
    or connectivity can invalidate quality values and periodic metadata.
    """

    points: npt.NDArray[np.float64]
    tetrahedra: npt.NDArray[np.int64]
    surface: npt.NDArray[np.int64]
    boundary_tags: npt.NDArray[np.uint8]
    periodic_pairs: tuple[npt.NDArray[np.int64], ...]
    diagnostics: dict
    parameters: npt.NDArray[np.float64] | None = None
    constraint_masks: npt.NDArray[np.uint8] | npt.NDArray[np.uint32] | None = None
    constraint_names: tuple[str, ...] = ()
    periodic_transforms: npt.NDArray[np.float64] | None = None

    def surface_for(self, name: str) -> npt.NDArray[np.int64]:
        """Return boundary triangles assigned to a named intersection constraint.

        Returns:
            Triangle connectivity referencing ``points``.

        Raises:
            ValueError: The mesh has no constraint with this name.
        """
        if name not in self.constraint_names or self.constraint_masks is None:
            raise ValueError(f"Unknown constraint: {name}")
        bit = 1 << self.constraint_names.index(name)
        masks = self.constraint_masks[self.surface]
        labels = masks[:, 0] & masks[:, 1] & masks[:, 2]
        return self.surface[(labels & bit) != 0]

    @property
    def feature_edges(self) -> npt.NDArray[np.int64]:
        """Boundary edges shared by at least two named constraints."""
        if self.constraint_masks is None:
            return np.empty((0, 2), dtype=np.int64)
        edges = np.unique(
            np.sort(
                np.concatenate([
                    self.surface[:, [0, 1]],
                    self.surface[:, [1, 2]],
                    self.surface[:, [2, 0]],
                ]),
                axis=1,
            ),
            axis=0,
        )
        common = self.constraint_masks[edges[:, 0]] & self.constraint_masks[edges[:, 1]]
        return edges[(common & (common - np.uint8(1))) != 0]

    def write_vtkhdf(self, path: str | Path, *, surface: bool = False) -> None:
        """Write a VTKHDF unstructured grid, with physical periodic shifts. Requires meshers[io].

        Quadratic cell Volume integrates the current curved mapping.
        CornerVolume and CornerMMGQuality identify the linear parent measures.
        Linear cells retain the Volume and MMGQuality arrays.
        """
        from ._vtkhdf import write_grid

        write_grid(self, path, surface=surface)
