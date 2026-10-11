import numpy as np
import pytest

from microgen_benchmark.cases import CASES
from microgen_benchmark.geometry import canonical_mesh, overlapping_pair, validate


def tet():
    return np.array([[0.0, 0, 0], [0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.3]])


def test_overlap_and_contact_are_distinct():
    a = tet()
    assert overlapping_pair(np.array([a, a + 0.02]))[0] == [0, 1]
    assert overlapping_pair(np.array([a, a + 0.3]))[0] is None


def test_incomplete_validation_cannot_pass():
    a = tet()
    with pytest.raises(RuntimeError, match="incomplete"):
        overlapping_pair(np.array([a, a + 0.01, a + 0.02]), max_pairs=1)


def test_duplicate_cells_rejected():
    p, t, _, _ = canonical_mesh(tet(), [[0, 1, 2, 3], [0, 1, 2, 3]])
    with pytest.raises(ValueError, match="Duplicate"):
        validate(p, t, CASES[0])


def test_outside_specimen_rejected():
    p, t, _, _ = canonical_mesh(tet() + 1, [[0, 1, 2, 3]])
    with pytest.raises(ValueError, match="specimen"):
        validate(p, t, CASES[0])


def test_vertex_only_material_contact_rejected():
    a = tet()
    points = np.vstack([a, -a[1:]])
    p, t, _, _ = canonical_mesh(points, [[0, 1, 2, 3], [0, 4, 5, 6]])
    with pytest.raises(ValueError, match="touch only"):
        validate(p, t, CASES[0])


def test_batched_overlap_matches_independent_face_and_edge_axes():
    # Vary scale and translation; include all 36 edge axes in the independent oracle.
    from microgen_benchmark.geometry import EDGES, FACES

    rng = np.random.default_rng(197)
    for _ in range(60):
        a, b = rng.normal(size=(2, 4, 3))
        ea = np.array([a[j] - a[i] for i, j in EDGES])
        eb = np.array([b[j] - b[i] for i, j in EDGES])
        axes = np.concatenate(
            [
                np.cross(a[FACES[:, 1]] - a[FACES[:, 0]], a[FACES[:, 2]] - a[FACES[:, 0]]),
                np.cross(b[FACES[:, 1]] - b[FACES[:, 0]], b[FACES[:, 2]] - b[FACES[:, 0]]),
                np.cross(ea[:, None], eb[None]).reshape(-1, 3),
            ]
        )
        axes = axes[np.linalg.norm(axes, axis=1) > 1e-15]
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        pa, pb = axes @ a.T, axes @ b.T
        hit = np.all(np.minimum(pa.max(1), pb.max(1)) - np.maximum(pa.min(1), pb.min(1)) > 1e-10)
        assert (overlapping_pair(np.array([a, b]))[0] is not None) == hit


def test_saved_vtk_gyroid_adjacent_overlap_is_rejected():
    # Two cells minimized from the n16 VTK failure, retaining the shared-face defect.
    points = np.array(
        [
            [-0.375, 0.4348050761781621, -0.5],
            [-0.3125, 0.4111906857729064, -0.5],
            [-0.3125, 0.375, -0.4375],
            [-0.375, 0.39330582617584076, -0.4375],
            [-0.3125, 0.375, -0.5],
        ]
    )
    p, t, _, _ = canonical_mesh(points, [[0, 1, 2, 3], [4, 0, 2, 1]])
    from microgen_benchmark.geometry import boundary

    with pytest.raises(ValueError, match="same side"):
        boundary(p, t)
    assert overlapping_pair(p[t])[0] == [0, 1]
