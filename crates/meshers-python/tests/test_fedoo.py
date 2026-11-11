import meshers
import numpy as np
import pytest

fd = pytest.importorskip("fedoo")


@pytest.mark.parametrize("order", [1, 2])
def test_fedoo_preserves_nodes_order_and_owns_arrays(order):
    mesh = meshers.generate(
        lambda x, y, z: x * 0 - 1,
        bounds=(-0.5, 0.5) * 3,
        cells=4,
        order=order,
        threads=1,
    )
    converted = mesh.to_fedoo()
    assert converted.elm_type == ("tet4" if order == 1 else "tet10")
    np.testing.assert_array_equal(converted.nodes, mesh.points)
    np.testing.assert_array_equal(converted.elements, mesh.tetrahedra)
    assert not np.shares_memory(converted.nodes, mesh.points)
    assert not np.shares_memory(converted.elements, mesh.tetrahedra)
    element = fd.lib_elements.element_list.get_element(converted.elm_type)()
    np.testing.assert_allclose(
        element.shape_function(element.xi_nd), np.eye(mesh.tetrahedra.shape[1]), atol=1e-14
    )


def test_fedoo_rejects_unknown_connectivity():
    mesh = meshers.Mesh(
        points=np.zeros((8, 3)),
        tetrahedra=np.arange(8)[None],
        surface=np.empty((0, 3), dtype=int),
        boundary_tags=np.empty(0, dtype=np.uint8),
        periodic_pairs=(),
        diagnostics={},
    )
    with pytest.raises(ValueError, match="Tet4 or Tet10"):
        mesh.to_fedoo()
