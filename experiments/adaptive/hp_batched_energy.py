"""Bounded cell batches for saved six-load h/p cross energies."""

import numpy as np


def integrate_batched(points, cells, h_fields, cubic_cells, p_fields, parents, maps, quadrature, *, batch_size=256, return_stress=False):
    """Return the full cross energy and per-cell p-minus-h Gram matrices."""
    from cubic_reference import register_cubic
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    assert batch_size > 0 and h_fields.shape == (6, len(points), 3)
    assert p_fields.shape[0] == 6 and p_fields.shape[2] == 3
    linear, quadratic, cubic = Tet4(), Tet10(), register_cubic()()
    locations, weights = quadrature
    assert np.isfinite(weights).all() and np.min(weights) > 0
    derivatives = quadratic.shape_function_derivative(locations)
    barycentric = linear.shape_function(locations)
    reference = linear.shape_function_derivative([[.25] * 3])[0]
    energy, minimum = np.zeros((12, 12)), np.inf
    mean_stress = np.zeros((6, 12))
    local_gram, volumes = np.zeros((len(cells), 6, 6)), np.zeros(len(cells))
    for start in range(0, len(cells), batch_size):
        stop = min(start + batch_size, len(cells))
        mapping = maps[start:stop]
        root_locations = np.einsum('pa,cab->cpb', barycentric, mapping) @ linear.xi_nd
        transform = np.einsum('ab,cbn->can', reference, mapping) @ linear.xi_nd
        root_derivatives = cubic.shape_function_derivative(root_locations.reshape(-1, 3)).reshape(stop-start, len(weights), 3, 20)
        p_derivatives = np.einsum('cab,cqbn->cqan', transform, root_derivatives)
        geometry = points[cells[start:stop]]
        h_values = h_fields[:, cells[start:stop]]
        p_values = p_fields[:, cubic_cells[parents[start:stop]]]
        for q, (derivative_h, weight) in enumerate(zip(derivatives, weights, strict=True)):
            jacobian = np.einsum('an,tnb->tab', derivative_h, geometry)
            determinant = np.linalg.det(jacobian)
            assert np.isfinite(determinant).all() and determinant.min() > 0
            minimum = min(minimum, float(determinant.min()))
            inverse = np.linalg.inv(jacobian)
            gradient_h = np.einsum('tba,an->tbn', inverse, derivative_h)
            gradient_p = np.einsum('tba,tan->tbn', inverse, p_derivatives[:, q])
            gradients = np.concatenate((np.einsum('ltni,tjn->ltij', h_values, gradient_h),
                np.einsum('ltni,tjn->ltij', p_values, gradient_p)))
            strain = (gradients + gradients.swapaxes(-1, -2)) / 2
            stress = 2 * strain + np.trace(strain, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)
            weighted = weight * determinant
            energy += np.einsum('t,itkl,jtkl->ij', weighted, strain, stress, optimize=True)
            local_gram[start:stop] += np.einsum('t,itkl,jtkl->tij', weighted,
                strain[6:] - strain[:6], stress[6:] - stress[:6], optimize=True)
            if return_stress:
                mean = np.einsum('t,ltij->lij', weighted, stress)
                mean_stress += mean[:, [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]].T
            volumes[start:stop] += weighted
    result = energy, float(volumes.sum()), minimum, local_gram, volumes
    return result + (mean_stress,) if return_stress else result
