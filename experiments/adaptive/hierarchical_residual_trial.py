"""Same-domain periodic cubic-bubble residual modes; never a native mesh candidate."""

from itertools import combinations

import numpy as np

EDGES = np.array(list(combinations(range(4), 2)))
FACES = np.array(list(combinations(range(4), 3)))


def values(barycentric):
    l = barycentric
    return np.stack([l[..., i]*l[..., j]*(l[..., i]-l[..., j]) for i, j in EDGES]
        + [l[..., i]*l[..., j]*l[..., k] for i, j, k in FACES], axis=-1)


def basis(locations):
    from fedoo.lib_elements.tetrahedron import Tet4, Tet10

    linear = Tet4()
    np.testing.assert_array_equal(linear.xi_nd, Tet10().xi_nd[:4])
    l = linear.shape_function(locations)
    dl = np.asarray(linear.shape_function_derivative(locations))
    derivative = []
    for i, j in EDGES:
        derivative.append((dl[:, :, i]*l[:, j, None]+dl[:, :, j]*l[:, i, None])
            *(l[:, i]-l[:, j])[:, None]+l[:, i, None]*l[:, j, None]*(dl[:, :, i]-dl[:, :, j]))
    for i, j, k in FACES:
        derivative.append(dl[:, :, i]*(l[:, j]*l[:, k])[:, None]
            +dl[:, :, j]*(l[:, i]*l[:, k])[:, None]+dl[:, :, k]*(l[:, i]*l[:, j])[:, None])
    return values(l), np.stack(derivative, axis=-1)


def topology(points, cells, surface, periodic=True):
    if periodic:
        from periodic import periodic_pairs
        _, classes, checked = periodic_pairs(points, surface)
    else:
        classes, checked = np.arange(len(points)), None
    corners = classes[cells[:, :4]]
    assert np.all(np.diff(np.sort(corners, axis=1), axis=1) > 0)
    physical = points[cells[:, :4]]
    assert not periodic or np.linalg.norm(physical[:, EDGES[:, 0]]-physical[:, EDGES[:, 1]], axis=-1).max() < .5
    edge_keys = np.sort(corners[:, EDGES], axis=-1).reshape(-1, 2)
    face_keys = np.sort(corners[:, FACES], axis=-1).reshape(-1, 3)
    edges, edge_ids = np.unique(edge_keys, axis=0, return_inverse=True)
    faces, face_ids, face_counts = np.unique(face_keys, axis=0, return_inverse=True, return_counts=True)
    assert face_counts.max() <= 2
    ids = np.column_stack((edge_ids.reshape(-1, 6), len(edges)+face_ids.reshape(-1, 4)))
    signs = np.column_stack((np.where(corners[:, EDGES[:, 0]] < corners[:, EDGES[:, 1]], 1., -1.),
        np.ones((len(cells), 4))))
    return dict(ids=ids, signs=signs, corner_classes=corners, face_ids=face_ids,
        edges=len(edges), faces=len(faces), features=len(edges)+len(faces), periodic=checked)


def stress(gradient):
    strain = (gradient+gradient.swapaxes(-1, -2))/2
    return strain, 2*strain+np.trace(strain, axis1=-2, axis2=-1)[..., None, None]*np.eye(3)


def build_modes(points, cells, fields, space, quadrature, batch=1024):
    from fedoo.lib_elements.tetrahedron import Tet10

    locations, weights = quadrature
    derivatives = Tet10().shape_function_derivative(locations)
    _, bubbles = basis(locations)
    diagonal = np.zeros((space['features'], 3, 3))
    residual = np.zeros((space['features'], 3, 6))
    minimum = np.inf
    for start in range(0, len(cells), batch):
        subset = slice(start, min(start+batch, len(cells)))
        geometry, old = points[cells[subset]], fields[:, cells[subset]]
        ids, signs = space['ids'][subset], space['signs'][subset]
        for d, db, weight in zip(derivatives, bubbles, weights, strict=True):
            jacobian = np.einsum('an,tnb->tab', d, geometry)
            determinant = np.linalg.det(jacobian)
            minimum = min(minimum, float(determinant.min()))
            assert determinant.min() > 0 and weight > 0
            inverse = np.linalg.inv(jacobian)
            grad = np.einsum('tba,an->tbn', inverse, d)
            bubble_grad = np.einsum('tba,ak->tkb', inverse, db)*signs[:, :, None]
            _, sigma = stress(np.einsum('ltni,tjn->ltij', old, grad))
            weighted = weight*determinant
            blocks = ((bubble_grad*bubble_grad).sum(axis=-1)[:, :, None, None]*np.eye(3)
                +2*np.einsum('tki,tkj->tkij', bubble_grad, bubble_grad))*weighted[:, None, None, None]
            force = np.einsum('ltij,tkj->tkil', sigma, bubble_grad)*weighted[:, None, None, None]
            np.add.at(diagonal, ids.ravel(), blocks.reshape(-1, 3, 3))
            np.add.at(residual, ids.ravel(), force.reshape(-1, 3, 6))
    assert np.isfinite(residual).all() and np.linalg.eigvalsh(diagonal).min() > 0
    coefficients = -np.linalg.solve(diagonal, residual)
    assert np.isfinite(coefficients).all()
    return coefficients, residual, dict(minimum_sampled_jacobian=minimum,
        residual_maximum=float(abs(residual).max()), diagonal_minimum_eigenvalue=float(np.linalg.eigvalsh(diagonal).min()))


def integrate(points, cells, fields, space, coefficients, quadrature, batch=1024):
    from fedoo.lib_elements.tetrahedron import Tet10

    locations, weights = quadrature
    derivatives = Tet10().shape_function_derivative(locations)
    _, bubbles = basis(locations)
    energy = np.zeros((12, 12))
    local_modes = np.zeros((len(cells), 6, 6))
    minimum, volume = np.inf, 0.
    for start in range(0, len(cells), batch):
        subset = slice(start, min(start+batch, len(cells)))
        geometry, old = points[cells[subset]], fields[:, cells[subset]]
        c = coefficients[space['ids'][subset]]
        signs = space['signs'][subset]
        for d, db, weight in zip(derivatives, bubbles, weights, strict=True):
            jacobian = np.einsum('an,tnb->tab', d, geometry)
            determinant = np.linalg.det(jacobian)
            minimum = min(minimum, float(determinant.min()))
            assert determinant.min() > 0 and weight > 0
            inverse = np.linalg.inv(jacobian)
            grad = np.einsum('tba,an->tbn', inverse, d)
            bubble_grad = np.einsum('tba,ak->tkb', inverse, db)*signs[:, :, None]
            parent_gradient = np.einsum('ltni,tjn->ltij', old, grad)
            mode_gradient = np.einsum('tkcm,tkj->mtcj', c, bubble_grad, optimize=True)
            strain, sigma = stress(np.concatenate((parent_gradient, mode_gradient)))
            weighted = weight*determinant
            energy += np.einsum('t,itkl,jtkl->ij', weighted, strain, sigma, optimize=True)
            local_modes[subset] += np.einsum('t,itkl,jtkl->tij', weighted, strain[6:], sigma[6:], optimize=True)
            volume += float(weighted.sum())
    np.testing.assert_allclose(local_modes.sum(axis=0), energy[6:, 6:], atol=1e-12, rtol=1e-10)
    return energy, local_modes, dict(minimum_sampled_jacobian=minimum, volume=volume)


def trace_check(space, coefficients, batch=8192):
    ids, signs, corners = space['ids'], space['signs'], space['corner_classes']
    face_ids = space['face_ids']
    # Six P2 trace nodes are zeros. These four points complete a unisolvent
    # ten-point test for the cubic trace, including every edge and face mode.
    probes = np.array(((2/3,1/3,0.), (2/3,0.,1/3), (0.,2/3,1/3), (1/3,1/3,1/3)))
    all_values = np.empty((len(face_ids), len(probes), 3, 6))
    for start in range(0, len(face_ids), batch):
        records = np.arange(start, min(start+batch, len(face_ids)))
        t, f = records//4, records%4
        local_vertices = FACES[f]
        labels = np.take_along_axis(corners[t], local_vertices, axis=1)
        ordered = np.take_along_axis(local_vertices, np.argsort(labels, axis=1), axis=1)
        barycentric = np.zeros((len(records), len(probes), 4))
        for i in range(3):
            barycentric[np.arange(len(records))[:, None], np.arange(len(probes))[None, :], ordered[:, i, None]] = probes[None, :, i]
        phi = values(barycentric)*signs[t, None, :]
        all_values[records] = np.einsum('rqk,rkcm->rqcm', phi, coefficients[ids[t]], optimize=True)
    first = np.full(space['faces'], len(face_ids), dtype=np.int64)
    np.minimum.at(first, face_ids, np.arange(len(face_ids)))
    maximum = 0.
    for start in range(0, len(face_ids), batch):
        records = np.arange(start, min(start+batch, len(face_ids)))
        maximum = max(maximum, float(abs(all_values[records]-all_values[first[face_ids[records]]]).max()))
    assert maximum <= 1e-12*max(1., float(abs(all_values).max()))
    return dict(shared_and_periodic_face_records=len(face_ids), probes_per_face=len(probes),
        zero_p2_trace_nodes=6, complete_cubic_trace_check=True, maximum_trace_jump=maximum)


def trial(energy):
    from hp_cross_energy import combine

    h, cross, modes = energy[:6, :6], energy[:6, 6:], energy[6:, 6:]
    # Reuse the independently checked h / (h+mode) six-mode energy optimizer.
    equivalent = np.block([[h, h+cross], [h+cross.T, h+cross+cross.T+modes]])
    return combine(equivalent)
