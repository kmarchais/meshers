"""One geometry-only periodic curve-distance grading hypothesis, not an error bound."""

import numpy as np
from scipy.spatial import cKDTree

from lattice_crease_mechanism import opening_angles
from periodic_crease_audit import features, gradient, project_box_curve


def make_metric(mesh, function, minimum, maximum, cells):
    """Use a fresh seed's confirmed re-entrant wall/box curves, with no FE data."""
    feature = features(dict(points=mesh.points, tetrahedra=mesh.tetrahedra, surface=mesh.surface))
    edges = np.asarray(feature['curve_edges'])
    mids = feature['curve_midpoints']
    probes, projection = project_box_curve(function, mesh.points[mids].copy())
    locked = np.isclose(abs(probes), .5, atol=1e-10, rtol=0)
    selected = np.zeros(len(edges), dtype=bool)
    axes = []
    for axis in range(3):
        eligible = (locked.sum(axis=1) == 1) & locked[:, axis]
        ids = np.flatnonzero(eligible)
        if not len(ids):
            continue
        low, high = probes[ids].copy(), probes[ids].copy()
        low[:, axis], high[:, axis] = -.5, .5
        mismatch = float(abs(function(*low.T)-function(*high.T)).max())
        assert mismatch < 1e-9
        estimates, valid = [], np.ones(len(ids), dtype=bool)
        for epsilon in (2e-6, 1e-6):
            first, second = low.copy(), high.copy()
            first[:, axis] += epsilon; second[:, axis] -= epsilon
            angle, good = opening_angles(gradient(function, first, epsilon/8),
                gradient(function, second, epsilon/8), axis)
            estimates.append(angle); valid &= good
        stable = abs(estimates[1]-estimates[0]) <= .1
        keep = valid & stable & (estimates[1] > 190)
        selected[ids[keep]] = True
        axes.append(dict(axis=axis, examined_edges=len(ids), selected_edges=int(keep.sum()),
            scalar_trace_mismatch=mismatch, maximum_angle_change_deg=float(abs(estimates[1]-estimates[0]).max())))
    assert selected.any(), 'No confirmed re-entrant single-plane curves'
    # Sample each quadratic curve more finely than the requested smallest edge.
    clouds = []
    for edge, mid in zip(edges[selected], mids[selected], strict=True):
        a, b, m = mesh.points[[*edge, mid]]
        length_bound = 4*(np.linalg.norm(a-m)+np.linalg.norm(m-b))
        count = max(8, int(np.ceil(length_bound/(minimum/4))))
        t = np.linspace(0, 1, count+1)
        clouds.append((1-t)[:, None]*(1-2*t)[:, None]*a
            + t[:, None]*(2*t-1)[:, None]*b+4*(t*(1-t))[:, None]*m)
    points, cloud_projection = project_box_curve(function, np.concatenate(clouds))
    start, maximum_spacing = 0, 0.
    for segment in clouds:
        length = len(segment)
        maximum_spacing = max(maximum_spacing, float(np.linalg.norm(
            np.diff(points[start:start+length], axis=0), axis=1).max()))
        start += length
    assert maximum_spacing <= minimum/4
    canonical = np.mod(points+.5, 1.)
    cloud = np.unique(canonical, axis=0)
    tree = cKDTree(cloud, boxsize=1.)
    radius = 1/cells
    exponent = .75

    def target(x, y, z):
        coordinates = np.stack(np.broadcast_arrays(x, y, z), axis=-1)
        shape = coordinates.shape[:-1]
        distance = tree.query(np.mod(coordinates.reshape(-1, 3)+.5, 1.), workers=1)[0]
        return np.maximum(minimum, maximum*np.minimum(distance/radius, 1.)**exponent).reshape(shape)

    queries = np.random.default_rng(20261025).uniform(-.5, .5, (4096, 3))
    values = target(*queries.T)
    assert np.isfinite(values).all() and (values >= minimum).all() and (values <= maximum).all()
    for axis in range(3):
        low, high = queries.copy(), queries.copy()
        low[:, axis], high[:, axis] = -.5, .5
        np.testing.assert_array_equal(target(*low.T), target(*high.T))
    np.testing.assert_allclose(target(*(cloud-.5).T), minimum, atol=1e-15, rtol=0)
    metadata = dict(law=dict(kind='periodic-reentrant-curve-distance', exponent=exponent,
        minimum_length=minimum, maximum_length=maximum, transition_distance=radius,
        curve_sample_spacing_target=minimum/4, minimum_reentrant_angle_deg=190),
        examined_curve_edges=len(edges), selected_curve_edges=int(selected.sum()),
        periodic_curve_samples=len(cloud), axes=axes, probe_projection=projection,
        measured_maximum_curve_sample_spacing=maximum_spacing,
        cloud_projection=cloud_projection, periodic_verification_points_per_axis=len(queries),
        fe_fields_used=False, saved_seed_replayed=False,
        scope='Geometry-only grading hypothesis from freshly generated curves; no 3D singularity theorem or accuracy guarantee')
    return target, metadata, cloud


def directional_metric(scalar, cloud, function, maximum=.1, ratio=5.):
    """Keep the fixed normal target and allow longer edges along each curve."""
    frame = np.zeros_like(cloud)
    physical = cloud-.5
    planes = np.isclose(abs(physical), .5, atol=1e-10, rtol=0)
    for axis in range(3):
        ids = np.flatnonzero((planes.sum(axis=1) == 1) & planes[:, axis])
        probes = physical[ids].copy()
        probes[:, axis] -= np.sign(probes[:, axis])*1e-6
        normal = gradient(function, probes, 1e-6/8)
        tangent = np.cross(np.eye(3)[axis], normal)
        tangent /= np.linalg.norm(tangent, axis=1)[:, None]
        frame[ids] = tangent
    tree = cKDTree(cloud, boxsize=1.)

    def evaluate(points):
        nearest = tree.query(np.mod(points+.5, 1.), workers=1)[1]
        tangent = frame[nearest]
        normal_size = scalar(*points.T)
        along_size = np.minimum(maximum, ratio*normal_size)
        matrix = np.eye(3)[None]/normal_size[:, None, None]**2
        return matrix+(along_size**-2-normal_size**-2)[:, None, None]*np.einsum('ni,nj->nij', tangent, tangent)

    return evaluate
