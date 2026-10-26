"""Small saved-field traction diagnostic; no equilibrated estimator or error bound."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np

import microgen_workflow as w
from nested_refinement import FACES


def side(points, cells, fields, owners, face, bary):
    from fedoo.lib_elements.tetrahedron import Tet4, Tet10
    reference = Tet4().xi_nd
    local = (cells[owners, :4, None] == face[:, None, :]).argmax(axis=1)
    assert np.array_equal(np.take_along_axis(cells[owners, :4], local, axis=1), face)
    xi = np.einsum('a,fad->fd', bary, reference[local])
    element = Tet10()
    derivative = element.shape_function_derivative(xi)
    geometry = points[cells[owners]]
    jacobian = np.einsum('fan,fnb->fab', derivative, geometry)
    assert np.linalg.det(jacobian).min() > 0
    gradients = np.einsum('fba,fan->fbn', np.linalg.inv(jacobian), derivative)
    gradient = np.einsum('lfni,fjn->flij', fields[:, cells[owners]], gradients)
    strain = (gradient + gradient.swapaxes(-1, -2)) / 2
    stress = 2 * strain + np.trace(strain, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)
    tangents = np.einsum('fqa,fab->fqb', reference[local[:, 1:]] - reference[local[:, :1]], jacobian)
    normal = np.cross(tangents[:, 0], tangents[:, 1])
    density = np.linalg.norm(normal, axis=1)
    assert density.min() > 0
    probe = np.einsum('fn,fnd->fd', element.shape_function(xi), geometry)
    return stress, normal / density[:, None], density, probe


def integrate(points, cells, fields, owner_a, face_a, owner_b=None, face_b=None, *, order=3):
    """h_F integral of six-load traction Gram matrices on canonical curved faces."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    matrix = np.zeros((len(face_a), 6, 6))
    area = np.zeros(len(face_a))
    h = np.max(np.linalg.norm(points[face_a[:, [1, 2, 0]]] - points[face_a], axis=2), axis=1)
    maximum_mapping_error = 0.
    for r, wr in zip(nodes, weights, strict=True):
        for s, ws in zip(nodes, weights, strict=True):
            bary = np.array([(1-r)*(1-s), r, (1-r)*s])
            stress, normal, density, probe = side(points, cells, fields, owner_a, face_a, bary)
            if owner_b is not None:
                other, other_normal, other_density, other_probe = side(points, cells, fields, owner_b, face_b, bary)
                translation = points[face_b[:, 0]] - points[face_a[:, 0]]
                error = float(abs(other_probe - probe - translation).max())
                maximum_mapping_error = max(maximum_mapping_error, error)
                assert error < 1e-10
                np.testing.assert_allclose(other_density, density, rtol=1e-9, atol=1e-12)
                np.testing.assert_allclose(other_normal, normal, rtol=1e-9, atol=1e-10)
                stress = stress - other
            traction = np.einsum('flij,fj->fli', stress, normal)
            measure = wr * ws * (1-r) * density
            matrix += np.einsum('f,fli,fmi->flm', h * measure, traction, traction)
            area += measure
    return matrix, area, maximum_mapping_error


def face_groups(points, cells, surface):
    from periodic import periodic_pairs
    from periodic_crease_audit import features
    raw = np.sort(cells[:, FACES].reshape(-1, 3), axis=1)
    order = np.lexsort(raw.T[::-1])
    sorted_faces = raw[order]
    starts = np.flatnonzero(np.r_[True, np.any(sorted_faces[1:] != sorted_faces[:-1], axis=1)])
    counts = np.diff(np.r_[starts, len(raw)])
    assert np.isin(counts, [1, 2]).all()
    keys = sorted_faces[starts]
    first = order[starts] // 4
    shared = counts == 2
    second = order[starts[shared] + 1] // 4
    assert {tuple(f) for f in keys[~shared]} == {tuple(f) for f in np.sort(surface[:, :3], axis=1)}
    feature = features(dict(points=points, tetrahedra=cells, surface=surface))
    curve = feature['groups']['curve_vertex_incident']
    bulk = feature['groups']['neither_wall_nor_cap_vertex']
    a, b, face = first[shared], second, keys[shared]
    curve_pair = curve[a] | curve[b]
    bulk_pair = bulk[a] & bulk[b]
    groups = {}
    for name, mask in dict(interior_curve=curve_pair, interior_bulk=bulk_pair,
                           interior_other_boundary=~curve_pair & ~bulk_pair).items():
        groups[name] = (a[mask], face[mask], b[mask], face[mask])
    a, face = first[~shared], keys[~shared]
    corners = points[face]
    cap = np.any(np.all(w.on_box_plane(corners, -.5), axis=1) |
                 np.all(w.on_box_plane(corners, .5), axis=1), axis=1)
    for name, mask in dict(wall_curve=~cap & curve[a], wall_other=~cap & ~curve[a]).items():
        groups[name] = (a[mask], face[mask], None, None)
    _, classes, _ = periodic_pairs(points, surface)
    cap_a, cap_face = a[cap], face[cap]
    cap_face = np.take_along_axis(cap_face, np.argsort(classes[cap_face], axis=1), axis=1)
    axis = np.argmax(np.all(w.on_box_plane(points[cap_face], -.5), axis=1) |
                     np.all(w.on_box_plane(points[cap_face], .5), axis=1), axis=1)
    keys = np.column_stack((axis, classes[cap_face]))
    ordered = np.lexsort(keys.T[::-1])
    assert len(ordered) % 2 == 0
    low, high = ordered[::2], ordered[1::2]
    np.testing.assert_array_equal(keys[low], keys[high])
    translation = points[cap_face[high]] - points[cap_face[low]]
    np.testing.assert_allclose(translation, np.broadcast_to(translation[:, :1], translation.shape), atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(translation[:, 0], axis=1), 1., atol=1e-12)
    groups['periodic_cap'] = (cap_a[low], cap_face[low], cap_a[high], cap_face[high])
    return groups


def controls():
    from fedoo.lib_elements.tetrahedron import Tet4, Tet10
    reference = Tet4().xi_nd
    maps = np.array([reference, reference.copy()])
    maps[1, 3, 0] = -1
    maps[1, [0, 1]] = maps[1, [1, 0]]
    points = (Tet4().shape_function(Tet10().xi_nd)[None] @ maps).reshape(-1, 3)
    cells = np.arange(20).reshape(2, 10)
    a, b = np.array([0]), np.array([1])
    fa, fb = np.array([[0, 1, 2]]), np.array([[11, 10, 12]])
    rng = np.random.default_rng(19)
    strains = rng.normal(size=(6, 3, 3)); strains = (strains + strains.transpose(0, 2, 1)) / 2
    fields = np.einsum('lij,nj->lni', strains, points)
    jump, area, _ = integrate(points, cells, fields, a, fa, b, fb, order=4)
    assert abs(jump).max() < 1e-24
    np.testing.assert_allclose(area, .5, atol=1e-14)
    # Exact nonconstant P2 displacement has continuous stress on both affine maps.
    quadratic = fields + rng.normal(size=(6, 1, 3)) * points[None, :, :1]**2
    jump, _, _ = integrate(points, cells, quadratic, a, fa, b, fb, order=4)
    assert abs(jump).max() < 1e-24
    # A continuous kink has a genuine traction jump despite matching interface values.
    kink = fields + np.abs(points[None, :, :1]) * np.eye(6, 3)[:, None]
    jump, _, _ = integrate(points, cells, kink, a, fa, b, fb, order=4)
    assert np.trace(jump[0]) > 1.
    reverse, _, _ = integrate(points, cells, kink, b, fb[:, ::-1], a, fa[:, ::-1], order=4)
    np.testing.assert_allclose(jump, reverse, atol=1e-12)
    # A uniaxial tangential stress has zero traction on the x-normal face.
    strain = np.diag([-1/15, 4/15, -1/15])
    tangential = np.broadcast_to(points @ strain.T, fields.shape).copy()
    wall, _, _ = integrate(points, cells, tangential, a, fa, order=4)
    assert abs(wall).max() < 1e-25
    active, _, _ = integrate(points, cells, fields, a, fa, order=4)
    assert np.trace(active[0]) > 1.
    return dict(affine_jump_zero=True, quadratic_jump_zero=True, kink_detected=True,
                orientation_invariant=True, tangential_wall_zero=True, nonzero_wall_detected=True)


def run(args):
    w.imports()
    result = dict(controls=controls(),seed=args.seed,sample_cap_per_group=args.samples,
                  surface_rules=[9,16],routes={},inputs={},new_solve=False,
                  reference_promoted=False,full_objective_complete=False,
                  scope='Stratified saved-field traction-jump and free-wall residual diagnostic. Not an equilibrated stress, stiffness correction, error bound or fresh workflow timing.')
    paths = [Path(__file__), w.HERE/'nested_refinement.py', w.HERE/'periodic_crease_audit.py',
             w.HERE/'microgen_workflow.py', w.ROOT/'handoff/macbook/traction-probe-selection-20261026.json']
    if not args.controls_only:
        evidence = json.loads((w.ROOT/'handoff/macbook/comparison-findings-20261026.json').read_text())
        for name, expected in evidence['lattice_stiffness']['inputs'].items():
            assert w.digest(Path(name)) == expected
        result['retained_models'] = 24
        for label, name in [('native','macbook-lattice-relocation-solve-v1'),
                            ('actual','macbook-lattice-actual-live-solve-v1')]:
            source = w.HERE/'runs'/name
            row = json.loads((source/'result.json').read_text())
            assert row['qualified'] and row['completed_loads'] == 6
            paths.extend(source/n for n in ('result.json','mesh.npz','fields.npz'))
            with np.load(source/'mesh.npz') as data:
                points,cells,surface = (data[k] for k in ('points','tetrahedra','surface'))
            fields = np.load(source/'fields.npz')['displacement']
            assert fields.shape == (6,len(points),3) and len(points)<=250000 and len(cells)<=250000
            groups = face_groups(points,cells,surface)
            inverse = np.linalg.inv(row['stiffness'])
            rng = np.random.default_rng(args.seed)
            rows = []
            for name, arrays in groups.items():
                population = len(arrays[0]);assert population
                selected = np.sort(rng.choice(population,min(args.samples,population),replace=False))
                arguments = [None if value is None else value[selected] for value in arrays]
                rules = []
                archive = dict(selected=selected,**{f'input_{i}':value for i,value in enumerate(arguments) if value is not None})
                for order in (3,4):
                    matrix, area, error = integrate(points,cells,fields,*arguments,order=order)
                    scalar = np.einsum('ij,fji->f',inverse,matrix)/6
                    assert scalar.min() >= -1e-12
                    archive[f'gram_{order}'] = matrix;archive[f'area_{order}'] = area
                    estimate = matrix.sum(axis=0)*population/len(selected)
                    rules.append(dict(points=order**2,normalized_trace=float(np.trace(inverse@estimate)/6),
                        estimated_gram=estimate.tolist(),mapping_error=error,
                        sampled_area=float(area.sum()),sampled_trace_quantiles=np.quantile(scalar,[0,.5,.9,.99,1]).tolist()))
                path = args.output.parent/f'{label}-{name}.npz';np.savez_compressed(path,**archive);paths.append(path)
                change = float(np.linalg.norm(archive['gram_4']-archive['gram_3'])/max(np.linalg.norm(archive['gram_4']),1e-30))
                rows.append(dict(group=name,population=population,sampled=len(selected),rules=rules,
                                 local_gram_relative_quadrature_change=change))
            result['routes'][label] = dict(source=str(source),nodes=len(points),elements=len(cells),groups=rows)
    for i,path in enumerate(paths):
        result['inputs'][str(path)] = w.digest(path)
        if path.suffix=='.py':shutil.copyfile(path,args.output.parent/f'input-{i}-{path.name}')
    w.save(args.output,result)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--controls-only',action='store_true');parser.add_argument('--samples',type=int,default=512)
    parser.add_argument('--seed',type=int,default=20261026);args=parser.parse_args()
    assert not args.output.exists() and 1<=args.samples<=2048
    run(args)
