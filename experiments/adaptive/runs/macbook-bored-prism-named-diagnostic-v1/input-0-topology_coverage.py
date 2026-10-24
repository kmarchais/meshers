"""Bounded analytic geometry controls for failure coverage, without FE solves."""

import argparse
from pathlib import Path
import shutil
import time
import traceback

import microgen_workflow as w


def boundary_topology(surface, np):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    edges = np.unique(np.sort(surface[:, ((0, 1), (1, 2), (2, 0))].reshape(-1, 2), axis=1), axis=0)
    vertices, inverse = np.unique(edges, return_inverse=True)
    inverse = inverse.ravel()
    graph = coo_matrix((np.ones(inverse.size//2), (inverse[::2], inverse[1::2])),
        shape=(len(vertices), len(vertices))).tocsr()
    return connected_components(graph, directed=False, return_labels=False), len(vertices)-len(edges)+len(surface)


def fixture(name, np):
    angles = 2*np.pi*np.arange(96)/96
    if name == 'bored-prism':
        half = (.35, .35, .3)
        center, radius = np.array((.03, -.02)), .12
        def field(x, y, z):
            box = np.maximum(np.maximum(np.abs(x)-half[0], np.abs(y)-half[1]), np.abs(z)-half[2])
            return np.maximum(box, radius-np.sqrt((x-center[0])**2+(y-center[1])**2))
        grid = np.linspace(-1., 1., 65)
        u, v = np.meshgrid(grid, grid, indexing='ij')
        sides = []
        for axis in (0, 1):
            other = 1-axis
            for sign in (-1., 1.):
                p = np.zeros((u.size, 3))
                p[:, axis] = sign*half[axis]
                p[:, other] = half[other]*u.ravel()
                p[:, 2] = half[2]*v.ravel()
                sides.append(p)
        xy = np.column_stack((half[0]*u.ravel(), half[1]*v.ravel()))
        xy = xy[np.linalg.norm(xy-center, axis=1) >= radius]
        for sign in (-1., 1.):
            sides.append(np.column_stack((xy, np.full(len(xy), sign*half[2]))))
        theta, z = np.meshgrid(2*np.pi*np.arange(256)/256, half[2]*grid, indexing='ij')
        sides.append(np.column_stack((center[0]+radius*np.cos(theta.ravel()),
            center[1]+radius*np.sin(theta.ravel()), z.ravel())))
        return field, np.concatenate(sides), (8*np.prod(half)-np.pi*radius**2*2*half[2]), 1, 0
    if name in ('octahedron', 'notched-prism'):
        vertices, faces = octahedron(np) if name == 'octahedron' else notched_prism(np)
        bary = np.array([(i/48, j/48, 1-(i+j)/48) for i in range(49) for j in range(49-i)])
        queries = np.einsum('sn,fnd->fsd', bary, vertices[faces]).reshape(-1, 3)
        if name == 'octahedron':
            def field(x, y, z):
                return (np.abs(x)+np.abs(y)+np.abs(z)-.35)/np.sqrt(3)
            return field, queries, 4*.35**3/3, 1, 2
        def field(x, y, z):
            return np.maximum(np.maximum(np.maximum(np.abs(x)-.35, np.abs(y)-.35),
                np.abs(z)-.3), np.minimum(x-.05, y-.05))
        return field, queries, (.7*.7-.3*.3)*.6, 1, 2
    if name == 'torus':
        major, minor = .25, .08
        def field(x, y, z):
            return np.sqrt((np.sqrt(x*x+y*y)-major)**2+z*z)-minor
        u, v = np.meshgrid(angles, angles, indexing='ij')
        queries = np.stack(((major+minor*np.cos(v))*np.cos(u),
            (major+minor*np.cos(v))*np.sin(u), minor*np.sin(v)), axis=-1).reshape(-1, 3)
        return field, queries, 2*np.pi**2*major*minor**2, 1, 0
    indices = np.arange(8192)
    z = 1-2*(indices+.5)/len(indices)
    angle = indices*np.pi*(3-np.sqrt(5))
    sphere = np.column_stack((np.sqrt(1-z*z)*np.cos(angle), np.sqrt(1-z*z)*np.sin(angle), z))
    if name == 'shell':
        inner, outer = .18, .32
        def field(x, y, z):
            radius = np.sqrt(x*x+y*y+z*z)
            return np.maximum(radius-outer, inner-radius)
        return field, np.concatenate((inner*sphere, outer*sphere)), 4*np.pi/3*(outer**3-inner**3), 2, 4
    assert name == 'disjoint'
    radius, offset = .15, .23
    def field(x, y, z):
        return np.minimum(np.sqrt((x-offset)**2+y*y+z*z)-radius,
            np.sqrt((x+offset)**2+y*y+z*z)-radius)
    return field, np.concatenate((radius*sphere+[offset, 0, 0], radius*sphere-[offset, 0, 0])), 8*np.pi/3*radius**3, 2, 4


def octahedron(np):
    from itertools import product

    vertices = .35*np.concatenate((np.eye(3), -np.eye(3)))
    faces = []
    for signs in product((1, -1), repeat=3):
        face = [axis+(3 if sign == -1 else 0) for axis, sign in enumerate(signs)]
        if np.prod(signs) < 0:
            face[1], face[2] = face[2], face[1]
        faces.append(face)
    return vertices, np.array(faces)


def notched_prism(np):
    """Exact oriented boundary of an L prism, independent of implicit meshing."""
    polygon = np.array(((-.35, -.35), (.35, -.35), (.35, .05),
        (.05, .05), (.05, .35), (-.35, .35)))
    vertices = np.concatenate((np.column_stack((polygon, np.full(6, -.3))),
        np.column_stack((polygon, np.full(6, .3)))))
    triangles = ((0, 1, 3), (1, 2, 3), (0, 3, 5), (3, 4, 5))
    faces = [list(reversed(t)) for t in triangles]+[[v+6 for v in t] for t in triangles]
    for a in range(6):
        b = (a+1)%6
        faces.extend(((a, b, b+6), (a, b+6, a+6)))
    return vertices, np.array(faces)


def bored_prism_distance(points, np):
    """Exact distance to four rectangles, two perforated caps and a finite cylinder."""
    half = np.array((.35, .35, .3))
    center, radius = np.array((.03, -.02)), .12
    result = np.full(len(points), np.inf)
    for axis in (0, 1):
        for sign in (-1., 1.):
            closest = np.clip(points, -half, half)
            closest[:, axis] = sign*half[axis]
            result = np.minimum(result, np.linalg.norm(points-closest, axis=1))
    rho = np.linalg.norm(points[:, :2]-center, axis=1)
    result = np.minimum(result, np.hypot(rho-radius, np.maximum(np.abs(points[:, 2])-half[2], 0.)))
    closest_xy = np.clip(points[:, :2], -half[:2], half[:2])
    in_hole = np.linalg.norm(closest_xy-center, axis=1) < radius
    cap_xy_distance = np.where(in_hole, np.abs(rho-radius), np.linalg.norm(points[:, :2]-closest_xy, axis=1))
    for sign in (-1., 1.):
        result = np.minimum(result, np.hypot(cap_xy_distance, points[:, 2]-sign*half[2]))
    return result


def bored_prism_witness(point, distance, np):
    half = np.array((.35, .35, .3))
    edge_distance = np.inf
    for free in range(3):
        fixed = [a for a in range(3) if a != free]
        for first in (-1., 1.):
            for second in (-1., 1.):
                q = np.clip(point, -half, half)
                q[fixed[0]], q[fixed[1]] = first*half[fixed[0]], second*half[fixed[1]]
                edge_distance = min(edge_distance, np.linalg.norm(point-q))
    radial = abs(np.linalg.norm(point[:2]-np.array((.03, -.02)))-.12)
    rim = np.hypot(radial, abs(abs(point[2])-.3))
    return dict(point=point.tolist(), distance=float(distance),
        distance_to_outer_edges=float(edge_distance), distance_to_bore_cap_rims=float(rim),
        distance_to_cylinder_wall=float(radial))


def bored_prism_constraints(np):
    constraints = {}
    for axis, half in enumerate((.35, .35, .3)):
        for sign in (-1., 1.):
            def plane(x, y, z, axis=axis, half=half, sign=sign):
                return sign*(x, y, z)[axis]-half
            constraints[f'box{axis}_{sign}'] = plane
    def bore(x, y, z):
        return .12-np.sqrt((x-.03)**2+(y+.02)**2)
    constraints['bore'] = bore
    return constraints


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--case', choices=('torus', 'shell', 'disjoint', 'octahedron', 'notched-prism', 'bored-prism'), required=True)
    parser.add_argument('--tolerance', type=float, default=.002)
    parser.add_argument('--named', action='store_true', help='Retain octahedron planes or bored-prism fields as named constraints')
    parser.add_argument('--no-adapt', action='store_true')
    parser.add_argument('--optimization-passes', type=int, default=1)
    parser.add_argument('--snap', type=float, default=.2)
    parser.add_argument('--rotation', type=float, default=0., help='Rotation about z, in radians')
    parser.add_argument('--translation', type=float, nargs=3, default=(0., 0., 0.))
    parser.add_argument('--diagnostic-unchecked', action='store_true',
        help='Save rejected geometry for diagnosis; external admission gates remain unchanged')
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, _, meshers, _, _, scorecard = w.imports()
    from facet_incidence_audit import audit
    import meshers._meshers
    paths = (Path(__file__), Path(w.__file__), Path(meshers._meshers.__file__),
        Path(meshers.__file__),Path(meshers.expression.__file__),Path(meshers.intersection.__file__),
        w.HERE/'scorecard.py', w.HERE/'facet_incidence_audit.py')
    hashes = {str(p): w.digest(p) for p in paths}
    w.save(args.run/'analysis-inputs.json', hashes)
    for index, path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path, args.run/f'input-{index}-{path.name}')
    shutil.copyfile(__file__, args.run/'analysis-driver.py')
    row = dict(case=args.case, qualified=False, no_fe_solve=True,
        options=dict(bounds=w.BOUNDS, cells=24, geometry_tolerance=args.tolerance,
            minimum_quality=.1, max_tetrahedra=w.CAP, threads=1,
            adapt=not args.no_adapt, min_length=.002, max_length=.1, gradation=.3,
            adapt_iterations=4, quality_guard=.3, optimize_passes=args.optimization_passes, order=1, snap=args.snap),
        gates=dict(distance=.002, relative_volume=.01, minimum_quality=.1,
            positive_jacobians=True, coherent_facets=True, exact_boundary_topology=True,
            maximum_nodes=250000, maximum_tetrahedra=w.CAP))
    row['diagnostic_unchecked'] = args.diagnostic_unchecked
    if args.diagnostic_unchecked:
        row['options'].update(geometry_tolerance=1., minimum_quality=0., adapt_tolerance=args.tolerance)
    row['representation'] = ('named-constraints' if args.case == 'bored-prism' else 'named-planes') if args.named else 'scalar-field'
    row['transform'] = dict(rotation_z=args.rotation, translation=args.translation)
    cosine, sine = np.cos(args.rotation), np.sin(args.rotation)
    rotation = np.array(((cosine, -sine, 0.), (sine, cosine, 0.), (0., 0., 1.)))
    transformed = bool(args.rotation or any(args.translation))
    assert args.case != 'bored-prism' or not transformed
    assert np.isfinite(rotation).all() and np.isfinite(args.translation).all()
    if args.named:
        assert not transformed
    if args.named:
        assert args.case in ('octahedron', 'bored-prism')
        row['options'] = {k: v for k, v in row['options'].items() if k in
            ('bounds', 'cells', 'geometry_tolerance', 'minimum_quality', 'max_tetrahedra', 'optimize_passes', 'snap')}
    started = time.perf_counter()
    try:
        field, queries, exact_volume, components, euler = fixture(args.case, np)
        if transformed:
            original = field
            cx, cy, cz = args.translation
            def field(x, y, z):
                return original(cosine*(x-cx)+sine*(y-cy), -sine*(x-cx)+cosine*(y-cy), z-cz)
            queries = queries@rotation.T+args.translation
        np.savez_compressed(args.run/'reference.npz', points=queries)
        row['reference_sha256'] = w.digest(args.run/'reference.npz')
        phase = time.perf_counter()
        if args.named:
            from itertools import product
            constraints = bored_prism_constraints(np) if args.case == 'bored-prism' else {}
            if args.case == 'octahedron':
                for index, (a, b, c) in enumerate(product((1, -1), repeat=3)):
                    def plane(x, y, z, a=a, b=b, c=c):
                        return (a*x+b*y+c*z-.35)/np.sqrt(3)
                    constraints[f'plane{index}'] = plane
            mesh = meshers.generate_intersection(constraints, **row['options'])
        else:
            mesh = meshers.generate(field, **row['options'])
        row.update(generation_s=time.perf_counter()-phase, nodes=len(mesh.points),
            elements=len(mesh.tetrahedra), diagnostics=mesh.diagnostics)
        assert row['nodes'] <= 250000 and row['elements'] <= w.CAP
        points, cells, surface = mesh.points, mesh.tetrahedra, mesh.surface
        np.savez_compressed(args.run/'mesh.npz', points=points, tetrahedra=cells, surface=surface)
        row['mesh_sha256'] = w.digest(args.run/'mesh.npz')
        probes = np.einsum('sn,fnd->fsd', scorecard.BARY, points[surface]).reshape(-1, 3)
        forward = float(np.max(np.abs(field(*probes.T))))
        if args.case == 'bored-prism':
            distances = bored_prism_distance(probes, np)
            forward = float(distances.max())
            row['forward_witness'] = bored_prism_witness(probes[distances.argmax()], forward, np)
        if args.case in ('octahedron', 'notched-prism'):
            target_points, target_surface = octahedron(np) if args.case == 'octahedron' else notched_prism(np)
            if transformed:
                target_points = target_points@rotation.T+args.translation
            forward = float(scorecard.distance_to_mesh(target_points, target_surface, probes).max())
        reverse_distances = scorecard.distance_to_mesh(points, surface, queries)
        reverse = float(reverse_distances.max())
        if args.case == 'bored-prism':
            row['reverse_witness'] = bored_prism_witness(queries[reverse_distances.argmax()], reverse, np)
        q, _, volumes = scorecard.tet_metrics(points, cells)
        observed_components, observed_euler = boundary_topology(surface, np)
        corners = points[surface]
        boundary_volume = float(np.einsum('ij,ij->i', corners[:, 0],
            np.cross(corners[:, 1], corners[:, 2])).sum()/6)
        facets = audit(points, cells, surface)
        row.update(sampled_distance=max(forward, reverse), forward_distance=forward,
            reverse_distance=reverse, relative_volume_error=float(volumes.sum()/exact_volume-1),
            exact_volume=exact_volume, minimum_quality=float(q.min()),
            minimum_jacobian=float(6*volumes.min()), facets=facets,
            expected_boundary_components=components, boundary_components=observed_components,
            expected_boundary_euler=euler, boundary_euler=observed_euler,
            boundary_volume=boundary_volume, volume_balance=abs(abs(boundary_volume)-float(volumes.sum())),
            periodic_requested=False, used_nodes=len(np.unique(cells)))
        row['qualified'] = bool(max(forward, reverse) <= .002 and abs(row['relative_volume_error']) <= .01
            and q.min() >= .1 and volumes.min() > 0 and facets['oriented_facet_consistency_passed']
            and observed_components == components and observed_euler == euler
            and row['volume_balance'] <= 1e-12
            and row['used_nodes'] == len(points))
        row['status'] = 'admitted' if row['qualified'] else 'rejected'
    except Exception as error:
        row.update(status='failed', error=str(error), traceback=traceback.format_exc())
    row['total_s'] = time.perf_counter()-started
    w.save(args.run/'result.json', row)
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    print(f'{args.case}: {row["status"]}; {row.get("elements")} tets; '
        f'distance {row.get("sampled_distance")}; volume error {row.get("relative_volume_error")}', flush=True)
    assert row['qualified'], 'Read result.json for the supported-input failure'


if __name__ == '__main__':
    main()
