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
    if name == 'octahedron':
        vertices, faces = octahedron(np)
        bary = np.array([(i/48, j/48, 1-(i+j)/48) for i in range(49) for j in range(49-i)])
        queries = np.einsum('sn,fnd->fsd', bary, vertices[faces]).reshape(-1, 3)
        def field(x, y, z):
            return (np.abs(x)+np.abs(y)+np.abs(z)-.35)/np.sqrt(3)
        return field, queries, 4*.35**3/3, 1, 2
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--case', choices=('torus', 'shell', 'disjoint', 'octahedron'), required=True)
    parser.add_argument('--tolerance', type=float, default=.002)
    parser.add_argument('--named', action='store_true', help='Preserve the octahedron planes as named constraints')
    parser.add_argument('--no-adapt', action='store_true')
    parser.add_argument('--optimization-passes', type=int, default=1)
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
            adapt_iterations=4, quality_guard=.3, optimize_passes=args.optimization_passes, order=1),
        gates=dict(distance=.002, relative_volume=.01, minimum_quality=.1,
            positive_jacobians=True, coherent_facets=True, exact_boundary_topology=True,
            maximum_nodes=250000, maximum_tetrahedra=w.CAP))
    row['representation'] = 'named-planes' if args.named else 'scalar-field'
    if args.named:
        assert args.case == 'octahedron'
        row['options'] = {k: v for k, v in row['options'].items() if k in
            ('bounds', 'cells', 'geometry_tolerance', 'minimum_quality', 'max_tetrahedra', 'optimize_passes')}
    started = time.perf_counter()
    try:
        field, queries, exact_volume, components, euler = fixture(args.case, np)
        np.savez_compressed(args.run/'reference.npz', points=queries)
        row['reference_sha256'] = w.digest(args.run/'reference.npz')
        phase = time.perf_counter()
        if args.named:
            from itertools import product
            constraints = {}
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
        if args.case == 'octahedron':
            target_points, target_surface = octahedron(np)
            forward = float(scorecard.distance_to_mesh(target_points, target_surface, probes).max())
        reverse = float(scorecard.distance_to_mesh(points, surface, queries).max())
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
