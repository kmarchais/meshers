"""MMG tensor control on a saved native seed; never an actual Microgen comparator."""

import argparse
import json
from pathlib import Path
import shutil
import time
import traceback
from types import SimpleNamespace

import microgen_workflow as w


def face_keys(points, faces):
    return sorted(tuple(sorted(tuple(p) for p in face)) for face in points[faces[:, :3]])


def make_backend(points, cells, surface, tensor):
    import mmgpy
    import numpy as np

    backend = mmgpy.MmgMesh3D()
    backend.set_mesh_size(vertices=len(points), tetrahedra=len(cells), triangles=len(surface))
    backend.set_vertices(points)
    backend.set_tetrahedra(cells.astype(np.int32))
    backend.set_triangles(surface.astype(np.int32))
    backend.set_field('tensor', tensor)
    np.testing.assert_array_equal(backend.get_field('tensor'), tensor)
    return backend


def adapter_control():
    import numpy as np

    points = .01*np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1]])
    cells = np.array([[0, 1, 2, 3]])
    surface = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]])
    tensor = np.tile([1., 0, 0, 4, 0, 9], (4, 1))
    backend = make_backend(points, cells, surface, tensor)
    assert len(backend.get_triangles()) == 4 and len(backend.get_tetrahedra()) == 1
    backend.remesh(nosurf=True, hmin=.002, hmax=.1, hgrad=1.3, anisosize=True, verbose=-1)
    assert len(backend.get_tetrahedra()) > 0
    return dict(passed=True, tensor_round_trip_exact=True, explicit_triangle_storage=True,
        public_native_adapter=True, anisotropic_remesh_options_accepted=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, _, _, _, _, _ = w.imports()
    import mmgpy
    import mmgpy._mmgpy
    from mmgpy.metrics import matrix_to_tensor, tensor_to_matrix
    from scipy.spatial import cKDTree
    from periodic import periodic_pairs
    from periodic_crease_audit import gradient
    from periodic_crease_grading import make_metric
    from microgen_tet10_control import straight_elevation
    from facet_incidence_audit import FACES

    row = dict(case='smooth-lattice', order=2, qualified=False, preflight_passed=False,
        native_only=False, actual_microgen_route=False, saved_native_seed_replayed=True,
        comparison_eligible=False, new_solve=False, reference_promoted=False,
        source=str(args.source), status='started', maximum_directional_ratio=5.,
        scope='MMG directional feasibility on a saved Meshers seed. Not native Meshers performance or an actual Microgen comparator.')
    paths = [Path(__file__), Path(w.__file__), Path(mmgpy._mmgpy.__file__),
        args.source/'result.json', args.source/'manifest.json', args.source/'crease-cloud.npz',
        args.source/'feature-seed/mesh.npz']
    paths.extend(w.HERE/n for n in ('periodic_crease_grading.py', 'lattice_crease_mechanism.py',
        'periodic_crease_audit.py', 'microgen_tet10_control.py', 'facet_incidence_audit.py'))
    for suffix in ('-samples.npz', '-specimen.json'):
        p = args.source/('smooth-lattice'+suffix)
        paths.append(p); shutil.copyfile(p, args.run/p.name)
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    shutil.copyfile(__file__, args.run/'analysis-driver.py')
    w.save(args.run/'analysis-inputs.json', hashes)
    start = time.perf_counter()
    try:
        row['adapter_control'] = adapter_control()
        original = json.loads((args.source/'result.json').read_text())
        assert not original['qualified'] and '250002 > 250000' in original['error']
        assert w.digest(args.source/'feature-seed/mesh.npz') == original['crease_grading']['seed_mesh_sha256']
        assert original['crease_grading']['seed_geometry']['geometry_admitted']
        source = np.load(args.source/'feature-seed/mesh.npz')
        seed = SimpleNamespace(**{k: source[k] for k in ('points', 'tetrahedra', 'surface')})
        corners = np.unique(seed.tetrahedra[:, :4])
        np.testing.assert_array_equal(corners, np.arange(len(corners)))
        points, cells = seed.points[:len(corners)], seed.tetrahedra[:, :4]
        case, _, _ = w.specimen('smooth-lattice')
        scalar, metadata, cloud = make_metric(seed, case.function, .002, .1, 26)
        np.testing.assert_array_equal(cloud, np.load(args.source/'crease-cloud.npz')['canonical_points'])
        frame = np.zeros_like(cloud)
        physical_cloud = cloud-.5
        planes = np.isclose(abs(physical_cloud), .5, atol=1e-10, rtol=0)
        for axis in range(3):
            ids = np.flatnonzero((planes.sum(axis=1) == 1) & planes[:, axis])
            p = physical_cloud[ids].copy()
            p[:, axis] -= np.sign(p[:, axis])*1e-6
            normal = gradient(case.function, p, 1e-6/8)
            tangent = np.cross(np.eye(3)[axis], normal)
            tangent /= np.linalg.norm(tangent, axis=1)[:, None]
            frame[ids] = tangent
        tree = cKDTree(cloud, boxsize=1.)
        nearest = tree.query(np.mod(points+.5, 1.), workers=1)[1]
        tangent = frame[nearest]
        normal_size = scalar(*points.T)
        along_size = np.minimum(.1, 5*normal_size)
        matrix = np.eye(3)[None]/normal_size[:, None, None]**2
        matrix = matrix+(along_size**-2-normal_size**-2)[:, None, None]*np.einsum('ni,nj->nij', tangent, tangent)
        eigenvalues = np.linalg.eigvalsh(matrix)
        assert np.isfinite(matrix).all() and eigenvalues.min() > 0
        assert (eigenvalues[:, -1]/eigenvalues[:, 0]).max() <= 25*(1+1e-12)
        packed = matrix_to_tensor(matrix)
        np.testing.assert_array_equal(tensor_to_matrix(packed), matrix)
        np.savez_compressed(args.run/'metric.npz', tensor=packed)
        backend = make_backend(points.copy(), cells, seed.surface[:, :3], packed)
        options = dict(nosurf=True, hmin=.002, hmax=.1, hgrad=1.3, anisosize=True, verbose=-1)
        row.update(metric=metadata, mmg_options=options, metric_vertex_count=len(points),
            metric_condition_max=float((eigenvalues[:, -1]/eigenvalues[:, 0]).max()),
            internal_transient_counts_available=False)
        w.save(args.run/'result.json', row)
        phase = time.perf_counter()
        backend.remesh(**options)
        p, t = backend.get_vertices(), backend.get_tetrahedra()
        t, reoriented = w.reorient_linear(p, t)
        row.update(remesh_s=time.perf_counter()-phase, raw_nodes=len(p), raw_tets=len(t), reoriented_tets=reoriented)
        assert max(len(p), len(t)) <= w.CAP, 'MMG output stage cap'
        used = np.unique(t)
        inverse = np.full(len(p), -1, dtype=np.int64); inverse[used] = np.arange(len(used))
        p, t = p[used], inverse[t]
        faces = t[:, FACES].reshape(-1, 3)
        _, indices, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
        assert counts.max() <= 2
        surface = faces[indices[counts == 1]]
        np.savez_compressed(args.run/'raw-mmg-mesh.npz', points=p, tetrahedra=t, surface=surface)
        old_boundary = np.unique(seed.surface[:, :3])
        new_boundary = np.unique(surface)
        row.update(old_boundary_corner_nodes=len(old_boundary), new_boundary_corner_nodes=len(new_boundary))
        assert len(old_boundary) == len(new_boundary), 'Boundary corner-node count changed'
        distances, nearest = cKDTree(seed.points[old_boundary]).query(p[new_boundary], workers=1)
        matched = old_boundary[nearest]
        limit = 64*np.finfo(float).eps*np.linalg.norm(np.ptp(seed.points, axis=0))
        row.update(boundary_corner_maximum_raw_difference=float(distances.max()),
            permitted_roundoff_restoration=limit, boundary_corner_mapping_one_to_one=len(np.unique(matched)) == len(matched))
        assert len(np.unique(matched)) == len(matched), 'Boundary correspondence is not one-to-one'
        assert distances.max() <= limit, 'Boundary coordinate change exceeds roundoff'
        mapping = np.full(len(p), -1, dtype=np.int64); mapping[new_boundary] = matched
        expected = np.sort(seed.surface[:, :3], axis=1)
        observed = np.sort(mapping[surface], axis=1)
        np.testing.assert_array_equal(expected[np.lexsort(expected.T[::-1])], observed[np.lexsort(observed.T[::-1])])
        p[new_boundary] = seed.points[matched]
        row['boundary_corner_roundoff_restored_exactly'] = True
        assert face_keys(p, surface) == face_keys(seed.points, seed.surface), 'Frozen surface changed'
        elevated = straight_elevation(p, t, surface)
        row.update(nodes=elevated['nodes'], elements=len(t))
        assert elevated['within_node_cap'], 'Quadratic node cap'
        p2, t2, s2 = (elevated[k] for k in ('points', 'tetrahedra', 'surface'))
        old_mids = {}
        for face in seed.surface:
            for k, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
                key = tuple(sorted((tuple(seed.points[face[a]]), tuple(seed.points[face[b]]))))
                value = seed.points[face[k+3]]
                if key in old_mids:
                    np.testing.assert_array_equal(old_mids[key], value)
                old_mids[key] = value
        for face in s2:
            for k, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
                key = tuple(sorted((tuple(p2[face[a]]), tuple(p2[face[b]]))))
                p2[face[k+3]] = old_mids[key]
        np.savez_compressed(args.run/'mesh.npz', points=p2, tetrahedra=t2, surface=s2)
        geometry = w.geometry(args.run, 'smooth-lattice', p2, t2, s2)
        _, _, periodic = periodic_pairs(p2, s2)
        volume_change = geometry['volume']/original['crease_grading']['seed_geometry']['volume']-1
        assert abs(volume_change) < 1e-10
        row.update(geometry=geometry, periodic=periodic, relative_volume_change=volume_change,
            boundary_corner_faces_exact=True, boundary_quadratic_nodes_restored_exactly=True,
            approximation_spaces_nested=False, preflight_passed=geometry['geometry_admitted'],
            status='diagnostic-preflight-passed' if geometry['geometry_admitted'] else 'diagnostic-geometry-rejected')
    except Exception as error:
        row.update(status='diagnostic-failed', error=str(error), traceback=traceback.format_exc())
    row.update(elapsed_s=time.perf_counter()-start, inputs=hashes)
    w.save(args.run/'result.json', row)
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    assert row['preflight_passed'], row.get('error', row['status'])
    print('Directional diagnostic passes geometry preflight; not a native or actual-Microgen competitor', flush=True)


if __name__ == '__main__':
    main()
