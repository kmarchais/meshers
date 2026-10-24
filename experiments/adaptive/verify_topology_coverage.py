"""Replay new analytic controls and the original four admitted Tet4 artifacts."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
from topology_coverage import boundary_topology, fixture, octahedron


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, _, _, _, scorecard = w.imports()
    import meshers._meshers
    from periodic import periodic_pairs
    from facet_incidence_audit import audit
    binary = Path(meshers._meshers.__file__)
    inputs = {str(p): w.digest(p) for p in (Path(__file__), w.HERE/'topology_coverage.py',
        Path(w.__file__), w.HERE/'scorecard.py', binary)}
    results = []
    analytic = ('macbook-topology-torus-volume-v1', 'macbook-topology-shell-v1',
        'macbook-topology-disjoint-volume-v1', 'macbook-topology-octahedron-named-v1')
    fixed = (
        ('gyroid', 'macbook-mesh-only-gyroid-repeated-v1/repeat-0-meshers-tet4-24'),
        ('schwarz-solid', 'macbook-mesh-only-schwarz-repeated-v1/repeat-0-meshers-tet4-24'),
        ('spinodoid', 'macbook-spinodoid-tet4-coverage-v1/repeat-0-meshers-tet4-historical32'),
        ('smooth-lattice', 'macbook-lattice-tet4-coverage-v1/repeat-0-meshers-tet4-historical32'))
    for name, relative in [(None, path) for path in analytic]+list(fixed):
        directory = w.HERE/'runs'/relative
        row = json.loads((directory/'result.json').read_text())
        assert row['qualified']
        run = directory if name is None else directory.parent
        measured = json.loads((run/'analysis-inputs.json').read_text())
        assert measured[str(binary)] == w.digest(binary)
        for source, expected in measured.items():
            path = Path(source)
            if w.digest(path) != expected and path.suffix == '.py':
                matches = [p for p in run.glob('*.py') if w.digest(p) == expected]
                assert matches, source
                path = matches[0]
            assert w.digest(path) == expected, source
        inputs[str(directory/'result.json')] = w.digest(directory/'result.json')
        assert w.digest(directory/'mesh.npz') == row['mesh_sha256']
        inputs[str(directory/'mesh.npz')] = row['mesh_sha256']
        arrays = np.load(directory/'mesh.npz')
        points, cells, surface = (arrays[k] for k in ('points', 'tetrahedra', 'surface'))
        assert cells.shape[1] == 4 and surface.shape[1] == 3
        assert len(points) <= 250000 and len(cells) <= w.CAP
        assert len(np.unique(cells)) == len(points)
        facets = audit(points, cells, surface)
        assert facets['oriented_facet_consistency_passed']
        quality, _, volumes = scorecard.tet_metrics(points, cells)
        assert quality.min() >= .1 and volumes.min() > 0
        if name is None:
            name = row['case']
            field, queries, exact_volume, components, euler = fixture(name, np)
            assert w.digest(directory/'reference.npz') == row['reference_sha256']
            np.testing.assert_array_equal(queries, np.load(directory/'reference.npz')['points'])
            assert boundary_topology(surface, np) == (components, euler)
            bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
            probes = np.einsum('sn,fnd->fsd', bary, points[surface]).reshape(-1, 3)
            forward = float(np.max(np.abs(field(*probes.T))))
            if name == 'octahedron':
                target_points, target_surface = octahedron(np)
                forward = float(scorecard.distance_to_mesh(target_points, target_surface, probes).max())
            distance = max(forward, float(scorecard.distance_to_mesh(points, surface, queries).max()))
            volume_error = float(volumes.sum()/exact_volume-1)
            assert distance <= .002 and abs(volume_error) <= .01
            faces = points[surface]
            boundary_volume = np.einsum('ij,ij->i', faces[:, 0], np.cross(faces[:, 1], faces[:, 2])).sum()/6
            assert abs(boundary_volume-volumes.sum()) < 1e-12
        else:
            geometry = w.geometry(directory, name, points, cells, surface, box_projection_fallback=True)
            assert geometry['geometry_admitted']
            _, _, periodic = periodic_pairs(points, surface)
            assert periodic == row['periodic']
            distance, volume_error = geometry['sampled_two_sided_max'], geometry['relative_volume_error']
            for suffix in ('-specimen.json', '-samples.npz'):
                assert w.digest(directory/(name+suffix)) == w.digest(w.HERE/'runs/g3-microgen-defaults-v4'/(name+suffix))
        results.append(dict(case=name, directory=str(directory), nodes=len(points), elements=len(cells),
            sampled_distance=distance, relative_volume_error=volume_error, minimum_quality=float(quality.min()), qualified=True))
    w.save(args.output, dict(cases=results, inputs=inputs, native_changed=False,
        new_equilibrium_solve=False, reference_promoted=False, comparative_timing_claim=False))
    print(f'All {len(results)} saved analytic and fixed-study Tet4 artifacts pass replay', flush=True)


if __name__ == '__main__':
    main()
