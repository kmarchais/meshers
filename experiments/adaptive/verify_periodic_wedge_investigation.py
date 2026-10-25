"""Replay periodic-seam rejections and the separate manufactured feasibility mesh."""

import argparse
import json
from pathlib import Path

import numpy as np
from vtkmodules.vtkCommonCore import vtkSMPTools

import microgen_workflow as w
import scorecard
from facet_incidence_audit import audit
from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_wedge


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    assert vtkSMPTools.GetBackend() == 'Sequential' and vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    paths = [Path(__file__), Path(w.__file__), Path(scorecard.__file__),
        w.HERE/'topology_coverage.py', w.HERE/'facet_incidence_audit.py',
        w.HERE/'periodic_wedge_feasibility.py', w.HERE/'test_topology_coverage.py']
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    field, queries, exact_volume, components, euler = fixture('periodic-wedge', np)
    failures = []
    expected_errors = (
        ('macbook-periodic-wedge-scalar-v1', 'sampled geometry error', 1),
        ('macbook-periodic-wedge-mapped-v1', 'sampled surface error', 1),
        ('macbook-periodic-wedge-mapped-seed-v1', 'minimum MMG quality', 0),
    )
    for name, error_kind, passes in expected_errors:
        run = w.HERE/'runs'/name
        record = json.loads((run/'result.json').read_text())
        worker = json.loads((run/'worker.json').read_text())
        assert record['case'] == 'periodic-wedge' and not record['qualified'] and record['status'] == 'failed'
        assert error_kind in record['error']
        assert worker['returncode'] == 1 and worker['limit_failure'] is None
        assert record['options']['cells'] == 24 and record['options']['optimize_passes'] == passes
        assert record['options']['periodic'] == [False, False, True]
        assert record['gates'] == dict(distance=.002, relative_volume=.01, minimum_quality=.1,
            positive_jacobians=True, coherent_facets=True, exact_boundary_topology=True,
            maximum_nodes=250000, maximum_tetrahedra=250000)
        assert record['options']['geometry_tolerance'] == .002 and record['options']['minimum_quality'] == .1
        np.testing.assert_array_equal(queries, np.load(run/'reference.npz')['points'])
        assert w.digest(run/'reference.npz') == record['reference_sha256']
        assert not (run/'mesh.npz').exists()
        for path, expected in json.loads((run/'analysis-inputs.json').read_text()).items():
            resolved = Path(path)
            if w.digest(resolved) != expected and resolved.suffix == '.py':
                matches = [p for p in run.glob('*.py') if w.digest(p) == expected]
                assert matches, path
                resolved = matches[0]
            assert w.digest(resolved) == expected, path
            hashes[str(resolved)] = expected
        for filename in ('result.json', 'worker.json', 'reference.npz'):
            hashes[str(run/filename)] = w.digest(run/filename)
        failures.append(dict(run=str(run), representation=record['representation'],
            options=record['options'], gates=record['gates'], error=record['error'],
            no_successful_mesh_returned=True, resource_limit_failure=None))
    run = w.HERE/'runs/macbook-periodic-wedge-feasibility-v1'
    record = json.loads((run/'result.json').read_text())
    assert record['qualified'] and record['manufactured_feasibility_control']
    assert not record['native_generator_called'] and record['no_fe_solve']
    for path, expected in json.loads((run/'analysis-inputs.json').read_text()).items():
        assert w.digest(Path(path)) == expected
        hashes[path] = expected
    arrays = np.load(run/'mesh.npz')
    p, t, f = (arrays[key] for key in ('points', 'tetrahedra', 'surface'))
    target, target_faces = periodic_wedge(np)
    np.testing.assert_array_equal(p[:len(target)], target)
    np.testing.assert_array_equal(f, target_faces)
    assert len(p) == 13 and len(t) == 20 and len(np.unique(t)) == len(p)
    assert audit(p, t, f)['oriented_facet_consistency_passed']
    assert boundary_topology(f, np) == (components, euler)
    periodic = audit_column_periodicity(p, f, np)
    assert periodic == record['periodic']
    quality, _, volumes = scorecard.tet_metrics(p, t)
    assert quality.min() >= .1 and volumes.min() > 0 and abs(volumes.sum()/exact_volume-1) < .01
    bary = np.array([(i/8, j/8, 1-(i+j)/8) for i in range(9) for j in range(9-i)])
    probes = np.einsum('sn,fnd->fsd', bary, p[f]).reshape(-1, 3)
    forward = float(scorecard.distance_to_mesh(target, target_faces, probes).max())
    reverse = float(scorecard.distance_to_mesh(p, f, queries).max())
    assert max(forward, reverse) <= .002
    boundary = p[f]
    boundary_volume = np.einsum('ij,ij->i', boundary[:, 0], np.cross(boundary[:, 1], boundary[:, 2])).sum()/6
    assert abs(boundary_volume-volumes.sum()) < 1e-12
    assert w.digest(run/'mesh.npz') == record['mesh_sha256']
    assert w.digest(run/'reference.npz') == record['reference_sha256']
    np.testing.assert_array_equal(queries, np.load(run/'reference.npz')['points'])
    for filename in ('result.json', 'worker.json', 'mesh.npz', 'reference.npz'):
        hashes[str(run/filename)] = w.digest(run/filename)
    w.save(args.output, dict(failures=failures, manufactured_feasibility=dict(run=str(run),
        nodes=len(p), elements=len(t), minimum_quality=float(quality.min()),
        minimum_jacobian=float(6*volumes.min()), forward_distance=forward, reverse_distance=reverse,
        relative_volume_error=float(volumes.sum()/exact_volume-1), periodic=periodic,
        native_generator_called=False), inputs=hashes, no_fe_solve=True, reference_promoted=False,
        comparative_timing_claim=False, native_periodic_seam_failure_resolved=False,
        scope='Exact-domain geometry/quality feasibility and preserved native/map-route rejections; root cause and native fix remain open'))
    print('Three genuine rejections and the separate 20-tet exact-domain feasibility mesh pass replay; native seam failure remains open', flush=True)


if __name__ == '__main__':
    main()
