"""Replay four fixed controls and a retained-feature specimen under unchanged gates."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
from mesh_comparison_report import verify_run
from verify_polyhedral_change import canonical_mesh


def replay_periodic_column(run, np, scorecard, hashes, baseline=None):
    from facet_incidence_audit import audit
    from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_column, periodic_wedge
    record = json.loads((run/'result.json').read_text())
    case = record['case']
    assert case in ('periodic-column', 'periodic-wedge') and record['periodic_requested']
    inputs = json.loads((run/'analysis-inputs.json').read_text())
    for source, expected in inputs.items():
        path = Path(source)
        if w.digest(path) != expected and path.suffix == '.so' and baseline is not None:
            path = baseline
        if w.digest(path) != expected and path.suffix == '.py':
            matches = [p for p in run.glob('*.py') if w.digest(p) == expected]
            assert matches, source
            path = matches[0]
        assert w.digest(path) == expected, source
        hashes[str(path)] = expected
    for name, key in (('mesh.npz', 'mesh_sha256'), ('reference.npz', 'reference_sha256')):
        assert w.digest(run/name) == record[key]
        hashes[str(run/name)] = record[key]
    hashes[str(run/'result.json')] = w.digest(run/'result.json')
    mesh = np.load(run/'mesh.npz')
    p, t, f = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
    assert len(p) <= 250000 and len(t) <= 250000 and len(np.unique(t)) == len(p)
    assert audit(p, t, f)['oriented_facet_consistency_passed']
    periodic = audit_column_periodicity(p, f, np)
    assert periodic == record['periodic']
    _, queries, exact_volume, components, euler = fixture(case, np)
    target, faces = periodic_column(np) if case == 'periodic-column' else periodic_wedge(np)
    pose = record.get('transform', dict(rotation_z=0., translation=(0., 0., 0.)))
    if pose['rotation_z'] or any(pose['translation']):
        assert case == 'periodic-wedge' and pose['translation'][2] == 0
        c, s = np.cos(pose['rotation_z']), np.sin(pose['rotation_z'])
        rotation = np.array(((c, -s, 0.), (s, c, 0.), (0., 0., 1.)))
        target = target@rotation.T+pose['translation']
        queries = queries@rotation.T+pose['translation']
        assert np.max(abs(target[:, :2])) < .5
    np.testing.assert_array_equal(queries, np.load(run/'reference.npz')['points'])
    bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
    probes = np.einsum('sn,fnd->fsd', bary, p[f]).reshape(-1, 3)
    forward = float(scorecard.distance_to_mesh(target, faces, probes).max())
    distances = scorecard.distance_to_mesh(p, f, queries)
    reverse = float(distances.max())
    q, _, volumes = scorecard.tet_metrics(p, t)
    assert q.min() >= .1 and volumes.min() > 0
    volume_error = float(volumes.sum()/exact_volume-1)
    assert abs(volume_error) <= .01
    assert boundary_topology(f, np) == (components, euler)
    boundary_volume = np.einsum('ij,ij->i', p[f[:, 0]], np.cross(p[f[:, 1]], p[f[:, 2]])).sum()/6
    assert abs(boundary_volume-volumes.sum()) < 1e-12
    qualified = max(forward, reverse) <= .002
    assert qualified == record['qualified']
    witness = queries[distances.argmax()]
    c, s = np.cos(.31), np.sin(.31)
    local_xy = np.array((c*witness[0]+s*witness[1], -s*witness[0]+c*witness[1]))
    return dict(run=str(run), qualified=qualified, nodes=len(p), elements=len(t),
        minimum_quality=float(q.min()), forward_distance=forward, reverse_distance=reverse,
        relative_volume_error=volume_error, periodic=periodic, options=record['options'],
        reverse_witness=witness.tolist(), reverse_witness_distance_to_corner_line=(float(np.linalg.norm(abs(local_xy)-.2))
            if case == 'periodic-column' else None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prefix', default='macbook-adaptive-feature')
    parser.add_argument('--previous-prefix', default='macbook-polyhedral')
    parser.add_argument('--baseline-receipt', default='native-baseline-before-adaptive-polyhedral-20261024.json')
    parser.add_argument('--feature-case', choices=('octahedron', 'notched-prism', 'bored-prism', 'periodic-column'), default='octahedron')
    parser.add_argument('--feature-run', default='macbook-adaptive-feature-octahedron-final-v1')
    parser.add_argument('--matched-feature-run')
    parser.add_argument('--shell-run')
    parser.add_argument('--failure-run')
    parser.add_argument('--named-periodic-run')
    parser.add_argument('--allow-changed-adaptive-controls', action='store_true')
    parser.add_argument('--control-run', nargs=2, action='append', default=[], metavar=('CASE', 'DIRECTORY'))
    parser.add_argument('--previous-control-run', nargs=2, action='append', default=[], metavar=('CASE', 'DIRECTORY'))
    args = parser.parse_args()
    assert not args.output.exists()
    overrides = dict(args.control_run)
    assert len(overrides) == len(args.control_run)
    assert overrides.keys() <= {'gyroid', 'schwarz-solid', 'spinodoid', 'smooth-lattice'}
    previous_overrides = dict(args.previous_control_run)
    assert len(previous_overrides) == len(args.previous_control_run)
    assert previous_overrides.keys() <= {'gyroid', 'schwarz-solid', 'spinodoid', 'smooth-lattice'}
    np, _, _, _, _, scorecard = w.imports()
    import meshers._meshers
    from facet_incidence_audit import audit
    binary = Path(meshers._meshers.__file__)
    receipt_path = w.ROOT/'handoff/macbook'/args.baseline_receipt
    receipt = json.loads(receipt_path.read_text())
    baseline = next(Path(row['copy']) for row in receipt['files'] if row['original'] == str(binary))
    assert w.digest(baseline) == receipt['binary_sha256'] != w.digest(binary)
    hashes = {str(p): w.digest(p) for p in (Path(__file__), binary, baseline, receipt_path)}
    controls = []
    for case, label, config in (
        ('gyroid', 'gyroid', 'meshers-tet4-24'),
        ('schwarz-solid', 'schwarz', 'meshers-tet4-24'),
        ('spinodoid', 'spinodoid', 'meshers-tet4-historical32'),
        ('smooth-lattice', 'lattice', 'meshers-tet4-historical32'),
    ):
        current = w.HERE/'runs'/overrides.get(case, f'{args.prefix}-{label}-control-v1')
        old = w.HERE/'runs'/previous_overrides.get(case, f'{args.previous_prefix}-{label}-control-v1')
        verification = verify_run(current, replay_geometry=True)
        previous_verification = verify_run(old, replay_geometry=True, native_baseline=baseline)
        directories = [run/f'repeat-0-{config}' for run in (current, old)]
        meshes = [np.load(directory/'mesh.npz') for directory in directories]
        canonical = [canonical_mesh(mesh, np) for mesh in meshes]
        identical = all(np.array_equal(canonical[0][name], canonical[1][name]) for name in canonical[0])
        assert identical or (args.allow_changed_adaptive_controls and config == 'meshers-tet4-historical32')
        for suffix in ('-samples.npz', '-specimen.json'):
            assert w.digest(directories[0]/(case+suffix)) == w.digest(directories[1]/(case+suffix))
        for directory in directories:
            row = json.loads((directory/'result.json').read_text())
            assert row['qualified'] and 'polyhedral_constraints' not in row['diagnostics']
            for name in ('mesh.npz', 'result.json'):
                hashes[str(directory/name)] = w.digest(directory/name)
        records = [json.loads((directory/'result.json').read_text()) for directory in directories]
        assert records[0]['config'] == records[1]['config']
        controls.append(dict(case=case, current=verification, baseline=previous_verification,
            canonical_arrays_identical=identical,
            current_mesh=dict(nodes=records[0]['nodes'], elements=records[0]['elements'], geometry=records[0]['geometry']),
            baseline_mesh=dict(nodes=records[1]['nodes'], elements=records[1]['elements'], geometry=records[1]['geometry'])))
    poly = w.HERE/'runs'/args.feature_run
    row = json.loads((poly/'result.json').read_text())
    assert row['qualified']
    if args.feature_case == 'bored-prism':
        assert row['diagnostics']['retained_constraints'] == 7
    else:
        assert row['diagnostics']['polyhedral_constraints'] == (4 if args.feature_case == 'periodic-column' else 8)
    if args.feature_case == 'notched-prism':
        assert row['diagnostics']['polyhedral_regions'] == 2
    assert w.digest(poly/'mesh.npz') == row['mesh_sha256']
    arrays = np.load(poly/'mesh.npz')
    points, cells, surface = (arrays[name] for name in ('points', 'tetrahedra', 'surface'))
    assert audit(points, cells, surface)['oriented_facet_consistency_passed']
    q, _, volumes = scorecard.tet_metrics(points, cells)
    reverse = float(scorecard.distance_to_mesh(points, surface, np.load(poly/'reference.npz')['points']).max())
    from topology_coverage import bored_prism_distance, boundary_topology, fixture, notched_prism, octahedron, periodic_column
    _, queries, exact_volume, _, _ = fixture(args.feature_case, np)
    target, faces = (periodic_column(np) if args.feature_case == 'periodic-column'
        else notched_prism(np) if args.feature_case == 'notched-prism' else octahedron(np))
    transform = row.get('transform', dict(rotation_z=0., translation=(0., 0., 0.)))
    angle, translation = transform['rotation_z'], transform['translation']
    if angle or any(translation):
        c, s = np.cos(angle), np.sin(angle)
        rotation = np.array(((c, -s, 0.), (s, c, 0.), (0., 0., 1.)))
        target, queries = target@rotation.T+translation, queries@rotation.T+translation
    np.testing.assert_array_equal(queries, np.load(poly/'reference.npz')['points'])
    probes = np.einsum('sn,fnd->fsd', scorecard.BARY, points[surface]).reshape(-1, 3)
    forward = float((bored_prism_distance(probes, np) if args.feature_case == 'bored-prism'
        else scorecard.distance_to_mesh(target, faces, probes)).max())
    _, _, _, components, euler = fixture(args.feature_case, np)
    assert boundary_topology(surface, np) == (components, euler)
    assert q.min() >= .1 and volumes.min() > 0 and abs(volumes.sum()/exact_volume-1) < .01
    assert forward <= .002
    assert reverse <= .002 and row['sampled_distance'] <= .002
    assert len(points) <= 250000 and len(cells) <= 250000
    for name in ('result.json', 'mesh.npz', 'reference.npz'):
        hashes[str(poly/name)] = w.digest(poly/name)
    periodic_result = None
    if args.feature_case == 'periodic-column':
        candidate = replay_periodic_column(poly, np, scorecard, hashes)
        assert candidate['qualified']
        periodic_result = dict(candidate=candidate, no_comparative_timing_claim=True)
        if args.failure_run or args.named_periodic_run:
            assert args.failure_run and args.named_periodic_run
            failed = replay_periodic_column(w.HERE/'runs'/args.failure_run, np, scorecard, hashes, baseline)
            named = replay_periodic_column(w.HERE/'runs'/args.named_periodic_run, np, scorecard, hashes, baseline)
            assert named['qualified'] and not failed['qualified'] and failed['reverse_distance'] > .002
            for other in (failed, named):
                for key in ('bounds', 'cells', 'periodic', 'geometry_tolerance', 'minimum_quality',
                    'max_tetrahedra', 'optimize_passes', 'snap'):
                    assert other['options'][key] == candidate['options'][key]
            assert failed['options'] == candidate['options']
            periodic_result.update(rejected_scalar=failed, named_seed=named, fixed_domain_and_gates=True)
    matched_result = None
    if args.matched_feature_run:
        matched = w.HERE/'runs'/args.matched_feature_run
        other = json.loads((matched/'result.json').read_text())
        assert other['qualified'] and other['case'] == row['case']
        for key in ('bounds', 'cells', 'geometry_tolerance', 'minimum_quality',
            'max_tetrahedra', 'optimize_passes', 'snap'):
            assert other['options'][key] == row['options'][key]
        assert w.digest(matched/'mesh.npz') == other['mesh_sha256']
        assert w.digest(matched/'reference.npz') == w.digest(poly/'reference.npz')
        canonical = canonical_mesh(np.load(matched/'mesh.npz'), np)
        current = canonical_mesh(arrays, np)
        for name in canonical:
            np.testing.assert_array_equal(canonical[name], current[name])
        for name in ('result.json', 'mesh.npz', 'reference.npz'):
            hashes[str(matched/name)] = w.digest(matched/name)
        matched_result = dict(run=str(matched), canonical_arrays_identical=True,
            matched_geometry_and_quality_settings=True)
    shell_result = None
    if args.shell_run:
        shell = w.HERE/'runs'/args.shell_run
        record = json.loads((shell/'result.json').read_text())
        assert record['qualified'] and record['case'] == 'shell'
        assert record['diagnostics']['retained_constraints'] == 2
        mesh = np.load(shell/'mesh.npz')
        p, t, f = (mesh[key] for key in ('points', 'tetrahedra', 'surface'))
        field, queries, volume, components, euler = fixture('shell', np)
        np.testing.assert_array_equal(queries, np.load(shell/'reference.npz')['points'])
        assert w.digest(shell/'mesh.npz') == record['mesh_sha256']
        assert w.digest(shell/'reference.npz') == record['reference_sha256']
        assert audit(p, t, f)['oriented_facet_consistency_passed']
        assert boundary_topology(f, np) == (components, euler)
        q, _, volumes = scorecard.tet_metrics(p, t)
        assert q.min() >= .1 and volumes.min() > 0 and abs(volumes.sum()/volume-1) < .01
        bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
        probes = np.einsum('sn,fnd->fsd', bary, p[f]).reshape(-1, 3)
        forward_shell = float(np.abs(field(*probes.T)).max())
        reverse_shell = float(scorecard.distance_to_mesh(p, f, queries).max())
        assert max(forward_shell, reverse_shell) <= .002
        assert len(p) <= 250000 and len(t) <= 250000 and len(np.unique(t)) == len(p)
        boundary_volume = np.einsum('ij,ij->i', p[f[:, 0]], np.cross(p[f[:, 1]], p[f[:, 2]])).sum()/6
        assert abs(boundary_volume-volumes.sum()) < 1e-12
        for name in ('mesh.npz', 'reference.npz', 'result.json'):
            hashes[str(shell/name)] = w.digest(shell/name)
        shell_result = dict(run=str(shell), admitted=True, sampled_distance=max(forward_shell, reverse_shell),
            minimum_quality=float(q.min()), relative_volume_error=float(volumes.sum()/volume-1))
    for path in (w.HERE/'topology_coverage.py', w.HERE/'scorecard.py',
        w.HERE/'facet_incidence_audit.py', w.HERE/'mesh_comparison_report.py'):
        hashes[str(path)] = w.digest(path)
    from shutil import copyfile
    archived_source = poly/'verification-source.py'
    assert not archived_source.exists()
    copyfile(__file__, archived_source)
    w.save(args.output, dict(study_controls=controls, inputs=hashes, feature_case=args.feature_case,
        matched_feature=matched_result, shell_regression=shell_result, periodic_feature=periodic_result,
        native_baseline_sha256=w.digest(baseline), candidate_sha256=w.digest(binary),
        polyhedral_reverse_distance=reverse, polyhedral_forward_distance=forward, polyhedral_result=str(poly/'result.json'),
        no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False))
    print(f'Four admitted study meshes and the {args.feature_case} pass independent replay; exact matches: {sum(row["canonical_arrays_identical"] for row in controls)}/4', flush=True)


if __name__ == '__main__':
    main()
