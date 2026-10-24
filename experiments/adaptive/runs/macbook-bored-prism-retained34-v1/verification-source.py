"""Replay four fixed controls and a retained-feature specimen under unchanged gates."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
from mesh_comparison_report import verify_run
from verify_polyhedral_change import canonical_mesh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prefix', default='macbook-adaptive-feature')
    parser.add_argument('--previous-prefix', default='macbook-polyhedral')
    parser.add_argument('--baseline-receipt', default='native-baseline-before-adaptive-polyhedral-20261024.json')
    parser.add_argument('--feature-case', choices=('octahedron', 'notched-prism', 'bored-prism'), default='octahedron')
    parser.add_argument('--feature-run', default='macbook-adaptive-feature-octahedron-final-v1')
    parser.add_argument('--matched-feature-run')
    args = parser.parse_args()
    assert not args.output.exists()
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
        current = w.HERE/'runs'/f'{args.prefix}-{label}-control-v1'
        old = w.HERE/'runs'/f'{args.previous_prefix}-{label}-control-v1'
        verification = verify_run(current, replay_geometry=True)
        previous_verification = verify_run(old, replay_geometry=True, native_baseline=baseline)
        directories = [run/f'repeat-0-{config}' for run in (current, old)]
        meshes = [np.load(directory/'mesh.npz') for directory in directories]
        canonical = [canonical_mesh(mesh, np) for mesh in meshes]
        for name in canonical[0]:
            np.testing.assert_array_equal(canonical[0][name], canonical[1][name])
        for suffix in ('-samples.npz', '-specimen.json'):
            assert w.digest(directories[0]/(case+suffix)) == w.digest(directories[1]/(case+suffix))
        for directory in directories:
            row = json.loads((directory/'result.json').read_text())
            assert row['qualified'] and 'polyhedral_constraints' not in row['diagnostics']
            for name in ('mesh.npz', 'result.json'):
                hashes[str(directory/name)] = w.digest(directory/name)
        controls.append(dict(case=case, current=verification, baseline=previous_verification,
            canonical_arrays_identical=True))
    poly = w.HERE/'runs'/args.feature_run
    row = json.loads((poly/'result.json').read_text())
    assert row['qualified']
    if args.feature_case == 'bored-prism':
        assert row['diagnostics']['retained_constraints'] == 7
    else:
        assert row['diagnostics']['polyhedral_constraints'] == 8
    if args.feature_case == 'notched-prism':
        assert row['diagnostics']['polyhedral_regions'] == 2
    assert w.digest(poly/'mesh.npz') == row['mesh_sha256']
    arrays = np.load(poly/'mesh.npz')
    points, cells, surface = (arrays[name] for name in ('points', 'tetrahedra', 'surface'))
    assert audit(points, cells, surface)['oriented_facet_consistency_passed']
    q, _, volumes = scorecard.tet_metrics(points, cells)
    reverse = float(scorecard.distance_to_mesh(points, surface, np.load(poly/'reference.npz')['points']).max())
    from topology_coverage import bored_prism_distance, boundary_topology, fixture, notched_prism, octahedron
    _, queries, exact_volume, _, _ = fixture(args.feature_case, np)
    target, faces = notched_prism(np) if args.feature_case == 'notched-prism' else octahedron(np)
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
    for path in (w.HERE/'topology_coverage.py', w.HERE/'scorecard.py',
        w.HERE/'facet_incidence_audit.py', w.HERE/'mesh_comparison_report.py'):
        hashes[str(path)] = w.digest(path)
    from shutil import copyfile
    archived_source = poly/'verification-source.py'
    assert not archived_source.exists()
    copyfile(__file__, archived_source)
    w.save(args.output, dict(study_controls=controls, inputs=hashes, feature_case=args.feature_case,
        matched_feature=matched_result,
        native_baseline_sha256=w.digest(baseline), candidate_sha256=w.digest(binary),
        polyhedral_reverse_distance=reverse, polyhedral_forward_distance=forward, polyhedral_result=str(poly/'result.json'),
        no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False))
    print(f'Four unchanged study meshes and the {args.feature_case} pass independent replay', flush=True)


if __name__ == '__main__':
    main()
