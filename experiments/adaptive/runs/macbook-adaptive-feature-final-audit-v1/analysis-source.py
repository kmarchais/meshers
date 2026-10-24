"""Replay four fixed controls and the adaptive polyhedron under unchanged gates."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
from mesh_comparison_report import verify_run
from verify_polyhedral_change import canonical_mesh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, _, _, _, scorecard = w.imports()
    import meshers._meshers
    from facet_incidence_audit import audit
    binary = Path(meshers._meshers.__file__)
    receipt_path = w.ROOT/'handoff/macbook/native-baseline-before-adaptive-polyhedral-20261024.json'
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
        current = w.HERE/'runs'/f'macbook-adaptive-feature-{label}-control-v1'
        old = w.HERE/'runs'/f'macbook-polyhedral-{label}-control-v1'
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
    poly = w.HERE/'runs/macbook-adaptive-feature-octahedron-final-v1'
    row = json.loads((poly/'result.json').read_text())
    assert row['qualified'] and row['diagnostics']['polyhedral_constraints'] == 8
    assert w.digest(poly/'mesh.npz') == row['mesh_sha256']
    arrays = np.load(poly/'mesh.npz')
    points, cells, surface = (arrays[name] for name in ('points', 'tetrahedra', 'surface'))
    assert audit(points, cells, surface)['oriented_facet_consistency_passed']
    q, _, volumes = scorecard.tet_metrics(points, cells)
    reverse = float(scorecard.distance_to_mesh(points, surface, np.load(poly/'reference.npz')['points']).max())
    assert q.min() >= .1 and volumes.min() > 0 and abs(volumes.sum()/(4*.35**3/3)-1) < .01
    assert reverse <= .002 and row['sampled_distance'] <= .002
    assert len(points) <= 250000 and len(cells) <= 250000
    for name in ('result.json', 'mesh.npz', 'reference.npz'):
        hashes[str(poly/name)] = w.digest(poly/name)
    w.save(args.output, dict(study_controls=controls, inputs=hashes,
        native_baseline_sha256=w.digest(baseline), candidate_sha256=w.digest(binary),
        polyhedral_reverse_distance=reverse, polyhedral_result=str(poly/'result.json'),
        no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False))
    print('Four unchanged study meshes and the adaptive exact polyhedron pass independent replay', flush=True)


if __name__ == '__main__':
    main()
