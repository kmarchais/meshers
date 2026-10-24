"""Audit exact polyhedral admission and byte equality of the four study controls."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import microgen_workflow as w
from mesh_comparison_report import verify_run


def canonical_mesh(arrays, np):
    points = arrays['points']
    assert len(np.unique(points, axis=0)) == len(points)
    order = np.lexsort(points[:, ::-1].T)
    ranks = np.empty(len(points), dtype=np.int64)
    ranks[order] = np.arange(len(points))
    result = dict(points=points[order])
    for name in ('tetrahedra', 'surface'):
        cells = np.sort(ranks[arrays[name]], axis=1)
        result[name] = cells[np.lexsort(cells[:, ::-1].T)]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, _, _, _, scorecard = w.imports()
    import meshers._meshers
    from facet_incidence_audit import audit
    binary = Path(meshers._meshers.__file__)
    baseline_receipt = w.ROOT/'handoff/macbook/native-baseline-before-polyhedral-20261024.json'
    receipt = json.loads(baseline_receipt.read_text())
    baseline = next(Path(row['copy']) for row in receipt['files'] if row['original'] == str(binary))
    assert w.digest(baseline) == receipt['binary_sha256'] != w.digest(binary)
    cases = (
        ('gyroid', 'gyroid', 'meshers-tet4-24', 'macbook-mesh-only-gyroid-repeated-v1'),
        ('schwarz-solid', 'schwarz', 'meshers-tet4-24', 'macbook-mesh-only-schwarz-repeated-v1'),
        ('spinodoid', 'spinodoid', 'meshers-tet4-historical32', 'macbook-spinodoid-tet4-coverage-v1'),
        ('smooth-lattice', 'lattice', 'meshers-tet4-historical32', 'macbook-lattice-tet4-coverage-v1'))
    controls = []
    hashes = {str(p): w.digest(p) for p in (Path(__file__), binary, baseline, baseline_receipt)}
    for case, label, config, old in cases:
        run = w.HERE/'runs'/f'macbook-polyhedral-{label}-control-v1'
        verification = verify_run(run, replay_geometry=True)
        current = run/f'repeat-0-{config}'
        original = w.HERE/'runs'/old/f'repeat-0-{config}'
        new_row = json.loads((current/'result.json').read_text())
        old_row = json.loads((original/'result.json').read_text())
        assert new_row['qualified'] and old_row['qualified']
        assert new_row['config'] == old_row['config']
        assert 'polyhedral_constraints' not in new_row['diagnostics']
        assert w.digest(original/'mesh.npz') == old_row['mesh_sha256']
        expected = json.loads((original.parent/'analysis-inputs.json').read_text())
        assert expected[str(binary)] == receipt['binary_sha256']
        new, previous = np.load(current/'mesh.npz'), np.load(original/'mesh.npz')
        canonical_new, canonical_old = canonical_mesh(new, np), canonical_mesh(previous, np)
        raw = {name: bool(np.array_equal(new[name], previous[name])) for name in canonical_new}
        for name in canonical_new:
            np.testing.assert_array_equal(canonical_new[name], canonical_old[name])
        for suffix in ('-samples.npz', '-specimen.json'):
            assert w.digest(current/(case+suffix)) == w.digest(original/(case+suffix))
        controls.append(dict(case=case, verified=verification, canonical_arrays_identical=True, raw_array_equality=raw,
            original=str(original), current=str(current), elements=new_row['elements'], nodes=new_row['nodes']))
        for directory in (current, original):
            hashes[str(directory/'result.json')] = w.digest(directory/'result.json')
            hashes[str(directory/'mesh.npz')] = w.digest(directory/'mesh.npz')
    poly = w.HERE/'runs/macbook-octahedron-auto-polyhedral-v1'
    row = json.loads((poly/'result.json').read_text())
    assert row['qualified'] and row['diagnostics']['polyhedral_constraints'] == 8
    assert w.digest(poly/'mesh.npz') == row['mesh_sha256']
    arrays = np.load(poly/'mesh.npz')
    points, cells, surface = (arrays[k] for k in ('points', 'tetrahedra', 'surface'))
    assert audit(points, cells, surface)['oriented_facet_consistency_passed']
    q, _, volumes = scorecard.tet_metrics(points, cells)
    assert q.min() >= .1 and volumes.min() > 0 and abs(volumes.sum()/(4*.35**3/3)-1) < .01
    queries = np.load(poly/'reference.npz')['points']
    reverse = float(scorecard.distance_to_mesh(points, surface, queries).max())
    assert reverse <= .002 and row['sampled_distance'] <= .002
    assert len(np.unique(cells)) == len(points)
    for name in ('result.json', 'mesh.npz', 'reference.npz'):
        hashes[str(poly/name)] = w.digest(poly/name)
    historical = args.output.with_name(args.output.stem+'-historical.json')
    subprocess.run([sys.executable, str(w.HERE/'verify_topology_coverage.py'),
        '--output', str(historical), '--native-baseline', str(baseline)], check=True)
    hashes[str(historical)] = w.digest(historical)
    w.save(args.output, dict(study_controls=controls, native_baseline_sha256=receipt['binary_sha256'],
        candidate_sha256=w.digest(binary), polyhedral_reverse_distance=reverse,
        polyhedral_result=str(poly/'result.json'), inputs=hashes, no_fe_solve=True,
        reference_promoted=False, comparative_timing_claim=False))
    print('All four study meshes are exactly identical after relabeling; exact polyhedron and historical replay pass', flush=True)


if __name__ == '__main__':
    main()
