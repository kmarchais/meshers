"""Replay the actual public normal-jump recipe and preserve its prior rejection."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
from verify_adaptive_polyhedral_change import replay_periodic_column


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, meshers, _, _, scorecard = w.imports()
    import meshers._meshers
    from topology_coverage import fixture
    binary = Path(meshers._meshers.__file__)
    paths = [Path(__file__), w.HERE/'verify_adaptive_polyhedral_change.py',
        w.HERE/'topology_coverage.py', w.HERE/'facet_incidence_audit.py', Path(scorecard.__file__), binary]
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    run = w.HERE/'runs/macbook-periodic-affine-traces-prism-v1'
    recovered = replay_periodic_column(run, np, scorecard, hashes)
    assert recovered['qualified'] and recovered['minimum_quality'] >= .3
    record = json.loads((run/'result.json').read_text())
    assert record['representation'] == 'scalar-field' and not record['diagnostic_unchecked']
    assert record['diagnostics']['polyhedral_constraints'] == 6
    assert not record.get('mapped', False)
    old = w.HERE/'runs/macbook-periodic-wedge-public-wrapped-v1'
    failed = json.loads((old/'result.json').read_text())
    assert not failed['qualified'] and failed['status'] == 'failed'
    assert 'sampled geometry error' in failed['error'] and not (old/'mesh.npz').exists()
    assert failed['options'] == record['options'] and failed['gates'] == record['gates']
    assert failed['reference_sha256'] == record['reference_sha256']
    for source, expected in json.loads((old/'analysis-inputs.json').read_text()).items():
        path = Path(source)
        if not path.exists() or w.digest(path) != expected:
            if path.suffix == '.so':
                path = w.ROOT/'.local-work/native-baselines'/expected/'meshers'/path.name
            else:
                matches = [p for p in old.glob('*.py') if w.digest(p) == expected]
                assert matches, source
                path = matches[0]
        assert w.digest(path) == expected, source
        hashes[str(path)] = expected
    function, _, _, _, _ = fixture('periodic-wedge', np)
    field = meshers.compile_field(function)
    planes = np.asarray(field._halfspaces)
    points = np.random.default_rng(42).uniform(-.5, .5, (1024, 3))
    np.testing.assert_allclose((points@planes[:, :3].T+planes[:, 3]).max(axis=1),
        field.values(points), atol=3e-16, rtol=0)
    planes = planes/np.linalg.norm(planes[:, :3], axis=1)[:, None]
    traces = [np.column_stack((planes[:, :2], planes[:, 3]+side*planes[:, 2])) for side in (-.5, .5)]
    assert sorted(map(tuple, traces[0])) == sorted(map(tuple, traces[1]))
    assert np.any(planes[:, 2] != 0)
    controls = []
    for label, code, text in [('red', 1, '1 failed'), ('green', 0, '2 passed')]:
        control = w.HERE/'runs'/f'macbook-periodic-affine-traces-{label}-v1'
        worker = json.loads((control/'worker.json').read_text())
        assert worker['returncode'] == code and text in (control/'worker.log').read_text()
        for name in ('worker.json', 'worker.log', 'manifest.json'):
            hashes[str(control/name)] = w.digest(control/name)
        controls.append(dict(run=str(control), returncode=code))
    for name in ('result.json', 'reference.npz', 'worker.json'):
        hashes[str(old/name)] = w.digest(old/name)
    w.save(args.output, dict(passed=True, public_recipe_recovered=True, recovered=recovered,
        previous_rejection=dict(run=str(old), error=failed['error'], qualified=False),
        fixed_options_and_gates=True, fixed_sample_hash=record['reference_sha256'],
        complete_normalized_affine_trace_multiset_equal=True, red_green_controls=controls,
        inputs=hashes, no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False,
        scope='One previously rejected public periodic affine prism recipe recovered; no universal coverage or original study FE convergence claim'))
    print('Actual public prism recovery passes independent replay at unchanged specimen, options and gates', flush=True)


if __name__ == '__main__':
    main()
