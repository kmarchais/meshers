"""Replay saved-seed affine mobility evidence without claiming a public route."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--candidate-run', default='macbook-periodic-wedge-fixed-seed-adapt-v3')
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, _, _, _, scorecard = w.imports()
    from facet_incidence_audit import audit
    from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_wedge
    runs = [w.HERE/'runs/macbook-periodic-wedge-fixed-seed-adapt-v2', w.HERE/'runs'/args.candidate_run]
    assert w.digest(runs[0]/'seed-bits.json') == w.digest(runs[1]/'seed-bits.json')
    _, queries, volume, components, euler = fixture('periodic-wedge', np)
    target, faces = periodic_wedge(np)
    hashes = {str(Path(__file__).resolve()): w.digest(Path(__file__))}
    sources, rows = [], []
    for run in runs:
        inputs = json.loads((run/'analysis-inputs.json').read_text())
        sources.append({p: h for p, h in inputs.items() if p.endswith('.rs')})
        for source, expected in inputs.items():
            path = Path(source)
            if not path.exists() or w.digest(path) != expected:
                matches = [p for p in run.glob('input-*') if p.is_file() and w.digest(p) == expected]
                assert matches, source
                path = matches[0]
            assert w.digest(path) == expected
            hashes[str(path)] = expected
        result = json.loads((run/'result.json').read_text())
        assert not result['qualified'] and result['fixed_seed_private_core_research']
        assert result['exact_seed_float_bits'] and not result['comparative_timing_claim']
        assert w.digest(run/'mesh.npz') == result['mesh_sha256']
        assert w.digest(run/'adapted-bits.json') == result['output_bits_sha256']
        mesh = np.load(run/'mesh.npz')
        p, t, f = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
        assert len(p) <= 250000 and len(t) <= 250000 and len(np.unique(t)) == len(p)
        assert audit(p, t, f)['oriented_facet_consistency_passed']
        assert boundary_topology(f, np) == (components, euler)
        periodic = audit_column_periodicity(p, f, np)
        q, _, volumes = scorecard.tet_metrics(p, t)
        assert volumes.min() > 0 and abs(volumes.sum()/volume-1) <= .01
        bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
        probes = np.einsum('sn,fnd->fsd', bary, p[f]).reshape(-1, 3)
        distance = max(float(scorecard.distance_to_mesh(target, faces, probes).max()),
            float(scorecard.distance_to_mesh(p, f, queries).max()))
        assert distance < 1e-12
        bad = int(np.sum(q < .1))
        admitted = bad == 0
        assert result['native_finish_admitted'] == admitted == result['external_geometry_admitted']
        rows.append(dict(run=str(run), nodes=len(p), elements=len(t), minimum_quality=float(q.min()),
            below_quality_gate=bad, sampled_distance=distance, relative_volume_error=float(volumes.sum()/volume-1),
            periodic=periodic, native_finish_admitted=admitted, adapt=result['adapt']))
        for name in ('result.json', 'mesh.npz', 'seed-bits.json', 'adapted-bits.json', 'analysis-inputs.json'):
            hashes[str(run/name)] = w.digest(run/name)
    assert sources[0].keys() == sources[1].keys()
    changed = [p for p in sources[0] if sources[0][p] != sources[1][p]]
    assert {Path(p).name for p in changed} == {'tet_adapt.rs', 'polyhedral.rs'}
    assert rows[0]['below_quality_gate'] == 10 and rows[1]['below_quality_gate'] == 0
    assert rows[1]['minimum_quality'] >= .3
    controls = []
    for name, expected in (('red', 101), ('green', 0)):
        run = w.HERE/'runs'/f'macbook-affine-curve-{name}-v1'
        worker = json.loads((run/'worker.json').read_text())
        assert worker['returncode'] == expected
        log = (run/'worker.log').read_text()
        assert 'running 1 test' in log and 'affine_curve_mobility_counts_independent_constraints' in log
        assert ('two planes leave a movable line' in log) if name == 'red' else ('1 passed' in log)
        for filename in ('worker.json', 'worker.log', 'manifest.json'):
            hashes[str(run/filename)] = w.digest(run/filename)
        controls.append(dict(run=str(run), returncode=expected))
    w.save(args.output, dict(passed=True, rows=rows, changed_rust_sources=changed,
        red_green_controls=controls, inputs=hashes, public_route_promoted=False,
        comparative_timing_claim=False, reference_promoted=False,
        scope='Exact saved-seed private-core recovery only; the public periodic prism route remains unresolved'))
    print('Affine mobility replay passes: identical seed, ten bad cells to zero; no public route promoted', flush=True)


if __name__ == '__main__':
    main()
