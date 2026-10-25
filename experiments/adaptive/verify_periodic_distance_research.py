"""Replay the controlled periodic-query correction without promoting the failed route."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, _, _, _, scorecard = w.imports()
    from facet_incidence_audit import audit
    from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_wedge
    baseline_receipt = w.ROOT/'handoff/macbook/native-baseline-before-periodic-distance-20261025.json'
    receipt = json.loads(baseline_receipt.read_text())
    baseline = next(Path(row['copy']) for row in receipt['files'] if row['original'].endswith('_meshers.abi3.so'))
    assert w.digest(baseline) == receipt['binary_sha256']
    hashes = {str(p.resolve()): w.digest(p) for p in (Path(__file__), baseline_receipt, baseline,
        w.HERE/'topology_coverage.py', w.HERE/'scorecard.py', w.HERE/'facet_incidence_audit.py')}
    def inputs(run):
        data = json.loads((run/'analysis-inputs.json').read_text())
        for source, expected in data.items():
            path = Path(source)
            if w.digest(path) != expected and path.suffix == '.so':
                path = baseline
            if w.digest(path) != expected and path.suffix in ('.py', '.rs'):
                matches = [p for p in run.glob('input-*') if p.is_file() and w.digest(p) == expected]
                assert matches, source
                path = matches[0]
            assert w.digest(path) == expected, source
            hashes[str(path)] = expected
        return data
    _, queries, volume, components, euler = fixture('periodic-wedge', np)
    target, faces = periodic_wedge(np)
    parent = w.HERE/'runs/macbook-periodic-wedge-physical-seed-construction-v1'
    planes = np.load(parent/'mesh.npz')['planes']
    rows = []
    runs = [w.HERE/'runs'/f'macbook-periodic-wedge-fixed-seed-adapt-v{i}' for i in (1, 2)]
    assert w.digest(runs[0]/'seed-bits.json') == w.digest(runs[1]/'seed-bits.json')
    source_records = []
    for run in runs:
        source_records.append(inputs(run))
        record = json.loads((run/'result.json').read_text())
        assert not record['qualified'] and not record['native_finish_admitted']
        assert 'minimum MMG quality' in record['native_finish_error']
        assert record['fixed_seed_private_core_research'] and record['exact_seed_float_bits']
        assert not record['external_geometry_admitted']
        assert w.digest(run/'mesh.npz') == record['mesh_sha256']
        assert w.digest(run/'adapted-bits.json') == record['output_bits_sha256']
        mesh = np.load(run/'mesh.npz')
        p, t, f = (mesh[key] for key in ('points', 'tetrahedra', 'surface'))
        assert audit(p, t, f)['oriented_facet_consistency_passed']
        periodic = audit_column_periodicity(p, f, np)
        assert periodic == record['periodic'] and boundary_topology(f, np) == (components, euler)
        q, _, volumes = scorecard.tet_metrics(p, t)
        assert len(p) <= 250000 and len(t) <= 250000 and len(np.unique(t)) == len(p)
        assert volumes.min() > 0 and abs(volumes.sum()/volume-1) <= .01
        bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
        probes = np.einsum('sn,fnd->fsd', bary, p[f]).reshape(-1, 3)
        distance = max(float(scorecard.distance_to_mesh(target, faces, probes).max()),
            float(scorecard.distance_to_mesh(p, f, queries).max()))
        assert distance < 1e-12
        bad = np.flatnonzero(q < .1)
        assert len(bad) == 10
        diagnostics = []
        for cell in bad:
            corners = p[t[cell]]
            on_planes = abs(corners@planes[:, :3].T+planes[:, 3])/np.linalg.norm(planes[:, :3], axis=1) < 1e-11
            diagnostics.append(dict(cell=int(cell), quality=float(q[cell]),
                points=corners.tolist(), vertices=t[cell].tolist(),
                on_low_cap=bool(np.any(abs(corners[:, 2]+.5) < 1e-10)),
                on_high_cap=bool(np.any(abs(corners[:, 2]-.5) < 1e-10)),
                on_midplane=bool(np.any(abs(corners[:, 2]) < 1e-10)),
                active_plane_counts=on_planes.sum(axis=1).tolist()))
        for filename in ('mesh.npz', 'result.json', 'adapted-bits.json', 'seed-bits.json', 'preflight.json'):
            hashes[str(run/filename)] = w.digest(run/filename)
        rows.append(dict(run=str(run), nodes=len(p), elements=len(t), minimum_quality=float(q.min()),
            sampled_distance=distance, periodic=periodic, adapt=record['adapt'],
            native_finish_error=record['native_finish_error'], bad_cells=diagnostics))
    source_records = [{path: value for path, value in record.items() if not path.endswith('/seed-bits.json')} for record in source_records]
    assert source_records[0].keys() == source_records[1].keys()
    changed_sources = [p for p in source_records[0] if source_records[0][p] != source_records[1][p]]
    assert len(changed_sources) == 1 and changed_sources[0].endswith('/tet_adapt.rs')
    assert rows[0]['adapt']['error_splits'] == 111 and rows[1]['adapt']['error_splits'] == 0
    assert rows[0]['adapt']['maximum_error'] > .002 and rows[1]['adapt']['maximum_error'] < 1e-12
    public_rows = []
    for name in ('macbook-periodic-wedge-scalar-v1', 'macbook-periodic-wedge-public-wrapped-v1'):
        run = w.HERE/'runs'/name
        inputs(run)
        record = json.loads((run/'result.json').read_text())
        assert not record['qualified'] and record['status'] == 'failed'
        assert 'sampled geometry error' in record['error'] and not (run/'mesh.npz').exists()
        np.testing.assert_array_equal(np.load(run/'reference.npz')['points'], queries)
        hashes[str(run/'result.json')] = w.digest(run/'result.json')
        public_rows.append(dict(run=str(run), options=record['options'], error=record['error']))
    assert public_rows[0]['options'] == public_rows[1]['options']
    for name, code in (('macbook-periodic-distance-red-v1', 101), ('macbook-periodic-distance-green-v1', 0)):
        run = w.HERE/'runs'/name
        worker = json.loads((run/'worker.json').read_text())
        assert worker['returncode'] == code and worker['limit_failure'] is None
        for filename in ('worker.json', 'worker.log', 'manifest.json'):
            hashes[str(run/filename)] = w.digest(run/filename)
    w.save(args.output, dict(fixed_seed_controls=rows, exact_seed_bits_identical=True,
        changed_input_sources=changed_sources, public_generator_controls=public_rows, inputs=hashes,
        no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False,
        public_periodic_wedge_admitted=False,
        scope='Verified periodic image-query correction and removal of false error splits; both saved-seed quality failures and the original public feature failure remain rejected'))
    print('Periodic query correction replays; false splits111to0; both ten-cell quality failures and public wedge rejection remain explicit', flush=True)


if __name__ == '__main__':
    main()
