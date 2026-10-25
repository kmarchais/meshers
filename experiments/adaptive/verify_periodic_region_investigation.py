"""Replay the re-entrant periodic failure, feasibility and private exact-seed controls."""

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
    from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_notched_wedge
    _, queries, volume, components, euler = fixture('periodic-notched-wedge', np)
    target, faces = periodic_notched_wedge(np)
    hashes = {str(Path(__file__).resolve()): w.digest(Path(__file__))}

    def sources(run):
        for source, expected in json.loads((run/'analysis-inputs.json').read_text()).items():
            path = Path(source)
            if not path.exists() or w.digest(path) != expected:
                matches = [p for p in run.glob('input-*') if p.is_file() and w.digest(p) == expected]
                if path.suffix == '.so':
                    matches.append(w.ROOT/'.local-work/native-baselines'/expected/'meshers'/path.name)
                assert matches, source
                path = matches[0]
            assert w.digest(path) == expected, source
            hashes[str(path)] = expected
        for name in ('result.json', 'worker.json', 'mesh.npz', 'reference.npz', 'analysis-inputs.json'):
            if (run/name).exists():
                hashes[str(run/name)] = w.digest(run/name)

    failed = w.HERE/'runs/macbook-periodic-notched-wedge-scalar-v1'
    sources(failed)
    failure = json.loads((failed/'result.json').read_text())
    assert not failure['qualified'] and failure['status'] == 'failed'
    assert 'sampled geometry error' in failure['error'] and not (failed/'mesh.npz').exists()
    np.testing.assert_array_equal(queries, np.load(failed/'reference.npz')['points'])
    rows = []
    runs = [w.HERE/'runs'/f'macbook-periodic-notched-wedge-{label}-v1'
        for label in ('feasibility', 'physical-seed', 'fixed-seed-adapt')]
    for run, floor in zip(runs, (.1, 0., .1), strict=True):
        sources(run)
        record = json.loads((run/'result.json').read_text())
        assert w.digest(run/'mesh.npz') == record['mesh_sha256']
        mesh = np.load(run/'mesh.npz')
        p, t, f = (mesh[key] for key in ('points', 'tetrahedra', 'surface'))
        assert len(p) <= 250000 and len(t) <= 250000 and len(np.unique(t)) == len(p)
        assert audit(p, t, f)['oriented_facet_consistency_passed']
        assert boundary_topology(f, np) == (components, euler)
        q, _, volumes = scorecard.tet_metrics(p, t)
        assert q.min() >= floor and volumes.min() > 0 and abs(volumes.sum()/volume-1) <= .01
        periodic = audit_column_periodicity(p, f, np)
        bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
        probes = np.einsum('sn,fnd->fsd', bary, p[f]).reshape(-1, 3)
        distance = max(float(scorecard.distance_to_mesh(target, faces, probes).max()),
            float(scorecard.distance_to_mesh(p, f, queries).max()))
        assert distance < 1e-12
        rows.append(dict(run=str(run), nodes=len(p), elements=len(t), minimum_quality=float(q.min()),
            positive_quadrature_corner_jacobians=True, sampled_distance=distance,
            relative_volume_error=float(volumes.sum()/volume-1), periodic=periodic,
            public_periodic_candidate=False))
    feasibility, seed, adapted = [json.loads((r/'result.json').read_text()) for r in runs]
    assert feasibility['manufactured_feasibility_control'] and not feasibility['native_generator_called']
    assert not seed['qualified'] and seed['region_masks'] == [127, 191]
    encoded = json.loads((runs[2]/'seed-bits.json').read_text())
    mesh = np.load(runs[1]/'mesh.npz')
    np.testing.assert_array_equal(np.asarray(encoded['point_bits'], dtype=np.uint64).view(np.float64), mesh['points'])
    np.testing.assert_array_equal(encoded['tetrahedra'], mesh['tetrahedra'])
    np.testing.assert_array_equal(encoded['surface'], mesh['surface'])
    assert encoded['region_masks'] == seed['region_masks']
    planes = mesh['planes']/np.linalg.norm(mesh['planes'][:, :3], axis=1)[:, None]
    traces = []
    for z in (-.5, .5):
        equations = np.column_stack((planes[:, :2], planes[:, 3]+z*planes[:, 2]))
        traces.append(sorted(tuple(sorted(tuple(equations[i]) for i in range(len(planes)) if mask & (1 << i)))
            for mask in seed['region_masks']))
    assert traces[0] == traces[1]
    assert not adapted['qualified'] and adapted['fixed_seed_private_core_research']
    assert adapted['native_finish_admitted'] and adapted['external_geometry_admitted']
    assert w.digest(runs[2]/'adapted-bits.json') == adapted['output_bits_sha256']
    for name in ('seed-bits.json', 'adapted-bits.json'):
        hashes[str(runs[2]/name)] = w.digest(runs[2]/name)
    w.save(args.output, dict(passed=True, public_failure=dict(run=str(failed), error=failure['error']),
        controls=rows, exact_seed_bits_preserved=True, complete_region_trace_multisets_equal=True,
        public_recipe_recovered=False, no_fe_solve=True, reference_promoted=False,
        comparative_timing_claim=False, inputs=hashes,
        scope='Manufactured feasibility and private exact-seed recovery; original public periodic region recipe remains rejected'))
    print('Periodic region investigation passes independent replay; public failure retained', flush=True)


if __name__ == '__main__':
    main()
