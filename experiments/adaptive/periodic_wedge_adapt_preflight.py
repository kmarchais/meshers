"""Run one private-core saved-seed experiment; never a public workflow timing."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, _, meshers, _, _, scorecard = w.imports()
    from facet_incidence_audit import audit
    from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_wedge
    parent = w.HERE/'runs/macbook-periodic-wedge-physical-seed-construction-v1'
    source = json.loads((parent/'result.json').read_text())
    assert not source['qualified'] and source['end_correspondence']['caps_match']
    assert w.digest(parent/'mesh.npz') == source['mesh_sha256']
    seed = np.load(parent/'mesh.npz')
    p, t, f, planes = (seed[key] for key in ('points', 'tetrahedra', 'surface', 'planes'))
    function, queries, exact_volume, components, euler = fixture('periodic-wedge', np)
    field = meshers.compile_field(function)
    probes = np.vstack((np.random.default_rng(42).uniform(-.5, .5, (1024, 3)),
        ((.35, .2, 0.), (.2, .2, -.5), (.2, .2, .5), (0., 0., 0.))))
    input_path = args.run.resolve()/'seed-bits.json'
    output_path = args.run.resolve()/'adapted-bits.json'
    w.save(input_path, dict(point_bits=p.view(np.uint64).tolist(), tetrahedra=t.tolist(),
        surface=f.tolist(), plane_bits=planes.view(np.uint64).tolist(), probe_bits=probes.view(np.uint64).tolist()))
    encoded = json.loads(input_path.read_text())
    np.testing.assert_array_equal(np.asarray(encoded['point_bits'], dtype=np.uint64).view(np.float64), p)
    command = ['cargo', '+1.95.0', 'test', '--release', '--locked', '-p', 'meshers-core',
        '--all-features', '--lib', 'implicit::periodic_feature_research::adapt_saved_wedge_seed',
        '--', '--ignored', '--exact', '--nocapture']
    paths = [Path(__file__), Path(w.__file__), w.HERE/'topology_coverage.py',
        w.HERE/'facet_incidence_audit.py', Path(scorecard.__file__), w.ROOT/'Cargo.lock',
        parent/'mesh.npz', parent/'result.json', input_path]
    paths.extend(sorted((w.ROOT/'crates/meshers-core/src').rglob('*.rs')))
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        if path.suffix in ('.py', '.rs'):
            shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    w.save(args.run/'preflight.json', dict(parent=str(parent), command=command,
        exact_seed_float_bits=True, minimum_quality=.1, distance=.002, periodic=(False, False, True),
        maximum_nodes=250000, maximum_tetrahedra=250000, outer_worker_timeout_s=900,
        no_public_periodic_gate_relaxed=True, no_live_generation_or_timing_claim=True))
    env = os.environ.copy()
    env.update(MESHERS_WEDGE_SEED_JSON=str(input_path), MESHERS_WEDGE_OUTPUT_JSON=str(output_path))
    subprocess.run(command, cwd=w.ROOT, env=env, check=True)
    raw = json.loads(output_path.read_text())
    actual_values = np.asarray(raw['probe_value_bits'], dtype=np.uint64).view(np.float64)
    actual_gradients = np.asarray(raw['probe_gradient_bits'], dtype=np.uint64).view(np.float64)
    np.testing.assert_allclose(actual_values, field.values(probes), atol=1e-15, rtol=0)
    expected_gradients = np.array([field.gradient(*point) for point in probes])
    np.testing.assert_allclose(actual_gradients, expected_gradients, atol=1e-14, rtol=0)
    p = np.asarray(raw['point_bits'], dtype=np.uint64).view(np.float64)
    t, f = (np.asarray(raw[key], dtype=np.int64) for key in ('tetrahedra', 'surface'))
    assert len(p) <= 250000 and len(t) <= 250000 and len(np.unique(t)) == len(p)
    facets = audit(p, t, f)
    quality, _, volumes = scorecard.tet_metrics(p, t)
    target, target_faces = periodic_wedge(np)
    samples = np.einsum('sn,fnd->fsd', scorecard.BARY, p[f]).reshape(-1, 3)
    forward = float(scorecard.distance_to_mesh(target, target_faces, samples).max())
    reverse = float(scorecard.distance_to_mesh(p, f, queries).max())
    periodic_error = None
    try:
        periodic = audit_column_periodicity(p, f, np)
    except AssertionError:
        periodic = None
        periodic_error = 'Independent periodic node/cap matching failed'
    volume_error = float(volumes.sum()/exact_volume-1)
    topology = boundary_topology(f, np)
    geometry_admitted = bool(quality.min() >= .1 and volumes.min() > 0
        and facets['oriented_facet_consistency_passed'] and max(forward, reverse) <= .002
        and abs(volume_error) <= .01 and periodic is not None and topology == (components, euler))
    np.savez_compressed(args.run/'mesh.npz', points=p, tetrahedra=t, surface=f)
    w.save(args.run/'result.json', dict(case='periodic-wedge', qualified=False,
        fixed_seed_private_core_research=True, parent=str(parent), native_finish_admitted=raw['native_finish_admitted'],
        native_finish_error=raw.get('native_finish_error'), external_geometry_admitted=geometry_admitted,
        nodes=len(p), elements=len(t), minimum_quality=float(quality.min()), minimum_jacobian=float(6*volumes.min()),
        forward_distance=forward, reverse_distance=reverse, relative_volume_error=volume_error,
        periodic=periodic, periodic_error=periodic_error, facets=facets, topology=topology,
        adapt=raw['adapt'], field_and_gradient_probe_count=len(probes), exact_seed_float_bits=True,
        mesh_sha256=w.digest(args.run/'mesh.npz'), output_bits_sha256=w.digest(output_path),
        no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False,
        scope='Saved-seed torus adaptation research with the actual core final gate; not a public generator workflow, admitted production route or timing comparison'))
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    print(f'Fixed-seed research: native final gate {raw["native_finish_admitted"]}, external geometry {geometry_admitted}, {len(t)} tets, qmin {quality.min():.6f}; no public workflow promoted', flush=True)


if __name__ == '__main__':
    main()
