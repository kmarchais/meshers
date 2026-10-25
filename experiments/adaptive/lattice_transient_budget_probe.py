"""Diagnose temporary adaptation growth; this probe is never comparison-admitted."""

import argparse
import json
from pathlib import Path
import shutil

import admission_pilot as pilot
import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--failed', type=Path, required=True)
    parser.add_argument('--density', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    import meshers
    import meshers._meshers

    failed = json.loads((args.failed/'result.json').read_text())
    density = json.loads(args.density.read_text())
    assert not failed['qualified'] and 'exceeded max_tetrahedra' in failed['error']
    assert failed['max_length'] == density['equal_regular_tet_edge_at_actual_count']
    assert w.CAP == 250_000
    manifest = json.loads((args.failed/'manifest.json').read_text())
    assert w.digest(Path(meshers._meshers.__file__)) == manifest['native_sha256']
    paths = [Path(__file__), Path(pilot.__file__), Path(w.__file__), args.density,
        args.failed/'result.json', args.failed/'manifest.json', Path(meshers._meshers.__file__),
        w.HERE/'gradient_tie_probe.py']
    for suffix in ('-samples.npz', '-specimen.json'):
        source = args.failed/('smooth-lattice'+suffix)
        shutil.copyfile(source, args.run/source.name)
        paths.append(source)
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    shutil.copyfile(__file__, args.run/'analysis-driver.py')
    parameters = {key: failed[key] for key in ('case', 'cells', 'min_length', 'max_length',
        'order', 'native_tolerance', 'optimization_passes', 'adapt_tolerance', 'solve_node_cap')}
    parameters.update(run=args.run, compiled_spinodoid=False, symmetric_union=True,
        interior_refine_nodes=failed.get('interior_refine_nodes', 0),
        clip_void=.02, no_adapt=False, cap_length=None, cap_band=.05,
        quadrature_points=15, amg=True, solve=False)
    assert parameters['case'] == 'smooth-lattice' and parameters['order'] == 2
    original = meshers.generate
    calls = []
    def generate(*positional, **keywords):
        assert keywords['max_tetrahedra'] == w.CAP
        calls.append(dict(original_native_limit=w.CAP, diagnostic_native_limit=500_000,
            unchanged_parameters={k: v for k, v in keywords.items() if k != 'size_field'}))
        keywords['max_tetrahedra'] = 500_000
        return original(*positional, **keywords)
    meshers.generate = generate
    try:
        pilot.worker(argparse.Namespace(**parameters))
    finally:
        meshers.generate = original
    assert len(calls) == 1
    observation_path = args.run/'output-observation.json'
    (args.run/'result.json').rename(observation_path)
    observation = json.loads(observation_path.read_text())
    counts = [observation.get(k) for k in ('nodes', 'elements')]
    fits = all(isinstance(n, int) and 0 < n <= w.CAP for n in counts)
    output_hashes = {str(observation_path.resolve()): w.digest(observation_path)}
    if (args.run/'mesh.npz').exists():
        output_hashes[str((args.run/'mesh.npz').resolve())] = w.digest(args.run/'mesh.npz')
    w.save(args.run/'analysis-inputs.json', hashes)
    w.save(args.run/'result.json', dict(qualified=False, comparison_eligible=False,
        scope='Temporary-capacity diagnostic only; normal intermediate/output admission remains 250000, no FE solve',
        normal_tet_and_node_admission_cap=w.CAP, diagnostic_native_limit=500_000,
        original_failure=failed['error'], source_status=observation['status'],
        nodes=observation.get('nodes'), elements=observation.get('elements'),
        final_counts_fit_normal_admission=fits, output_geometry_passed=observation.get('qualified', False),
        output_error=observation.get('error'), calls=calls, inputs=hashes,
        output_artifact_sha256=output_hashes, new_solve=False, reference_promoted=False))
    assert all(w.digest(Path(p)) == value for p, value in hashes.items())
    print(f'Diagnostic output: nodes={counts[0]}, tets={counts[1]}, counts fit normal cap={fits}; comparison eligible=False', flush=True)


if __name__ == '__main__':
    main()
