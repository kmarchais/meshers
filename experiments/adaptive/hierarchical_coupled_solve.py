"""Six zero-macro responses and a small Schur solve on an immutable Tet10 domain."""

import argparse
import json
from pathlib import Path
import shutil
import signal
import time

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--preflight', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    from amg_reference import AuditedAMG
    from capture_reference_system import CapturedSystem
    from hp_coupled_preflight import inherited_inputs
    from hp_coupled_solve import condense
    from scipy.linalg import eigvalsh
    from stress_recovery import solve
    from periodic import periodic_pairs
    checked = json.loads((args.preflight/'result.json').read_text())
    assert checked['passed'] and checked['consistent_rule_passed']
    assert checked['resources']['amg_setup_performed'] and not checked['resources']['cg_started']
    assert not checked['resources']['dense_mode_rows_added_to_sparse_matrix']
    source, trial = Path(checked['source']), Path(checked['trial'])
    paths = [Path(__file__), args.preflight/'result.json', args.preflight/'coupling.npz',
        args.preflight/'consistent.json', args.preflight/'worker.json']
    paths.extend(w.HERE/name for name in ('hp_coupled_solve.py', 'amg_reference.py', 'stress_recovery.py',
        'hierarchical_residual_trial.py', 'test_hierarchical_residual_trial.py', 'test_hp_coupled_solve.py',
        'hp_coupled_preflight.py', 'microgen_workflow.py'))
    hashes = {**inherited_inputs(args.preflight), **{str(p.resolve()): w.digest(p) for p in paths}}
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    mesh = np.load(source/'mesh.npz')
    points, cells, surface = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
    fields = np.load(source/'fields.npz')['displacement']
    assert len(points) <= 250000 and len(cells) <= 250000
    consistent = json.loads((args.preflight/'consistent.json').read_text())
    assert consistent['coupling_sha256'] == w.digest(args.preflight/'coupling.npz')
    data = np.load(args.preflight/'coupling.npz')
    forces, energy = data['forces'], data['energy']

    class Correction:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr()
            transform = self.problem._MatCB
            full_fields = np.zeros((transform.shape[0],6))
            full_fields[:3*len(points)] = fields.transpose(0,2,1).reshape(6,-1).T
            internal = self.problem.get_A()@full_fields
            residuals = np.linalg.norm(transform.T@internal,axis=0)/np.maximum(np.linalg.norm(internal,axis=0),1e-30)
            assert residuals.max() < 2e-9
            padded = np.zeros_like(full_fields); padded[:len(forces)] = forces
            reduced = transform.T@padded
            responses = np.empty_like(reduced)
            amg = AuditedAMG(elastic_blocks=True); amg.problem = self.problem
            started = time.monotonic()

            def expired(*_):
                raise TimeoutError('Auxiliary response exceeded the unchanged 240-second cap')

            old_handler = signal.signal(signal.SIGALRM, expired)
            try:
                for column in range(6):
                    signal.setitimer(signal.ITIMER_REAL,240)
                    try:
                        responses[:,column] = amg(matrix,reduced[:,column])
                    finally:
                        signal.setitimer(signal.ITIMER_REAL,0)
                    assert amg.key[1] == checked['resources']['matrix_sha256']
                    full = transform@responses[:,column]
                    physical = full[:3*len(points)].reshape(3,len(points)).T
                    np.savez_compressed(args.run/f'response-{column+1}.npz', displacement=physical)
                    w.save(args.run/'checkpoint.json', dict(completed_auxiliary_responses=column+1,
                        solver_records=amg.records, elapsed_s=time.monotonic()-started))
                    print(f'Auxiliary response {column+1}/6 completed', flush=True)
            finally:
                signal.signal(signal.SIGALRM,old_handler)
            alpha, tensor, gain, schur = condense(energy,reduced,responses)
            physical = (transform@responses)[:3*len(points)].T.reshape(6,3,len(points)).transpose(0,2,1)
            pairs, _, periodic = periodic_pairs(points,surface)
            jump = max(float(abs(physical[:,p[:,1]]-physical[:,p[:,0]]).max(initial=0)) for p in pairs)
            pin = int(np.linalg.norm(points,axis=1).argmin())
            assert jump < 1e-9 and abs(physical[:,pin]).max() < 1e-12
            modal_residual = energy[6:,:6]-energy[:6,:6]+schur@alpha
            mode_relative = float(np.linalg.norm(modal_residual)/np.linalg.norm(energy[6:,:6]-energy[:6,:6]))
            assert mode_relative < 1e-10
            archive = args.run/'coupled-fields.npz'
            np.savez_compressed(archive,response_displacement=physical,mode_coefficients=alpha)
            w.save(args.run/'result.json', dict(status='coupled-solved-needs-independent-audit',
                source=str(source), trial=str(trial), preflight=str(args.preflight),
                geometry_unchanged=True, old_approximation_space_preserved=True,
                physical_nodes=len(points), physical_tets=len(cells), extra_scalar_mode_unknowns=6,
                reduced_parent_unknowns=len(rhs), augmented_unknowns=len(rhs)+6,
                fresh_auxiliary_responses=6, saved_parent_reduced_residuals=residuals.tolist(),
                solver_records=amg.records, predicted_energy_tensor=tensor.tolist(),
                energy_gain_tensor=gain.tolist(), maximum_relative_gain=float(eigvalsh(gain,energy[:6,:6]).max()),
                schur_eigenvalues=eigvalsh(schur).tolist(), modal_relative_residual=mode_relative,
                zero_macro_response_jump_max=jump, response_pin_max=float(abs(physical[:,pin]).max()),
                periodic=periodic, auxiliary_solve_s=time.monotonic()-started,
                field_archive_sha256=w.digest(archive), field_representation=
                'Immutable parent Tet10 fields plus (analytic cubic bubbles minus quadratic responses) times a 6x6 coefficient matrix; nodal values alone omit the bubble field',
                reference_promoted=False, native_candidate_promoted=False, comparative_timing_claim=False,
                scope='Coupled same-domain six-mode equilibrium awaiting independent complete-field and positive-rule audits; not native workflow timing or continuum/local stress certification'))
            raise CapturedSystem()

    try:
        solve(points,cells,n_elm_gp=15,solver=Correction(),release_reduced_cache=True)
    except CapturedSystem:
        assert (args.run/'result.json').exists()
    else:
        raise AssertionError('Do not return a fake macro solution from the capture')
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())


if __name__ == '__main__':
    main()
