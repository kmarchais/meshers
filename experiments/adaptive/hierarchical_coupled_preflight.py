"""Consistent-rule coupling and sparse-parent resource/setup gates; no CG solve."""

import argparse
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--trial', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    import hierarchical_residual_trial as ht
    from hp_coupled_preflight import inherited_inputs
    from quadrature_audit import duffy_rule
    from capture_reference_system import CapturedSystem
    from stress_recovery import solve
    from amg_reference import AuditedAMG
    checked = json.loads((args.trial/'result.json').read_text())
    assert checked['passed'] and checked['gain_signal_passed'] and checked['positive_rule_precision_passed']
    source = Path(checked['source'])
    row = json.loads((source/'result.json').read_text())
    paths = [Path(__file__), Path(ht.__file__), args.trial/'result.json', args.trial/'modes.npz']
    paths.extend(w.HERE/name for name in ('hp_coupled_preflight.py', 'hp_cross_energy.py',
        'hp_coupled_solve.py', 'amg_reference.py', 'stress_recovery.py', 'capture_reference_system.py',
        'quadrature_audit.py', 'test_hierarchical_residual_trial.py', 'microgen_workflow.py'))
    hashes = {**inherited_inputs(args.trial), **{str(p.resolve()): w.digest(p) for p in paths}}
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    mesh = np.load(source/'mesh.npz')
    points, cells, surface = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
    fields = np.load(source/'fields.npz')['displacement']
    assert len(points) <= 250000 and len(cells) <= 250000
    space = ht.topology(points, cells, surface)
    modes = np.load(args.trial/'modes.npz')
    for key, stored in (('ids','feature_ids'), ('signs','signs'), ('corner_classes','corner_classes'), ('face_ids','face_ids')):
        np.testing.assert_array_equal(space[key], modes[stored])
    coefficients = modes['coefficients']
    element = Tet10(15)
    print('Assemble consistent 15-point mode forces and energy', flush=True)
    energy, _, geometry, coupling = ht.integrate(points, cells, fields, space, coefficients,
        (element.xi_pg, element.w_pg), return_coupling=True)
    combined = ht.trial(energy)
    work = fields.transpose(0,2,1).reshape(6,-1)@coupling
    np.testing.assert_allclose(work, energy[:6,6:], atol=1e-12, rtol=1e-9)
    assert np.linalg.norm(coupling.reshape(3,len(points),6).sum(axis=1)) < 1e-10
    print('Check consistent mode Gram/gain against positive 125-point integration', flush=True)
    high, _, _ = ht.integrate(points, cells, fields, space, coefficients, duffy_rule(5))
    high_trial = ht.trial(high)
    np.testing.assert_allclose(high_trial['trial_energy_tensor'], checked['rules'][-1]['trial']['trial_energy_tensor'],
        atol=1e-12, rtol=1e-9)
    gram_difference = float(np.linalg.norm(energy[6:,6:]-high[6:,6:])/np.linalg.norm(high[6:,6:]))
    gain_difference = float(np.linalg.norm(np.asarray(combined['gain_tensor'])-high_trial['gain_tensor'])
        /np.linalg.norm(high_trial['gain_tensor']))
    assert max(gram_difference, gain_difference) <= .01
    h, cross, gram = energy[:6,:6], energy[:6,6:], energy[6:,6:]
    equivalent = np.block([[h,h+cross],[h+cross.T,h+cross+cross.T+gram]])
    archive = args.run/'coupling.npz'
    np.savez_compressed(archive, forces=coupling, energy=equivalent, parent_mode_energy=energy)
    w.save(args.run/'consistent.json', dict(passed=True, source=str(source), trial_source=str(args.trial),
        quadrature_points=15, positive_rule_points=125, mode_gram_relative_difference=gram_difference,
        gain_relative_difference=gain_difference, precision_limit=.01, geometry=geometry,
        full_coupling_shape=list(coupling.shape), coupling_bytes=coupling.nbytes,
        coupling_sha256=w.digest(archive), trial=combined,
        extra_scalar_mode_unknowns=6, equilibrium_solve_performed=False))
    captured = {}

    class Capture:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr()
            assert np.isfinite(matrix.data).all() and (matrix.diagonal() > 0).all()
            transform = self.problem._MatCB
            assert np.array_equal(self.problem.space.get_rank_vector('Disp'), np.arange(3))
            full_fields = np.zeros((transform.shape[0],6))
            full_fields[:3*len(points)] = fields.transpose(0,2,1).reshape(6,-1).T
            internal = self.problem.get_A()@full_fields
            residuals = np.linalg.norm(transform.T@internal,axis=0)/np.maximum(np.linalg.norm(internal,axis=0),1e-30)
            assert residuals.max() < 2e-9
            padded = np.zeros_like(full_fields); padded[:len(coupling)] = coupling
            reduced = transform.T@padded
            rng = np.random.default_rng(20261025)
            vector = rng.normal(size=len(rhs))
            np.testing.assert_allclose(vector@reduced, (transform@vector)@padded, atol=1e-10)
            first = fields[0].T.ravel()[self.problem._dof_free]
            first_residual = float(np.linalg.norm(matrix@first-rhs)/np.linalg.norm(rhs))
            assert first_residual < 2e-9
            amg = AuditedAMG(elastic_blocks=True); amg.problem = self.problem
            column = int(np.linalg.norm(reduced,axis=0).argmax())
            with patch('scipy.sparse.linalg.cg', side_effect=CapturedSystem) as intercepted:
                try:
                    amg(matrix,reduced[:,column])
                except CapturedSystem:
                    assert intercepted.call_count == 1 and not amg.records
                else:
                    raise AssertionError('Setup must stop at the CG boundary')
            assert amg.symmetry_key == amg.key
            assert amg.key[1] == row['solver_records'][0]['reduced_matrix_sha256']
            captured.update(reduced_shape=list(matrix.shape), reduced_nnz=int(matrix.nnz),
                csr_bytes=matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes,
                reduced_coupling_shape=list(reduced.shape), reduced_coupling_bytes=reduced.nbytes,
                matrix_sha256=amg.key[1], saved_field_relative_residuals=residuals.tolist(),
                first_saved_field_rhs_residual=first_residual, amg_setup_performed=True,
                amg_levels=len(amg.hierarchy.levels), amg_operator_complexity=amg.hierarchy.operator_complexity(),
                cg_started=False, dense_mode_rows_added_to_sparse_matrix=False)
            raise CapturedSystem()

    print('Measure actual sparse parent reduction and AMG setup; stop before CG', flush=True)
    try:
        solve(points,cells,n_elm_gp=15,solver=Capture(),release_reduced_cache=True)
    except CapturedSystem:
        pass
    else:
        raise AssertionError('Assembly capture must not return a solution')
    assert captured and not captured['cg_started']
    hashes[str(archive.resolve())] = w.digest(archive)
    hashes[str((args.run/'consistent.json').resolve())] = w.digest(args.run/'consistent.json')
    w.save(args.run/'result.json', dict(passed=True, source=str(source), trial=str(args.trial),
        physical_nodes=len(points), physical_tets=len(cells), extra_scalar_mode_unknowns=6,
        consistent_rule_passed=True, mode_gram_relative_difference=gram_difference,
        gain_relative_difference=gain_difference, resources=captured, inputs=hashes,
        new_equilibrium_solve=False, reference_promoted=False, native_candidate=False,
        scope='Consistent-rule and actual sparse-parent/setup feasibility only; no CG response or fake solution'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(f"Consistent coupling and setup pass, {captured['reduced_shape'][0]} free rows, "
        f"{captured['amg_levels']} AMG levels; no CG solve started", flush=True)


if __name__ == '__main__':
    main()
