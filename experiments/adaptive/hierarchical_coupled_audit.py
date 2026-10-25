"""Audit complete bubble/P2 coupled fields with fresh quadrature and assembly."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--solve', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    from amg_reference import reduced_matrix_key
    from capture_reference_system import CapturedSystem
    from homogenization import strain_tensor
    from hp_coupled_preflight import inherited_inputs
    import hierarchical_residual_trial as ht
    from periodic import periodic_pairs, tensor_difference
    from quadrature_audit import duffy_rule
    from stress_recovery import solve

    solved = json.loads((args.solve/'result.json').read_text())
    assert solved['status'] == 'coupled-solved-needs-independent-audit'
    source, trial, preflight = (Path(solved[k]) for k in ('source','trial','preflight'))
    archive = args.solve/'coupled-fields.npz'
    assert w.digest(archive) == solved['field_archive_sha256']
    paths = [Path(__file__), Path(ht.__file__), args.solve/'result.json', archive,
             preflight/'coupling.npz', preflight/'consistent.json', trial/'result.json']
    paths.extend(w.HERE/name for name in ('quadrature_audit.py','stress_recovery.py',
        'amg_reference.py','hp_coupled_preflight.py','microgen_workflow.py'))
    hashes = {**inherited_inputs(args.solve), **{str(p.resolve()):w.digest(p) for p in paths}}
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path,args.run/f'input-{i}-{path.name}')
    mesh = np.load(source/'mesh.npz')
    points, cells, surface = (mesh[k] for k in ('points','tetrahedra','surface'))
    fields = np.load(source/'fields.npz')['displacement']
    corrected = np.load(archive)
    responses, alpha = corrected['response_displacement'], corrected['mode_coefficients']
    assert responses.shape == fields.shape == (6,len(points),3)
    assert len(points) <= 250000 and len(cells) <= 250000
    space = ht.topology(points,cells,surface)
    modes = np.load(trial/'modes.npz')
    for key, stored in (('ids','feature_ids'),('signs','signs'),('corner_classes','corner_classes'),('face_ids','face_ids')):
        np.testing.assert_array_equal(space[key],modes[stored])
    traces = ht.trace_check(space,modes['coefficients'])
    pairs, _, periodic = periodic_pairs(points,surface)
    jump = max(float(abs(responses[:,p[:,1]]-responses[:,p[:,0]]).max(initial=0)) for p in pairs)
    pin = int(np.linalg.norm(points,axis=1).argmin())
    assert jump < 1e-9 and abs(responses[:,pin]).max() < 1e-12
    predicted = np.asarray(solved['predicted_energy_tensor'])
    basis = np.vstack((np.eye(6),alpha))
    rules = []
    previous = None
    element = Tet10(15)
    for count, rule in ((15,(element.xi_pg,element.w_pg)),(64,duffy_rule(4)),(125,duffy_rule(5))):
        print(f'Integrate complete coupled fields with {count} positive points',flush=True)
        energy, local, geometry = ht.integrate(points,cells,fields,space,modes['coefficients'],rule,
            responses=responses,return_stress=True)
        tensor = basis.T@energy@basis
        mean_stress = np.asarray(geometry.pop('mean_stress_matrix'))@basis
        difference = tensor_difference(tensor,predicted)
        stress_error = float(np.linalg.norm(mean_stress-tensor)/np.linalg.norm(tensor))
        assert difference['worst_relative_strain_energy_difference'] < (1e-7 if count == 15 else 1e-4)
        assert stress_error < (1e-7 if count == 15 else 1e-4)
        correction = np.einsum('mi,tmn,nj->tij',alpha,local,alpha,optimize=True)
        if count == 15:
            balance = float(np.linalg.norm(correction.sum(axis=0)-(energy[:6,:6]-tensor))
                            /np.linalg.norm(correction.sum(axis=0)))
            assert balance < 1e-7
        change = None if previous is None else tensor_difference(tensor,previous)
        if change is not None:
            assert change['worst_relative_strain_energy_difference'] < 1e-4
        np.savez_compressed(args.run/f'local-{count}.npz',correction_energy=correction)
        rules.append(dict(points=count,geometry=geometry,energy_tensor=tensor.tolist(),
            energy_vs_solved=difference,mean_stress_relative_error=stress_error,vs_previous=change))
        previous = tensor
    data = np.load(preflight/'coupling.npz')
    forces, original_energy = data['forces'], data['parent_mode_energy']
    captured = {}

    class Capture:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self,matrix,rhs,**kwargs):
            matrix = matrix.tocsr()
            assert reduced_matrix_key(matrix)[1] == solved['solver_records'][0]['reduced_matrix_sha256']
            transform = self.problem._MatCB
            assert np.array_equal(self.problem.space.get_rank_vector('Disp'),np.arange(3))
            flat = fields.transpose(0,2,1).reshape(6,-1).T
            z = responses.transpose(0,2,1).reshape(6,-1).T
            padded = np.zeros((transform.shape[0],6)); padded[:len(forces)] = forces
            reduced = transform.T@padded
            free = z[self.problem._dof_free]
            residual = np.linalg.norm(matrix@free-reduced,axis=0)/np.linalg.norm(reduced,axis=0)
            assert residual.max() < 2e-9
            np.testing.assert_allclose((transform@free)[:len(z)],z,atol=1e-12,rtol=1e-10)
            nodal = flat-z@alpha
            full = np.zeros_like(padded); full[:len(nodal)] = nodal
            internal = self.problem.get_A()@full+padded@alpha
            coarse = np.linalg.norm(transform.T@internal,axis=0)/np.linalg.norm(internal,axis=0)
            assert coarse.max() < 2e-9
            modal = forces.T@nodal+original_energy[6:,6:]@alpha
            modal_error = float(np.linalg.norm(modal)/np.linalg.norm(original_energy[6:,:6]))
            assert modal_error < 1e-7
            affine = np.array([points@strain_tensor(e).T for e in np.eye(6)]).transpose(0,2,1).reshape(6,-1)
            reaction = affine@internal[:len(nodal)]
            reaction_error = float(np.linalg.norm(reaction-predicted)/np.linalg.norm(predicted))
            assert reaction_error < 1e-7
            captured.update(matrix_sha256=reduced_matrix_key(matrix)[1],nnz=int(matrix.nnz),
                response_residuals=residual.tolist(),coarse_equilibrium_residuals=coarse.tolist(),
                modal_relative_residual=modal_error,reaction_relative_error=reaction_error,
                fresh_operator=True,cg_started=False)
            raise CapturedSystem()

    print('Freshly assemble and audit all coarse and modal equilibrium rows',flush=True)
    try:
        solve(points,cells,n_elm_gp=15,solver=Capture(),release_reduced_cache=True)
    except CapturedSystem:
        pass
    else:
        raise AssertionError('Audit must not return a fake solution')
    assert captured
    models = []
    for model in json.loads((trial/'result.json').read_text())['models']:
        path = Path(model['source'])
        reference = np.asarray(json.loads(path.read_text())['models'][model['index']]['stiffness'])
        models.append(dict(source=str(path),index=model['index'],difference=tensor_difference(predicted,reference)))
    assert len(models) == 24
    worst = max(m['difference']['worst_relative_strain_energy_difference'] for m in models)
    for count in (15,64,125):
        path = args.run/f'local-{count}.npz'
        hashes[str(path.resolve())] = w.digest(path)
    w.save(args.run/'result.json',dict(passed=True,source=str(source),solve=str(args.solve),
        geometry_unchanged=True,old_approximation_space_preserved=True,
        physical_nodes=len(points),physical_tets=len(cells),complete_field_representation=solved['field_representation'],
        periodic=periodic,response_jump_max=jump,cubic_trace_check=traces,rules=rules,
        galerkin_energy_balance_relative_error=balance,equilibrium=captured,
        retained_model_count=24,models=models,worst_model_disagreement=worst,
        passes_one_percent_all_models=worst <= .01,reference_promoted=False,
        native_candidate_promoted=False,comparative_timing_claim=False,inputs=hashes,
        scope='Same-domain six-mode enriched FE diagnostic; not a native mesh workflow or continuum/local stress certificate'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(f'Complete field audit passes; all-model disagreement {worst:.9%}; no reference promoted',flush=True)


if __name__ == '__main__':
    main()
