"""Independent algebra, localization and sampled grouped-integration replay."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh
import microgen_workflow as workflow


def verify(run, result):
    workflow.imports()
    from hp_parent25_probe import PREPARED, H_SOLVE, P_MESH, P_SOLVE
    from hp_cross_energy import integrate
    from hp_batched_energy import integrate_batched
    from quadrature_audit import duffy_rule
    from cubic_refinement import MODEL_REPORTS
    from periodic_crease_audit import features

    assert result['passed'] and not result['new_equilibrium_solve']
    assert not result['reference_promoted'] and not result['native_candidate_promoted']
    benchmark = json.loads((Path(result['benchmark'])/'result.json').read_text())
    shared = json.loads((Path(benchmark['preflight'])/'result.json').read_text())
    assert shared['passed'] and shared['parent25_shared_domain_preflight'] and benchmark['passed']
    mesh, cubic = (dict(np.load(path/'mesh.npz')) for path in (H_SOLVE, P_MESH))
    h, p = (np.load(path/'fields.npz')['displacement'] for path in (H_SOLVE, P_SOLVE))
    maps = dict(np.load(PREPARED/'prolongation.npz'))
    groups = features(mesh)['groups']
    selected = np.linspace(0, len(mesh['tetrahedra'])-1, 64, dtype=int)
    canonical = [(str(workflow.HERE/'runs'/name), index, np.asarray(model['stiffness']))
        for name in MODEL_REPORTS for index, model in enumerate(json.loads((workflow.HERE/'runs'/name).read_text())['models'])]
    assert len(canonical) == 24
    for record, order in zip(result['rules'], (5, 6), strict=True):
        assert record['positive_rule_points'] == order**3
        energy = np.asarray(record['energy_matrix'])
        coefficients = np.asarray(record['coefficients'])
        base, cross, cubic_energy = energy[:6,:6], energy[:6,6:], energy[6:,6:]
        gram = base+cubic_energy-cross-cross.T
        residual = cross.T-base
        np.testing.assert_allclose(gram @ coefficients, -residual, atol=1e-12)
        gain = coefficients.T @ gram @ coefficients
        trial = np.asarray(record['trial_energy_tensor'])
        np.testing.assert_allclose(base-trial, gain, atol=1e-12)
        basis = np.vstack((np.eye(6)-coefficients, coefficients))
        np.testing.assert_allclose(basis.T @ energy @ basis, trial, atol=1e-12)
        assert eigvalsh(trial).min() > 0 and eigvalsh(gain).min() > -1e-12
        archive = run/f'local-rule-{order**3}.npz'
        assert workflow.digest(archive) == record['local_archive_sha256']
        local = np.load(archive)
        np.testing.assert_allclose(local['mode_gram'].sum(axis=0), gram, atol=1e-12)
        local_gain = np.einsum('ij,tik,kl->tjl', coefficients, local['mode_gram'], coefficients, optimize=True)
        np.testing.assert_allclose(local_gain.sum(axis=0), gain, atol=1e-12)
        scores = np.einsum('ij,tji->t', np.linalg.inv(base), local_gain)
        np.testing.assert_allclose(scores, local['normalized_gain'], atol=1e-14)
        np.testing.assert_allclose(local['volumes'].sum(), record['volume'], atol=1e-12)
        for name, mask in groups.items():
            saved = record['localization'][name]
            assert saved['tets'] == int(mask.sum())
            np.testing.assert_allclose(saved['normalized_gain_fraction'], scores[mask].sum()/scores.sum(), atol=1e-12)
            np.testing.assert_allclose(saved['volume_fraction'], local['volumes'][mask].sum()/record['volume'], atol=1e-12)
        worst = 0.
        for saved, (source, index, tensor) in zip(record['models'], canonical, strict=True):
            assert (saved['source'], saved['index']) == (source, index)
            eigenvalues = eigvalsh(trial-tensor, tensor)
            np.testing.assert_allclose(eigenvalues, saved['signed_relative_energy_eigenvalues'], atol=1e-12)
            worst = max(worst, float(abs(eigenvalues).max()))
        np.testing.assert_allclose(worst, record['worst_retained_model_disagreement'], atol=1e-12)
        arguments = (mesh['points'], mesh['tetrahedra'][selected], h, cubic['tetrahedra'], p,
            maps['parents'][selected], maps['parent_barycentric_maps'][selected], duffy_rule(order))
        grouped = integrate(*arguments)
        batched = integrate_batched(*arguments, batch_size=17)
        np.testing.assert_allclose(grouped[0], batched[0], atol=1e-12, rtol=1e-10)
        np.testing.assert_allclose(grouped[1:3], batched[1:3], atol=1e-12, rtol=1e-10)
        np.testing.assert_allclose(batched[3], local['mode_gram'][selected], atol=1e-12, rtol=1e-10)
    precision = np.linalg.norm(np.asarray(result['rules'][1]['gain_tensor'])-result['rules'][0]['gain_tensor'])/np.linalg.norm(result['rules'][1]['gain_tensor'])
    np.testing.assert_allclose(precision, result['gain_rule_relative_change'], atol=1e-12)
    assert precision <= .01
    return dict(run=run.name, model_comparisons=48, grouped_full_rule_sample_cells=64,
        worst_disagreement=result['rules'][-1]['worst_retained_model_disagreement'], reference_promoted=False)


def verify_coupling_resource(run, result):
    from hp_parent25_probe import H_SOLVE
    assert result['parent25'] and not result['new_equilibrium_solve'] and not result['reference_promoted']
    assert result['nested_source'] == str(H_SOLVE)
    mesh = np.load(H_SOLVE/'mesh.npz')
    assert result['physical_nodes'] == len(mesh['points']) <= 250000
    assert result['physical_tets'] == len(mesh['tetrahedra']) <= 250000
    worker = json.loads((run/'worker.json').read_text())
    assert worker['returncode'] == 0 and worker['peak_rss_bytes'] < workflow.MEMORY
    if 'coupling_sha256' in result:
        assert workflow.digest(run/'coupling.npz') == result['coupling_sha256']
        data = np.load(run/'coupling.npz')
        energy, coupling = data['energy'], data['forces']
        fields = np.load(H_SOLVE/'fields.npz')['displacement']
        work = fields.transpose(0,2,1).reshape(6,-1) @ coupling
        np.testing.assert_allclose(work, energy[:6,6:]-energy[:6,:6], atol=1e-12, rtol=1e-9)
        assert np.linalg.norm(coupling.reshape(3,len(mesh['points']),6).sum(axis=1)) < 1e-10
        gram = energy[:6,:6]+energy[6:,6:]-energy[:6,6:]-energy[6:,:6]
        np.testing.assert_allclose(gram, data['mode_matrix'], atol=1e-12)
        np.testing.assert_allclose(eigvalsh(gram), result['mode_matrix_eigenvalues'], atol=1e-12)
        trial = json.loads((Path(result['source'])/'result.json').read_text())['rules'][-1]
        old = np.asarray(trial['energy_matrix'])
        old_gram = old[:6,:6]+old[6:,6:]-old[:6,6:]-old[6:,:6]
        difference = np.linalg.norm(gram-old_gram)/np.linalg.norm(old_gram)
        np.testing.assert_allclose(difference,result['mode_gram_relative_difference'],atol=1e-12)
        assert difference <= .01 and result['gain_relative_difference'] <= .01
        return dict(run=run.name, coupling_virtual_work_replayed=True, mode_gram_relative_difference=float(difference))
    assert result['no_solve'] and result['complete_symmetry_check_passed'] and result['constraint_virtual_work_passed']
    assert result['saved_first_load_reduced_relative_residual'] < 2e-9
    if result['amg_setup_performed']:
        assert not result['cg_started'] and result['amg_levels'] > 1
    return dict(run=run.name, resource_capture_record_checked=True, independent_matrix_reassembly=False,
        peak_rss_bytes=worker['peak_rss_bytes'])


def verify_modal_fields(run, result):
    workflow.imports()
    from hp_cross_energy import integrate
    from hp_batched_energy import integrate_batched
    from hp_coupled_audit import cubic_nodal_values
    from hp_shared_domain import periodic_jumps
    from quadrature_audit import duffy_rule

    source = Path(result['source'])
    solved = json.loads((source/'result.json').read_text())
    assert workflow.digest(source/'hp-fields.npz') == solved['field_archive_sha256']
    mesh, cubic = np.load(Path(solved['nested_source'])/'mesh.npz'), np.load(solved['cubic_mesh'])
    h = np.load(Path(solved['nested_source'])/'fields.npz')['displacement']
    p = np.load(Path(solved['cubic_source'])/'fields.npz')['displacement']
    fields = np.load(source/'hp-fields.npz')
    response, coefficients = fields['response_displacement'], fields['mode_coefficients']
    maps = np.load(Path(solved['prepared'])/'prolongation.npz')
    nodal_p = cubic_nodal_values(mesh['points'], mesh['tetrahedra'], cubic['tetrahedra'], p,
        maps['parents'], maps['parent_barycentric_maps'])
    pin = int(np.linalg.norm(mesh['points'],axis=1).argmin())
    gauge = nodal_p[:,pin]-h[:,pin]
    nodal = h+np.einsum('ij,ind->jnd',coefficients,nodal_p-h-response-gauge[:,None])
    np.testing.assert_allclose(nodal,np.load(run/'nodal-values.npz')['displacement'],atol=1e-12)
    assert abs(nodal[:,pin]).max() < 1e-12
    assert periodic_jumps(mesh['points'],mesh['surface'],nodal) == result['periodic']
    selected = np.linspace(0,len(mesh['tetrahedra'])-1,64,dtype=int)
    basis = np.vstack((np.eye(6)-coefficients,coefficients,-coefficients))
    summed = np.vstack((np.eye(6),np.eye(6)))
    for order in (5,6):
        args = (mesh['points'],mesh['tetrahedra'][selected],h,cubic['tetrahedra'],p,
            maps['parents'][selected],maps['parent_barycentric_maps'][selected],duffy_rule(order))
        energy, volume, minimum, mean = integrate(*args, additional_fields=response, return_stress=True)
        quadratic = h-np.einsum('ij,ind->jnd',coefficients,h+response)
        cubic_part = np.einsum('ij,ind->jnd',coefficients,p)
        folded, measured_volume, measured_minimum, _, _, stress = integrate_batched(
            args[0],args[1],quadratic,args[3],cubic_part,*args[5:],batch_size=17,return_stress=True)
        np.testing.assert_allclose(summed.T@folded@summed,basis.T@energy@basis,atol=1e-12,rtol=1e-10)
        np.testing.assert_allclose(stress@summed,mean@basis,atol=1e-12,rtol=1e-10)
        np.testing.assert_allclose([volume,minimum],[measured_volume,measured_minimum],atol=1e-12)
    return dict(run=run.name, full_nodal_field_reconstructed=True, positive_rule_grouped_sample_cells=64,
        unfolded_energy_and_stress_replayed=True)
