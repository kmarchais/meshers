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
