"""Inspect signed periodic junctions and the current remaining trial, without meshing."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def opening_angles(low, high, axis):
    """Material is below each wall branch; low is the positive-axis branch."""
    import numpy as np

    tangent_low, tangent_high = low.copy(), high.copy()
    tangent_low[:, axis] = tangent_high[:, axis] = 0
    a, b = np.linalg.norm(tangent_low, axis=1), np.linalg.norm(tangent_high, axis=1)
    valid = (a > 1e-8) & (b > 1e-8)
    alignment = np.einsum('ij,ij->i', tangent_low, tangent_high)/np.maximum(a*b, 1e-30)
    valid &= alignment > 1-1e-5
    angles = np.degrees(np.pi+np.arctan2(high[:, axis], b)-np.arctan2(low[:, axis], a))
    return angles, valid


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--remaining', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    from periodic_crease_audit import features, gradient, project_box_curve
    from hp_coupled_preflight import inherited_inputs

    low = np.array([[-1., 1, 0], [1., 1, 0], [0., 1, 0]])
    high = np.array([[1., 1, 0], [-1., 1, 0], [0., 1, 0]])
    angles, valid = opening_angles(low, high, 0)
    np.testing.assert_allclose(angles, [270, 90, 180], atol=1e-12)
    assert valid.all()
    np.testing.assert_allclose(opening_angles(3*low, 2*high, 0)[0], angles, atol=1e-12)
    remaining = json.loads((args.remaining/'result.json').read_text())
    assert remaining['passed'] and remaining['parent_field_includes_all_bubbles']
    assert Path(remaining['source']).resolve() == args.source.resolve()
    paths = [Path(__file__), Path(w.__file__), w.HERE/'periodic_crease_audit.py',
        w.HERE/'scorecard.py', args.remaining/'result.json', args.source/'mesh.npz',
        args.source/'smooth-lattice-specimen.json', args.source/'smooth-lattice-samples.npz']
    hashes = {**inherited_inputs(args.remaining), **{str(p.resolve()): w.digest(p) for p in paths}}
    mesh = np.load(args.source/'mesh.npz')
    feature = features(mesh)
    case, _, _ = w.specimen('smooth-lattice')
    points, projection = project_box_curve(case.function, mesh['points'][feature['curve_nodes']].copy())
    on_plane = np.isclose(abs(points), .5, atol=1e-10, rtol=0)
    rows = []
    for axis in range(3):
        selected = (on_plane.sum(axis=1) == 1) & np.isclose(points[:, axis], -.5, atol=1e-10, rtol=0)
        low_points = points[selected]
        assert len(low_points)
        high_points = low_points.copy(); high_points[:, axis] = .5
        mismatch = float(abs(case.function(*low_points.T)-case.function(*high_points.T)).max())
        assert mismatch < 1e-9
        estimates, masks = [], []
        for epsilon in (2e-6, 1e-6):
            first, second = low_points.copy(), high_points.copy()
            first[:, axis] += epsilon; second[:, axis] -= epsilon
            angles, valid = opening_angles(gradient(case.function, first, epsilon/8),
                gradient(case.function, second, epsilon/8), axis)
            estimates.append(angles); masks.append(valid)
        valid = masks[0] & masks[1]
        assert valid.any()
        values = estimates[-1][valid]
        rows.append(dict(axis=axis, examined_single_plane_points=len(low_points), valid_points=int(valid.sum()),
            excluded_tangentially_inconsistent_points=int((~valid).sum()), scalar_trace_mismatch=mismatch,
            opening_angle_quantiles_deg=np.quantile(values, [0, .1, .5, .9, 1]).tolist(),
            reentrant_above_190_deg=int((values > 190).sum()), convex_below_170_deg=int((values < 170).sum()),
            maximum_two_epsilon_change_deg=float(abs(estimates[1][valid]-estimates[0][valid]).max())))
    local = np.load(args.remaining/'local-125.npz')['correction_energy']
    gain = np.asarray(remaining['rules'][-1]['trial']['gain_tensor'])
    np.testing.assert_allclose(local.sum(axis=0), gain, atol=1e-12, rtol=1e-8)
    assert np.linalg.eigvalsh(local).min() >= -1e-12
    parent = np.asarray(remaining['rules'][-1]['trial']['trial_energy_tensor'])+gain
    scalar = np.einsum('ij,tji->t', np.linalg.inv(parent), local)
    localization = {name: dict(elements=int(mask.sum()), energy_fraction=float(scalar[mask].sum()/scalar.sum()))
        for name, mask in feature['groups'].items()}
    for p in (args.remaining/'local-125.npz',):
        hashes[str(p.resolve())] = w.digest(p)
    shutil.copyfile(__file__, args.run/'analysis-driver.py')
    w.save(args.run/'result.json', dict(passed=True, analytic_roof_valley_plane_controls_passed=True,
        source=str(args.source), remaining=str(args.remaining), projection=projection, axes=rows,
        remaining_trial_localization=localization,
        additional_trial_gain=remaining['maximum_additional_trial_gain'],
        existing_gain_signal_passed=remaining['gain_signal_passed'], new_meshing_pilot_selected=False,
        reference_promoted=False, new_solve=False, meshes_changed=False, inputs=hashes,
        scope='Signed local wall-junction geometry and resolved remaining-trial energy only; no 3D singularity proof, error bound, revived sizing recipe or native winner'))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    print('Signed junction and remaining-trial diagnostics complete; no meshing pilot selected', flush=True)


if __name__ == '__main__':
    main()
