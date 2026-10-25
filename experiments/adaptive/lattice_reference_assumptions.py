"""Replay all retained lattice extrapolations without changing their admission role."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import microgen_workflow as w

PAIRS = (
    ('g3-lattice-curve-grading-reference-comparison-v1.json',
     'g3-lattice-adapt24-tol003-v1', 'g3-lattice-adapt24-tol0028-v1'),
    ('g3-lattice-reference-sensitivity-v2.json',
     'g3-lattice-curved28-symmetric-v1', 'g3-lattice-curved32-reference-v1'),
    ('g3-lattice-microgen-reference-comparison-v1.json',
     'g3-lattice-adapt24-solve-v1', 'g3-lattice-adapt24-tol003-v1'),
)


def resolve(name):
    path = Path(name.replace('/home/kmarchais/dev/',str(w.ROOT.parent)+'/'))
    return path if path.is_absolute() else w.ROOT/path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_ = w.imports()
    from periodic import tensor_difference
    native_path = w.HERE/'runs/macbook-lattice26-current-native-amg-v1/result.json'
    native = np.asarray(json.loads(native_path.read_text())['stiffness'])
    inputs = {str(p.resolve()):w.digest(p) for p in (Path(__file__),native_path,w.HERE/'adaptive_reference.py')}
    shutil.copyfile(__file__,args.run/'analysis-driver.py')
    families, all_models = [], []
    for filename,coarse_name,fine_name in PAIRS:
        path = w.HERE/'runs'/filename
        report = json.loads(path.read_text())
        assert report['case'] == 'smooth-lattice' and report['stiffness_target'] == .01
        assert len(report['models']) == 8
        inputs[str(path.resolve())] = w.digest(path)
        checked, missing = [], []
        for name,expected in report['inputs'].items():
            p = resolve(name)
            if p.exists():
                assert w.digest(p) == expected, str(p)
                inputs[str(p.resolve())] = expected
                checked.append(name)
            else:
                missing.append(name)
        parents = [w.HERE/'runs'/name for name in (coarse_name,fine_name)]
        rows = [json.loads((p/'result.json').read_text()) for p in parents]
        for directory,row in zip(parents,rows,strict=True):
            assert row['qualified'] and row['geometry']['geometry_admitted'] and len(row['solver_records']) == 6
            assert max(r['relative_residual'] for r in row['solver_records']) < 2e-9
            for name in ('result.json','mesh.npz','fields.npz'):
                p = directory/name
                assert str(p.resolve()) in inputs, f'Required pair evidence unavailable: {p}'
        coarse,fine = (np.asarray(row['stiffness']) for row in rows)
        np.testing.assert_allclose(fine,report['primary_reference_stiffness'],atol=1e-14,rtol=1e-13)
        assert resolve(report['primary_reference_source']).resolve() == (parents[1]/'result.json').resolve()
        models = []
        for index,model in enumerate(report['models']):
            ratio,order = model['length_ratio'],model['assumed_order']
            assert 0 < ratio < 1 and order in (1,2,3,4)
            if model['length_proxy'] == 'element_density':
                np.testing.assert_allclose(ratio,(rows[0]['elements']/rows[1]['elements'])**(1/3),atol=1e-14)
            elif model['length_proxy'] == 'uniform_seed':
                np.testing.assert_allclose(ratio,rows[0]['cells']/rows[1]['cells'],atol=1e-14)
            else:
                assert model['length_proxy'] == 'target_length'
                np.testing.assert_allclose(ratio,(rows[1]['adapt_tolerance']/rows[0]['adapt_tolerance'])**.5,atol=1e-14)
            weight = ratio**order
            replay = (fine-weight*coarse)/(1-weight)
            retained = np.asarray(model['stiffness'],dtype=np.float64)
            np.testing.assert_allclose(replay,retained,atol=1e-14,rtol=1e-13)
            assert np.linalg.eigvalsh(retained).min() > 0
            item = dict(index=index,length_proxy=model['length_proxy'],length_ratio=ratio,assumed_order=order,
                matrix_sha256=hashlib.sha256(retained.tobytes()).hexdigest(),
                fine_coefficient=1/(1-weight),coarse_coefficient=-weight/(1-weight),
                extrapolation_beyond_fine=tensor_difference(retained,fine),
                native_difference=tensor_difference(native,retained))
            models.append(item);all_models.append(dict(source=str(path),**item))
        change = rows[1]['geometry']['volume']/rows[0]['geometry']['volume']-1
        families.append(dict(source=str(path),coarse=str(parents[0]),fine=str(parents[1]),
            provenance_hashes_verified=len(checked),unavailable_optional_provenance=missing,
            rows=[dict(source=str(p),nodes=r['nodes'],tets=r['elements'],geometry=r['geometry'])
                  for p,r in zip(parents,rows,strict=True)],relative_volume_change=change,
            identical_physical_domain_proven=False,exact_domain_difference_detected=abs(change)>1e-10,
            nested_approximation_spaces_proven=False,observed_convergence_order_proven=False,
            pair_difference=tensor_difference(coarse,fine),models=models))
    assert len(all_models) == 24
    worst = max(all_models,key=lambda m:m['native_difference']['worst_relative_strain_energy_difference'])
    w.save(args.run/'analysis-inputs.json',inputs)
    w.save(args.run/'result.json',dict(passed=True,families=families,retained_model_count=24,
        retained_model_tensor_hashes=[m['matrix_sha256'] for m in all_models],worst_native_model=worst,
        all_retained_formulas_reproduced=True,models_removed=0,target_changed=False,
        reference_promoted=False,native_candidate_promoted=False,new_solve=False,inputs=inputs,
        scope='Provenance and extrapolation-assumption audit; all24 retained as sensitivity tests, no continuum certification or relaxation of admission'))
    assert all(w.digest(Path(p)) == h for p,h in inputs.items())
    print('All24 retained extrapolations reproduce; no order/domain-nesting certificate inferred',flush=True)


if __name__ == '__main__':
    main()
