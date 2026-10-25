"""One measured-edge nested preflight on the current immutable native lattice."""

import argparse
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_ = w.imports()
    import hierarchical_residual_trial as ht
    from amg_reference import AuditedAMG
    from capture_reference_system import CapturedSystem
    from hp_coupled_preflight import inherited_inputs
    from hp_shared_domain import periodic_jumps
    from nested_preflight import fixture_inputs
    from nested_refinement import EDGES, refine, interpolation_check
    from periodic import periodic_pairs
    from quadrature_audit import duffy_rule
    from stress_recovery import solve
    source = w.HERE/'runs/macbook-lattice26-current-native-amg-v1'
    trial = w.HERE/'runs/macbook-lattice26-hierarchical-residual-v1'
    checked = json.loads((trial/'result.json').read_text())
    assert checked['passed'] and checked['gain_signal_passed']
    fixtures = fixture_inputs(source)
    paths = [Path(__file__),Path(ht.__file__),trial/'modes.npz',trial/'result.json']
    paths.extend(w.HERE/name for name in ('nested_refinement.py','test_nested_refinement.py',
        'nested_preflight.py','amg_reference.py','stress_recovery.py','hp_shared_domain.py'))
    hashes = {**inherited_inputs(trial), **{str(p.resolve()):w.digest(p) for p in paths}}
    for i,path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path,args.run/f'input-{i}-{path.name}')
    for path in fixtures:
        target = args.run/path.name
        shutil.copyfile(path,target)
        assert w.digest(target) == w.digest(path)
        hashes[str(target.resolve())] = w.digest(target)
    manifest = json.loads((args.run/'manifest.json').read_text())
    manifest['inputs'].update({str((args.run/p.name).resolve()):w.digest(p) for p in fixtures})
    w.save(args.run/'manifest.json',manifest)
    row = json.loads((source/'result.json').read_text())
    mesh = np.load(source/'mesh.npz')
    points,cells,surface = (mesh[k] for k in ('points','tetrahedra','surface'))
    fields = np.load(source/'fields.npz')['displacement']
    space = ht.topology(points,cells,surface)
    baseline = np.asarray(row['stiffness'])
    original = np.load(trial/'modes.npz')
    np.testing.assert_array_equal(space['ids'],original['feature_ids'])
    scores, precision = [], None
    for order in (4,5):
        print(f'Measure shared edge residual potentials with positive{order**3}',flush=True)
        coefficients,residual,construction = ht.build_modes(points,cells,fields,space,duffy_rule(order))
        if order == 4:
            np.testing.assert_allclose(coefficients,original['coefficients'],atol=1e-12,rtol=1e-9)
        scalar = -np.einsum('fcl,fcm,lm->f',residual,coefficients,np.linalg.inv(baseline),optimize=True)
        assert scalar.min() >= -1e-12 and scalar.sum() > 0
        scores.append(np.maximum(scalar,0))
        del coefficients,residual
    precision = float(np.linalg.norm(scores[1]-scores[0])/np.linalg.norm(scores[1]))
    assert precision <= .01
    _,classes,_ = periodic_pairs(points,surface)
    edges = np.unique(np.sort(cells[:,EDGES].reshape(-1,2),axis=1),axis=0)
    canonical,ids = np.unique(np.sort(classes[edges],axis=1),axis=0,return_inverse=True)
    expected = np.unique(np.sort(space['corner_classes'][:,ht.EDGES],axis=-1).reshape(-1,2),axis=0)
    np.testing.assert_array_equal(canonical,expected)
    assert len(canonical) == space['edges']
    edge_potential = scores[-1][:space['edges']]
    physical = edge_potential[ids]/np.bincount(ids)[ids]
    np.testing.assert_allclose(physical.sum(),edge_potential.sum(),rtol=1e-12)
    edge_scores = {tuple(edge):float(value) for edge,value in zip(edges,physical,strict=True)}
    np.savez_compressed(args.run/'edge-indicators.npz',edges=edges,scores=physical,
                        canonical_edges=canonical,canonical_scores=edge_potential)
    print('Bisect complete periodic edge stars ranked by residual potential per added node',flush=True)
    refined = refine(points,cells,surface,fields,edge_scores=edge_scores)
    reproduction = interpolation_check(points,cells,fields,refined)
    periodic = periodic_jumps(refined['points'],refined['surface'],refined['displacement'])
    np.savez_compressed(args.run/'mesh.npz',**{k:refined[k] for k in ('points','tetrahedra','surface')})
    np.savez_compressed(args.run/'prolongation.npz',**{k:refined[k] for k in
        ('displacement','parents','parent_barycentric_maps')})
    geometry = w.geometry(args.run,'smooth-lattice',refined['points'],refined['tetrahedra'],refined['surface'])
    volume_change = geometry['volume']/row['geometry']['volume']-1
    assert geometry['geometry_admitted'] and abs(volume_change) < 1e-10
    assert len(refined['points']) <= 250000 and len(refined['tetrahedra']) <= 250000
    resource = {}

    class Capture:
        def __init__(self):
            self.problem,self.records = None,[]

        def __call__(self,matrix,rhs,**kwargs):
            matrix = matrix.tocsr()
            amg = AuditedAMG(elastic_blocks=True);amg.problem = self.problem
            with patch('scipy.sparse.linalg.cg',side_effect=CapturedSystem) as intercepted:
                try:
                    amg(matrix,rhs)
                except CapturedSystem:
                    assert intercepted.call_count == 1 and not amg.records
                else:
                    raise AssertionError('Stop at CG boundary')
            resource.update(reduced_shape=list(matrix.shape),nnz=int(matrix.nnz),matrix_sha256=amg.key[1],
                levels=len(amg.hierarchy.levels),operator_complexity=amg.hierarchy.operator_complexity(),
                cg_started=False,amg_setup_passed=True)
            raise CapturedSystem()

    print('Measure actual refined sparse operator and AMG setup without solving',flush=True)
    try:
        solve(refined['points'],refined['tetrahedra'],n_elm_gp=15,solver=Capture(),release_reduced_cache=True)
    except CapturedSystem:
        pass
    else:
        raise AssertionError('Capture must not return a fake solution')
    assert resource
    for name in ('mesh.npz','prolongation.npz','edge-indicators.npz'):
        path = args.run/name
        hashes[str(path.resolve())] = w.digest(path)
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(status='nested-preflight-passed',qualified=True,case='smooth-lattice',
        source=str(source),nodes=len(refined['points']),elements=len(refined['tetrahedra']),order=2,
        geometry=geometry,relative_volume_change=volume_change,reproduction=reproduction,periodic=periodic,
        selection=refined['diagnostics'],resources=resource,indicator_relative_rule_change=precision,
        edge_potential_fraction=float(edge_potential.sum()/scores[-1].sum()),
        indicator_scope='Block diagonal cubic edge residual potential, not an error bound or predicted h gain',
        old_quadratic_space_preserved=True,original_curved_domain_preserved=True,
        equilibrium_solve_performed=False,reference_promoted=False,native_candidate_promoted=False,
        comparative_timing_claim=False,inputs=hashes,
        scope='One data-directed fixed-domain Tet10 reference preflight; no new sizing field or native workflow claim'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(f'Preflight passes: {len(refined["points"])} nodes, {len(refined["tetrahedra"])} tets; CG not started',flush=True)


if __name__ == '__main__':
    main()
