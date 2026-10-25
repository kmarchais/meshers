"""One exact variational split ranking, fixed-domain trial and resource preflight."""

import argparse
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import microgen_workflow as w


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path,
        default=w.HERE/'runs/macbook-periodic-edge-normal-lattice26-amg-v1')
    parser.add_argument('--bulk-fraction',type=float,default=1.,
        help='Retain this fraction of the admissible disjoint-star normalized gain')
    args=parser.parse_args()
    args.source=args.source.resolve()
    assert 0<args.bulk_fraction<=1
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_=w.imports()
    from fedoo.lib_elements.tetrahedron import Tet4,Tet10
    from amg_reference import AuditedAMG
    from capture_reference_system import CapturedSystem
    from hp_shared_domain import periodic_jumps
    from nested_preflight import fixture_inputs
    from nested_refinement import EDGES,child_maps,refine,interpolation_check
    from periodic import periodic_pairs,tensor_difference
    from quadrature_audit import duffy_rule
    from split_patch_energy import local_blocks,patch_structure,assemble_patches,solve_patches
    from stress_recovery import solve
    from tet10_integration import integrated_tensors
    source=args.source
    fixtures=fixture_inputs(source)
    paths=[Path(__file__),*(source/name for name in ('result.json','mesh.npz','fields.npz','manifest.json'))]
    paths.extend(w.HERE/name for name in ('split_patch_energy.py','test_split_patch_energy.py',
        'nested_refinement.py','nested_preflight.py','amg_reference.py','stress_recovery.py',
        'tet10_integration.py','quadrature_audit.py','hp_shared_domain.py','hierarchical_residual_trial.py'))
    hashes={str(p.resolve()):w.digest(p) for p in paths}
    for i,path in enumerate(paths):
        if path.suffix=='.py':shutil.copyfile(path,args.run/f'input-{i}-{path.name}')
    for path in fixtures:
        target=args.run/path.name;shutil.copyfile(path,target)
        assert w.digest(path)==w.digest(target)
        hashes[str(path.resolve())]=w.digest(path)
        hashes[str(target.resolve())]=w.digest(target)
    manifest=json.loads((args.run/'manifest.json').read_text())
    manifest['inputs'].update(hashes)
    w.save(args.run/'manifest.json',manifest)
    row=json.loads((source/'result.json').read_text())
    assert row['qualified'] and row['completed_loads']==6
    mesh=np.load(source/'mesh.npz')
    points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
    fields=np.load(source/'fields.npz')['displacement']
    baseline=np.asarray(row['stiffness'])
    _,classes,_=periodic_pairs(points,surface)
    structure=patch_structure(cells,classes)
    element=Tet10(15);rule=element.xi_pg,element.w_pg
    construction=[]

    def blocks():
        for offset in range(6):
            print(f'Integrate exact split blocks {offset+1}/6',flush=True)
            gram,residual,record=local_blocks(points,cells,fields,offset,rule)
            construction.append(record)
            yield offset,gram,residual

    grams,residuals=assemble_patches(structure,blocks())
    coefficients,gains=solve_patches(grams,residuals)
    del grams,residuals
    scalar=np.einsum('ij,fji->f',np.linalg.inv(baseline),gains)
    assert scalar.min()>=-1e-12 and np.isfinite(scalar).all()
    scalar=np.maximum(scalar,0)
    edges=np.unique(np.sort(cells[:,EDGES].reshape(-1,2),axis=1),axis=0)
    canonical,ids=np.unique(np.sort(classes[edges],axis=1),axis=0,return_inverse=True)
    np.testing.assert_array_equal(canonical,structure['edges'])
    physical=scalar[ids]/np.bincount(ids)[ids]
    edge_scores={tuple(edge):float(value) for edge,value in zip(edges,physical,strict=True)}
    print('Select disjoint complete periodic stars by exact gain per added node',flush=True)
    refined=refine(points,cells,surface,fields,edge_scores=edge_scores)
    mapping_offsets={tuple(mapping.ravel()):offset for offset in range(6) for mapping in child_maps(offset)}

    def selected_groups(candidate):
        groups=set()
        for parent,mapping in zip(candidate['parents'],candidate['parent_barycentric_maps'],strict=True):
            if np.array_equal(mapping,np.eye(4)):continue
            offset=mapping_offsets[tuple(mapping.ravel())]
            groups.add(int(structure['ids'][parent,offset]))
        return groups

    selected=selected_groups(refined)
    bulk=dict(requested_fraction=args.bulk_fraction,unpruned_nodes=len(refined['points']),
        unpruned_elements=len(refined['tetrahedra']),unpruned_orbits=len(selected),
        unpruned_normalized_gain=float(scalar[sorted(selected)].sum()))
    if args.bulk_fraction<1:
        costs={}
        for group in selected:
            copies={}
            for parent,offset in structure['incidence'][group]:
                edge=tuple(sorted(cells[parent,EDGES[offset]]))
                copies.setdefault(edge,[]).append(parent)
            costs[group]=sum(2+len(set(cells[parents,:4].ravel())-set(edge)) for edge,parents in copies.items())
        ranked=sorted(selected,key=lambda group:(-scalar[group]/costs[group],group))
        cumulative=np.cumsum(scalar[ranked])
        count=int(np.searchsorted(cumulative,args.bulk_fraction*cumulative[-1]))+1
        retained=set(ranked[:count])
        filtered={tuple(edge):float(value) for edge,value,group in zip(edges,physical,ids,strict=True) if group in retained}
        print(f'Retain {count}/{len(selected)} admissible stars for fixed gain fraction {args.bulk_fraction}',flush=True)
        refined=refine(points,cells,surface,fields,edge_scores=filtered)
        selected=selected_groups(refined)
        assert selected==retained
    bulk.update(retained_orbits=len(selected),
        retained_normalized_gain_fraction=float(scalar[sorted(selected)].sum()/bulk['unpruned_normalized_gain']))
    reproduction=interpolation_check(points,cells,fields,refined)
    assert len(selected)==refined['diagnostics']['accepted_edge_orbits']
    incident=[parent for group in selected for parent,_ in structure['incidence'][group]]
    assert len(set(incident))==len(incident)==refined['diagnostics']['bisected_parent_elements']
    gain=gains[sorted(selected)].sum(axis=0)
    updates={}
    for group in selected:
        a,b=structure['edges'][group]
        keys=[((a,3),(b,1)),((a,1),(b,3))]
        keys.extend(tuple(sorted(((a,1),(b,1),(c,2)))) for c in structure['rings'][group])
        for key,value in zip(keys,coefficients[group].reshape(-1,3,6),strict=True):
            assert key not in updates
            updates[key]=value.T
    delta=np.zeros_like(refined['displacement'])
    nodes,first=np.unique(refined['tetrahedra'].ravel(),return_index=True)
    new=nodes>=len(points);nodes,first=nodes[new],first[new]
    children,local=np.divmod(first,10)
    bary=Tet4().shape_function(Tet10().xi_nd)
    for node,child,local_node in zip(nodes,children,local,strict=True):
        parent=refined['parents'][child]
        key_weights=np.rint(4*(bary[local_node]@refined['parent_barycentric_maps'][child])).astype(int)
        assert key_weights.sum()==4
        key=tuple(sorted((int(c),int(v)) for c,v in zip(structure['corners'][parent],key_weights,strict=True) if v))
        delta[:,node]=updates[key]
    assert abs(delta[:,:len(points)]).max()==0
    trial=refined['displacement']+delta
    pairs,_,_=periodic_pairs(refined['points'],refined['surface'])
    delta_jump=max(float(abs(delta[:,p[:,1]]-delta[:,p[:,0]]).max(initial=0)) for p in pairs)
    assert delta_jump<1e-10
    periodic=periodic_jumps(refined['points'],refined['surface'],trial)
    print('Verify full trial energy against the sum of independent split gains',flush=True)
    args_energy=(refined['points'],refined['tetrahedra'])
    _,old_energy,_=integrated_tensors(*args_energy,refined['displacement'])
    _,trial_energy,_=integrated_tensors(*args_energy,trial)
    np.testing.assert_allclose(old_energy,baseline,atol=1e-9,rtol=1e-7)
    np.testing.assert_allclose(old_energy-trial_energy,gain,atol=1e-10,rtol=1e-7)
    balance=float(np.linalg.norm(old_energy-trial_energy-gain)/np.linalg.norm(gain))
    precision=[];previous=None
    for order in (4,5):
        _,energy,volume=integrated_tensors(*args_energy,trial,duffy_rule(order))
        change=None if previous is None else tensor_difference(energy,previous)
        if change is not None:assert change['worst_relative_strain_energy_difference']<1e-4
        _,delta_energy,_=integrated_tensors(*args_energy,delta,duffy_rule(order))
        gain_error=float(np.linalg.norm(delta_energy-gain)/np.linalg.norm(gain))
        assert gain_error<.01
        precision.append(dict(points=order**3,volume=volume,vs_previous=change,gain_relative_error=gain_error))
        previous=energy
    for name,data in (('mesh',dict((k,refined[k]) for k in ('points','tetrahedra','surface'))),
        ('prolongation',dict((k,refined[k]) for k in ('displacement','parents','parent_barycentric_maps'))),
        ('split-gains',dict(canonical_edges=canonical,gains=gains,scalar=scalar,selected=np.asarray(sorted(selected)))),
        ('trial-fields',dict(displacement=trial,delta=delta))):
        np.savez_compressed(args.run/f'{name}.npz',**data)
    geometry=w.geometry(args.run,'smooth-lattice',*args_energy,refined['surface'])
    volume_change=geometry['volume']/row['geometry']['volume']-1
    assert geometry['geometry_admitted'] and abs(volume_change)<1e-10
    gain_eigenvalues=np.linalg.eigvalsh(np.linalg.solve(np.linalg.cholesky(baseline),
        np.linalg.solve(np.linalg.cholesky(baseline),gain).T).T)
    assert gain_eigenvalues.min()>-1e-10
    assert gain_eigenvalues.max()>=.0005, 'Read trial before spending on an unresolved solve'
    resources={}

    class Capture:
        def __init__(self):self.problem,self.records=None,[]
        def __call__(self,matrix,rhs,**kwargs):
            matrix=matrix.tocsr();amg=AuditedAMG(elastic_blocks=True);amg.problem=self.problem
            with patch('scipy.sparse.linalg.cg',side_effect=CapturedSystem) as intercepted:
                try:amg(matrix,rhs)
                except CapturedSystem:assert intercepted.call_count==1 and not amg.records
                else:raise AssertionError('Stop at CG boundary')
            resources.update(reduced_shape=list(matrix.shape),nnz=int(matrix.nnz),matrix_sha256=amg.key[1],
                levels=len(amg.hierarchy.levels),operator_complexity=amg.hierarchy.operator_complexity(),
                cg_started=False,amg_setup_passed=True)
            raise CapturedSystem()

    print('Measure reduced quadratic operator and AMG setup, stopping before CG',flush=True)
    try:solve(*args_energy,n_elm_gp=15,solver=Capture(),release_reduced_cache=True)
    except CapturedSystem:pass
    else:raise AssertionError('Capture must not return a fake solution')
    assert resources
    for name in ('mesh','prolongation','split-gains','trial-fields'):
        path=args.run/f'{name}.npz';hashes[str(path.resolve())]=w.digest(path)
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(status='nested-preflight-passed',qualified=True,case='smooth-lattice',
        source=str(source),nodes=len(refined['points']),elements=len(refined['tetrahedra']),order=2,
        geometry=geometry,relative_volume_change=volume_change,reproduction=reproduction,periodic=periodic,
        selection=refined['diagnostics'],bulk=bulk,resources=resources,construction=construction,
        trial_gain=gain.tolist(),trial_gain_eigenvalues=gain_eigenvalues.tolist(),trial_energy_balance=balance,
        trial_delta_periodic_jump=delta_jump,trial_precision=precision,trial_globally_equilibrated=False,
        selection_scope='Exact fixed-old-node six-load split gains per physical added-node upper bound; all regions eligible',
        old_quadratic_space_preserved=True,original_curved_domain_preserved=True,
        equilibrium_solve_performed=False,reference_promoted=False,native_candidate_promoted=False,
        comparative_timing_claim=False,inputs=hashes))
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    print(f'Preflight passes: {len(refined["points"])} nodes, {len(refined["tetrahedra"])} tets; maximum trial gain {gain_eigenvalues.max():.6%}',flush=True)


if __name__=='__main__':main()
