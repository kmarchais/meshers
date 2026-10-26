"""One reversible periodic leaf coarsening and metric-only refinement exchange."""
import argparse
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

import microgen_workflow as w
from nested_refinement import EDGES,FACES,refine,interpolation_check,positive_quality
from lattice_improvement_audit import cut_edges
from crease_directional_refinement import edge_metric_squared
from directional_nesting_audit import prolong


def edges_of(mesh):
    cells=mesh['tetrahedra']
    pairs=np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges,first=np.unique(pairs,axis=0,return_index=True)
    mids=cells[:,4:].ravel()[first]
    return edges,mids


def compact(points,cells,parents,maps):
    edges=np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    unique,first,inverse=np.unique(edges,axis=0,return_index=True,return_inverse=True)
    mids=cells[:,4:].ravel()[first]
    assert np.array_equal(mids[inverse],cells[:,4:].ravel())
    assert len(np.unique(mids))==len(mids)
    assert not np.isin(mids,np.unique(cells[:,:4])).any()
    mid_by_edge={tuple(e):int(m) for e,m in zip(unique,mids,strict=True)}
    faces=cells[:,FACES].reshape(-1,3)
    _,first,count=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_counts=True)
    assert count.max()<=2
    boundary=faces[first[count==1]]
    surface=np.array([list(f)+[mid_by_edge[tuple(sorted((int(f[a]),int(f[b]))))]
        for a,b in ((0,1),(1,2),(2,0))] for f in boundary])
    used=np.unique(cells);ranks=np.full(len(points),-1,dtype=int);ranks[used]=np.arange(len(used))
    result=dict(points=points[used],tetrahedra=ranks[cells],surface=ranks[surface],
        parents=parents,parent_barycentric_maps=maps,source_node_ids=used)
    from facet_incidence_audit import audit
    assert audit(result['points'],result['tetrahedra'],result['surface'])['oriented_facet_consistency_passed']
    return result


def coarsen(old,fine,parents,maps,tensor):
    from periodic import periodic_pairs
    np.testing.assert_array_equal(old['points'],fine['points'][:len(old['points'])])
    refined,cut,mids=cut_edges(old['tetrahedra'],parents,maps)
    ratio=np.sqrt(edge_metric_squared(old['points'],cut,mids,tensor))
    old_edges,_=edges_of(old)
    _,classes,_=periodic_pairs(old['points'],old['surface'])
    groups={}
    for edge in old_edges:
        groups.setdefault(tuple(sorted(classes[edge])),[]).append(tuple(edge))
    cut_parents={}
    for i,(parent,edge) in enumerate(zip(refined,cut,strict=True)):
        cut_parents.setdefault(tuple(edge),[]).append((int(parent),float(ratio[i])))
    old_q=positive_quality(old['points'][old['tetrahedra'][:,:4]])
    fine_q=positive_quality(fine['points'][fine['tetrahedra'][:,:4]])
    child_min=np.full(len(old['tetrahedra']),np.inf);np.minimum.at(child_min,parents,fine_q)
    restore=[];orbits=0;boundary_orbits=0
    surface_nodes=set(old['surface'].ravel())
    for group in groups.values():
        if not all(e in cut_parents for e in group):continue
        records=[record for e in group for record in cut_parents[e]]
        ids=np.array([r[0] for r in records]);values=np.array([r[1] for r in records])
        if not (values<.8).all() or not (old_q[ids]>=child_min[ids]-1e-12).all():continue
        # Recover the complete star, including all periodic images.
        for edge in group:
            u,v=edge
            actual=np.flatnonzero((old['tetrahedra'][:,:4]==u).any(axis=1)&
                (old['tetrahedra'][:,:4]==v).any(axis=1))
            np.testing.assert_array_equal(np.sort([r[0] for r in cut_parents[edge]]),actual)
        restore.extend(ids);orbits+=1
        boundary_orbits+=int(any(u in surface_nodes and v in surface_nodes for u,v in group))
    restore=np.array(sorted(restore),dtype=int)
    assert len(np.unique(restore))==len(restore)
    keep=~np.isin(parents,restore)
    out=compact(fine['points'],np.vstack((fine['tetrahedra'][keep],old['tetrahedra'][restore])),
        np.r_[parents[keep],restore],np.concatenate((maps[keep],np.repeat(np.eye(4)[None],len(restore),axis=0))))
    new_edges,new_mids=edges_of(out)
    old_fine_edges={tuple(e) for e in edges_of(fine)[0]}
    introduced=np.array([tuple(e) not in old_fine_edges for e in out['source_node_ids'][new_edges]])
    restored_ratio=np.sqrt(edge_metric_squared(out['points'],new_edges[introduced],new_mids[introduced],tensor))
    assert (restored_ratio < .8*(1+1e-12)).all(), 'Coarsening must not introduce overlong edges'
    return out,dict(coarsened_periodic_orbits=orbits,coarsened_boundary_orbits=boundary_orbits,
        restored_parent_tets=len(restore),introduced_coarse_edges=int(introduced.sum()),
        maximum_restored_edge_metric_ratio=float(restored_ratio.max()) if len(restored_ratio) else None,nodes_released=len(fine['points'])-len(out['points']),
        tets_removed=len(fine['tetrahedra'])-len(out['tetrahedra']))


def synthetic_fields(points):
    from homogenization import strain_tensor
    fields=np.array([points@strain_tensor(load).T for load in np.eye(6)])
    for load in range(6):fields[load,:,load%3]+=.01*np.sin(2*np.pi*points[:,0])*np.cos(2*np.pi*points[:,1])
    return fields


def exchange(old,fine,parents,maps,tensor,node_cap,tet_cap):
    from hp_shared_domain import compose,periodic_jumps
    coarsened,stats=coarsen(old,fine,parents,maps,tensor)
    fields=synthetic_fields(old['points'])
    coarsened['displacement']=prolong(coarsened['points'],coarsened['tetrahedra'],old['tetrahedra'],
        fields,coarsened['parents'],coarsened['parent_barycentric_maps'])
    check0=interpolation_check(old['points'],old['tetrahedra'],fields,coarsened)
    periodic_jumps(coarsened['points'],coarsened['surface'],coarsened['displacement'])
    edges,mids=edges_of(coarsened)
    squared=edge_metric_squared(coarsened['points'],edges,mids,tensor)
    scores={tuple(map(int,e)):float(v-(4/3)**2) for e,v in zip(edges,squared,strict=True) if v>(4/3)**2}
    result=refine(coarsened['points'],coarsened['tetrahedra'],coarsened['surface'],coarsened['displacement'],
        maximum_nodes=node_cap,maximum_tets=tet_cap,minimum_quality=.1,edge_scores=scores)
    check1=interpolation_check(coarsened['points'],coarsened['tetrahedra'],coarsened['displacement'],result)
    periodic_jumps(result['points'],result['surface'],result['displacement'])
    parents2,maps2=compose(coarsened['parents'],coarsened['parent_barycentric_maps'],result['parents'],result['parent_barycentric_maps'])
    common=dict(result,parents=parents2,parent_barycentric_maps=maps2)
    check2=interpolation_check(old['points'],old['tetrahedra'],fields,common)
    return coarsened,result,dict(stats,coarsen_mapping=check0,refine_mapping=check1,common_parent_mapping=check2,
        refinement=result['diagnostics'],old_fine_space_preserved=False,common_coarser_space_preserved=True)


def controls():
    _,_,meshers,*_=w.imports()
    mesh=meshers.generate(lambda x,y,z: -1.+0*x,bounds=(-.5,.5)*3,cells=4,
        periodic=(True,)*3,order=2,optimize_passes=0,threads=1)
    old={k:np.array(getattr(mesh,{'tetrahedra':'tetrahedra'}.get(k,k))) for k in ('points','tetrahedra','surface')}
    old['points'][:,2]+=.002*np.prod(np.sin(2*np.pi*old['points']),axis=1)
    edges,_=edges_of(old);delta=old['points'][edges[:,1]]-old['points'][edges[:,0]]
    selected=(abs(delta[:,0])>.2)&(np.linalg.norm(delta[:,1:],axis=1)<1e-10)
    scores={tuple(map(int,e)):1. for e in edges[selected]}
    fields=synthetic_fields(old['points'])
    fine=refine(old['points'],old['tetrahedra'],old['surface'],fields,edge_scores=scores)
    tensor=lambda p:np.broadcast_to(np.diag([1.44,36.,36.]),(len(p),3,3))
    coarse,result,stats=exchange(old,fine,fine['parents'],fine['parent_barycentric_maps'],tensor,
        len(fine['points']),len(fine['tetrahedra']))
    assert stats['coarsened_periodic_orbits']>0 and stats['coarsened_boundary_orbits']>0
    assert stats['nodes_released']>0 and result['diagnostics']['accepted_edge_orbits']>0
    assert len(result['points'])<=len(fine['points']) and len(result['tetrahedra'])<=len(fine['tetrahedra'])
    from quadratic_edge_removal_probe import measure
    volumes=[]
    for m in (old,coarse,result):
        minimum,volume=measure(m['points'][m['tetrahedra']]);assert minimum.min()>0
        volumes.append(float(volume.sum()))
    assert np.max(abs(np.array(volumes)-volumes[0]))<1e-12
    return stats


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path);args=parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    shutil.copyfile(__file__,args.run/'exchange-driver.py')
    checks=controls()
    if args.source is None:
        w.save(args.run/'result.json',dict(controls_passed=True,periodic_exchange=checks));return
    from periodic_crease_grading import make_metric,directional_metric
    from periodic import periodic_pairs
    from quadratic_edge_removal_probe import measure
    with np.load(args.source/'feature-seed/mesh.npz') as mesh:seed=SimpleNamespace(**{k:mesh[k] for k in mesh.files})
    case,_,_=w.specimen('smooth-lattice');scalar,metadata,cloud=make_metric(seed,case.function,.002,.1,26)
    source_row=json.loads((args.source/'result.json').read_text())
    assert metadata['law']==source_row['crease_bisection']['law']
    with np.load(w.HERE/'runs/macbook-lattice-crease-distance-preflight-v1/crease-cloud.npz') as fixed:
        np.testing.assert_array_equal(cloud,fixed['canonical_points'])
    tensor=directional_metric(scalar,cloud,case.function,.1,5.)
    base=args.source/'directional-refinement';old_dir=base/'level-32';fine_dir=base/'level-33'
    def load(directory):
        with np.load(directory/'mesh.npz') as f:return {k:f[k] for k in f.files}
    old,fine=load(old_dir),load(fine_dir)
    with np.load(fine_dir/'parent-map.npz') as f:parents,maps=f['parents'],f['parent_barycentric_maps']
    coarsened,result,stats=exchange(old,fine,parents,maps,tensor,w.CAP,w.CAP)
    rows={}
    for label,mesh in [('coarsened',coarsened),('exchanged',result)]:
        directory=args.run/label;directory.mkdir()
        for suffix in ('-samples.npz','-specimen.json'):
            shutil.copyfile(args.source/('smooth-lattice'+suffix),directory/('smooth-lattice'+suffix))
        geometry=w.geometry(directory,'smooth-lattice',mesh['points'],mesh['tetrahedra'],mesh['surface'])
        assert geometry['geometry_admitted']
        minimum,volume=measure(mesh['points'][mesh['tetrahedra']]);assert minimum.min()>0
        expected=json.loads((fine_dir/'result.json').read_text())['geometry']['volume']
        assert abs(volume.sum()/expected-1)<1e-10
        _,_,periodic=periodic_pairs(mesh['points'],mesh['surface'])
        edges,mids=edges_of(mesh);ratios=np.sqrt(edge_metric_squared(mesh['points'],edges,mids,tensor))
        np.savez_compressed(directory/'mesh.npz',**{k:mesh[k] for k in ('points','tetrahedra','surface')})
        q=positive_quality(mesh['points'][mesh['tetrahedra'][:,:4]])
        rows[label]=dict(qualified=True,nodes=len(mesh['points']),elements=len(mesh['tetrahedra']),geometry=geometry,
            periodic=periodic,remaining_metric_edges=int((ratios>4/3).sum()),maximum_metric_ratio=float(ratios.max()),
            quality_quantiles=np.quantile(q,[0,.01,.1,.5,.9,.99,1]).tolist())
        w.save(directory/'result.json',rows[label])
    paths=[Path(__file__),w.HERE/'nested_refinement.py',w.HERE/'lattice_improvement_audit.py',w.HERE/'directional_nesting_audit.py',
        old_dir/'mesh.npz',fine_dir/'mesh.npz',fine_dir/'parent-map.npz']
    w.save(args.run/'result.json',dict(passed=True,controls=checks,source=str(args.source),stages=rows,exchange=stats,
        policy=dict(coarsen_only_last_level_siblings=True,collapse_ratio=.8,refine_ratio=4/3,exchange_passes=1,
            directional_ratio=5,minimum_quality=.1,field_selection=False),new_equilibrium_fields=False,
        scope='One saved-hierarchy metric exchange at unchanged domain/caps; not live generation, original-fine nested space, accuracy or timing result',
        inputs={str(p):w.digest(p) for p in paths}))


if __name__=='__main__':main()
