"""Complete curved edge-star replacements with joint metric/quality acceptance."""
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

import microgen_workflow as w
from nested_refinement import EDGES,FACES,positive_quality
from quadratic_edge_removal_probe import boundary,measure,oriented,ring_of,controls as old_controls

THRESHOLD=(4/3)**2
POLICY=dict(valences=[3,4,5,6,7],maximum_examined_orbits=256,
    selection='Descending original metric excess, then edge IDs; no FE fields',
    local_quality='Empirical quality CDF must not worsen, with strict minimum gain1e-6',
    metric='Strict decrease of complete-cavity squared-length excess; no new overlong edge; maximum non-increasing',
    cavity_Tri6_boundary_fixed=True,minimum_quality=.1,relative_volume_tolerance=1e-10)


@lru_cache(None)
def triangulations(vertices):
    if len(vertices)<3:return ((),)
    result=[]
    for split in range(1,len(vertices)-1):
        for left in triangulations(vertices[:split+1]):
            for right in triangulations(vertices[split:]):
                result.append(left+right+((vertices[0],vertices[split],vertices[-1]),))
    return tuple(result)


def dominates(before,after):
    a,b=np.sort(before),np.sort(after)
    shifted=b+1e-12
    thresholds=np.unique(np.r_[a,shifted])
    old=np.searchsorted(a,thresholds,side='right')
    new=np.searchsorted(shifted,thresholds,side='right')
    return bool(np.all(new*len(a)<=old*len(b)))


def variants(points,cells,u,v,mids,existing_faces):
    ring=ring_of(cells,u,v)
    if ring is None:return
    old_boundary=boundary(cells)
    volume0=np.linalg.det(points[cells[:,1:4]]-points[cells[:,:1]]).sum()
    for triangles in triangulations(tuple(ring)):
        corners=[]
        for a,b,c in triangles:
            if np.linalg.det(points[[a,b,c]]-points[u])*np.linalg.det(points[[a,b,c]]-points[v])>=0:break
            corners.extend((oriented(points,[a,b,c,u]),oriented(points,[a,b,c,v])))
        if len(corners)!=2*(len(ring)-2):continue
        corners=np.array(corners)
        if abs(np.linalg.det(points[corners[:,1:]]-points[corners[:,:1]]).sum()/volume0-1)>1e-9:continue
        missing=sorted({tuple(sorted((int(c[a]),int(c[b])))) for c in corners for a,b in EDGES
            if tuple(sorted((int(c[a]),int(c[b])))) not in mids})
        if len(missing)!=len(ring)-3:continue
        additional=np.array([points[list(e)].mean(axis=0) for e in missing]).reshape(-1,3)
        lookup={e:len(points)+i for i,e in enumerate(missing)}
        full=np.array([list(c)+[mids.get(e,lookup.get(e)) for e in
            (tuple(sorted((int(c[a]),int(c[b])))) for a,b in EDGES)] for c in corners],dtype=int)
        if boundary(full)!=old_boundary:continue
        counts=Counter(tuple(sorted(int(c[i]) for i in f)) for c in full for f in FACES)
        if any(n==2 and f in existing_faces for f,n in counts.items()):continue
        coordinates=points[np.minimum(full,len(points)-1)].copy()
        for i,p in enumerate(additional):coordinates[full==len(points)+i]=p
        yield dict(cells=full,coordinates=coordinates,additional_points=additional,new_edges=missing,
            removed_mid=int(mids[tuple(sorted((u,v)))]),old_tets=len(cells),new_tets=len(full))


def metric(cells,coordinates,tensor):
    pairs=np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges,first=np.unique(pairs,axis=0,return_index=True)
    corner_xyz=coordinates[:,EDGES].reshape(-1,2,3)[first]
    mid_xyz=coordinates[:,4:].reshape(-1,3)[first]
    delta=corner_xyz[:,1]-corner_xyz[:,0]
    matrices=tensor(np.vstack((corner_xyz[:,0],mid_xyz,corner_xyz[:,1])))
    eig=np.linalg.eigvalsh(matrices)
    assert eig.min()>0 and (eig[:,-1]/eig[:,0]).max()<=25*(1+1e-12)
    squared=np.maximum.reduce([np.einsum('ni,nij,nj->n',delta,m,delta)
        for m in np.split(matrices,3)])
    return dict(edges=edges,squared=squared,excess=float(np.maximum(squared-THRESHOLD,0).sum()),
        maximum=float(squared.max()))


def joint_guard(points,old_cells,proposal,tensor):
    q0=positive_quality(points[old_cells[:,:4]])
    q1=positive_quality(proposal['coordinates'][:,:4])
    if q1.min()<.1 or q1.min()<=q0.min()+1e-6 or not dominates(q0,q1):return None,'quality'
    before=metric(old_cells,points[old_cells],tensor)
    after=metric(proposal['cells'],proposal['coordinates'],tensor)
    old_edges=set(map(tuple,before['edges']))
    new=np.array([tuple(e) not in old_edges for e in after['edges']])
    if (after['squared'][new]>THRESHOLD*(1+1e-12)).any():return None,'new-overlong-edge'
    if after['maximum']>before['maximum']*(1+1e-12) or after['excess']>=before['excess']-1e-10:
        return None,'metric'
    return dict(before_excess=before['excess'],after_excess=after['excess'],
        before_maximum_ratio=float(np.sqrt(before['maximum'])),after_maximum_ratio=float(np.sqrt(after['maximum'])),
        before_minimum_quality=float(q0.min()),after_minimum_quality=float(q1.min()),
        before_quality=q0.tolist(),after_quality=q1.tolist()),None


def controls():
    assert [len(triangulations(tuple(range(n)))) for n in range(3,8)]==[1,2,5,14,42]
    assert dominates(np.array([.2,.2,.2]),np.array([.3,.3]))
    assert not dominates(np.array([.2,.5]),np.array([.3,.4,.6]))
    tensor=lambda p:np.broadcast_to(.2*np.eye(3),(len(p),3,3))
    accepted=[]
    for points,cells,mids,_ in old_controls():
        n=len(cells);faces={tuple(sorted(int(c[i]) for i in f)) for c in cells for f in FACES}
        proposals=list(variants(points,cells,n,n+1,mids,faces));assert proposals
        valid=[]
        for proposal in proposals:
            row,reason=joint_guard(points,cells,proposal,tensor)
            if row is None:continue
            minimum,volume=measure(proposal['coordinates']);_,old=measure(points[cells])
            assert minimum.min()>0 and abs(volume.sum()/old.sum()-1)<1e-10
            valid.append(row)
        assert valid
        accepted.append(dict(valence=n,accepted=len(valid)))
        if n==4:
            bad=lambda p:np.broadcast_to(np.diag([1.,1.,.2]),(len(p),3,3))
            assert all(joint_guard(points,cells,proposal,bad)[0] is None for proposal in proposals)
    larger=[]
    for n in (5,6,7):
        angle=2*np.pi*np.arange(n)/n
        base=list(np.column_stack((np.cos(angle),np.sin(angle),np.zeros(n))))
        base.extend((np.array([0.,0.,2.]),np.array([0.,0.,-2.])))
        corners=[oriented(np.array(base),[n,n+1,i,(i+1)%n]) for i in range(n)]
        mids={};cells=[]
        for c in corners:
            row=list(c)
            for a,b in EDGES:
                edge=tuple(sorted((int(c[a]),int(c[b]))))
                if edge not in mids:
                    mids[edge]=len(base);base.append((base[edge[0]]+base[edge[1]])/2)
                row.append(mids[edge])
            cells.append(row)
        points,cells=np.array(base),np.array(cells)
        points[:,2]+=.01*points[:,0]**2
        faces={tuple(sorted(int(c[i]) for i in f)) for c in cells for f in FACES}
        proposals=list(variants(points,cells,n,n+1,mids,faces))
        assert len(proposals)==len(triangulations(tuple(range(n))))
        _,old_volume=measure(points[cells])
        for proposal in proposals:
            minimum,volume=measure(proposal['coordinates'])
            assert minimum.min()>0 and abs(volume.sum()/old_volume.sum()-1)<1e-10
        larger.append(dict(valence=n,curved_partitions_checked=len(proposals)))
    return dict(catalan_counts_passed=True,variable_count_distribution_controls_passed=True,
        metric_and_curved_cavity_controls=accepted,larger_curved_partition_controls=larger,new_overlong_diagonal_rejected=True)


def search(source,metric_source,maximum):
    from periodic import periodic_pairs
    from periodic_crease_grading import make_metric,directional_metric
    from crease_directional_refinement import edge_metric_squared
    from scipy.spatial import cKDTree
    source_row=json.loads((source/'result.json').read_text());assert source_row['qualified']
    assert source_row['case']=='smooth-lattice' and source_row['order']==2
    with np.load(source/'mesh.npz') as mesh:points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
    with np.load(metric_source/'feature-seed/mesh.npz') as mesh:seed=SimpleNamespace(**{k:mesh[k] for k in mesh.files})
    case,_,_=w.specimen('smooth-lattice');scalar,metadata,cloud=make_metric(seed,case.function,.002,.1,26)
    assert metadata['law']==json.loads((metric_source/'result.json').read_text())['crease_bisection']['law']
    with np.load(w.HERE/'runs/macbook-lattice-crease-distance-preflight-v1/crease-cloud.npz') as fixed:
        np.testing.assert_array_equal(cloud,fixed['canonical_points'])
    tensor=directional_metric(scalar,cloud,case.function,.1,5.)
    pairs=np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges,first,inverse,counts=np.unique(pairs,axis=0,return_index=True,return_inverse=True,return_counts=True)
    mid=cells[:,4:].ravel()[first]
    assert np.array_equal(mid[inverse],cells[:,4:].ravel()) and len(np.unique(mid))==len(mid)
    assert not np.isin(mid,np.unique(cells[:,:4])).any()
    mids={tuple(e):int(m) for e,m in zip(edges,mid,strict=True)}
    order=np.argsort(inverse,kind='stable');bounds=np.r_[0,np.cumsum(counts)]
    faces={tuple(sorted(int(c[i]) for i in f)) for c in cells for f in FACES}
    _,classes,_=periodic_pairs(points,surface);sizes=np.bincount(classes,minlength=len(points))
    boundary_nodes=set(surface.ravel());tree=cKDTree(points)
    ratio2=edge_metric_squared(points,edges,mid,tensor)
    supported=(counts>=3)&(counts<=7)
    eligible=[int(i) for i in np.flatnonzero((ratio2>THRESHOLD)&supported)
        if int(mid[i]) not in boundary_nodes and sizes[classes[mid[i]]]==1]
    eligible.sort(key=lambda i:(-float(ratio2[i]-THRESHOLD),tuple(edges[i])))
    take=eligible[:maximum]
    stats=dict(requested_edges=int((ratio2>THRESHOLD).sum()),eligible_complete_interior_orbits=len(eligible),
        examined_orbits=len(take),examined_valences={str(n):int(sum(counts[i]==n for i in take)) for n in range(3,8)},
        linear_trace_candidates=0,quality_rejected=0,new_overlong_edge_rejected=0,metric_rejected=0,
        collision_rejected=0,resource_rejected=0,jacobian_or_volume_rejected=0)
    _,old_volumes=measure(points[cells],positive_rule=False)
    assert abs(old_volumes.sum()/source_row['geometry']['volume']-1)<1e-10
    accepted=[]
    for i in take:
        parents=order[bounds[i]:bounds[i+1]]//6;u,v=map(int,edges[i])
        for proposal in variants(points,cells[parents],u,v,mids,faces):
            stats['linear_trace_candidates']+=1
            row,reason=joint_guard(points,cells[parents],proposal,tensor)
            if row is None:
                stats[{'quality':'quality_rejected','new-overlong-edge':'new_overlong_edge_rejected','metric':'metric_rejected'}[reason]]+=1
                continue
            additional=proposal['additional_points']
            if any(np.isclose(abs(p),.5,atol=1e-10,rtol=0).any() or
                any(d<1e-12 and node!=proposal['removed_mid'] for d,node in zip(*tree.query(p,k=2),strict=True))
                for p in additional) or (len(additional)>1 and cKDTree(additional).query_pairs(1e-12)):
                stats['collision_rejected']+=1;continue
            net=proposal['new_tets']-proposal['old_tets']
            if len(points)+net>w.CAP or len(cells)+net>w.CAP:stats['resource_rejected']+=1;continue
            minimum,volume=measure(proposal['coordinates'])
            error=abs(float(volume.sum()/old_volumes[parents].sum()-1))
            if minimum.min()<=0 or error>=1e-10:stats['jacobian_or_volume_rejected']+=1;continue
            accepted.append(dict(proposal,**row,edge=[u,v],parents=parents,
                minimum_sampled_jacobian=float(minimum.min()),relative_volume_difference=error))
    return points,cells,surface,tensor,accepted,stats


def summary(patch):
    return {k:patch[k] for k in ('edge','old_tets','new_tets','before_excess','after_excess',
        'before_maximum_ratio','after_maximum_ratio','before_minimum_quality','after_minimum_quality',
        'minimum_sampled_jacobian','relative_volume_difference')}


def apply_pass(points,cells,surface,patches):
    q0=positive_quality(points[cells[:,:4]])
    probabilities=[0,.01,.1,.5,.9,.99,1]
    initial=np.quantile(q0,probabilities);current=initial.copy()
    removed=np.zeros(len(cells),dtype=bool);additions=[];qualities=[]
    moved=list(points.copy());diagonals=set();new_positions=set();operations=[]
    projected_nodes,projected_tets=len(points),len(cells)
    counts=dict(overlap_rejected=0,budget_rejected=0,collision_rejected=0,global_quality_rejected=0,
        rejected_quality_quantiles={k:0 for k in ('min','p01','p10','median','p90','p99','max')})
    ordered=sorted(patches,key=lambda p:(-(p['before_excess']-p['after_excess'])/
        max(1,p['new_tets']-p['old_tets']),p['new_tets']-p['old_tets'],tuple(p['edge']),tuple(p['new_edges'])))
    for patch in ordered:
        ids=patch['parents'];net=patch['new_tets']-patch['old_tets']
        if removed[ids].any():counts['overlap_rejected']+=1;continue
        if projected_nodes+net>w.CAP or projected_tets+net>w.CAP:
            counts['budget_rejected']+=1;continue
        if any(tuple(e) in diagonals for e in patch['new_edges']) or any(
            tuple(p) in new_positions for p in patch['additional_points']):
            counts['collision_rejected']+=1;continue
        mask=removed.copy();mask[ids]=True
        q1=positive_quality(patch['coordinates'][:,:4])
        trial=np.r_[q0[~mask],*(qualities+[q1])]
        after=np.quantile(trial,probabilities)
        bad=after<current-1e-12
        if bad.any():
            counts['global_quality_rejected']+=1
            for key,rejected in zip(counts['rejected_quality_quantiles'],bad,strict=True):
                counts['rejected_quality_quantiles'][key]+=int(rejected)
            continue
        template=patch['cells'];new=template.copy()
        for i,position in enumerate(patch['additional_points']):
            if i==0:
                index=patch['removed_mid'];moved[index]=position
            else:
                index=len(moved);moved.append(position)
            new[template==len(points)+i]=index
        removed[ids]=True;additions.append(new);qualities.append(q1);current=after
        projected_nodes+=net;projected_tets+=net
        diagonals.update(map(tuple,patch['new_edges']));new_positions.update(map(tuple,patch['additional_points']))
        operations.append(dict(summary(patch),parents=ids.tolist(),net_nodes=net,
            additional_points=patch['additional_points'].tolist()))
    out=np.concatenate([cells[~removed],*additions]) if additions else cells.copy()
    used=np.unique(out);ranks=np.full(len(moved),-1,dtype=int);ranks[used]=np.arange(len(used))
    assert (ranks[surface]>=0).all()
    result=dict(points=np.asarray(moved)[used],tetrahedra=ranks[out],surface=ranks[surface])
    assert len(result['points'])==projected_nodes and len(out)==projected_tets
    np.testing.assert_array_equal(points[surface],result['points'][result['surface']])
    q1=positive_quality(result['points'][result['tetrahedra'][:,:4]])
    assert (np.quantile(q1,probabilities)>=initial-1e-12).all()
    return result,operations,counts


def admit_pass(run,source,points,cells,surface,tensor,patches):
    from periodic import periodic_pairs
    from scipy.spatial import cKDTree
    from crease_directional_refinement import edge_metric_squared
    from hierarchical_metric_exchange import edges_of
    out,operations,counts=apply_pass(points,cells,surface,patches)
    row=json.loads((source/'result.json').read_text());case=row['case']
    for suffix in ('-samples.npz','-specimen.json'):
        path=source/(case+suffix);shutil.copyfile(path,run/path.name)
        assert w.digest(path)==w.digest(run/path.name)
    assert len(np.unique(out['points'],axis=0))==len(out['points'])
    assert not cKDTree(out['points']).query_pairs(1e-12)
    geometry=w.geometry(run,case,out['points'],out['tetrahedra'],out['surface'])
    assert geometry['geometry_admitted']
    minimum,volume=measure(out['points'][out['tetrahedra']])
    assert minimum.min()>0 and abs(volume.sum()/row['geometry']['volume']-1)<1e-10
    _,_,periodic=periodic_pairs(out['points'],out['surface'])
    _,_,old_periodic=periodic_pairs(points,surface)
    assert periodic['periodic_pair_counts']==old_periodic['periodic_pair_counts']
    metrics=[]
    for mesh in (dict(points=points,tetrahedra=cells),out):
        edges,mids=edges_of(mesh)
        squared=edge_metric_squared(mesh['points'],edges,mids,tensor)
        metrics.append(dict(violating_edges=int((squared>THRESHOLD).sum()),
            excess=float(np.maximum(squared-THRESHOLD,0).sum()),maximum_ratio=float(np.sqrt(squared.max()))))
    assert metrics[1]['violating_edges']==metrics[0]['violating_edges']-len(operations)
    assert metrics[1]['maximum_ratio']<=metrics[0]['maximum_ratio']*(1+1e-12)
    delta=sum(p['before_excess']-p['after_excess'] for p in operations)
    assert abs(metrics[0]['excess']-metrics[1]['excess']-delta)<1e-8*max(1,delta)
    np.savez_compressed(run/'mesh.npz',**out);w.save(run/'operations.json',operations)
    q0=positive_quality(points[cells[:,:4]]);q1=positive_quality(out['points'][out['tetrahedra'][:,:4]])
    keys=('min','p01','p10','median','p90','p99','max')
    return dict(qualified=True,status='admitted-needs-reference' if operations else 'no-joint-global-operation',
        case=case,order=2,nodes=len(out['points']),elements=len(out['tetrahedra']),geometry=geometry,periodic=periodic,
        accepted_operations=len(operations),global_mesh_changed=bool(operations),application=counts,metric_before=metrics[0],metric_after=metrics[1],
        global_quality_guard='All seven reported quantiles non-decreasing per accepted operation',
        old_quality_quantiles=dict(zip(keys,map(float,np.quantile(q0,[0,.01,.1,.5,.9,.99,1])))),
        quality_quantiles=dict(zip(keys,map(float,np.quantile(q1,[0,.01,.1,.5,.9,.99,1])))),
        all_reported_quality_quantiles_non_decreasing=True,full_quality_cdf_non_decreasing=dominates(q0,q1),
        physical_boundary_Tri6_coordinates_identical=True,old_approximation_space_preserved=False,
        saved_seed_replayed=True,native_only=False,new_equilibrium_fields=False,
        comparative_timing_claim=False,reference_promoted=False,
        scope='One disjoint saved-mesh joint cavity pass; full admission and metric/quality improvement, not accuracy or live timing')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path);parser.add_argument('--metric-source',type=Path)
    parser.add_argument('--whole-pass',action='store_true')
    args=parser.parse_args();assert args.run.is_dir() and not (args.run/'result.json').exists()
    w.imports();shutil.copyfile(__file__,args.run/'joint-driver.py');checks=controls()
    if args.source is None:
        w.save(args.run/'result.json',dict(passed=True,controls=checks));return
    assert args.metric_source is not None
    limit=10000 if args.whole_pass else POLICY['maximum_examined_orbits']
    points,cells,surface,tensor,patches,statistics=search(args.source,args.metric_source,limit)
    paths=[Path(__file__),w.HERE/'quadratic_edge_removal_probe.py',w.HERE/'nested_refinement.py',
        w.HERE/'periodic_crease_grading.py',args.source/'mesh.npz',args.source/'result.json',
        args.metric_source/'feature-seed/mesh.npz',args.metric_source/'result.json',
        w.HERE/'crease_directional_refinement.py',w.HERE/'hierarchical_metric_exchange.py',
        w.HERE/'quadrature_audit.py',w.HERE/'microgen_workflow.py']
    for suffix in ('-samples.npz','-specimen.json'):
        paths.append(args.source/('smooth-lattice'+suffix))
    result=dict(passed=True,source=str(args.source),policy={**POLICY,'maximum_examined_orbits':limit},controls=checks,
        statistics=statistics,admitted_local_variants=len(patches),
        admitted_local_edges=len({tuple(p['edge']) for p in patches}),
        examples=[summary(p) for p in sorted(patches,key=lambda p:-(p['before_excess']-p['after_excess']))[:10]],
        global_mesh_changed=False,new_fields=False,reference_promoted=False,
        scope='Bounded complete interior cavities with unchanged outer and periodic traces. Local variants can overlap; no global mesh, stiffness or timing claim.',
        inputs={str(p):w.digest(p) for p in paths})
    if args.whole_pass:
        result.update(admit_pass(args.run,args.source,points,cells,surface,tensor,patches))
    w.save(args.run/'result.json',result)
    print('Joint cavity results:',statistics,'admitted',len(patches),flush=True)


if __name__=='__main__':main()
