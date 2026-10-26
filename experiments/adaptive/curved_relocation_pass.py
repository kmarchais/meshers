"""One disjoint curved interior quality pass on all four saved specimens."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np

import microgen_workflow as w
from nested_refinement import positive_quality
from curved_relocation_probe import propose,controls,POLICY
from quadratic_edge_removal_probe import measure


def quantiles(q):
    return dict(zip(('min','p01','p10','median','p90','p99','max'),
        map(float,np.quantile(q,[0,.01,.1,.5,.9,.99,1]))))


def quality_dominates(before,after):
    assert len(before)==len(after)
    return bool(np.all(np.sort(after)>=np.sort(before)-1e-12))


def distribution_controls():
    # Better local minimum and median can still lose the upper tail.
    assert not quality_dominates(np.array([.2,.6,.8]),np.array([.3,.61,.79]))
    assert quality_dominates(np.array([.2,.6,.8]),np.array([.21,.61,.81]))
    assert quality_dominates(np.array([.2,.6,.8]),np.array([.2,.6,.8]))


def apply_pass(points,cells,surface,dominance=False):
    q=positive_quality(points[cells[:,:4]])
    flat=cells[:,:4].ravel();order=np.argsort(flat,kind='stable')
    bounds=np.r_[0,np.cumsum(np.bincount(flat,minlength=len(points)))]
    rank=np.full(len(points),np.inf);np.minimum.at(rank,flat,np.repeat(q,4))
    boundary_nodes=set(surface.ravel())
    eligible=[int(u) for u in np.argsort(rank,kind='stable') if np.isfinite(rank[u]) and u not in boundary_nodes]
    touched=np.zeros(len(cells),dtype=bool);moved=points.copy();operations=[]
    counts=dict(eligible_vertices=len(eligible),overlap_skipped=0,corner_rejected=0,distribution_rejected=0,curved_path_rejected=0)
    for node in eligible:
        ids=order[bounds[node]:bounds[node+1]]//4
        if touched[ids].any():counts['overlap_skipped']+=1;continue
        row=propose(points,cells[ids],node,boundary_nodes)
        if row is None:counts['corner_rejected']+=1;continue
        if dominance and not quality_dominates(positive_quality(row['before'][:,:4]),positive_quality(row['after'][:,:4])):
            counts['distribution_rejected']+=1;continue
        _,old_volume=measure(row['before'],positive_rule=False)
        minimum,error=np.inf,0.
        for fraction in POLICY['sampled_motion_fractions']:
            determinant,volume=measure(row['before']+fraction*(row['after']-row['before']))
            minimum=min(minimum,float(determinant.min()))
            error=max(error,abs(float(volume.sum()/old_volume.sum()-1)))
        if minimum<=0 or error>=1e-10:counts['curved_path_rejected']+=1;continue
        moved[cells[ids]]=row['after']
        assert np.max(abs(moved[cells[ids]]-row['after']))<1e-12
        touched[ids]=True
        operations.append(dict(node=node,parents=ids.tolist(),delta=row['delta'].tolist(),
            before_min=row['initial_minimum_quality'],after_min=row['final_minimum_quality'],
            before_median=row['initial_median_quality'],after_median=row['final_median_quality'],
            minimum_path_jacobian=minimum,relative_volume_difference=error))
    np.testing.assert_array_equal(points[surface],moved[surface])
    return moved,operations,counts


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--dominance',action='store_true');args=parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    w.imports();controls();distribution_controls()
    from periodic import periodic_pairs
    from scipy.spatial import cKDTree
    receipt_path=w.ROOT/'handoff/macbook/adaptive-interior-validation-20261025.json'
    receipt=json.loads(receipt_path.read_text())
    sources={name:w.ROOT/receipt['cases'][name]['run'] for name in ('gyroid','schwarz','spinodoid')}
    sources['lattice']=args.source
    paths=[Path(__file__),w.HERE/'curved_relocation_probe.py',w.HERE/'quadratic_edge_removal_probe.py',
        w.HERE/'nested_refinement.py',w.HERE/'microgen_workflow.py',receipt_path]
    for i,p in enumerate(paths):
        if p.suffix=='.py':shutil.copyfile(p,args.run/f'input-{i}-{p.name}')
    results={}
    for name,source in sources.items():
        row=json.loads((source/'result.json').read_text());assert row['qualified']
        case=row['case'];directory=args.run/name;directory.mkdir()
        for suffix in ('-samples.npz','-specimen.json'):
            p=source/(case+suffix);shutil.copyfile(p,directory/p.name)
            assert w.digest(p)==w.digest(directory/p.name);paths.append(p)
        with np.load(source/'mesh.npz') as mesh:
            points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
        moved,operations,counts=apply_pass(points,cells,surface,args.dominance)
        assert not cKDTree(moved).query_pairs(1e-12),'Relocation created coincident nodes'
        geometry=w.geometry(directory,case,moved,cells,surface)
        assert geometry['geometry_admitted']
        minimum,volume=measure(moved[cells]);assert minimum.min()>0
        assert abs(volume.sum()/row['geometry']['volume']-1)<1e-10
        _,_,periodic=periodic_pairs(moved,surface)
        _,_,before_periodic=periodic_pairs(points,surface)
        assert periodic['periodic_pair_counts']==before_periodic['periodic_pair_counts']
        old_q,new_q=(positive_quality(p[cells[:,:4]]) for p in (points,moved))
        before,after=quantiles(old_q),quantiles(new_q)
        assert after['min']>=before['min']-1e-12
        global_dominance=quality_dominates(old_q,new_q)
        if args.dominance:assert global_dominance
        np.savez_compressed(directory/'mesh.npz',points=moved,tetrahedra=cells,surface=surface)
        w.save(directory/'operations.json',operations)
        result=dict(case=case,order=2,source=str(source),qualified=True,status='admitted-needs-reference',
            nodes=len(moved),elements=len(cells),geometry=geometry,periodic=periodic,
            accepted_operations=len(operations),counts=counts,full_sorted_quality_non_decreasing=global_dominance,strict_local_distribution_guard=args.dominance,old_quality_quantiles=before,quality_quantiles=after,
            all_recorded_quality_quantiles_no_worse=all(after[k]>=v-1e-12 for k,v in before.items()),
            curved_boundary_coordinates_identical=True,connectivity_identical=True,old_approximation_space_preserved=False,
            saved_seed_replayed=True,native_only=False,fe_fields_used=False,new_equilibrium_fields=False,
            timing_comparable=False,scope='One saved-mesh disjoint relocation pass; full geometry/quality/PBC admission only, no FE transfer or accuracy/performance claim')
        w.save(directory/'result.json',result);results[name]=result
        w.save(args.run/'partial-results.json',results)
        paths.extend(source/k for k in ('mesh.npz','result.json'))
    w.save(args.run/'result.json',dict(passed=True,cases=results,passes=1,disjoint_old_cells=True,
        policy={**POLICY,'strict_local_distribution_guard':args.dominance,'selection':'All interior corners, initial incident minimum quality then node ID; skip overlapping stars'},
        new_fields=False,reference_promoted=False,objective_complete=False,inputs={str(p):w.digest(p) for p in paths}))
    print('Accepted global moves:',{k:v['accepted_operations'] for k,v in results.items()},flush=True)
    print('Median quality:',{k:[v['old_quality_quantiles']['median'],v['quality_quantiles']['median']] for k,v in results.items()},flush=True)


if __name__=='__main__':main()
