"""One disjoint quality pass with unchanged curved boundary, never field transfer."""

import argparse
import json
from pathlib import Path
import shutil
import numpy as np

import microgen_workflow as w
from nested_refinement import EDGES,FACES,positive_quality
from quadratic_edge_removal_probe import controls,feasible_patches,measure,summary

QUALITY_PASS_POLICY = dict(passes=1,valences=[3,4],minimum_corner_quality_gain=1e-6,
    selection='Descending local minimum-quality gain, then edge and diagonal IDs',
    disjoint_old_cells=True,node_and_tet_increase_allowed=False,
    boundary_Tri6_coordinates_unchanged=True,old_approximation_space_preserved=False)


def apply_patches(points,cells,surface,patches):
    assert len(np.unique(cells))==len(points)
    changed = points.copy();removed = np.zeros(len(cells),dtype=bool)
    additions,operations = [],[]
    diagonals,new_coordinates = set(),set()
    ordered = sorted(patches,key=lambda c:(-(c['after']-c['before']),c['before'],
        tuple(c['edge']),tuple(c['diagonal'] or (-1,-1))))
    for patch in ordered:
        parents = patch['parents']
        if removed[parents].any():continue
        diagonal = patch['diagonal']
        if diagonal and (diagonal in diagonals or tuple(patch['new_point']) in new_coordinates):continue
        new = patch['cells'].copy();new[new==len(points)] = patch['removed_mid']
        if diagonal:
            changed[patch['removed_mid']] = patch['new_point']
            diagonals.add(diagonal);new_coordinates.add(tuple(patch['new_point']))
        removed[parents] = True;additions.append(new)
        operations.append(dict(summary(patch),parents=parents.tolist(),
            removed_mid=patch['removed_mid'],diagonal=diagonal,
            new_point=None if patch['new_point'] is None else patch['new_point'].tolist()))
    if not operations:
        return dict(points=points.copy(),tetrahedra=cells.copy(),surface=surface.copy()),[]
    out = np.concatenate([cells[~removed],*additions])
    assert len(out)<=len(cells)
    used = np.unique(out.ravel());ranks = np.full(len(points),-1,dtype=np.int64)
    ranks[used] = np.arange(len(used))
    assert (ranks[surface]>=0).all()
    final_points,final_cells,final_surface = changed[used],ranks[out],ranks[surface]
    assert len(final_points)<=len(points)
    assert np.array_equal(points[surface],final_points[final_surface])
    assert len(np.unique(final_points,axis=0))==len(final_points)
    pairs = np.sort(final_cells[:,EDGES].reshape(-1,2),axis=1)
    _,first,inverse = np.unique(pairs,axis=0,return_index=True,return_inverse=True)
    mids = final_cells[:,4:].ravel()[first]
    assert np.array_equal(mids[inverse],final_cells[:,4:].ravel())
    assert len(np.unique(mids))==len(mids)
    assert not np.isin(mids,np.unique(final_cells[:,:4])).any()
    assert positive_quality(final_points[final_cells[:,:4]]).min() >= positive_quality(points[cells[:,:4]]).min()-1e-12
    return dict(points=final_points,tetrahedra=final_cells,surface=final_surface),operations


def check_application():
    from facet_incidence_audit import audit
    for points,cells,mids,patches in controls():
        faces = cells[:,FACES].reshape(-1,3)
        _,first,count = np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_counts=True)
        faces = faces[first[count==1]]
        surface = np.array([list(f)+[mids[tuple(sorted((int(f[a]),int(f[b]))))]
            for a,b in ((0,1),(1,2),(2,0))] for f in faces])
        out,operations = apply_patches(points,cells,surface,patches+patches)
        assert len(operations)==1
        expected = 1 if len(cells)==3 else 0
        assert len(out['points'])==len(points)-expected
        assert len(out['tetrahedra'])==len(cells)-expected
        assert audit(out['points'],out['tetrahedra'],out['surface'])['oriented_facet_consistency_passed']
        minimum,volume = measure(out['points'][out['tetrahedra']])
        _,old = measure(points[cells])
        assert minimum.min()>0 and abs(volume.sum()/old.sum()-1)<1e-10
        same,empty = apply_patches(points,cells,surface,[])
        assert not empty and all(np.array_equal(same[k],v) for k,v in
            zip(('points','tetrahedra','surface'),(points,cells,surface),strict=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--run',required=True,type=Path)
    args = parser.parse_args();assert args.run.is_dir() and not (args.run/'result.json').exists()
    w.imports();check_application()
    from periodic import periodic_pairs
    from scipy.spatial import cKDTree
    source = json.loads((args.source/'result.json').read_text())
    assert source['qualified'] and source['geometry']['geometry_admitted'] and source['order']==2
    case = source['case']
    paths = [Path(__file__),w.HERE/'quadratic_edge_removal_probe.py',w.HERE/'nested_refinement.py',
        w.HERE/'facet_incidence_audit.py',w.HERE/'microgen_workflow.py',w.HERE/'scorecard.py',
        w.HERE/'quadrature_audit.py',w.ROOT/'crates/meshers-core/src/tet_adapt.rs']
    paths.extend(args.source/n for n in ('mesh.npz','result.json',case+'-samples.npz',case+'-specimen.json'))
    hashes = {str(p):w.digest(p) for p in paths};w.save(args.run/'analysis-inputs.json',hashes)
    for i,p in enumerate(paths):
        if p.suffix=='.py':shutil.copyfile(p,args.run/f'input-{i}-{p.name}')
    for suffix in ('-samples.npz','-specimen.json'):
        shutil.copyfile(args.source/(case+suffix),args.run/(case+suffix))
    mesh = np.load(args.source/'mesh.npz');points,cells,surface = (mesh[k] for k in ('points','tetrahedra','surface'))
    patches,statistics = feasible_patches(points,cells,surface,source['geometry']['volume'])
    out,operations = apply_patches(points,cells,surface,patches)
    assert max(len(out['points']),len(out['tetrahedra']))<=w.CAP
    assert not cKDTree(out['points']).query_pairs(1e-12),'New near-coincident coordinates need investigation'
    geometry = w.geometry(args.run,case,out['points'],out['tetrahedra'],out['surface'])
    assert geometry['geometry_admitted']
    assert abs(geometry['volume']/source['geometry']['volume']-1)<1e-10
    _,_,before = periodic_pairs(points,surface);_,_,periodic = periodic_pairs(out['points'],out['surface'])
    assert before['periodic_pair_counts']==periodic['periodic_pair_counts']
    old_q,new_q = (positive_quality(p[t[:,:4]]) for p,t in ((points,cells),(out['points'],out['tetrahedra'])))
    quantiles = lambda q:dict(zip(('min','p01','p10','median','p90','p99','max'),map(float,np.quantile(q,[0,.01,.1,.5,.9,.99,1]))))
    np.savez_compressed(args.run/'mesh.npz',**out)
    w.save(args.run/'operations.json',operations)
    w.save(args.run/'result.json',dict(status='admitted-needs-reference',qualified=True,case=case,order=2,
        nodes=len(out['points']),elements=len(out['tetrahedra']),geometry=geometry,periodic=periodic,
        source=str(args.source),source_nodes=len(points),source_elements=len(cells),local_statistics=statistics,
        accepted_disjoint_operations=len(operations),valence3_operations=sum(o['old_tets']==3 for o in operations),
        valence4_operations=sum(o['old_tets']==4 for o in operations),passes=1,
        old_quality_quantiles=quantiles(old_q),quality_quantiles=quantiles(new_q),
        curved_boundary_coordinates_identical=True,all_boundary_Tri6_associations_preserved=True,
        old_approximation_space_preserved=False,fe_fields_used=False,saved_seed_replayed=True,native_only=False,
        comparative_timing_claim=False,reference_promoted=False,stiffness_assessed=False,
        mesh_sha256=w.digest(args.run/'mesh.npz'),operations_sha256=w.digest(args.run/'operations.json'),
        quality_pass_policy=QUALITY_PASS_POLICY,
        inputs=hashes,scope='Saved Meshers prototype with geometry-only curved quality retopology; not live generation, a nested-space correction or an accuracy/performance result'))
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    print(f'{len(operations)} disjoint curved quality operations; median {np.median(old_q):.6f} to {np.median(new_q):.6f}; '
        'global geometry/periodicity admitted, fresh fields required',flush=True)


if __name__=='__main__':main()
