"""Extract an exact periodic cap-edge obstruction and its complete local stars."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--audit-run',type=Path,
        default=w.HERE/'runs/macbook-spinodoid-linear-wall-audit-v1')
    args=parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_,scorecard=w.imports()
    from periodic import periodic_pairs
    from nested_refinement import EDGES
    audit=args.audit_run
    measured=json.loads((audit/'result.json').read_text())['cases'][0]
    source=Path(measured['source']);mesh=np.load(source/'mesh.npz')
    points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
    worst=measured['geometry']['worst_boundary_samples'][0]
    weights=scorecard.BARY[worst['probe_index']]
    assert np.count_nonzero(weights)==2 and np.all(weights[weights!=0]==.5)
    edge=np.asarray(worst['node_ids'])[weights!=0]
    midpoint=points[edge].mean(axis=0);foot=np.asarray(worst['foot'])
    np.testing.assert_allclose(midpoint,worst['probe'],atol=1e-15)
    _,classes,_=periodic_pairs(points,surface)
    pairs=np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    physical_edges=np.unique(pairs,axis=0)
    key=np.sort(classes[edge])
    orbit=physical_edges[np.all(np.sort(classes[physical_edges],axis=1)==key,axis=1)]
    assert len(orbit)>=1
    records=[];core=[]
    for endpoints in orbit:
        parents=np.flatnonzero(np.isin(cells,endpoints).sum(axis=1)==2)
        assert len(parents)>0
        core.extend(cells[parents].ravel())
        center=points[endpoints].mean(axis=0)
        target=foot+(center-midpoint)
        before,_,volume=scorecard.tet_metrics(points,cells[parents])
        assert volume.min()>0
        attempts=[]
        for fraction in (0.,.25,.5,.75,1.):
            candidate=center+fraction*(target-center)
            coordinates=points[cells[parents]].copy()
            children=[]
            for endpoint in endpoints:
                child=coordinates.copy();child[cells[parents]==endpoint]=candidate
                children.append(child)
            children=np.concatenate(children)
            flat=children.reshape(-1,3);indices=np.arange(len(flat)).reshape(-1,4)
            quality,_,vol=scorecard.tet_metrics(flat,indices)
            admissible=vol.min()>0 and quality.min()>=.1
            attempts.append(dict(fraction=fraction,minimum_quality=float(quality.min()),
                minimum_determinant=float(6*vol.min()),passes_hard_quality=bool(admissible)))
        records.append(dict(edge=endpoints.tolist(),coordinates=points[endpoints].tolist(),
            parent_count=len(parents),minimum_parent_quality=float(before.min()),
            original_guard_floor=max(.5*min(float(before.min()),.3),.1),
            midpoint=center.tolist(),target=target.tolist(),attempts=attempts))
    core_classes=np.unique(classes[np.asarray(core,dtype=int)])
    pack=np.flatnonzero(np.isin(classes[cells],core_classes).any(axis=1))
    used=np.unique(cells[pack]);mapping=np.full(len(points),-1,dtype=int);mapping[used]=np.arange(len(used))
    boundary=surface[np.isin(surface,used).all(axis=1)]
    case,shape,_=w.specimen('spinodoid');f=shape._frep
    coefficients=[[*map(float,mode/f.cell_size),float(amplitude),float(phase)]
        for mode,amplitude,phase in zip(f.modes,f.amplitudes,f.phases,strict=True)]
    fixture=dict(points=points[used].tolist(),tets=mapping[cells[pack]].tolist(),surface=mapping[boundary].tolist(),
        bad_coordinates=points[edge].tolist(),coefficients=coefficients,threshold=float(f.threshold),
        periodic=[True]*3,minimum_quality=.1,source=str(source),source_mesh_sha256=w.digest(source/'mesh.npz'))
    w.save(args.run/'cap-star.json',fixture)
    paths=[Path(__file__),audit/'result.json',source/'result.json',source/'mesh.npz',args.run/'cap-star.json']
    paths.extend(w.HERE/name for name in ('scorecard.py','microgen_workflow.py','nested_refinement.py'))
    hashes={str(p.resolve()):w.digest(p) for p in paths}
    shutil.copyfile(__file__,args.run/'analysis-driver.py')
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(passed=True,records=records,
        fixture_nodes=len(used),fixture_tets=len(pack),coefficients=len(coefficients),
        inputs=hashes,mesh_changed=False,native_change_promoted=False,
        scope='Local split quality obstruction and exact fixed-field periodic star fixture; no mesh admission or FE solve'))
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    print(f'Extracted complete periodic cap-star neighborhood: {len(used)} nodes, {len(pack)} tets',flush=True)


if __name__=='__main__':main()
