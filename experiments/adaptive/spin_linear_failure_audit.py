"""Inspect frozen linear spinodoid wall failures without remeshing or changing gates."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--source',action='append',type=Path,
        help='Frozen mesh result directory; repeat to compare sources without remeshing')
    args=parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_,scorecard=w.imports()
    case,_,_=w.specimen('spinodoid')
    sources=args.source or [w.HERE/'runs'/name/leaf for name,leaf in (
        ('macbook-spinodoid-linear-repair-red-v1','repeat-0-meshers-tet4-adaptive'),
        ('macbook-spinodoid-linear-guard01-ablation-v1','repeat-0-meshers-tet4-adaptive-guard0.1'))]
    paths=[Path(__file__),w.HERE/'microgen_workflow.py',w.HERE/'scorecard.py']
    paths.extend(p/name for p in sources for name in ('result.json','mesh.npz','spinodoid-samples.npz','spinodoid-specimen.json'))
    hashes={str(p.resolve()):w.digest(p) for p in paths}
    for p in sources:
        for suffix in ('-samples.npz','-specimen.json'):
            origin=p/('spinodoid'+suffix);target=args.run/origin.name
            if target.exists():assert w.digest(target)==w.digest(origin)
            else:shutil.copyfile(origin,target)
    shutil.copyfile(__file__,args.run/'analysis-driver.py')
    manifest=json.loads((args.run/'manifest.json').read_text());manifest['inputs'].update(hashes)
    w.save(args.run/'manifest.json',manifest)
    results=[]
    for index,source in enumerate(sources):
        mesh=np.load(source/'mesh.npz');points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
        print(f'Inspect immutable mesh {index+1}/{len(sources)}',flush=True)
        geometry=w.geometry(args.run,'spinodoid',points,cells,surface,detail_count=12,box_projection_fallback=True)
        corners=points[surface]
        caps=np.any(np.all(w.on_box_plane(corners,-.5),axis=1)|np.all(w.on_box_plane(corners,.5),axis=1),axis=1)
        wall=surface[~caps];vertices=np.unique(wall)
        values=case.function(*points[vertices].T)
        gradient=scorecard.gradient(case.function,points[vertices])
        distance=abs(values)/np.maximum(np.linalg.norm(gradient,axis=1),1e-30)
        q,_,_=scorecard.tet_metrics(points,cells)
        details=[]
        for offset in np.argsort(distance)[-12:][::-1]:
            node=int(vertices[offset]);incident=np.flatnonzero(np.any(cells==node,axis=1))
            faces=wall[np.any(wall==node,axis=1)]
            used=np.unique(cells[incident]);mapping=np.full(len(points),-1,dtype=int);mapping[used]=np.arange(len(used))
            witness=dict(global_node=node,node=int(mapping[node]),points=points[used].tolist(),
                tetrahedra=mapping[cells[incident]].tolist(),wall_faces=mapping[faces].tolist(),
                coordinates=points[node].tolist(),field_value=float(values[offset]),
                first_order_wall_distance=float(distance[offset]),minimum_incident_quality=float(q[incident].min()))
            details.append(witness)
        result=dict(source=str(source),geometry=geometry,wall_vertices=len(vertices),
            first_order_wall_distance_quantiles=np.quantile(distance,[0,.5,.9,.99,1]).tolist(),
            wall_vertices_over_1e6=int((distance>1e-6).sum()),wall_vertices_over_gate=int((distance>.002).sum()),
            worst_wall_vertex_stars=details)
        w.save(args.run/f'case-{index}.json',result);results.append(result)
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(passed=True,cases=results,inputs=hashes,
        mesh_changed=False,gate_changed=False,new_fe_solve=False,native_change_promoted=False,
        scope='Frozen wall and incident-star failure evidence, not mesh admission or a parameter selection'))
    assert all(w.digest(Path(p))==h for p,h in hashes.items())


if __name__=='__main__':main()
