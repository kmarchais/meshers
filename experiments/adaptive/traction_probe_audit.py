"""Denser face-rule replay plus an independent finite-difference stress check."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np

import microgen_workflow as w
from traction_probe import integrate, side


def finite_difference_side(points,cells,fields,owner,face,bary):
    from fedoo.lib_elements.tetrahedron import Tet4,Tet10
    local=(cells[owner,:4,None]==face[:,None,:]).argmax(axis=1)
    reference=Tet4().xi_nd
    xi=np.einsum('a,fad->fd',bary,reference[local])
    geometry=points[cells[owner]];values=fields[:,cells[owner]]
    jacobian=np.empty((len(face),3,3));field_derivative=np.empty((len(face),6,3,3))
    for axis in range(3):
        shift=np.eye(3)[axis]*1e-5
        weights=(Tet10().shape_function(xi+shift)-Tet10().shape_function(xi-shift))/(2e-5)
        jacobian[:,axis]=np.einsum('fn,fnd->fd',weights,geometry)
        field_derivative[:,:,:,axis]=np.einsum('fn,lfni->fli',weights,values)
    gradient=np.einsum('flia,fba->flib',field_derivative,np.linalg.inv(jacobian))
    strain=(gradient+gradient.swapaxes(-1,-2))/2
    return 2*strain+np.trace(strain,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)


def curved_control():
    from fedoo.lib_elements.tetrahedron import Tet4,Tet10
    reference=Tet4().xi_nd
    maps=np.array([reference,reference.copy()]);maps[1,3,0]=-1;maps[1,[0,1]]=maps[1,[1,0]]
    points=(Tet4().shape_function(Tet10().xi_nd)[None]@maps).reshape(-1,3)
    points[[4,14],0]+=.04
    cells=np.arange(20).reshape(2,10)
    fields=np.einsum('lij,nj->lni',np.random.default_rng(18).normal(size=(6,3,3)),points)
    matrix,_,error=integrate(points,cells,fields,np.array([0]),np.array([[0,1,2]]),
                             np.array([1]),np.array([[11,10,12]]),order=6)
    assert abs(matrix).max()<1e-23 and error<1e-14
    return dict(affine_physical_field_on_curved_maps_zero_jump=True)


def run(args):
    w.imports();control=curved_control()
    source=json.loads((args.source/'result.json').read_text())
    paths=[Path(__file__),w.HERE/'traction_probe.py',args.source/'result.json',
           w.ROOT/'handoff/macbook/traction-quadrature-selection-20261026.json']
    hashes={str(p):w.digest(p) for p in paths}
    rows={};maximum_fd=0.
    for label,data in source['routes'].items():
        directory=Path(data['source'])
        for name in ('mesh.npz','fields.npz','result.json'):
            path=directory/name;assert w.digest(path)==source['inputs'][str(path)];hashes[str(path)]=w.digest(path)
        with np.load(directory/'mesh.npz') as mesh:points,cells=(mesh[k] for k in ('points','tetrahedra'))
        fields=np.load(directory/'fields.npz')['displacement']
        inverse=np.linalg.inv(json.loads((directory/'result.json').read_text())['stiffness'])
        groups=[]
        for group in data['groups']:
            path=args.source/f'{label}-{group["group"]}.npz'
            assert w.digest(path)==source['inputs'][str(path)]
            hashes[str(path)]=w.digest(path)
            archive=np.load(path)
            arguments=[archive[f'input_{i}'] if f'input_{i}' in archive else None for i in range(4)]
            # Independent derivatives use shape values, not the analytic derivative implementation.
            bary=np.array([.23,.31,.46])
            for owner,face in (arguments[:2],arguments[2:]):
                if owner is None:continue
                analytic=side(points,cells,fields,owner[:2],face[:2],bary)[0]
                numeric=finite_difference_side(points,cells,fields,owner[:2],face[:2],bary)
                error=float(np.linalg.norm(analytic-numeric)/max(np.linalg.norm(analytic),1e-30))
                assert error<1e-7;maximum_fd=max(maximum_fd,error)
            rules=[];previous=archive['gram_4']
            for order in args.orders:
                gram,area,mapping=integrate(points,cells,fields,*arguments,order=order)
                estimate=gram.sum(axis=0)*group['population']/group['sampled']
                change=float(np.linalg.norm(gram-previous)/max(np.linalg.norm(gram),1e-30))
                rules.append(dict(points=order**2,normalized_trace=float(np.trace(inverse@estimate)/6),
                                  gram_relative_change=change,mapping_error=mapping))
                previous=gram
            groups.append(dict(group=group['group'],population=group['population'],sampled=group['sampled'],
                               rules=rules,quadrature_stable=rules[-1]['gram_relative_change']<.01))
        rows[label]=groups
    assert all(w.digest(Path(path))==expected for path,expected in hashes.items())
    for i,path in enumerate((Path(__file__),w.HERE/'traction_probe.py')):
        shutil.copyfile(path,args.output.parent/f'input-{i}-{path.name}')
    w.save(args.output,dict(controls=control,source=str(args.source),routes=rows,
        maximum_independent_stress_relative_error=maximum_fd,
        all_groups_quadrature_stable=all(g['quadrature_stable'] for groups in rows.values() for g in groups),
        unchanged_inputs=True,inputs=hashes,new_solve=False,reference_promoted=False,
        scope='Same selected faces, denser positive quadrature and independent finite-difference stress. Not a stiffness estimate or certified error bound.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--orders',type=int,nargs='+',default=[6,8])
    args=p.parse_args();assert not args.output.exists() and all(5<=n<=12 for n in args.orders);run(args)
