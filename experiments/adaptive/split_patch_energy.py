"""Exact fixed-old-node variational gains of conforming quadratic edge-star splits."""

from functools import lru_cache

import numpy as np

from nested_refinement import EDGES, child_maps


@lru_cache(maxsize=6)
def template(offset):
    from fedoo.lib_elements.tetrahedron import Tet4, Tet10
    linear,quadratic=Tet4(),Tet10()
    bary=linear.shape_function(quadratic.xi_nd)
    old={tuple(row) for row in np.rint(4*bary).astype(int)}
    maps=child_maps(offset)
    children=[bary@mapping for mapping in maps]
    keys=sorted({tuple(row) for child in children for row in np.rint(4*child).astype(int)}-old)
    assert len(keys)==4
    selectors=[]
    for child in children:
        selected=[(node,keys.index(tuple(key))) for node,key in enumerate(np.rint(4*child).astype(int)) if tuple(key) in keys]
        assert len(selected)==3
        selectors.append(tuple(zip(*selected,strict=True)))
    return np.asarray(keys),maps,selectors


def local_blocks(points,cells,fields,offset,quadrature,batch=1024):
    from fedoo.lib_elements.tetrahedron import Tet4, Tet10
    from hierarchical_residual_trial import stress
    linear,quadratic=Tet4(),Tet10()
    locations,weights=quadrature
    keys,maps,selectors=template(offset)
    child_derivatives=np.asarray(quadratic.shape_function_derivative(locations))
    root_bary=linear.shape_function(locations)
    dl=np.asarray(linear.shape_function_derivative([[.25]*3]))[0]
    gram=np.zeros((len(cells),4,3,4,3));residual=np.zeros((len(cells),4,3,6))
    minimum,volume=np.inf,0.
    for start in range(0,len(cells),batch):
        subset=slice(start,min(start+batch,len(cells)))
        geometry=points[cells[subset]];values=fields[:,cells[subset]]
        for mapping,(nodes,slots) in zip(maps,selectors,strict=True):
            root_locations=root_bary@mapping@linear.xi_nd
            transform=dl@mapping@linear.xi_nd
            mapped=np.einsum('ab,qbn->qan',transform,quadratic.shape_function_derivative(root_locations))
            for dnew,dold,weight in zip(child_derivatives,mapped,weights,strict=True):
                jacobian=np.einsum('an,tnb->tab',dold,geometry)
                determinant=np.linalg.det(jacobian)
                assert determinant.min()>0 and weight>0
                minimum=min(minimum,float(determinant.min()))
                inverse=np.linalg.inv(jacobian)
                old_gradient=np.einsum('tba,an->tbn',inverse,dold)
                gradient=np.zeros((len(geometry),4,3))
                gradient[:,slots]=np.einsum('tba,an->tnb',inverse,dnew[:,nodes])
                weighted=weight*determinant
                dot=np.einsum('tai,tbi->tab',gradient,gradient)
                blocks=(dot[:,:,None,:,None]*np.eye(3)[None,None,:,None,:]
                    +np.einsum('taj,tbi->taibj',gradient,gradient)
                    +np.einsum('tai,tbj->taibj',gradient,gradient))
                gram[subset]+=weighted[:,None,None,None,None]*blocks
                _,sigma=stress(np.einsum('ltni,tjn->ltij',values,old_gradient))
                residual[subset]+=np.einsum('t,taj,ltij->tail',weighted,gradient,sigma,optimize=True)
                volume+=float(weighted.sum())
    return gram.reshape(len(cells),12,12),residual.reshape(len(cells),12,6),dict(minimum_jacobian=minimum,volume=volume)


def patch_structure(cells,classes):
    corners=classes[cells[:,:4]]
    assert np.all(np.diff(np.sort(corners,axis=1),axis=1)>0)
    edges,ids=np.unique(np.sort(corners[:,EDGES],axis=-1).reshape(-1,2),axis=0,return_inverse=True)
    ids=ids.reshape(-1,6)
    incidence=[[] for _ in edges]
    for parent,row in enumerate(ids):
        for offset,group in enumerate(row):incidence[group].append((parent,offset))
    rings=[];slots=[]
    for edge,items in zip(edges,incidence,strict=True):
        ring=sorted(set(corners[[p for p,_ in items]].ravel())-set(edge))
        rings.append(np.asarray(ring,dtype=np.int64));slots.append({node:2+i for i,node in enumerate(ring)})
    return dict(corners=corners,edges=edges,ids=ids,incidence=incidence,rings=rings,slots=slots)


def assemble_patches(structure,blocks):
    grams=[np.zeros((3*(2+len(r)),)*2) for r in structure['rings']]
    residuals=[np.zeros((len(g),6)) for g in grams]
    for offset,local_gram,local_residual in blocks:
        keys,_,_=template(offset)
        for parent,group in enumerate(structure['ids'][:,offset]):
            corner_classes=structure['corners'][parent]
            local_slots=[]
            for key in keys:
                corner=int(key.argmax());node=corner_classes[corner]
                local_slots.append((0 if node==structure['edges'][group,0] else 1) if key[corner]==3
                    else structure['slots'][group][node])
            dofs=(3*np.asarray(local_slots)[:,None]+np.arange(3)).ravel()
            grams[group][np.ix_(dofs,dofs)]+=local_gram[parent]
            residuals[group][dofs]+=local_residual[parent]
    return grams,residuals


def solve_patches(grams,residuals):
    from scipy.linalg import cho_factor,cho_solve
    coefficients=[];gains=[]
    for gram,residual in zip(grams,residuals,strict=True):
        np.testing.assert_allclose(gram,gram.T,atol=1e-12,rtol=1e-10)
        solution=cho_solve(cho_factor((gram+gram.T)/2,lower=True),residual)
        coefficients.append(-solution)
        gain=residual.T@solution
        gains.append((gain+gain.T)/2)
    return coefficients,np.asarray(gains)
