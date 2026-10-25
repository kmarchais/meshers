"""Compare patch predictions with a separately reassembled periodic FE operator."""

import numpy as np
import pytest

import microgen_workflow as w

w.imports()
from split_patch_energy import local_blocks,patch_structure,assemble_patches,solve_patches


@pytest.fixture(scope='module')
def patch_data():
    from fedoo.lib_elements.tetrahedron import Tet10
    from homogenization import strain_tensor
    from periodic import periodic_pairs
    np,_,meshers,*_=w.imports()
    mesh=meshers.generate(lambda x,*_:x*0-1,bounds=(-.5,.5)*3,cells=4,order=2,
        geometry_tolerance=1.,minimum_quality=.1,periodic=(True,True,True),threads=1)
    points=mesh.points.copy();cells=mesh.tetrahedra;surface=mesh.surface
    _,classes,_=periodic_pairs(points,surface)
    rng=np.random.default_rng(71)
    displacement=rng.normal(scale=.001,size=points.shape)[classes]
    displacement[np.isclose(abs(points),.5,atol=1e-12)]=0
    displacement[np.unique(cells[:,:4])]=0
    points+=displacement
    _,classes,_=periodic_pairs(points,surface)
    affine=np.array([points@strain_tensor(e).T for e in np.eye(6)])
    fields=affine+.01*rng.normal(size=(6,len(points),3))[:,classes]
    structure=patch_structure(cells,classes)
    element=Tet10(15);rule=element.xi_pg,element.w_pg
    blocks=[]
    for offset in range(6):
        g,r,_=local_blocks(points,cells,fields,offset,rule)
        blocks.append((offset,g,r))
    grams,residuals=assemble_patches(structure,blocks)
    coefficients,gains=solve_patches(grams,residuals)
    return points,cells,surface,classes,affine,fields,structure,rule,grams,residuals,coefficients,gains


def test_affine_periodic_fields_have_zero_new_node_residual_on_curved_maps(patch_data):
    points,cells,_,_,affine,_,structure,rule,*_=patch_data
    blocks=[]
    for offset in range(6):
        g,r,_=local_blocks(points,cells,affine,offset,rule)
        blocks.append((offset,g,r))
    _,residuals=assemble_patches(structure,blocks)
    assert max(abs(r).max() for r in residuals)<1e-11


@pytest.mark.parametrize('multiple_copies',[False,True])
def test_patch_gain_matches_independent_reduced_operator_and_full_energy(patch_data,multiple_copies):
    from capture_reference_system import CapturedSystem
    from nested_refinement import EDGES,refine,interpolation_check
    from stress_recovery import solve
    points,cells,surface,_,_,fields,structure,_,_,_,_,gains=patch_data
    def copies(group):
        return {tuple(sorted(cells[parent,EDGES[offset]])) for parent,offset in structure['incidence'][group]}
    group=next(i for i in range(len(gains)) if (len(copies(i))>1)==multiple_copies and np.trace(gains[i])>1e-10)
    preferred=[next(iter(copies(group)))]
    refined=refine(points,cells,surface,fields,preferred_edges=preferred,preferred_only=True)
    assert refined['diagnostics']['accepted_edge_orbits']==1
    interpolation_check(points,cells,fields,refined)
    checked={}

    class Capture:
        def __init__(self):self.problem,self.records=None,[]
        def __call__(self,matrix,rhs,**kwargs):
            n=len(refined['points']);free=self.problem._dof_free
            new=np.flatnonzero((free<3*n)&(free%n>=len(points)))
            assert len(new)>0
            transform=self.problem._MatCB
            full=np.zeros((transform.shape[0],6))
            full[:3*n]=refined['displacement'].transpose(0,2,1).reshape(6,-1).T
            operator=self.problem.get_A()
            residual=np.asarray(transform.T@(operator@full))[new]
            gram=matrix.tocsr()[new][:,new].toarray()
            update=-np.linalg.solve(gram,residual)
            predicted=-residual.T@update
            np.testing.assert_allclose(predicted,gains[group],atol=1e-11,rtol=1e-8)
            reduced_delta=np.zeros((len(rhs),6));reduced_delta[new]=update
            delta=transform@reduced_delta
            physical=delta[:3*n].T.reshape(6,3,n)
            assert abs(physical[:,:,:len(points)]).max()<1e-12
            trial=full+delta
            actual=full.T@(operator@full)-trial.T@(operator@trial)
            np.testing.assert_allclose(actual,predicted,atol=1e-11,rtol=1e-8)
            checked['passed']=True
            raise CapturedSystem()

    with pytest.raises(CapturedSystem):
        solve(refined['points'],refined['tetrahedra'],n_elm_gp=15,solver=Capture(),release_reduced_cache=True)
    assert checked['passed']
