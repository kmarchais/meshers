"""One bounded geometry-only tangent fit on a localized curved cap edge."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.optimize import minimize_scalar

import microgen_workflow as w
import comparison_fixtures as fixtures
from comparison_mesh_audit import observe
import scorecard


def feet(function, points, level, locked):
    original=np.asarray(points,dtype=float)
    foot=original.copy();mask=np.broadcast_to(locked,foot.shape)
    for _ in range(40):
        gradient=scorecard.gradient(function,foot);gradient[mask]=0
        residual=function(*foot.T)-level
        step=residual[:,None]*gradient/np.maximum((gradient*gradient).sum(1),1e-30)[:,None]
        foot-=step
        if np.linalg.norm(step,axis=1).max()<1e-12:break
    residual=abs(function(*foot.T)-level)/np.maximum(np.linalg.norm(scorecard.gradient(function,foot),axis=1),1e-30)
    good=np.isfinite(foot).all(1)&(residual<1e-9)&(abs(foot)<=.5+1e-12).all(1)
    return foot,good


def curve(ends,midpoint,t):
    return ((1-t)*(1-2*t))[:,None]*ends[0]+(t*(2*t-1))[:,None]*ends[1]+(4*t*(1-t))[:,None]*midpoint


def fit(function,ends,midpoint,level,locked,tangent,valid):
    samples=(np.polynomial.legendre.leggauss(16)[0]+1)/2
    dense=np.linspace(0.,1.,65)
    bound=.1*float(np.linalg.norm(ends[1]-ends[0]))
    def error(candidate,positions):
        probes=curve(ends,candidate,positions)
        projected,good=feet(function,probes,level,locked)
        return float(np.linalg.norm(projected-probes,axis=1).max()) if good.all() else 10.
    trials=[]
    def objective(displacement):
        candidate,good=feet(function,(midpoint+displacement*tangent)[None],level,locked)
        candidate=candidate[0];candidate[locked]=midpoint[locked]
        value=error(candidate,samples) if good.all() and valid(candidate) else 10.
        trials.append((value,float(displacement),candidate.copy()))
        return value
    before=error(midpoint,samples);dense_before=error(midpoint,dense)
    optimum=minimize_scalar(objective,bounds=(-bound,bound),method='bounded',options=dict(maxiter=48,xatol=1e-10))
    value,displacement,candidate=min(trials,key=lambda row:row[0])
    dense_after=error(candidate,dense)
    accepted=bool(value<before-1e-8 and dense_after<=dense_before+1e-10 and valid(candidate))
    return dict(accepted=accepted,midpoint=candidate.tolist(),scalar_displacement=displacement,
        gauss16_before=before,gauss16_after=value,dense65_before=dense_before,dense65_after=dense_after,
        optimizer_success=bool(optimum.success),evaluations=len(trials),maximum_scalar_displacement=bound)


def controls():
    radius=.3;theta=.4
    ends=radius*np.array([[np.cos(theta),-np.sin(theta),0],[np.cos(theta),np.sin(theta),0]])
    function=lambda x,y,z:x*x+y*y+z*z-radius*radius
    locked=np.array([False,False,True])
    shifted=radius*np.array([np.cos(.03),np.sin(.03),0])
    tangent=np.array([-np.sin(.03),np.cos(.03),0])
    result=fit(function,ends,shifted,0.,locked,tangent,lambda _:True)
    assert result['accepted'] and result['dense65_after']<result['dense65_before']
    symmetric=fit(function,ends,np.array([radius,0.,0.]),0.,locked,np.array([0.,1.,0.]),lambda _:True)
    assert not symmetric['accepted']
    return dict(shifted_circle_improves=True,symmetric_circle_stays_fixed=True)


def apply_edges(target,points,cells,edges):
    from quadratic_edge_removal_probe import measure
    points=points.copy();operations=[];touched=set()
    for edge in edges:
        mid=edge['midpoint'];incident=np.flatnonzero((cells==mid).any(1))
        if edge['tangent_dimensions']!=1 or any(i in touched for i in incident):continue
        before=points[cells[incident]].copy();mask=cells[incident]==mid
        def valid(candidate):
            coordinates=before.copy();coordinates[mask]=candidate
            determinant,_=measure(coordinates)
            return bool(determinant.min()>0)
        levels=np.array([v for v in target.levels if np.isfinite(v)])
        level=float(levels[np.argmin(abs(levels-target.function(*points[mid:mid+1].T)[0]))])
        result=fit(target.function,points[edge['endpoints']],points[mid],level,np.array(edge['locked_axes']),
                   np.asarray(edge['tangent_basis'][0]),valid)
        result.update(midpoint_node=mid,incident_elements=incident.tolist(),original_midpoint=points[mid].tolist())
        if result['accepted']:
            candidate=np.asarray(result['midpoint']);minimum=np.inf
            for fraction in (.25,.5,.75,1.):
                coordinates=before.copy();coordinates[mask]=points[mid]+fraction*(candidate-points[mid])
                determinant,_=measure(coordinates);minimum=min(minimum,float(determinant.min()))
            assert minimum>0
            result['minimum_sampled_path_jacobian']=minimum
            points[mid]=candidate;touched.update(map(int,incident))
        operations.append(result)
    return points,operations


def mark_cap_edges(target,points,surface):
    corners=points[surface[:,:3]]
    caps=np.any(np.all(w.on_box_plane(corners,-.5),axis=1)|np.all(w.on_box_plane(corners,.5),axis=1),axis=1)
    wall=surface[~caps]
    pairs=np.sort(wall[:,np.array([[0,1],[1,2],[2,0]])].reshape(-1,2),axis=1)
    ends,first=np.unique(pairs,axis=0,return_index=True);mids=wall[:,3:].ravel()[first]
    locked=np.isclose(points[ends[:,0]],points[ends[:,1]],atol=1e-12,rtol=0)&np.isclose(abs(points[ends[:,0]]),.5,atol=1e-12,rtol=0)
    selected=locked.any(1);ends,mids,locked=ends[selected],mids[selected],locked[selected]
    t=(np.polynomial.legendre.leggauss(16)[0]+1)/2
    weights=np.stack(((1-t)*(1-2*t),t*(2*t-1),4*t*(1-t)),axis=1)
    coordinates=np.stack((points[ends[:,0]],points[ends[:,1]],points[mids]),axis=1)
    probes=np.einsum('sn,end->esd',weights,coordinates).reshape(-1,3)
    levels=np.array([v for v in target.levels if np.isfinite(v)])
    nearest=levels[np.argmin(abs(target.function(*points[mids].T)[:,None]-levels),axis=1)]
    projected,good=feet(target.function,probes,np.repeat(nearest,len(t)),np.repeat(locked,len(t),axis=0))
    errors=np.linalg.norm(projected-probes,axis=1).reshape(-1,len(t)).max(1)
    feasible=good.reshape(-1,len(t)).all(1)
    flagged=np.flatnonzero((errors>.002)&feasible)
    assert len(flagged)<=8,'One bounded pass supports at most eight marked edges'
    rows=[]
    for i in sorted(flagged,key=lambda i:(-errors[i],tuple(ends[i]))):
        normal=scorecard.gradient(target.function,points[mids[i]:mids[i]+1])[0]
        constraints=np.vstack((normal/max(np.linalg.norm(normal),1e-30),np.eye(3)[locked[i]]))
        _,singular,vectors=np.linalg.svd(constraints,full_matrices=True)
        rank=int((singular>1e-12).sum())
        rows.append(dict(endpoints=ends[i].tolist(),midpoint=int(mids[i]),locked_axes=locked[i].tolist(),
                         tangent_dimensions=3-rank,tangent_basis=vectors[rank:].tolist(),gauss16_error=float(errors[i])))
    return rows,dict(cap_edges=len(ends),marked_edges=len(rows),maximum_gauss16_error=float(errors.max(initial=0)),
                     unprojectable_edges=int((~feasible).sum()),detection_uses_fixed_benchmark_samples=False)


def fit_cap_edges(target,points,cells,surface):
    marked,diagnostics=mark_cap_edges(target,points,surface)
    moved,operations=apply_edges(target,points,cells,marked)
    _,after=mark_cap_edges(target,moved,surface)
    assert np.array_equal(points[np.unique(cells[:,:4])],moved[np.unique(cells[:,:4])])
    return moved,dict(before=diagnostics,after=after,operations=operations,native_only=False,
                      source_seed_replayed=False,scope='Fresh native seed plus a measured Python cap-edge fitting pass')


def run(directory,source,localization,fixture_root):
    w.imports();control=controls()
    saved=json.loads((source/'result.json').read_text());localized=json.loads(localization.read_text())
    name=saved['case'];assert name=='graded-gyroid-density-linear' and not fixtures.periodic(name)
    assert localized['complete_above_gate_sample_set'] and not saved['qualified']
    assert len(localized['unique_defective_edges'])<=8
    inputs=[Path(__file__),source/'result.json',source/'mesh.npz',localization,
        fixture_root/'linear-density-calibration.json',w.HERE/'comparison_mesh_audit.py',
        w.HERE/'quadratic_edge_removal_probe.py',w.HERE/'scorecard.py',
        w.ROOT/'handoff/macbook/curved-density-fit-selection-20261026.json']
    inputs.extend(source/(name+suffix) for suffix in ('-samples.npz','-specimen.json'))
    hashes={str(p):w.digest(p) for p in inputs}
    assert hashes[str(source/'mesh.npz')]==saved['mesh_sha256']
    with np.load(source/'mesh.npz') as mesh:points,cells,surface=(mesh[k].copy() for k in ('points','tetrahedra','surface'))
    original=points.copy();target=fixtures.case(name,fixture_root)
    points,operations=apply_edges(target,points,cells,localized['unique_defective_edges'])
    assert np.array_equal(points[np.unique(cells[:,:4])],original[np.unique(cells[:,:4])])
    for suffix in ('-samples.npz','-specimen.json'):shutil.copyfile(source/(name+suffix),directory/(name+suffix))
    geometry,periodic,quality=observe(directory,name,target,points,cells,surface,False)
    np.testing.assert_array_equal(list(quality.values()),list(saved['quality_quantiles'].values()))
    np.savez_compressed(directory/'mesh.npz',points=points,tetrahedra=cells,surface=surface)
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    w.save(directory/'result.json',dict(case=name,source=str(source),operations=operations,controls=control,
        admitted=geometry['geometry_admitted'],geometry=geometry,periodic=periodic,quality_quantiles=quality,
        nodes=len(points),elements=len(cells),changed_nodes=np.flatnonzero((points!=original).any(1)).tolist(),
        original_geometry=saved['geometry'],inputs=hashes,mesh_sha256=w.digest(directory/'mesh.npz'),
        source_files_unchanged=True,physical_field_unchanged=True,old_curved_domain_preserved=False,
        old_approximation_space_preserved=False,new_FE_solve=False,native_only=False,
        comparative_timing_claim=False,reference_promoted=False,
        scope='Saved-mesh cap-curve fitting experiment with unchanged counts and corner quality. New curved approximation; no native workflow speed or stiffness claim.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--localization',type=Path,required=True)
    parser.add_argument('--fixtures',type=Path,required=True);args=parser.parse_args()
    assert not (args.run/'result.json').exists()
    run(args.run,args.source,args.localization,args.fixtures)
