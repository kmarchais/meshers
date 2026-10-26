"""Locate all above-gate samples in a frozen curved-density mesh, without changing it."""
import argparse
import json
from pathlib import Path

import numpy as np

import microgen_workflow as w
import comparison_fixtures as fixtures


def investigate(source, fixture_root):
    w.imports()
    from box_projection import constrained_foot
    from facet_incidence_audit import FACES
    import scorecard
    saved=json.loads((source/'result.json').read_text())
    name=saved['case'];assert name=='graded-gyroid-density-linear'
    assert not saved['qualified'] and saved['geometry']['failed_gates']==['geometry']
    inputs=[source/'result.json',source/'mesh.npz',source/(name+'-samples.npz'),source/(name+'-specimen.json'),
            fixture_root/'linear-density-calibration.json',Path(__file__),Path(w.__file__),
            w.HERE/'scorecard.py',w.HERE/'box_projection.py',w.HERE/'comparison_fixtures.py']
    hashes={str(p):w.digest(p) for p in inputs}
    assert hashes[str(source/'mesh.npz')]==saved['mesh_sha256']
    mesh=np.load(source/'mesh.npz');points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
    assert cells.shape[1]==10 and surface.shape[1]==6
    target=fixtures.case(name,fixture_root)
    original=w.specimen
    w.specimen=lambda requested:(target,None,None) if requested==name else original(requested)
    try:
        geometry=w.geometry(source,name,points,cells,surface,detail_count=64,
                            reverse_subdivisions=(8,16),box_projection_fallback=True)
    finally:w.specimen=original
    np.testing.assert_allclose(geometry['sampled_two_sided_max'],saved['geometry']['sampled_two_sided_max'],atol=1e-13,rtol=0)
    worst=geometry['worst_boundary_samples']
    assert len(worst)==64 and worst[-1]['distance']<.002,'64 witnesses do not cover all above-gate samples'
    bad=[r for r in worst if r['distance']>.002]
    assert bad,'The frozen rejection must still reproduce'
    bad_faces={r['face'] for r in bad}
    wanted={tuple(sorted(surface[f,:3])):f for f in bad_faces}
    owners={f:[] for f in bad_faces}
    for i,cell in enumerate(cells[:,:4]):
        for face in cell[FACES]:
            key=tuple(sorted(face))
            if key in wanted:owners[wanted[key]].append(i)
    assert all(len(v)==1 for v in owners.values())
    knots=np.asarray(json.loads((fixture_root/'linear-density-calibration.json').read_text())['z'])
    finite=np.array([v for v in target.levels if np.isfinite(v)])
    edges={};rows=[]
    for witness in bad:
        row=dict(witness,owner_elements=owners[witness['face']])
        bary=scorecard.BARY[witness['probe_index']]
        zero=np.flatnonzero(abs(bary)<1e-14)
        row['barycentric']=bary.tolist()
        row['on_quadratic_edge']=len(zero)==1
        if len(zero)==1:
            active=[i for i in range(3) if i!=zero[0]]
            local_edge=next(k for k,pair in enumerate(((0,1),(1,2),(2,0))) if set(pair)==set(active))
            face=surface[witness['face']]
            ends=sorted(int(face[i]) for i in active);mid=int(face[3+local_edge])
            locked=np.isclose(abs(points[ends]),.5,atol=1e-12,rtol=0).all(axis=0)
            # Both endpoints must lie on the same cap, not opposite caps.
            locked&=np.isclose(points[ends[0]],points[ends[1]],atol=1e-12,rtol=0)
            normal=scorecard.gradient(target.function,points[mid:mid+1])[0]
            constraints=np.vstack((normal/max(np.linalg.norm(normal),1e-30),np.eye(3)[locked]))
            _,singular,vectors=np.linalg.svd(constraints,full_matrices=True)
            rank=int((singular>1e-12).sum())
            zmin,zmax=sorted(points[ends,2]);length=float(np.linalg.norm(points[ends[1]]-points[ends[0]]))
            info=dict(endpoints=ends,midpoint=mid,locked_axes=locked.tolist(),length=length,
                      midpoint_displacement=float(np.linalg.norm(points[mid]-points[ends].mean(axis=0))),
                      tangent_dimensions=3-rank,tangent_basis=vectors[rank:].tolist(),
                      lookup_knots_strictly_inside_edge_z_range=int(((knots>zmin)&(knots<zmax)).sum()))
            edges[tuple(ends)]=info;row['edge']=info
        probe=np.array(witness['probe']);initial=np.array(witness['foot'])
        value=target.function(*probe[None].T)[0];level=finite[np.argmin(abs(value-finite))]
        foot,control=constrained_foot(target.function,probe,level,np.array(witness['locked']),initial)
        row['independent_feasible_foot']=dict(foot=foot.tolist(),**control)
        rows.append(row)
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    return dict(case=name,source=str(source),source_nodes=len(points),source_elements=len(cells),
                source_geometry=geometry,above_gate_sample_count=len(bad),above_gate_face_count=len(bad_faces),
                complete_above_gate_sample_set=True,unique_defective_edges=list(edges.values()),witnesses=rows,
                source_arrays_unchanged=True,source_hashes=hashes,mesh_changed=False,
                all_four_original_cases_still_in_scope=True,reference_promoted=False,
                scope='Localization of the unchanged finite sampled rejection, not a repaired mesh, proof of globally closest distances or a new stiffness result')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--fixtures',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    w.save(args.output,investigate(args.source,args.fixtures))
