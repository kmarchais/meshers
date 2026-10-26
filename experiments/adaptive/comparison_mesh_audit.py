"""Shared comparison admission, including quadratic edge/face conformity."""
import json
import numpy as np
import microgen_workflow as w
from nested_refinement import EDGES
from facet_incidence_audit import audit as facet_audit


def quadratic_conformity(points,cells,surface):
    assert len(np.unique(cells))==len(points),'Unused nodes cannot be counted as useful mesh nodes'
    assert len(np.unique(points,axis=0))==len(points),'Duplicate coordinates'
    if cells.shape[1]==4:return dict(passed=True,order=1)
    pairs=np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges,first,inverse=np.unique(pairs,axis=0,return_index=True,return_inverse=True)
    mids=cells[:,4:].ravel()[first]
    assert np.array_equal(mids[inverse],cells[:,4:].ravel()),'Quadratic crack: edge has different midpoint IDs'
    assert len(np.unique(mids))==len(mids) and not np.isin(mids,np.unique(cells[:,:4])).any(),'Midpoint role collision'
    lookup={tuple(e):int(m) for e,m in zip(edges,mids,strict=True)}
    expected=np.array([[lookup[tuple(sorted((int(f[a]),int(f[b]))))] for a,b in ((0,1),(1,2),(2,0))] for f in surface])
    assert np.array_equal(expected,surface[:,3:]),'Surface Tri6 association differs from volume'
    return dict(passed=True,order=2,unique_edges=len(edges))


def observe(run,name,case,points,cells,surface,periodic_required):
    import scorecard
    from periodic import periodic_pairs
    from quadratic_edge_removal_probe import measure
    from microgen_workflow import tessellate
    try:topology=quadratic_conformity(points,cells,surface)
    except AssertionError as error:topology=dict(passed=False,error=str(error))
    facets=facet_audit(points,cells,surface)
    q,_,volume_corner=scorecard.tet_metrics(points,cells[:,:4])
    straight=cells.shape[1]==10 and np.allclose(points[cells[:,4:]],points[cells[:,EDGES]].mean(axis=2),atol=1e-14,rtol=0)
    if cells.shape[1]==4 or straight:
        volume=float(volume_corner.sum());minimum=float(volume_corner.min()*6)
    else:
        jacobian,volumes=measure(points[cells]);volume=float(volumes.sum());minimum=float(jacobian.min())
    reference=json.loads((run/(name+'-specimen.json')).read_text())
    row=dict(volume=volume,relative_volume_error=volume/reference['volume']-1,q_min_corner=float(q.min()),q_p01_corner=float(np.quantile(q,.01)),
        minimum_sampled_jacobian=minimum,facet_consistency=facets,oriented_facet_consistency_passed=facets['oriented_facet_consistency_passed'],
        quadratic_conformity=topology,exact_straight_reduction_verified=bool(straight),geometry_admitted=False)
    failures=[]
    if not topology['passed']:failures.append('quadratic conformity')
    if q.min()<.1:failures.append('quality')
    if minimum<=0:failures.append('Jacobian')
    if abs(row['relative_volume_error'])>.01:failures.append('volume')
    if not facets['oriented_facet_consistency_passed']:failures.append('facets')
    periodic={}
    if periodic_required:
        try:_,_,periodic=periodic_pairs(points,surface)
        except ValueError as error:failures.append('periodicity');periodic=dict(error=str(error))
    else:periodic=dict(required=False,reason='Linear thickness grading has unequal opposite-z material traces; bounded nonperiodic case for both arms')
    if failures:
        samples=np.load(run/(name+'-samples.npz'));queries=np.concatenate([samples['walls'],samples['caps']])
        if cells.shape[1]==4 or straight:vertices,faces=points,surface[:,:3]
        else:vertices,faces=tessellate(points,surface,8)
        row.update(sampled_reverse_max=float(scorecard.distance_to_mesh(vertices,faces,queries).max()),
            distance_audit_complete=False,projection_failed=None,
            audit_scope='Fail-fast after required quality/volume/topology/periodicity rejection; no forward-distance or successful speed claim')
    else:
        original=w.specimen
        w.specimen=lambda requested:(case,None,None) if requested==name else original(requested)
        try:
            # Straight Tet10 is exactly affine; reducing only its observer geometry
            # is exact. The complete original Tet10 topology and periodicity were checked above.
            full=w.geometry(run,name,points,cells[:,:4] if straight else cells,surface[:,:3] if straight else surface,
                reverse_subdivisions=(8,16),box_projection_fallback=True)
        finally:w.specimen=original
        row.update(full);row.update(quadratic_conformity=topology,exact_straight_reduction_verified=bool(straight),distance_audit_complete=True)
        # Preserve the independent full-order positive125 check.
        row['minimum_sampled_jacobian']=min(row['minimum_sampled_jacobian'],minimum)
        if not row['geometry_admitted']:failures.append('geometry')
    row['failed_gates']=failures
    row['geometry_admitted']=not failures
    return row,periodic,dict(zip(('min','p01','p10','median','p90','p99','max'),map(float,np.quantile(q,[0,.01,.1,.5,.9,.99,1]))))


def controls():
    import scorecard
    from microgen_tet10_control import straight_elevation
    from quadratic_edge_removal_probe import measure
    points=np.array([[0.,0.,0.],[1.,0.,0.],[.5,np.sqrt(3)/2,0.],[.5,np.sqrt(3)/6,np.sqrt(2/3)]])
    q,_,volume=scorecard.tet_metrics(points,np.arange(4)[None])
    np.testing.assert_allclose(q,1,atol=1e-14);np.testing.assert_allclose(volume,np.sqrt(2)/12,atol=1e-14)
    p=np.vstack((points,points[EDGES].mean(axis=1)))
    minimum,v=measure(p[np.arange(10)[None]])
    assert minimum.min()>0;np.testing.assert_allclose(v,volume,atol=1e-14)
    p2=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
    f=np.array([[0,1,2],[0,2,3]]);queries=np.array([[.3,.4,.002],[.3,.4,-.003],[1.1,.3,0.]])
    np.testing.assert_allclose(scorecard.distance_to_mesh(p2,f,queries),[.002,.003,.1],atol=1e-14)
    sphere=scorecard.CASES['sphere'];probes=np.eye(3)*.36
    foot,good=scorecard.project(sphere,probes)
    assert good.all();np.testing.assert_allclose(np.linalg.norm(foot-probes,axis=1),.01,atol=1e-12)
    # Two tets share a face. A second midpoint ID on its edge is a real crack,
    # even though the old corner-only facet audit still passes.
    p3=np.vstack((points,points[3]*np.array([1,1,-1])))
    cells=np.array([[0,1,2,3],[0,2,1,4]])
    from facet_incidence_audit import FACES
    allfaces=cells[:,FACES].reshape(-1,3)
    _,first,counts=np.unique(np.sort(allfaces,axis=1),axis=0,return_index=True,return_counts=True)
    faces=allfaces[first[counts==1]]
    elevated=straight_elevation(p3,cells,faces)
    ep,ec,es=(elevated[k] for k in ('points','tetrahedra','surface'))
    assert quadratic_conformity(ep,ec,es)['passed']
    cracked=ec.copy();old=int(cracked[1,6]);cracked[1,6]=len(ep)
    badpoints=np.vstack((ep,ep[old]+[0,0,.001]))
    assert facet_audit(badpoints,cracked,es)['oriented_facet_consistency_passed']
    try:quadratic_conformity(badpoints,cracked,es)
    except AssertionError as error:assert 'Quadratic crack' in str(error)
    else:raise AssertionError('Deliberate quadratic crack was accepted')
    return dict(regular_quality_and_volume=True,affine_Tet10_measure=True,analytic_plane_distance=True,analytic_sphere_distance=True,
        shared_quadratic_edge=True,deliberate_quadratic_crack_rejected=True)
