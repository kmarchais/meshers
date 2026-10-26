"""Real zero-contact regression and negative controls for the independent gate."""
import json
from pathlib import Path
import numpy as np
import pytest
import microgen_workflow as w
from periodic_trace import periodic_pairs


def box():
    points=np.array([(x,y,z) for x in (-.5,.5) for y in (-.5,.5) for z in (-.5,.5)])
    lookup={tuple(p):i for i,p in enumerate(points)};faces=[]
    for a in range(3):
        other=[i for i in range(3) if i!=a]
        for side in (-.5,.5):
            nodes=[]
            for u,v in [(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]:
                p=np.zeros(3);p[a]=side;p[other]=[u,v];nodes.append(lookup[tuple(p)])
            faces.extend([[nodes[0],nodes[1],nodes[2]],[nodes[0],nodes[2],nodes[3]]])
    return points,np.array(faces)


@pytest.mark.parametrize('quadratic',[False,True])
def test_real_zero_contact_keeps_complete_caps_and_node_classes(quadratic):
    data=json.loads((Path(__file__).parent/'fixtures/schwarz-d-zero4-periodic-contact.json').read_text())
    points,cells,surface=(np.array(data[k]) for k in ('points','tetrahedra','surface'))
    if quadratic:
        from microgen_tet10_control import straight_elevation
        lifted=straight_elevation(points,cells,surface)
        points,surface=lifted['points'],lifted['surface']
    pairs,classes,check=periodic_pairs(points,surface)
    assert check['complete_existing_node_classes'] and check['periodic_triangles_match']
    assert sum(len(ids) for sides in check['unmatched_noncap_contact_nodes'] for ids in sides)==3
    canonical=points.copy();canonical[canonical==.5]=-.5
    _,groups=np.unique(np.round(canonical,10),axis=0,return_inverse=True)
    assert len(np.unique(groups))==len(np.unique(classes))
    for ids in [np.flatnonzero(groups==g) for g in np.unique(groups)]:
        assert len(np.unique(classes[ids]))==1
    # All six affine macro modes have the exact required jump, also over
    # transitive corner classes without a nonexistent eighth image.
    for i,j in [(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]:
        strain=np.zeros((3,3));strain[i,j]=strain[j,i]=1.
        displacement=points@strain.T
        for axis,pair in enumerate(pairs):
            np.testing.assert_allclose(displacement[pair[:,1]]-displacement[pair[:,0]],
                np.broadcast_to(strain[:,axis],(len(pair),3)),atol=1e-12)


def test_unmatched_cap_node_and_wrong_cap_triangulations_are_rejected():
    points,surface=box();points[4,1]+=.001
    with pytest.raises(ValueError,match='cap node'):periodic_pairs(points,surface)
    points,surface=box();a,b,c=surface[2];_,_,d=surface[3]
    surface[2]=[a,b,d];surface[3]=[b,c,d]
    with pytest.raises(ValueError,match='triangulations'):periodic_pairs(points,surface)
    seed_control=periodic_pairs(points,surface,require_facets=False)[2]
    assert not seed_control['periodic_triangles_match'] and not seed_control['full_facet_check_required']


def test_unmatched_noncap_nodes_cannot_hide_disconnected_classes():
    points,surface=box();points=np.vstack((points,[[-.5,-.5,.1],[.5,.5,.1]]))
    with pytest.raises(ValueError,match='Disconnected periodic'):periodic_pairs(points,surface)
    points=np.vstack((points,[[.5,-.5,.1]]))
    assert periodic_pairs(points,surface)[2]['complete_existing_node_classes']


def test_previous_four_constraint_spaces_are_exactly_unchanged():
    w.imports()
    from periodic import periodic_pairs as historical
    receipt=json.loads((w.ROOT/'handoff/macbook/adaptive-interior-validation-20261025.json').read_text())
    for case in receipt['cases'].values():
        mesh=np.load(w.ROOT/case['run']/'mesh.npz')
        old,old_classes,_=historical(mesh['points'],mesh['surface'])
        new,new_classes,check=periodic_pairs(mesh['points'],mesh['surface'])
        for a,b in zip(old,new,strict=True):np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(old_classes,new_classes)
        assert all(not ids for sides in check['unmatched_noncap_contact_nodes'] for ids in sides)
