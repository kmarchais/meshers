"""Local curved Tet10 feasibility of existing native valence-3/4 edge removal."""

import argparse
from collections import Counter
import json
from pathlib import Path
import shutil

import numpy as np
import microgen_workflow as w
from nested_refinement import EDGES,FACES,positive_quality


def oriented(points, corners):
    corners = np.array(corners,dtype=np.int64)
    xyz = points[corners]
    if np.linalg.det(xyz[1:]-xyz[:1]) < 0:
        corners[[0,1]] = corners[[1,0]]
    return corners


def boundary(cells):
    faces = Counter()
    for cell in cells:
        for face in FACES:
            corners = sorted(int(cell[i]) for i in face)
            mid = {tuple(sorted((int(cell[a]),int(cell[b])))):int(cell[4+i])
                for i,(a,b) in enumerate(EDGES)}
            key = tuple(corners)+tuple(mid[tuple(sorted((corners[a],corners[b])))]
                for a,b in ((0,1),(1,2),(2,0)))
            faces[key] += 1
    assert max(faces.values()) <= 2
    return {f for f,count in faces.items() if count == 1}


def measure(coordinates, *, positive_rule=True):
    from fedoo.lib_elements.tetrahedron import Tet10
    from quadrature_audit import duffy_rule
    element = Tet10(15)
    volume = np.zeros(len(coordinates))
    minimum = np.full(len(coordinates),np.inf)
    for d,weight in zip(element.shape_function_derivative(element.xi_pg),element.w_pg,strict=True):
        det = np.linalg.det(np.einsum('an,tnb->tab',d,coordinates))
        volume += weight*det
        minimum = np.minimum(minimum,det)
    probes = np.vstack((element.xi_nd,[[.25]*3],duffy_rule(5)[0])) if positive_rule else element.xi_nd
    for d in element.shape_function_derivative(probes):
        minimum = np.minimum(minimum,np.linalg.det(np.einsum('an,tnb->tab',d,coordinates)))
    return minimum,volume


def ring_of(cells,u,v):
    adjacency = {}
    for cell in cells:
        rest = [int(x) for x in cell[:4] if x not in (u,v)]
        if len(rest)!=2:return None
        a,b = rest
        adjacency.setdefault(a,set()).add(b)
        adjacency.setdefault(b,set()).add(a)
    if len(adjacency)!=len(cells) or any(len(x)!=2 for x in adjacency.values()):return None
    start = min(adjacency);ring = [start,min(adjacency[start])]
    while len(ring)<len(cells):
        other = adjacency[ring[-1]]-{ring[-2]}
        if len(other)!=1:return None
        node = next(iter(other))
        if node in ring:return None
        ring.append(node)
    return ring if start in adjacency[ring[-1]] else None


def variants(points,cells,u,v,mids,existing_faces):
    ring = ring_of(cells,u,v)
    if ring is None:return
    triangles = [[ring]] if len(ring)==3 else [
        [[ring[0],ring[1],ring[2]],[ring[0],ring[2],ring[3]]],
        [[ring[1],ring[2],ring[3]],[ring[1],ring[3],ring[0]]]]
    old_boundary = boundary(cells)
    before = float(positive_quality(points[cells[:,:4]]).min())
    old_linear_volume = sum(np.linalg.det(points[c[1:4]]-points[c[:1]]) for c in cells)
    for triangulation in triangles:
        new = []
        for a,b,c in triangulation:
            du = np.linalg.det(points[[a,b,c]]-points[u])
            dv = np.linalg.det(points[[a,b,c]]-points[v])
            if du*dv >= 0:break
            new.extend((oriented(points,[a,b,c,u]),oriented(points,[a,b,c,v])))
        if len(new)!=2*len(triangulation):continue
        new = np.array(new)
        after = float(positive_quality(points[new]).min())
        if after<=before+1e-6:continue
        linear_volume = sum(np.linalg.det(points[c[1:]]-points[c[:1]]) for c in new)
        if abs(linear_volume/old_linear_volume-1)>1e-9:continue
        missing = {tuple(sorted((int(c[a]),int(c[b])))) for c in new for a,b in EDGES
            if tuple(sorted((int(c[a]),int(c[b])))) not in mids}
        if len(missing)!=(len(ring)-3):continue
        diagonal = next(iter(missing)) if missing else None
        new_point = points[list(diagonal)].mean(axis=0) if diagonal else None
        full = []
        for c in new:
            full.append(list(c)+[mids.get(tuple(sorted((int(c[a]),int(c[b])))),len(points)) for a,b in EDGES])
        full = np.array(full,dtype=np.int64)
        if boundary(full)!=old_boundary:continue
        # A new interior face cannot already belong to an external tetrahedron.
        face_counts = Counter(tuple(sorted(int(c[i]) for i in f)) for c in full for f in FACES)
        if any(count==2 and f in existing_faces for f,count in face_counts.items()):continue
        xyz = points[np.minimum(full,len(points)-1)].copy()
        if diagonal:xyz[full==len(points)] = new_point
        yield dict(cells=full,coordinates=xyz,before=before,after=after,diagonal=diagonal,
            new_point=new_point,removed_mid=int(mids[tuple(sorted((u,v)))]))


def controls():
    for n in (3,4):
        angle = 2*np.pi*np.arange(n)/n
        base = list(np.column_stack((np.cos(angle),np.sin(angle),np.zeros(n))))
        base += [np.array([0.,0.,2.]),np.array([0.,0.,-2.])]
        corners = [oriented(np.array(base),[n,n+1,i,(i+1)%n]) for i in range(n)]
        mids = {};cells = []
        for c in corners:
            row = list(c)
            for a,b in EDGES:
                edge = tuple(sorted((int(c[a]),int(c[b]))))
                if edge not in mids:
                    mids[edge] = len(base);base.append((base[edge[0]]+base[edge[1]])/2)
                row.append(mids[edge])
            cells.append(row)
        points,cells = np.array(base),np.array(cells)
        points[:,2] += .01*points[:,0]**2
        old_faces = {tuple(sorted(int(c[i]) for i in f)) for c in cells for f in FACES}
        results = list(variants(points,cells,n,n+1,mids,old_faces))
        assert results
        _,old_volume = measure(points[cells])
        for row in results:
            minimum,volume = measure(row['coordinates'])
            assert minimum.min()>0 and abs(volume.sum()/old_volume.sum()-1)<1e-10
            assert row['removed_mid'] not in row['cells']
            assert row['after']>row['before']+1e-6


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args();assert not args.output.exists()
    shutil.copyfile(__file__,args.output.parent/'analysis-driver.py')
    w.imports();controls()
    from periodic import periodic_pairs
    from scipy.spatial import cKDTree
    mesh = np.load(args.source/'mesh.npz');points,cells,surface = (mesh[k] for k in ('points','tetrahedra','surface'))
    row = json.loads((args.source/'result.json').read_text())
    assert row['qualified'] and row['geometry']['geometry_admitted'] and row['order']==2
    assert len(points)<=w.CAP and len(cells)<=w.CAP
    pairs = np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges,first,inverse,counts = np.unique(pairs,axis=0,return_index=True,return_inverse=True,return_counts=True)
    midpoint = cells[:,4:].ravel()[first]
    assert np.array_equal(midpoint[inverse],cells[:,4:].ravel())
    assert not np.isin(midpoint,np.unique(cells[:,:4])).any()
    mids = {tuple(e):int(m) for e,m in zip(edges,midpoint,strict=True)}
    order = np.argsort(inverse,kind='stable');bounds = np.r_[0,np.cumsum(counts)]
    faces = {tuple(sorted(int(c[i]) for i in f)) for c in cells for f in FACES}
    boundary_nodes = set(surface.ravel());_,classes,_ = periodic_pairs(points,surface)
    class_count = np.bincount(classes,minlength=len(points));tree = cKDTree(points)
    old_minimum,volumes = measure(points[cells],positive_rule=False)
    assert old_minimum.min()>0 and abs(volumes.sum()/row['geometry']['volume']-1)<1e-10
    candidates = [];examined = dict(valence3=0,valence4=0)
    for edge_index in np.flatnonzero((counts==3)|(counts==4)):
        n = int(counts[edge_index]);examined['valence'+str(n)] += 1
        mid = int(midpoint[edge_index])
        if mid in boundary_nodes or class_count[classes[mid]]!=1:continue
        parents = order[bounds[edge_index]:bounds[edge_index+1]]//6
        u,v = map(int,edges[edge_index])
        for proposal in variants(points,cells[parents],u,v,mids,faces):
            if proposal['new_point'] is not None:
                distances,ids = tree.query(proposal['new_point'],k=2)
                if any(distance<1e-12 and node!=mid for distance,node in zip(distances,ids,strict=True)):continue
            proposal.update(edge=[u,v],parents=parents,old_volume=float(volumes[parents].sum()))
            candidates.append(proposal)
    admissible = []
    if candidates:
        sizes = np.array([len(c['cells']) for c in candidates]);starts = np.r_[0,np.cumsum(sizes)]
        coordinates = np.concatenate([c['coordinates'] for c in candidates])
        minimum,volume = measure(coordinates)
        for i,c in enumerate(candidates):
            sl = slice(starts[i],starts[i+1]);difference = abs(float(volume[sl].sum())/c['old_volume']-1)
            if minimum[sl].min()>0 and difference<1e-10:
                admissible.append(dict(edge=c['edge'],old_tets=len(c['parents']),new_tets=len(c['cells']),
                    before=c['before'],after=c['after'],minimum_sampled_jacobian=float(minimum[sl].min()),
                    relative_volume_difference=difference))
    inputs = [Path(__file__),w.HERE/'nested_refinement.py',w.HERE/'quadrature_audit.py',
        w.ROOT/'crates/meshers-core/src/tet_adapt.rs',args.source/'mesh.npz',args.source/'result.json']
    w.save(args.output,dict(source=str(args.source),controls_passed=True,examined_edges=examined,
        corner_improving_variants=len(candidates),locally_curved_admitted_variants=len(admissible),
        locally_curved_admitted_edges=len({tuple(c['edge']) for c in admissible}),
        valence3_admitted=sum(c['old_tets']==3 for c in admissible),
        valence4_admitted=sum(c['old_tets']==4 for c in admissible),
        examples=sorted(admissible,key=lambda c:c['after']-c['before'],reverse=True)[:10],
        physical_mesh_changed=False,global_mesh_admitted=False,approximation_space_preserved=False,
        fe_fields_used=False,new_equilibrium_solve=False,reference_promoted=False,
        scope='Independent local edge-removal feasibility only. Cavity boundary Tet10 traces retained; no disjoint global mesh, field transfer, nested-space or performance claim.',
        inputs={str(p):w.digest(p) for p in inputs}))
    print(f'{len(admissible)} curved local variants on {len({tuple(c["edge"]) for c in admissible})} edges; '
        'no global mesh changed or admitted',flush=True)


if __name__=='__main__':main()
