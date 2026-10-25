"""Replay fixed-case mesh admission and the actual periodic-normal before/after."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_,scorecard = w.imports()
    from facet_incidence_audit import audit
    from periodic import periodic_pairs
    comparison = w.ROOT/'handoff/macbook/mesh-comparison-20261024.json'
    old = json.loads(comparison.read_text())
    hashes = {str(p.resolve()):w.digest(p) for p in (Path(__file__),comparison,
        w.HERE/'facet_incidence_audit.py',w.HERE/'scorecard.py')}
    shutil.copyfile(__file__,args.run/'analysis-driver.py')
    pairs = []
    for case,tag in (('gyroid','gyroid'),('schwarz-solid','schwarz'),('spinodoid','spinodoid')):
        old_row = next(r for c in old['cases'] if c['case']==case for r in c['rows'] if r['name']=='meshers-curved-known')
        before = Path(old_row['run'])/'repeat-0-meshers-curved-known'
        after = w.HERE/'runs'/f'macbook-periodic-edge-normal-{tag}-curved-v1/repeat-0-meshers-curved-known'
        pairs.append((case,before,after))
    pairs.append(('smooth-lattice',w.HERE/'runs/macbook-lattice26-current-native-amg-v1',
                  w.HERE/'runs/macbook-periodic-edge-normal-lattice26-amg-v1'))
    results = []
    for case,before,after in pairs:
        rows = [json.loads((p/'result.json').read_text()) for p in (before,after)]
        meshes = [np.load(p/'mesh.npz') for p in (before,after)]
        for p,row,mesh in zip((before,after),rows,meshes,strict=True):
            assert row['qualified'] and row['geometry']['geometry_admitted']
            points,cells,surface = (mesh[k] for k in ('points','tetrahedra','surface'))
            assert len(points) == row['nodes'] <= 250000 and len(cells) == row['elements'] <= 250000
            q,_,volume = scorecard.tet_metrics(points,cells)
            assert q.min() >= .1 and volume.min() > 0
            assert audit(points,cells,surface)['oriented_facet_consistency_passed']
            _,_,periodic = periodic_pairs(points,surface)
            assert periodic['periodic_triangles_match']
            for name in ('result.json','mesh.npz'):
                path = p/name;hashes[str(path.resolve())] = w.digest(path)
        identical_connectivity = all(np.array_equal(meshes[0][k],meshes[1][k]) for k in ('tetrahedra','surface'))
        displacement = None
        if meshes[0]['points'].shape == meshes[1]['points'].shape:
            displacement = float(abs(meshes[0]['points']-meshes[1]['points']).max())
        if case != 'smooth-lattice':
            assert identical_connectivity and displacement < 1e-12
        results.append(dict(case=case,before=str(before),after=str(after),
            before_nodes=rows[0]['nodes'],after_nodes=rows[1]['nodes'],
            before_tets=rows[0]['elements'],after_tets=rows[1]['elements'],
            before_geometry=rows[0]['geometry'],after_geometry=rows[1]['geometry'],
            identical_connectivity=identical_connectivity,maximum_same_index_coordinate_change=displacement))
    for before,after in ((p[1],p[2]) for p in pairs):
        for suffix in ('-samples.npz','-specimen.json'):
            # Case directory and frozen specimen names are explicit in their result.
            case=json.loads((after/'result.json').read_text())['case']
            a,b=before/(case+suffix),after/(case+suffix)
            assert w.digest(a)==w.digest(b)
            hashes[str(a.resolve())]=w.digest(a);hashes[str(b.resolve())]=w.digest(b)
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(passed=True,cases=results,unchanged_fixed_specimen_hashes=True,
        all_four_geometry_quality_periodic_gates_pass=True,no_new_reference=True,
        comparative_timing_claim=False,inputs=hashes,
        scope='Actual saved mesh and specimen replay; lattice accuracy is assessed separately against all24, no all-metric victory inferred'))
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    print('Four fixed specimen/mesh gates replay; lattice numerical accuracy remains a separate gate',flush=True)


if __name__ == '__main__':
    main()
