"""Compare saved fixed native/actual domains without perturbing either mesh."""

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
    np, *_ = w.imports()
    from hp_ranked_localization import wall_residuals
    from periodic_crease_audit import features
    case,_,_ = w.specimen('smooth-lattice')
    sources = dict(native=w.HERE/'runs/macbook-lattice26-current-native-amg-v1',
                   actual=w.HERE/'runs/macbook-lattice-actual-live-solve-v1')
    paths = [Path(__file__),w.HERE/'hp_ranked_localization.py',w.HERE/'periodic_crease_audit.py',Path(w.__file__)]
    paths.extend(p/name for p in sources.values() for name in
                 ('mesh.npz','result.json','smooth-lattice-specimen.json','smooth-lattice-samples.npz'))
    hashes = {str(p.resolve()):w.digest(p) for p in paths}
    assert len({w.digest(p/'smooth-lattice-specimen.json') for p in sources.values()}) == 1
    assert len({w.digest(p/'smooth-lattice-samples.npz') for p in sources.values()}) == 1
    w.save(args.run/'analysis-inputs.json',hashes)
    for i,path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path,args.run/f'input-{i}-{path.name}')
    results = {}
    for label,source in sources.items():
        print(f'Inspect {label} saved wall geometry without displacement or regeneration',flush=True)
        mesh = np.load(source/'mesh.npz')
        points,cells,surface = (mesh[k] for k in ('points','tetrahedra','surface'))
        assert len(points) <= 250000 and len(cells) <= 250000
        row = json.loads((source/'result.json').read_text())
        assert row['completed_loads'] == 6 and row['geometry']['geometry_admitted']
        corners = points[surface[:,:3]]
        cap = np.any(np.all(w.on_box_plane(corners,-.5),axis=1)|np.all(w.on_box_plane(corners,.5),axis=1),axis=1)
        wall = surface[~cap]
        curve_nodes = features(mesh)['curve_endpoints']
        curve = np.isin(wall[:,:3],curve_nodes).any(axis=1)
        rules = []
        for order in (3,4,6,8):
            arrays = wall_residuals(points,wall,case.function,order)
            np.savez_compressed(args.run/f'{label}-wall-{order**2}.npz',**arrays)
            groups = {}
            for name,mask in dict(all=np.ones(len(wall),dtype=bool),curve_vertex_incident=curve,
                                  without_curve_vertex=~curve).items():
                area = float(arrays['area'][mask].sum())
                groups[name] = dict(faces=int(mask.sum()),area=area,
                    rms_projection=float(np.sqrt(arrays['squared_distance_integral'][mask].sum()/area)),
                    signed_normal_integral=float(arrays['signed_normal_integral'][mask].sum()),
                    maximum_projection=float(arrays['maximum_sampled_projection'][mask].max(initial=0)),
                    failed_projection_points=int(arrays['failed_projections'][mask].sum()),
                    maximum_box_overshoot=float(arrays['maximum_box_overshoot'][mask].max(initial=0)),
                    maximum_remaining_field_distance=float(arrays['maximum_remaining_field_distance'][mask].max(initial=0)))
            rules.append(dict(surface_points_per_face=order**2,groups=groups))
        results[label] = dict(source=str(source),nodes=len(points),tets=len(cells),
                             geometry=row['geometry'],wall_rules=rules)
    for path in args.run.glob('*-wall-*.npz'):
        hashes[str(path.resolve())] = w.digest(path)
    w.save(args.run/'result.json',dict(passed=True,domains=results,identical_fixed_specimen_inputs=True,
        meshes_unchanged=True,no_solve=True,no_geometry_offset=True,no_reference_promoted=True,inputs=hashes,
        scope='Sampled wall projection diagnostics on saved domains, not stiffness shape derivatives, continuum bounds or fresh workflow timing'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(json.dumps({k:v['wall_rules'][-1]['groups'] for k,v in results.items()}),flush=True)


if __name__ == '__main__':
    main()
