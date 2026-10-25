"""Bounded exact extruded nonlinear feature control, no changed main specimens."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args = parser.parse_args()
    np,_,meshers,*_ = w.imports()
    def field(x,y,z):
        return np.maximum(np.maximum(abs(x)-.35,abs(y)-.35),.12-np.sqrt((x-.03)**2+(y+.02)**2))
    inputs = {str(Path(__file__).resolve()):w.digest(Path(__file__))}
    import shutil
    shutil.copyfile(__file__,args.run/'analysis-driver.py')
    result = dict(specimen='Fixed square prism with offset cylindrical bore, unit z extent',
        exact_volume=.7*.7-np.pi*.12**2,periodic_axes=[False,False,True],inputs=inputs)
    try:
        mesh = meshers.generate(field,bounds=(-.5,.5)*3,cells=24,periodic=(False,False,True),
            adapt=True,min_length=.002,max_length=.1,gradation=.3,adapt_iterations=4,
            geometry_tolerance=.002,minimum_quality=.1,max_tetrahedra=250000,
            optimize_passes=1,threads=1,order=1)
    except meshers.MeshingError as error:
        result.update(status='rejected',error=str(error),admitted=False)
    else:
        low = np.flatnonzero(abs(mesh.points[:,2]+.5) < 1e-10)
        high = np.flatnonzero(abs(mesh.points[:,2]-.5) < 1e-10)
        assert len(low) == len(high) > 0
        keys = lambda nodes:{tuple(np.round(mesh.points[i,:2],11)):int(i) for i in nodes}
        lower,upper = keys(low),keys(high)
        assert len(lower) == len(low) and len(upper) == len(high) and lower.keys() == upper.keys()
        expected = {(lower[key],upper[key]) for key in lower}
        assert set(map(tuple,mesh.periodic_pairs[2])) == expected
        assert not len(mesh.periodic_pairs[0]) and not len(mesh.periodic_pairs[1])
        def caps(z):
            faces = mesh.surface[np.all(abs(mesh.points[mesh.surface,2]-z) < 1e-10,axis=1)]
            return {tuple(sorted(tuple(np.round(mesh.points[i,:2],11)) for i in f)) for f in faces}
        assert caps(-.5) == caps(.5)
        periodic = dict(paired_z_nodes=len(low),matching_cap_triangles=len(caps(-.5)),
                        coordinate_error=max(float(abs(mesh.points[b,:2]-mesh.points[a,:2]).max()) for a,b in expected))
        volume_error = mesh.diagnostics['volume']/result['exact_volume']-1
        assert abs(volume_error) <= .01
        assert len(mesh.points) <= 250000 and len(mesh.tetrahedra) <= 250000
        assert mesh.diagnostics['minimum_mmg_quality'] >= .1
        np.savez_compressed(args.run/'mesh.npz',points=mesh.points,tetrahedra=mesh.tetrahedra,surface=mesh.surface)
        result.update(status='admitted',admitted=True,nodes=len(mesh.points),tets=len(mesh.tetrahedra),
            volume_error=volume_error,diagnostics=mesh.diagnostics,periodic=periodic,
            paired_z_nodes=len(low))
    w.save(args.run/'result.json',result)
    assert all(w.digest(Path(p)) == h for p,h in inputs.items())
    print(result['status'],result.get('error',''),flush=True)


if __name__ == '__main__':
    main()
