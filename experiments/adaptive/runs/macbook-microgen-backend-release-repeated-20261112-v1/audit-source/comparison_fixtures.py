"""Independent fixed fields and actual Microgen constructors for the comparison."""
from pathlib import Path
import json
import numpy as np
import microgen_workflow as w
from scorecard import Case

CORE=('gyroid','schwarz-solid','spinodoid','smooth-lattice')
EXTRA=('schwarz-d-solid','graded-gyroid-periodic','graded-gyroid-linear','spinodoid-anisotropic')
CASES=CORE+EXTRA+('graded-gyroid-density-linear',)
ORIGINAL_SPECIMEN=w.specimen


def periodic(name):
    return name not in ('graded-gyroid-linear','graded-gyroid-density-linear')


def raw_gyroid(x,y,z):
    x,y,z=(2*np.pi*a for a in (x,y,z))
    return np.sin(x)*np.cos(y)+np.sin(y)*np.cos(z)+np.sin(z)*np.cos(x)


def raw_diamond(x,y,z):
    x,y,z=(2*np.pi*a for a in (x,y,z))
    return np.sin(x)*np.sin(y)*np.sin(z)+np.sin(x)*np.cos(y)*np.cos(z)+np.cos(x)*np.sin(y)*np.cos(z)+np.cos(x)*np.cos(y)*np.sin(z)


def half_width(name,x,y,z):
    return .55*(1+.25*np.cos(2*np.pi*z)) if name.endswith('periodic') else .55*(1+.5*z)


def coefficients_path(name,fixtures):
    return (w.HERE/'runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json' if name=='spinodoid'
        else fixtures/'spinodoid-anisotropic-coefficients.json')


def spectral_field(path):
    data=json.loads(path.read_text());coefficients=[tuple(map(float,row)) for row in data['coefficients']];threshold=data['threshold']
    def field(x,y,z):
        value=threshold
        for a,b,c,amplitude,phase in coefficients:
            value=value-amplitude*np.cos(2*np.pi*(a*(x+.5)+b*(y+.5)+c*(z+.5))+phase)
        return value
    return field


def case(name,fixtures):
    if name=='graded-gyroid-density-linear':
        from linear_density_fixture import field
        return Case(name,field(fixtures),(-1.,1.))
    if name in ('spinodoid','spinodoid-anisotropic'):
        return Case(name,spectral_field(coefficients_path(name,fixtures)),(-np.inf,0.))
    if name in CORE:return ORIGINAL_SPECIMEN(name)[0]
    if name=='schwarz-d-solid':return Case(name,raw_diamond,(-np.inf,0.),.5)
    if name.startswith('graded-gyroid'):
        def field(x,y,z):return raw_gyroid(x,y,z)/half_width(name,x,y,z)
        return Case(name,field,(-1.,1.))
    raise ValueError(name)


def actual_shape(name,fixtures,resolution=None):
    import microgen
    from microgen.shape import surface_functions
    if name in CORE:
        _,shape,part=ORIGINAL_SPECIMEN(name)
        if resolution is not None:
            if name in ('gyroid','schwarz-solid'):
                shape=microgen.Tpms(surface_functions.gyroid if name=='gyroid' else surface_functions.schwarz_p,
                    offset=1.1 if name=='gyroid' else 0.,resolution=resolution)
            elif name=='spinodoid':shape.resolution=resolution
        return shape,part
    if name=='schwarz-d-solid':return microgen.Tpms(surface_functions.schwarz_d,offset=0.,**({} if resolution is None else dict(resolution=resolution))),'lower skeletal'
    if name=='graded-gyroid-density-linear':
        from linear_density_fixture import width
        thickness=width(fixtures)
        def offset(x,y,z):return 2*thickness(x,y,z)
        return microgen.Tpms(surface_functions.gyroid,offset=offset,**({} if resolution is None else dict(resolution=resolution))),'sheet'
    if name.startswith('graded-gyroid'):
        def offset(x,y,z):return 2*half_width(name,x,y,z)
        return microgen.Tpms(surface_functions.gyroid,offset=offset,**({} if resolution is None else dict(resolution=resolution))),'sheet'
    if name=='spinodoid-anisotropic':
        shape=microgen.Spinodoid(k0=4*np.pi,bandwidth=3,offset=0.,anisotropy=(1.,1.,2.),
            cell_size=1.,resolution=16,seed=20261026,mode_amplitude_threshold=.05)
        frozen=json.loads(coefficients_path(name,fixtures).read_text());f=shape._frep
        expected=np.array(frozen['coefficients'])
        np.testing.assert_array_equal(f.modes/f.cell_size,expected[:,:3])
        raw=np.column_stack((f.modes/f.cell_size,f.amplitudes,f.phases))
        assert abs(raw-expected).max()<1e-12 and f.threshold==frozen['threshold']
        f.amplitudes,f.phases=expected[:,3].copy(),expected[:,4].copy()
        if resolution is not None:shape.resolution=resolution
        return shape,'solid'
    raise ValueError(name)


def prepare(root):
    from scipy.stats import qmc
    import microgen
    import scorecard
    shape=microgen.Spinodoid(k0=4*np.pi,bandwidth=3,offset=0.,anisotropy=(1.,1.,2.),
        cell_size=1.,resolution=16,seed=20261026,mode_amplitude_threshold=.05)
    f=shape._frep
    w.save(root/'spinodoid-anisotropic-coefficients.json',dict(coefficients=np.column_stack((f.modes/f.cell_size,f.amplitudes,f.phases)).tolist(),
        threshold=float(f.threshold),seed=20261026,anisotropy=[1.,1.,2.],fft_resolution=16))
    rows={}
    # The separately calibrated density ramp has its own immutable preparation.
    for name in CORE+EXTRA:
        directory=root/name;directory.mkdir()
        target=case(name,root)
        if name in CORE:
            import shutil
            for suffix in ('-samples.npz','-specimen.json'):
                p=w.HERE/'runs/g3-microgen-defaults-v4'/(name+suffix);shutil.copyfile(p,directory/p.name)
        else:
            estimates=[]
            for seed in range(4):
                count=0
                for block in range(4):
                    # Disjoint blocks from one Sobol sequence per scramble.
                    if block==0:sobol=qmc.Sobol(3,scramble=True,seed=seed)
                    probes=sobol.random_base2(18) if block==0 else sobol.random(2**18)
                    value=target.function(*(probes-.5).T)
                    count+=int(((value>=target.levels[0])&(value<=target.levels[1])).sum())
                estimates.append(count/2**20)
            walls=scorecard.wall_samples(target,count=2**15)
            assert len(walls)>=4096
            uv=qmc.Sobol(2,scramble=True,seed=119).random_base2(12)-.5;caps=[]
            for axis in range(3):
                for side in (-.5,.5):
                    p=np.empty((len(uv),3));p[:,axis]=side;p[:,[a for a in range(3) if a!=axis]]=uv
                    value=target.function(*p.T);caps.append(p[(value>=target.levels[0])&(value<=target.levels[1])])
            np.savez_compressed(directory/(name+'-samples.npz'),walls=walls,caps=np.concatenate(caps))
            w.save(directory/(name+'-specimen.json'),dict(volume=target.volume if target.volume is not None else float(np.mean(estimates)),
                volume_sobol_estimates=estimates,volume_scramble_spread=float(np.ptp(estimates)),volume_samples_per_scramble=2**20,
                wall_samples=len(walls),definition=name,periodic_axes=[periodic(name)]*3))
        probes=np.random.default_rng(20261026).uniform(-.5,.5,(4096,3));periodic_errors=[]
        for axis in range(3):
            a,b=probes.copy(),probes.copy();a[:,axis]=-.5;b[:,axis]=.5
            periodic_errors.append(float(abs(target.function(*a.T)-target.function(*b.T)).max()))
        if periodic(name):
            if name=='smooth-lattice':
                for axis in range(3):
                    a,b=probes.copy(),probes.copy();a[:,axis]=-.5;b[:,axis]=.5
                    va,vb=target.function(*a.T),target.function(*b.T)
                    np.testing.assert_array_equal(va<=0,vb<=0)
                    assert abs(np.minimum(va,.02)-np.minimum(vb,.02)).max()<1e-10
            else:assert max(periodic_errors)<1e-10
        else:assert periodic_errors[2]>1e-3
        shape,part=actual_shape(name,root,20)
        if name.startswith('graded-gyroid'):
            raw=shape.raw_field(*probes.T);offset=shape._offset_func(*probes.T)
            np.testing.assert_allclose(raw/(offset/2),target.function(*probes.T),atol=1e-12,rtol=1e-12)
        elif name=='schwarz-d-solid':np.testing.assert_allclose(shape.raw_field(*probes.T),target.function(*probes.T),atol=1e-12)
        elif name in ('spinodoid','spinodoid-anisotropic'):
            np.testing.assert_allclose(shape.evaluate(*(probes+.5).T),target.function(*probes.T),atol=1e-12,rtol=1e-12)
        else:
            legacy=ORIGINAL_SPECIMEN(name)[0]
            np.testing.assert_array_equal(legacy.function(*probes.T),target.function(*probes.T))
        rows[name]=dict(periodic_axes=[periodic(name)]*3,field_identity_points=4096,periodic_field_errors=periodic_errors,
            fixture_hashes={p.name:w.digest(p) for p in directory.iterdir()})
    w.save(root/'result.json',dict(passed=True,cases=rows,fixed_original_four_preserved=True,
        added_cases=list(EXTRA),linear_grading_is_nonperiodic=True,full_mechanical_accuracy_not_inferred=True))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args()
    w.imports()
    from comparison_mesh_audit import controls
    w.save(args.run/'controls.json',controls())
    prepare(args.run)
