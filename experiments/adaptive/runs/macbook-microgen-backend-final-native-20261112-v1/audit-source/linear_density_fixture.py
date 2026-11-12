"""Freeze a Gyroid sheet with a calibrated linear cross-sectional density ramp."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from scipy.stats import qmc
import microgen_workflow as w

NAME='graded-gyroid-density-linear'


def width(fixtures):
    table=json.loads((fixtures/'linear-density-calibration.json').read_text())
    z=np.asarray(table['z']);half=np.asarray(table['half_width'])
    def evaluate(x,y,position):return np.interp(position,z,half)
    return evaluate


def field(fixtures):
    from comparison_fixtures import raw_gyroid
    thickness=width(fixtures)
    def evaluate(x,y,z):return raw_gyroid(x,y,z)/thickness(x,y,z)
    return evaluate


def prepare(root):
    import scorecard
    from comparison_fixtures import raw_gyroid
    w.imports()
    assert root.is_dir() and not (root/'result.json').exists()
    z=np.linspace(-.5,.5,257);target=.27+.18*(z+.5)
    uv=qmc.Sobol(2,scramble=True,seed=20261026).random_base2(18)-.5
    x,y=uv.T;xy=np.sin(2*np.pi*x)*np.cos(2*np.pi*y)
    sy,cx=np.sin(2*np.pi*y),np.cos(2*np.pi*x)
    half=[float(np.quantile(abs(xy+sy*np.cos(2*np.pi*p)+np.sin(2*np.pi*p)*cx),rho))
        for p,rho in zip(z,target,strict=True)]
    w.save(root/'linear-density-calibration.json',dict(z=z.tolist(),half_width=half,
        target_density=[.27,.45],definition='Cross-sectional area fraction of material at each physical z',
        calibration_samples_per_plane=2**18,interpolation='Piecewise linear, shared exactly by both arms',
        boundary_policy='Bounded nonperiodic for both arms'))
    scalar=field(root)
    probes=np.r_[-.5,qmc.Sobol(1,scramble=True,seed=731).random_base2(6).ravel()-.5,.5]
    validation=[]
    for seed in range(4):
        uv=qmc.Sobol(2,scramble=True,seed=1000+seed).random_base2(17)-.5
        measured=[float(np.mean(abs(raw_gyroid(uv[:,0],uv[:,1],p))<=np.interp(p,z,half))) for p in probes]
        validation.append(measured)
    validation=np.asarray(validation);desired=.27+.18*(probes+.5)
    maximum=float(abs(validation-desired).max())
    w.save(root/'density-controls.json',dict(z=probes.tolist(),target=desired.tolist(),measured=validation.tolist(),
        maximum_absolute_density_difference=maximum,admission=.001,
        independent_scrambles=4,samples_per_plane_per_scramble=2**17,
        scope='Finite independent cross-section checks, not a certified pointwise density bound'))
    assert maximum<=.001,'Density calibration control failed; do not mesh'
    case=scorecard.Case(NAME,scalar,(-1.,1.))
    estimates=[]
    for seed in range(4):
        points=qmc.Sobol(3,scramble=True,seed=seed).random_base2(20)-.5
        estimates.append(float(np.mean(abs(scalar(*points.T))<=1)))
    directory=root/NAME;directory.mkdir()
    walls=scorecard.wall_samples(case,count=2**15);assert len(walls)>=4096
    uv=qmc.Sobol(2,scramble=True,seed=119).random_base2(12)-.5;caps=[]
    for axis in range(3):
        for side in (-.5,.5):
            points=np.empty((len(uv),3));points[:,axis]=side;points[:,[a for a in range(3) if a!=axis]]=uv
            caps.append(points[abs(scalar(*points.T))<=1])
    np.savez_compressed(directory/(NAME+'-samples.npz'),walls=walls,caps=np.concatenate(caps))
    w.save(directory/(NAME+'-specimen.json'),dict(volume=float(np.mean(estimates)),
        volume_sobol_estimates=estimates,volume_scramble_spread=float(np.ptp(estimates)),
        volume_samples_per_scramble=2**20,target_density=[.27,.45],target_integrated_volume=.36,
        definition='Calibrated linear cross-sectional material fraction; immutable interpolated width',
        periodic_axes=[False]*3))
    # The benchmark manifest always retains both fixed spectral-field hashes.
    source=w.HERE/'runs/macbook-comparison-fixtures-v1/spinodoid-anisotropic-coefficients.json'
    shutil.copyfile(source,root/source.name)
    import microgen
    from microgen.shape import surface_functions
    thickness=width(root)
    shape=microgen.Tpms(surface_functions.gyroid,offset=lambda x,y,z:2*thickness(x,y,z),resolution=20)
    points=np.random.default_rng(20261026).uniform(-.5,.5,(4096,3))
    np.testing.assert_allclose(shape.raw_field(*points.T)/(shape._offset_func(*points.T)/2),scalar(*points.T),atol=1e-12,rtol=1e-12)
    w.save(root/'result.json',dict(passed=True,case=NAME,existing_eight_cases_unchanged=True,
        density_difference_maximum=maximum,field_identity_probes=4096,geometry_is_not_the_old_linear_halfwidth_case=True,
        hashes={str(p.relative_to(root)):w.digest(p) for p in [root/'linear-density-calibration.json',root/'density-controls.json',*directory.iterdir()]}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    prepare(parser.parse_args().run)
