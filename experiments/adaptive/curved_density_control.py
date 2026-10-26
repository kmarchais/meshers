"""Regenerate the exact rejected uniform recipe with the new automatic evaluator."""
import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w
from comparison_benchmark import worker


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args();w.imports()
    source=w.HERE/'runs/macbook-linear-density-comparison-pilot-v1/repeat-0-meshers-curved-tet10'
    fixtures=w.HERE/'runs/macbook-linear-density-fixtures-v1'
    previous=json.loads((source/'result.json').read_text())
    assert not previous['qualified'] and previous['geometry']['failed_gates']==['geometry']
    assert w.digest(source/'mesh.npz')==previous['mesh_sha256']
    name=previous['case'];config=dict(previous['config'],name='meshers-uniform16-auto-interp',compile=True,native_interp_auto=True)
    directory=args.run/'mesh';directory.mkdir()
    for path in (fixtures/name).iterdir():shutil.copyfile(path,directory/path.name)
    import meshers
    inputs={str(p):w.digest(p) for p in (Path(__file__),Path(worker.__code__.co_filename),Path(w.__file__),
        Path(meshers.expression.__file__),source/'result.json',source/'mesh.npz',
        fixtures/'linear-density-calibration.json',w.ROOT/'handoff/macbook/curved-density-control-selection-20261026.json')}
    w.save(args.run/'inputs.json',inputs)
    worker(directory,name,config,fixtures)
    result=json.loads((directory/'result.json').read_text())
    assert result['status'] in ('admitted','admission-rejected'),result
    assert result['diagnostics']['evaluator']=='compiled' and result['diagnostics']['callback_calls']==0
    assert all(w.digest(Path(p))==h for p,h in inputs.items())
    w.save(args.run/'result.json',dict(source=str(source),source_config=previous['config'],config=config,
        generated_directory=str(directory),admitted=result['qualified'],nodes=result['nodes'],elements=result['elements'],
        geometry=result['geometry'],quality_quantiles=result['quality_quantiles'],generation_pilot_s=result['generation_s'],
        inputs=inputs,new_FE_solve=False,reference_promoted=False,
        scope='One fixed uniform16 recipe with public automatic interpolation; original callback failure remains frozen. No repeated timing or inherited FE result.'))
