import sys,json,shutil,subprocess,os
from pathlib import Path
sys.path[:0]=[str(Path.home()/'bin'),'/Users/kmarchais/dev/meshers/experiments/adaptive']
from pueue_progress import Job
from job_parts import worker_part
import microgen_workflow as w
root=Path('/Users/kmarchais/dev/meshers/.local-work/overnight-release-20261112')
run=root/'clean-install-afc-v1';run.mkdir(exist_ok=False)
env=run/'venv';python=env/'bin/python';uv='/Users/kmarchais/.local/bin/uv'
microgen=Path('/Users/kmarchais/dev/microgen')
wheel=root/'ci-artifacts-afc/meshers-0.2.0-cp310-abi3-macosx_11_0_arm64.whl'
assert wheel.exists()
constraints=run/'constraints.txt';constraints.write_text('numpy==2.5.4\nscipy==1.18.1\npyvista==0.49.1\nvtk==9.7.1\ngmsh==4.15.2\nmeshio==5.3.5\nautograd==1.9.1\npytest==9.1.1\nh5py==3.16.0\n')
for n in ['test_meshers_backend.py','test_mmg_backend.py','test_cad_optional.py']:shutil.copyfile(microgen/'tests'/n,run/n)
shutil.copyfile(root/'clean_vtk_control.py',run/'control.py')
os.chdir(run)
os.environ['PATH']='/Users/kmarchais/dev/meshers/.local-work/macbook-runtime/venv/lib/python3.12/site-packages/mmgpy/bin'+os.pathsep+os.environ['PATH']
os.environ.pop('PYTHONPATH',None)
commands=[
 ('create',[uv,'venv','--python','3.12.9',str(env)]),
 ('core-install',[uv,'pip','install','--python',str(python),'-c',str(constraints),str(microgen),'pytest','h5py']),
 ('core-control',[str(python),'control.py','core','core-vtk-control.json']),
 ('optional-install',[uv,'pip','install','--python',str(python),'-c',str(constraints),str(microgen)+'[meshers]',str(wheel)]),
 ('optional-control',[str(python),'control.py','optional','core-vtk-control.json']),
 ('integration-tests',[str(python),'-m','pytest','-q','test_meshers_backend.py','test_mmg_backend.py','test_cad_optional.py','--junitxml=tests.xml']),
 ('dependencies',[uv,'pip','check','--python',str(python)])]
results={}
with Job(parts=[n for n,_ in commands]) as job:
 for name,cmd in commands:
  with worker_part(job,name,memory_bytes=12*1024**3) as part:
   result=w.bounded(cmd,run/(name+'.log'),part,timeout=600,memory_bytes=12*1024**3)
   results[name]=result
   (run/'result.json').write_text(json.dumps({'complete':False,'steps':results},indent=2)+'\n')
   job.result(name,returncode=result['returncode'],seconds=result['wall_s'],peak_gib=result['peak_rss_bytes']/1024**3)
   assert result['returncode']==0,(name,(run/(name+'.log')).read_text()[-9000:])
 (run/'result.json').write_text(json.dumps({'complete':True,'steps':results,'wheel_sha256':w.digest(wheel),'microgen_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=microgen,text=True).strip(),'cad_installed':False},indent=2)+'\n')
