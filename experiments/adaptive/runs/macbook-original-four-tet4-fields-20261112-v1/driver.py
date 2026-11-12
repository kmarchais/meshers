"""Use the validated saved Tet4 worker on the remaining original specimens."""
import hashlib,json,os,shutil,sys
from pathlib import Path
sys.path[:0]=['/Users/kmarchais/dev/meshers/experiments/adaptive',str(Path.home()/'bin')]
import microgen_workflow as w
from job_parts import worker_part
from pueue_progress import Job
root=Path('/Users/kmarchais/dev/meshers');run=root/'experiments/adaptive/runs/macbook-original-four-tet4-fields-20261112-v1';run.mkdir(exist_ok=False)
shutil.copyfile(__file__,run/'driver.py')
old=root/'experiments/adaptive/runs/macbook-gyroid-tet4-stiffness-tradeoff-20261112-v2/driver.py';text=old.read_text();needle='assert original["case"] == "gyroid" and original["periodic_required"]';assert needle in text;text=text.replace(needle,'assert original["case"] in w.CASES and original["periodic_required"]');(run/'worker.py').write_text(text)
os.environ['PYTHONPATH']=str(w.HERE)
names=('schwarz-solid','spinodoid','smooth-lattice');source=root/'experiments/adaptive/runs/macbook-microgen-backend-release-final-20261112-v1'
inputs={str(source/f'repeat-0-{name}-native'/file):w.digest(source/f'repeat-0-{name}-native'/file) for name in names for file in ('result.json','mesh.npz')}
helpers=('microgen_workflow.py','job_parts.py','amg_reference.py','stress_recovery.py','tet10_integration.py','quadrature_audit.py','scorecard.py');snap=run/'audit-source';snap.mkdir()
for name in helpers:shutil.copyfile(w.HERE/name,snap/name)
for name in ('periodic.py','homogenization.py'):shutil.copyfile(root.parent/'periodic-microgen/experiments/periodic-fea'/name,snap/name)
w.save(run/'manifest.json',dict(inputs=inputs,source_sha256={p.name:w.digest(p) for p in snap.iterdir()},worker_sha256=w.digest(run/'worker.py'),numerical_threads=1,memory_bytes=12*1024**3,timeout_s=600,protocol='Six independent engineering strain loads; exact fixed meshes; reaction/energy/periodic and positive quadrature checks. No reference or workflow timing promotion.',existing_gyroid='experiments/adaptive/runs/macbook-gyroid-tet4-stiffness-tradeoff-20261112-v2/uniform',new_fea=True,reference_promoted=False))
records=[]
with Job(parts=list(names)) as job:
 for name in names:
  with worker_part(job,name,memory_bytes=12*1024**3) as part:
   resource=w.bounded([sys.executable,str(run/'worker.py'),'--run',str(run/name),'--worker',str(source/f'repeat-0-{name}-native')],run/(name+'.log'),part,timeout=600,memory_bytes=12*1024**3)
   w.save(run/(name+'-resource.json'),resource);assert resource['returncode']==0,(name,'Read saved log')
   row=json.loads((run/name/'result.json').read_text());assert row['field_audit_passed'] and row['positive_quadrature_passed'];records.append(dict(case=name,directory=str(run/name),result_sha256=w.digest(run/name/'result.json'),resource=resource))
   job.result(name,loads=6,seconds=resource['wall_s'],peak_gib=resource['peak_rss_bytes']/1024**3)
assert all(w.digest(Path(p))==digest for p,digest in inputs.items())
w.save(run/'assessment.json',dict(complete=True,records=records,new_fea=True,reference_promoted=False,scope='Original fixed specimens saved-mesh FE consistency only; no continuum/stress certification, no live workflow timing or matched accuracy comparator claim'))
