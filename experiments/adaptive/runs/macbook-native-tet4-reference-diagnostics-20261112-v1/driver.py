"""Assess fresh native tensors against the retained historical reference evidence."""
import json,shutil,sys
from pathlib import Path
sys.path[:0]=['/Users/kmarchais/dev/meshers/experiments/adaptive',str(Path.home()/'bin')]
import microgen_workflow as w
from pueue_progress import Job
from job_parts import worker_part
w.imports()
import numpy as np
from periodic import tensor_difference
from cubic_refinement import MODEL_REPORTS
root=Path('/Users/kmarchais/dev/meshers');run=root/'experiments/adaptive/runs/macbook-native-tet4-reference-diagnostics-20261112-v1';run.mkdir(exist_ok=False);shutil.copyfile(__file__,run/'driver.py')
sources={
'gyroid-uniform':root/'experiments/adaptive/runs/macbook-gyroid-tet4-stiffness-tradeoff-20261112-v2/uniform/result.json',
'gyroid-adaptive':root/'experiments/adaptive/runs/macbook-gyroid-tet4-stiffness-tradeoff-20261112-v2/adaptive/result.json',
'schwarz-solid':root/'experiments/adaptive/runs/macbook-original-four-tet4-fields-20261112-v1/schwarz-solid/result.json',
'smooth-lattice':root/'experiments/adaptive/runs/macbook-original-four-tet4-fields-20261112-v1/smooth-lattice/result.json'}
reference=w.HERE/'runs/g3-microgen-defaults-v4/summary.json';inputs={str(p):w.digest(p) for p in [reference,*sources.values(),w.HERE/'stiffness_pair.py',w.HERE/'cubic_refinement.py']}
orders=[2,3,4,5,6]
w.save(run/'manifest.json',dict(inputs=inputs,declared_orders=orders,reference_justification='Existing qualified curved Tet10 16/20 extrapolations and documented Gyroid linear64/80 cross-check in RESULTS.md. Orders 2 through 6 are declared before this assessment. Lattice retains all 24 historical model tensors; no lattice reference is promoted.',new_reference=False,target=.01,spinodoid='No justified reference selected here; FE consistency only'))
records={}
with Job(parts=list(sources)) as job:
 for name,p in sources.items():
  with worker_part(job,name,memory_bytes=12*1024**3) as part:
   if name!='smooth-lattice':
    resource=w.bounded([sys.executable,str(w.HERE/'stiffness_pair.py'),'--first',str(p),'--reference-summary',str(reference),'--orders',*map(str,orders),'--output',str(run/(name+'.json'))],run/(name+'.log'),part,timeout=120,memory_bytes=12*1024**3)
    assert resource['returncode']==0,(name,'Read worker log');r=json.loads((run/(name+'.json')).read_text());records[name]=dict(base_p4=r['difference']['worst_relative_strain_energy_difference'],declared_model_range=r['declared_order_difference_range'],target_passes_all=r['passes_target_under_all_declared_models'])
   else:
    row=json.loads(p.read_text());assert row['qualified'] and row['field_audit_passed'] and row['positive_quadrature_passed'];tensor=np.asarray(row['stiffness']);models=[]
    for n in MODEL_REPORTS:
     path=w.HERE/'runs'/n;report=json.loads(path.read_text());assert report['case']=='smooth-lattice';inputs[str(path)]=w.digest(path)
     for i,m in enumerate(report['models']):models.append(dict(source=n,index=i,difference=tensor_difference(tensor,np.asarray(m['stiffness']))))
    assert len(models)==24;errors=[m['difference']['worst_relative_strain_energy_difference'] for m in models];records[name]=dict(retained_models=24,declared_model_range=[min(errors),max(errors)],target_passes_all=max(errors)<=.01);w.save(run/(name+'.json'),dict(models=models,source_sha256=inputs,new_reference=False,reference_promoted=False))
   job.result(name,worst_disagreement=records[name]['declared_model_range'][1],target_passes_all=records[name]['target_passes_all'])
assert all(w.digest(Path(p))==digest for p,digest in inputs.items())
w.save(run/'assessment.json',dict(complete=True,records=records,inputs=inputs,spinodoid_reference_assessed=False,reference_promoted=False,scope='Previously retained empirical reference comparisons under declared extrapolation assumptions, not continuum bounds or new actual workflow comparisons'))
