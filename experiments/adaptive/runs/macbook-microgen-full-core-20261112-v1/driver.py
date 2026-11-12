"""Exercise the full existing Microgen suite against the clean installed optional backend."""
import json,os,shutil,sys,subprocess
from pathlib import Path
sys.path[:0]=['/Users/kmarchais/dev/meshers/experiments/adaptive',str(Path.home()/'bin')]
import microgen_workflow as w
from pueue_progress import Job
from job_parts import worker_part
root=Path('/Users/kmarchais/dev/meshers');work=root/'.local-work/overnight-release-20261112';run=root/'experiments/adaptive/runs/macbook-microgen-full-core-20261112-v1';run.mkdir(exist_ok=False);shutil.copyfile(__file__,run/'driver.py')
repo=Path('/Users/kmarchais/dev/microgen');tests=run/'tests';shutil.copytree(repo/'tests',tests,ignore=shutil.ignore_patterns('__pycache__','.pytest_cache','data'));shutil.copyfile(repo/'pyproject.toml',run/'pyproject.toml')
python=work/'clean-install-afc-v1/venv/bin/python';os.environ['PATH']='/Users/kmarchais/dev/meshers/.local-work/macbook-runtime/venv/lib/python3.12/site-packages/mmgpy/bin'+os.pathsep+os.environ['PATH'];os.environ.pop('PYTHONPATH',None);os.chdir(run)
inputs={str(p.relative_to(run)):w.digest(p) for p in tests.rglob('*.py')}
w.save(run/'manifest.json',dict(microgen_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),release_commit='afc49e583a0450454b02440557f39b89bbc39ca6',python=str(python),tests_sha256=inputs,cad_available=False,workers=1,numerical_threads=1,timeout_s=1200,memory_bytes=12*1024**3,source_tree_on_sys_path=False))
with Job(parts=['full-suite']) as job:
 with worker_part(job,'full-suite',memory_bytes=12*1024**3) as part:
  r=w.bounded([str(python),'-m','pytest','-q','-ra','tests','--junitxml=tests.xml'],run/'tests.log',part,timeout=1200,memory_bytes=12*1024**3)
  w.save(run/'result.json',dict(complete=True,resource=r));job.result('full-suite',returncode=r['returncode'],seconds=r['wall_s'],peak_gib=r['peak_rss_bytes']/1024**3)
  assert all(w.digest(run/p)==d for p,d in inputs.items())
  if r['returncode']:print((run/'tests.log').read_text()[-12000:]);raise SystemExit(r['returncode'] if r['returncode']>0 else 1)
