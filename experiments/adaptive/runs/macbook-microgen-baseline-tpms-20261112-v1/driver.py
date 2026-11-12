"""Reproduce the full-suite TPMS failures on exact pre-integration source."""
import io,json,os,shutil,subprocess,sys,tarfile
from pathlib import Path
sys.path[:0]=['/Users/kmarchais/dev/meshers/experiments/adaptive',str(Path.home()/'bin')]
import microgen_workflow as w
from job_parts import worker_part
from pueue_progress import Job
root=Path('/Users/kmarchais/dev/meshers');run=root/'experiments/adaptive/runs/macbook-microgen-baseline-tpms-20261112-v1';run.mkdir(exist_ok=False);shutil.copyfile(__file__,run/'driver.py')
repo=Path('/Users/kmarchais/dev/microgen');base=subprocess.check_output(['git','rev-parse','38e1912'],cwd=repo,text=True).strip();archive=subprocess.check_output(['git','archive',base,'microgen'],cwd=repo);src=run/'baseline-source';src.mkdir()
with tarfile.open(fileobj=io.BytesIO(archive)) as t:t.extractall(src,filter='data')
old=root/'experiments/adaptive/runs/macbook-microgen-full-core-20261112-v2';shutil.copytree(old/'tests',run/'tests',ignore=shutil.ignore_patterns('__pycache__'));shutil.copyfile(old/'pyproject.toml',run/'pyproject.toml')
python=root/'.local-work/overnight-release-20261112/clean-install-afc-v1/venv/bin/python';os.environ['PYTHONPATH']=str(src);os.chdir(run)
script="import microgen,pytest; from pathlib import Path; assert Path(microgen.__file__).is_relative_to(Path('baseline-source').resolve()); print(microgen.__file__); raise SystemExit(pytest.main(['-q','-ra','tests/shapes/test_tpms.py','--junitxml=tests.xml']))"
w.save(run/'manifest.json',dict(baseline_commit=base,changed_integration_commit='7fb26cd',source_mode='Exact git archive of pre-integration Microgen, selected with explicit PYTHONPATH; identical installed dependencies',source_sha256={str(p.relative_to(src)):w.digest(p) for p in src.rglob('*.py')},test_sha256=w.digest(run/'tests/shapes/test_tpms.py'),python=str(python)))
with Job(parts=['baseline-tpms']) as job:
 with worker_part(job,'baseline-tpms',memory_bytes=12*1024**3) as part:
  r=w.bounded([str(python),'-c',script],run/'tests.log',part,timeout=600,memory_bytes=12*1024**3);w.save(run/'result.json',dict(complete=True,resource=r));job.result('baseline-tpms',returncode=r['returncode'],seconds=r['wall_s'],peak_gib=r['peak_rss_bytes']/1024**3)
  assert r['returncode']==1,'Expected baseline assertion failures must be inspected, not treated as a passing suite'
