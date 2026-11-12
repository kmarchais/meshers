"""Verify downloaded CI artifact digests, contents, and source against the release checkout."""
import argparse,email,hashlib,io,json,tarfile,zipfile
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--metadata',type=Path,required=True);p.add_argument('--directory',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();data=json.loads(a.metadata.read_text());records=[]
sha=lambda b:hashlib.sha256(b).hexdigest()
for artifact in data['artifacts']:
 assert artifact['workflow_run']['head_sha']==data['commit']
 archive=a.directory/(str(artifact['id'])+'.zip');content=archive.read_bytes();assert sha(content)==artifact['digest'].removeprefix('sha256:')
 with zipfile.ZipFile(io.BytesIO(content)) as z:
  assert len(z.namelist())==1
  name=z.namelist()[0];assert '/' not in name and '\\' not in name
  payload=z.read(name);(a.directory/name).write_bytes(payload)
 record=dict(artifact_id=artifact['id'],name=name,archive_sha256=sha(content),file_sha256=sha(payload))
 if name.endswith('.whl'):
  with zipfile.ZipFile(io.BytesIO(payload)) as wheel:
   names=wheel.namelist();record['files']=len(names)
   metadata=email.message_from_bytes(wheel.read(next(n for n in names if n.endswith('/METADATA'))));assert metadata['Version']=='0.2.0';assert metadata['Requires-Python']=='>=3.10';assert 'numpy>=1.23' in metadata.get_all('Requires-Dist')
   record['tags']=wheel.read(next(n for n in names if n.endswith('/WHEEL'))).decode().splitlines()
   assert any(n.endswith('/THIRDPARTY.json') for n in names)
   matched={}
   for n in names:
    if n.startswith('meshers/') and n.endswith(('.py','.pyi')):
     source=a.release/'crates/meshers-python/python'/n;b=wheel.read(n)
     assert b.replace(b'\r\n',b'\n')==source.read_bytes().replace(b'\r\n',b'\n'),(name,n)
     matched[n]=sha(b)
   assert len(matched)==6
   record['python_source_sha256']=matched
 else:
  assert name.endswith('.tar.gz')
  with tarfile.open(fileobj=io.BytesIO(payload)) as tar:
   names=tar.getnames();record['files']=len(names);matched={}
   assert not any(any(s in n for s in ['/experiments/','/.audit/','/handoff/','volume_error_probe','profile_pmy','meshers-integration.md']) for n in names)
   for member in tar.getmembers():
    relative=Path(member.name).relative_to('meshers-0.2.0')
    if relative.suffix not in ('.rs','.py','.pyi'):continue
    if relative.parts[0] in ('meshers-core','meshers-python'):source=a.release/'crates'/relative
    elif relative.parts[0]=='python':source=a.release/'crates/meshers-python'/relative
    else:raise AssertionError(relative)
    b=tar.extractfile(member).read();assert b==source.read_bytes(),relative;matched[str(relative)]=sha(b)
   assert len(matched)>60
   record['source_sha256']=matched;record['study_and_editorial_excluded']=True
 records.append(record)
assert len(records)==5
assert sum(r['name'].endswith('.whl') for r in records)==4
result=dict(commit=data['commit'],passed=True,artifacts=records)
a.output.write_text(json.dumps(result,indent=2)+'\n');print('Verified four platform wheels and isolated source distribution; all code matches release checkout')
