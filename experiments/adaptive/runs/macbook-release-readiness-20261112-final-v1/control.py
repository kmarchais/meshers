import hashlib,importlib.util,json,sys
from pathlib import Path
from vtkmodules.vtkCommonCore import vtkSMPTools
vtkSMPTools.SetBackend('Sequential')
import microgen
from microgen import Sphere,Tpms,Spinodoid
from microgen.shape.surface_functions import gyroid
assert Path(microgen.__file__).is_relative_to(Path(sys.prefix))
assert importlib.util.find_spec('OCP') is None
optional=importlib.util.find_spec('meshers') is not None
assert optional == (sys.argv[1]=='optional')
rows={}
for name,shape in [('sphere',Sphere(radius=.35)),('gyroid',Tpms(gyroid,offset=1.1,resolution=8)),('spinodoid',Spinodoid(offset=0,resolution=8,seed=4))]:
 a=shape.generate_volume_mesh(resolution=8);b=shape.generate_volume_mesh(backend='vtk',resolution=8)
 def digest(g):
  return hashlib.sha256(g.points.tobytes()+g.cells.tobytes()).hexdigest()
 assert digest(a)==digest(b)
 rows[name]={'sha256':digest(a),'points':a.n_points,'cells':a.n_cells}
p=Path(sys.argv[2])
if optional:assert rows==json.loads(p.read_text())['meshes']
record={'microgen':microgen.__file__,'meshers_available':optional,'OCP_available':False,'meshes':rows}
(p if not optional else p.with_name('optional-vtk-control.json')).write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
