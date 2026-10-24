import argparse,json,shutil,sys
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'experiments/adaptive'))
import microgen_workflow as w
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
np,_,meshers,_,_,scorecard=w.imports()
shutil.copyfile(__file__,a.run/'analysis-driver.py')
c,s=np.cos(.37),np.sin(.37)
def field(x,y,z):
 u=c*(x-.02)+s*(y+.03);v=-s*(x-.02)+c*(y+.03)
 return np.maximum(np.maximum(np.maximum(abs(u)-.35,abs(v)-.35),abs(z-.01)-.3),np.minimum(u-.05,v-.05))
f=meshers.compile_field(field);planes,regions=f._polyhedral_regions
m=meshers.generate(f,bounds=(-.5,.5)*3,cells=24,snap=0,optimize_passes=1,minimum_quality=0,geometry_tolerance=.002,max_tetrahedra=250000)
planes=np.array(planes);res=m.points@planes[:,:3].T+planes[:,3];active=np.abs(res)<1e-11
support=np.zeros(len(m.points),dtype=np.uint32)
for face in m.surface:
 mask=sum((1<<j) for j in range(8) if active[face,j].all())
 support[face]|=mask
raw=(active.astype(np.uint32)*(1<<np.arange(8,dtype=np.uint32))).sum(axis=1)
false=(np.array([int(v).bit_count() for v in raw])>=2)&(np.array([int(v).bit_count() for v in support])<2)
quality,_,_=scorecard.tet_metrics(m.points,m.tetrahedra);bad=quality<.1
np.savez_compressed(a.run/'mesh.npz',points=m.points,tetrahedra=m.tetrahedra,surface=m.surface)
row=dict(nodes=len(m.points),elements=len(m.tetrahedra),minimum_quality=float(quality.min()),bad_cells=int(bad.sum()),artificially_frozen_vertices=int(false.sum()),bad_cells_touching_artificial_features=int(np.any(false[m.tetrahedra[bad]],axis=1).sum()),interior_artificial_features=int(np.sum(false&(field(*m.points.T)<-1e-10))),mesh_sha256=w.digest(a.run/'mesh.npz'),native_sha256=w.digest(Path(meshers._meshers.__file__)),no_admission_claim=True)
w.save(a.run/'result.json',row);print(json.dumps(row))
