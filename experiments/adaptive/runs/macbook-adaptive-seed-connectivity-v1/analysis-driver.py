from pathlib import Path
import argparse,json,math,shutil,sys
sys.path.insert(0,str(Path.cwd()/'experiments/adaptive'))
import microgen_workflow as w
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);args=p.parse_args()
np,_,meshers,_,_,scorecard=w.imports()
shutil.copyfile(__file__,args.run/'analysis-driver.py')
cosine,sine=np.cos(.37),np.sin(.37)
def field(x,y,z):
 u=cosine*(x-.02)+sine*(y+.03);v=-sine*(x-.02)+cosine*(y+.03)
 return (abs(u)+abs(v)+abs(z-.01)-.3)/np.sqrt(3)
raw=meshers.compile_field(field)._halfspaces
unit=[]
for row in raw:
 scale=max(abs(v) for v in row[:3]);q=[v/scale for v in row]
 length=math.sqrt(sum(v*v for v in q[:3]));unit.append([v/length for v in q])
results={}
for label,planes in [('raw',raw),('unit',unit)]:
 funcs={}
 for index,row in enumerate(planes):
  def plane(x,y,z,row=row):return row[0]*x+row[1]*y+row[2]*z+row[3]
  funcs[f'plane{index}']=plane
 mesh=meshers.generate_intersection(funcs,bounds=(-.5,.5)*3,cells=24,minimum_quality=0,optimize_passes=1,max_tetrahedra=250000)
 points,cells=mesh.points,mesh.tetrahedra
 keys=np.floor((points+.5)*1e11+.5).astype(np.int64)
 unique,inverse=np.unique(keys,axis=0,return_inverse=True);inverse=inverse.ravel()
 sorted_cells=np.sort(inverse[cells],axis=1)
 bad=np.any(np.diff(sorted_cells,axis=1)==0,axis=1)
 quality,_,_=scorecard.tet_metrics(points,cells)
 row=dict(nodes=len(points),elements=len(cells),coordinate_buckets=len(unique),merged_nodes=len(points)-len(unique),collapsed_tetrahedra=int(bad.sum()),minimum_quality=float(quality.min()))
 if bad.any():
  t=int(np.flatnonzero(bad)[0]);row['first_collapsed_cell']=cells[t].tolist();row['coordinates']=points[cells[t]].tolist()
 np.savez_compressed(args.run/f'{label}-seed.npz',points=points,tetrahedra=cells,surface=mesh.surface)
 row['mesh_sha256']=w.digest(args.run/f'{label}-seed.npz');results[label]=row
w.save(args.run/'result.json',dict(cases=results,raw_planes=raw,unit_planes=unit));print(json.dumps(results))
