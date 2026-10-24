from pathlib import Path
import argparse,shutil,sys
sys.path.insert(0,str(Path.cwd()/'experiments/adaptive'))
import microgen_workflow as w
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
np,_,meshers,_,_,_=w.imports();c,s=np.cos(.37),np.sin(.37)
def field(x,y,z):
 u=c*(x-.02)+s*(y+.03);v=-s*(x-.02)+c*(y+.03)
 return np.maximum(np.maximum(np.maximum(abs(u)-.35,abs(v)-.35),abs(z-.01)-.3),np.minimum(u-.05,v-.05))
shutil.copyfile(__file__,a.run/'analysis-driver.py')
m=meshers.generate(field,bounds=(-.5,.5)*3,cells=24,snap=0,optimize_passes=2,minimum_quality=.1,geometry_tolerance=.002,max_tetrahedra=250000)
w.save(a.run/'result.json',dict(nodes=len(m.points),elements=len(m.tetrahedra),diagnostics=m.diagnostics,options=dict(snap=0,optimize_passes=2,minimum_quality=.1,geometry_tolerance=.002),no_fe_solve=True))
print('Two-pass rotated control passes:',len(m.tetrahedra),m.diagnostics['minimum_mmg_quality'])
