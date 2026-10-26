"""Independently minimize distances for an unchanged near-gate mesh witness."""
import argparse
import json
from pathlib import Path
import numpy as np
import microgen_workflow as w
import comparison_fixtures as fixtures


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists();w.imports()
    from box_projection import constrained_foot
    old=json.loads((a.source/'result.json').read_text());name=old['case'];target=fixtures.case(name,a.fixtures)
    mesh=np.load(a.source/'mesh.npz');points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
    original=w.specimen;w.specimen=lambda requested:(target,None,None) if requested==name else original(requested)
    try:geometry=w.geometry(a.source,name,points,cells,surface,detail_count=8,box_projection_fallback=True)
    finally:w.specimen=original
    np.testing.assert_allclose(geometry['sampled_two_sided_max'],old['geometry']['sampled_two_sided_max'],atol=1e-13)
    rows=[];finite=np.array([v for v in target.levels if np.isfinite(v)])
    for witness in geometry['worst_boundary_samples']:
        probe=np.array(witness['probe']);initial=np.array(witness['foot']);value=target.function(*probe[None].T)[0]
        level=finite[np.argmin(abs(value-finite))]
        foot,info=constrained_foot(target.function,probe,level,np.array(witness['locked']),initial)
        rows.append(dict(probe=probe.tolist(),initial_distance=witness['distance'],optimized_foot=foot.tolist(),**info))
    w.save(a.output,dict(case=name,source=str(a.source),mesh_sha256=w.digest(a.source/'mesh.npz'),
        original_sampled_distance=geometry['sampled_two_sided_max'],gate=.002,witnesses=rows,
        closest_distance_certified=False,source_result_unchanged=True,
        scope='Independent feasible constrained minimization of eight worst samples only; no reclassification or proof of global nearest distance.'))
