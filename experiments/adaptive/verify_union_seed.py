"""Independently admit the saved native nonconvex seed, without a solve."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w
from topology_coverage import boundary_topology, fixture, notched_prism


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, _, _, _, _, scorecard = w.imports()
    from facet_incidence_audit import audit
    source = args.run/'mesh.json'
    raw = json.loads(source.read_text())
    assert raw['no_fe_solve'] and raw['admission_not_claimed']
    points, cells, surface = (np.asarray(raw[name]) for name in ('points', 'tetrahedra', 'surface'))
    _, queries, exact_volume, components, euler = fixture('notched-prism', np)
    target, faces = notched_prism(np)
    probes = np.einsum('sn,fnd->fsd', scorecard.BARY, points[surface]).reshape(-1, 3)
    forward = float(scorecard.distance_to_mesh(target, faces, probes).max())
    reverse = float(scorecard.distance_to_mesh(points, surface, queries).max())
    q, _, volumes = scorecard.tet_metrics(points, cells)
    facets = audit(points, cells, surface)
    observed_topology = boundary_topology(surface, np)
    corners = points[surface]
    boundary_volume = float(np.einsum('ij,ij->i', corners[:, 0],
        np.cross(corners[:, 1], corners[:, 2])).sum()/6)
    volume_error = float(volumes.sum()/exact_volume-1)
    qualified = bool(forward <= .002 and reverse <= .002 and abs(volume_error) <= .01
        and q.min() >= .1 and volumes.min() > 0 and facets['oriented_facet_consistency_passed']
        and observed_topology == (components, euler) and abs(boundary_volume-volumes.sum()) <= 1e-12
        and len(points) <= 250000 and len(cells) <= 250000 and len(np.unique(cells)) == len(points))
    paths = [Path(__file__), w.HERE/'topology_coverage.py', w.HERE/'scorecard.py',
        w.HERE/'facet_incidence_audit.py', source]
    for path in paths[:-1]:
        destination = args.output.parent/('union-audit-input-'+path.name)
        assert not destination.exists()
        shutil.copyfile(path, destination)
    w.save(args.output, dict(qualified=qualified, nodes=len(points), elements=len(cells),
        forward_distance=forward, reverse_distance=reverse, relative_volume_error=volume_error,
        minimum_quality=float(q.min()), minimum_determinant=float(6*volumes.min()),
        facets=facets, topology=observed_topology, boundary_volume=boundary_volume,
        inputs={str(path): w.digest(path) for path in paths}, no_fe_solve=True,
        gates=dict(distance=.002, relative_volume=.01, minimum_quality=.1, node_cap=250000, tet_cap=250000),
        production_dispatch_fixed=False, comparative_timing_claim=False))
    print(f'Nonconvex seed admitted={qualified}: {len(cells)} tets, qmin={q.min():.8g}, distance={max(forward, reverse):.8g}', flush=True)
    assert qualified


if __name__ == '__main__':
    main()
