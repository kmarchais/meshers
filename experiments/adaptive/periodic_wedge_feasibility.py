"""Audit a manufactured exact-domain mesh; this is not a Meshers candidate."""

import argparse
from pathlib import Path
import shutil

import numpy as np
from vtkmodules.vtkCommonCore import vtkSMPTools

import microgen_workflow as w
import scorecard
from facet_incidence_audit import audit
from topology_coverage import audit_column_periodicity, boundary_topology, fixture, periodic_wedge, periodic_notched_wedge


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--case', choices=('periodic-wedge', 'periodic-notched-wedge'), default='periodic-wedge')
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    assert vtkSMPTools.GetBackend() == 'Sequential' and vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    paths = [Path(__file__), Path(w.__file__), Path(scorecard.__file__),
        w.HERE/'topology_coverage.py', w.HERE/'facet_incidence_audit.py']
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    points, surface = (periodic_notched_wedge(np) if args.case == 'periodic-notched-wedge' else periodic_wedge(np))
    center = len(points)
    interior = (-.075, -.075, 0.) if args.case == 'periodic-notched-wedge' else (0., 0., 0.)
    points = np.vstack((points, interior))
    cells = np.column_stack((np.full(len(surface), center), surface))
    quality, _, volumes = scorecard.tet_metrics(points, cells)
    facets = audit(points, cells, surface)
    periodic = audit_column_periodicity(points, surface, np)
    field, queries, exact_volume, components, euler = fixture(args.case, np)
    bary = np.array([(i/4, j/4, 1-(i+j)/4) for i in range(5) for j in range(5-i)])
    probes = np.einsum('sn,fnd->fsd', bary, points[surface]).reshape(-1, 3)
    closed_residual = np.maximum(field(*probes.T), abs(probes[:, 2])-.5)
    assert abs(closed_residual).max() < 1e-15
    reverse = float(scorecard.distance_to_mesh(points, surface, queries).max())
    assert reverse <= .002 and quality.min() >= .1 and volumes.min() > 0
    assert abs(volumes.sum()/exact_volume-1) <= .01
    assert facets['oriented_facet_consistency_passed']
    assert boundary_topology(surface, np) == (components, euler)
    faces = points[surface]
    boundary_volume = np.einsum('ij,ij->i', faces[:, 0], np.cross(faces[:, 1], faces[:, 2])).sum()/6
    assert abs(boundary_volume-volumes.sum()) < 1e-12
    assert len(np.unique(cells)) == len(points) <= 250000 and len(cells) <= 250000
    np.savez_compressed(args.run/'mesh.npz', points=points, tetrahedra=cells, surface=surface)
    np.savez_compressed(args.run/'reference.npz', points=queries)
    w.save(args.run/'result.json', dict(case=args.case, qualified=True,
        manufactured_feasibility_control=True, native_generator_called=False,
        no_fe_solve=True, reference_promoted=False, comparative_timing_claim=False,
        nodes=len(points), elements=len(cells), minimum_quality=float(quality.min()),
        minimum_jacobian=float(6*volumes.min()), volume=float(volumes.sum()),
        exact_volume=exact_volume, relative_volume_error=float(volumes.sum()/exact_volume-1),
        sampled_reverse_distance=reverse, maximum_exact_surface_residual=float(abs(closed_residual).max()),
        periodic=periodic, facets=facets, boundary_components=components, boundary_euler=euler,
        volume_balance=float(abs(boundary_volume-volumes.sum())),
        mesh_sha256=w.digest(args.run/'mesh.npz'), reference_sha256=w.digest(args.run/'reference.npz'),
        scope='Exact star-shaped-cell center-star feasibility with matching periodic caps; not a native meshing result, workflow or stress/continuum certificate'))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    print(f'Manufactured {args.case} passes: {len(cells)} tets, minimum quality {quality.min():.6f}; native failure remains open', flush=True)


if __name__ == '__main__':
    main()
