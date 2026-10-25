"""Inspect a nonperiodic physical-plane seed; never admit it as a periodic route."""

import argparse
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--construction-only', action='store_true')
    parser.add_argument('--case', choices=('periodic-wedge', 'periodic-notched-wedge'), default='periodic-wedge')
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, _, meshers, _, _, scorecard = w.imports()
    import meshers._meshers
    from facet_incidence_audit import audit
    from topology_coverage import fixture, periodic_wedge, periodic_notched_wedge
    paths = [Path(__file__), Path(w.__file__), Path(meshers.__file__),
        Path(meshers.expression.__file__), Path(meshers.intersection.__file__),
        Path(meshers._meshers.__file__), Path(scorecard.__file__),
        w.HERE/'topology_coverage.py', w.HERE/'facet_incidence_audit.py']
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    options = dict(bounds=w.BOUNDS, cells=24, periodic=(False, False, False),
        adapt=False, geometry_tolerance=.002, minimum_quality=0., max_tetrahedra=250000,
        optimize_passes=0 if args.construction_only else 1, snap=.2, threads=1)
    w.save(args.run/'preflight.json', dict(options=options,
        final_admission_gates=dict(distance=.002, absolute_relative_volume=.01,
            minimum_quality=.1, positive_jacobians=True, coherent_facets=True,
            matching_periodic_nodes_and_caps=True, maximum_nodes=250000, maximum_tetrahedra=250000),
        scope='Diagnostic nonperiodic retained-plane seed, generation quality floor0 to expose the actual seed; no periodic candidate is admitted'))
    function, queries, exact_volume, _, _ = fixture(args.case, np)
    field = meshers.compile_field(function)
    regions = None
    if args.case == 'periodic-notched-wedge':
        raw_planes, regions = field._polyhedral_regions
    else:
        raw_planes = field._halfspaces
    planes = np.asarray(raw_planes)
    assert planes.shape == (8 if regions else 6, 4)
    probes = np.random.default_rng(42).uniform(-.5, .5, (1024, 3))
    values = probes@planes[:, :3].T+planes[:, 3]
    reproduced = (np.min([values[:, [i for i in range(len(planes)) if mask & (1 << i)]].max(axis=1)
        for mask in regions], axis=0) if regions else values.max(axis=1))
    np.testing.assert_allclose(reproduced, field.values(probes), atol=3e-16)
    mesh = meshers.generate(field, **options)
    assert mesh.diagnostics['polyhedral_constraints'] == len(planes)
    if regions:
        assert mesh.diagnostics['polyhedral_regions'] == len(regions)
    p, t, f, masks = mesh.points, mesh.tetrahedra, mesh.surface, mesh.constraint_masks
    assert len(p) <= 250000 and len(t) <= 250000
    assert all(len(pair) == 0 for pair in mesh.periodic_pairs)
    facets = audit(p, t, f)
    quality, _, volumes = scorecard.tet_metrics(p, t)
    assert facets['oriented_facet_consistency_passed'] and volumes.min() > 0
    assert len(np.unique(t)) == len(p)
    target, target_faces = (periodic_notched_wedge(np) if regions else periodic_wedge(np))
    samples = np.einsum('sn,fnd->fsd', scorecard.BARY, p[f]).reshape(-1, 3)
    forward = float(scorecard.distance_to_mesh(target, target_faces, samples).max())
    reverse = float(scorecard.distance_to_mesh(p, f, queries).max())
    low = np.flatnonzero(abs(p[:, 2]+.5) <= 1e-10)
    high = np.flatnonzero(abs(p[:, 2]-.5) <= 1e-10)
    low_keys = [tuple(x) for x in np.rint(p[low, :2]/1e-10).astype(np.int64)]
    high_keys = [tuple(x) for x in np.rint(p[high, :2]/1e-10).astype(np.int64)]
    low_map = dict(zip(low_keys, low, strict=True))
    nodes_match = len(low) > 0 and len(low_map) == len(low) == len(set(high_keys)) == len(high) and set(low_keys) == set(high_keys)
    low_faces = f[np.all(abs(p[f, 2]+.5) <= 1e-10, axis=1)]
    high_faces = f[np.all(abs(p[f, 2]-.5) <= 1e-10, axis=1)]
    correspondence = dict(low_nodes=len(low), high_nodes=len(high), nodes_match=nodes_match,
        low_caps=len(low_faces), high_caps=len(high_faces), metadata_pairs_declared=False)
    signatures = [np.column_stack((planes[:, :2], planes[:, 3]+z*planes[:, 2])) for z in (-.5, .5)]
    permutation = []
    for signature in signatures[0]:
        matches = np.flatnonzero(np.max(abs(signatures[1]-signature), axis=1) <= 1e-13)
        permutation.append(int(matches[0]) if len(matches) == 1 else None)
    unique_permutation = None not in permutation and len(set(permutation)) == len(planes)
    if nodes_match:
        paired = np.array([low_map[key] for key in high_keys])
        mapping = np.arange(len(p)); mapping[high] = paired
        a = {tuple(sorted(face)) for face in low_faces}
        b = {tuple(sorted(face)) for face in mapping[high_faces]}
        correspondence.update(caps_match=bool(a) and a == b,
            low_only_caps=len(a-b), high_only_caps=len(b-a),
            maximum_coordinate_error=float(abs(p[high]-p[paired]-(0, 0, 1)).max()),
            raw_mask_mismatch_count=(int(np.sum(masks[high] != masks[paired])) if masks is not None else None))
        if unique_permutation and masks is not None:
            converted = np.array([sum(1 << j for i, j in enumerate(permutation) if int(mask) & (1 << i)) for mask in masks[paired]], dtype=np.uint64)
            correspondence['permuted_mask_mismatch_count'] = int(np.sum(converted != masks[high]))
    if unique_permutation and regions:
        mapped_regions = [sum(1 << j for i,j in enumerate(permutation) if mask & (1 << i)) for mask in regions]
        correspondence.update(permuted_region_masks=mapped_regions,
            region_trace_set_matches=sorted(mapped_regions) == sorted(regions),
            exact_plane_trace_match=bool(np.array_equal(signatures[0], signatures[1][permutation])))
    np.savez_compressed(args.run/'mesh.npz', points=p, tetrahedra=t, surface=f,
        constraint_masks=masks if masks is not None else np.empty(0, dtype=np.uint64), planes=planes)
    np.savez_compressed(args.run/'reference.npz', points=queries)
    w.save(args.run/'result.json', dict(case=args.case, qualified=False,
        region_masks=regions, constraint_masks_available=masks is not None,
        diagnostic_nonperiodic_seed=True, native_periodic_request=False, no_fe_solve=True,
        reference_promoted=False, comparative_timing_claim=False, options=options,
        nodes=len(p), elements=len(t), minimum_quality=float(quality.min()),
        minimum_jacobian=float(6*volumes.min()), forward_distance=forward, reverse_distance=reverse,
        relative_volume_error=float(volumes.sum()/exact_volume-1), facets=facets,
        end_correspondence=correspondence, low_to_high_plane_trace_permutation=permutation,
        unique_trace_permutation=unique_permutation, trace_signature_tolerance=1e-13,
        mesh_sha256=w.digest(args.run/'mesh.npz'), reference_sha256=w.digest(args.run/'reference.npz'),
        scope='Physical affine seed with generator periodicity off; independently inspect end geometry and labels before designing a periodic native change; never a periodic admission'))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    print(f'Physical seed: {len(t)} tets, qmin {quality.min():.6f}, end nodes match {nodes_match}, caps match {correspondence.get("caps_match")}; no periodic candidate admitted', flush=True)


if __name__ == '__main__':
    main()
