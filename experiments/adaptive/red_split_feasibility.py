"""Local full-edge split feasibility, without claiming a conforming mesh."""

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import microgen_workflow as w
from nested_refinement import EDGES, positive_quality


def red_templates():
    reference = np.vstack((np.eye(4), np.eye(4)[EDGES].mean(axis=1)))
    corners = [[0, 4, 6, 7], [1, 5, 4, 8], [2, 6, 5, 9], [3, 7, 8, 9]]
    templates = []
    for (a, b), cycle in (((4, 9), (5, 6, 7, 8)),
            ((5, 7), (4, 6, 9, 8)), ((6, 8), (4, 5, 9, 7))):
        cells = np.array(corners+[[a, b, cycle[i], cycle[(i+1)%4]] for i in range(4)])
        xyz = reference[cells, 1:]
        negative = np.linalg.det(xyz[:, 1:]-xyz[:, :1]) < 0
        cells[negative, :2] = cells[negative, 1::-1]
        maps = reference[cells]
        np.testing.assert_allclose(np.linalg.det(maps[:, 1:, 1:]-maps[:, :1, 1:]), 1/8,
            rtol=0, atol=1e-15)
        # All choices must induce the same four triangles on every parent face.
        from collections import Counter
        faces = Counter(tuple(sorted(face)) for c in cells for face in
            ((c[1], c[2], c[3]), (c[0], c[2], c[3]), (c[0], c[1], c[3]), (c[0], c[1], c[2])))
        boundary = {face for face, count in faces.items() if count == 1}
        assert len(boundary) == 16
        if templates:
            assert boundary == previous_boundary
        previous_boundary = boundary
        templates.append(cells)
    return np.array(templates)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    w.imports()
    templates = red_templates()
    paths = [Path(__file__), w.HERE/'nested_refinement.py', args.source/'mesh.npz',
        args.source/'result.json']
    inputs = {str(p): w.digest(p) for p in paths}
    shutil.copyfile(__file__, args.output.parent/'analysis-driver.py')
    w.save(args.output.parent/'analysis-inputs.json', inputs)
    mesh = np.load(args.source/'mesh.npz')
    points, cells = mesh['points'], mesh['tetrahedra']
    allowed = np.ones((len(cells), 6), dtype=bool)
    for offset, endpoints in enumerate(EDGES):
        for endpoint in endpoints:
            xyz = points[cells[:, :4]].copy()
            xyz[:, endpoint] = points[cells[:, offset+4]]
            allowed[:, offset] &= positive_quality(xyz) >= .1+1e-12
    no_single = ~allowed.any(axis=1)
    minimum = np.empty((len(cells), 3))
    for variant, template in enumerate(templates):
        minimum[:, variant] = positive_quality(points[cells[:, template]].reshape(-1, 4, 3)).reshape(-1, 8).min(axis=1)
    best = minimum.max(axis=1)
    feasible = best >= .1+1e-12
    report = dict(source=str(args.source), inputs=inputs, elements=len(cells),
        no_legal_single_edge=int(no_single.sum()),
        no_single_but_legal_red=int((no_single & feasible).sum()),
        no_single_and_red_fails=int((no_single & ~feasible).sum()),
        all_cells_red_feasible=int(feasible.sum()),
        best_red_quality_quantiles=np.quantile(best, [0, .01, .5, 1]).tolist(),
        no_single_best_red_quality_quantiles=np.quantile(best[no_single], [0, .01, .5, 1]).tolist(),
        template_controls_passed=True, conforming_mesh_generated=False,
        admission_passed=False, reference_promoted=False,
        scope='Independent per-parent full-edge template feasibility only. No green closure, periodic closure, Tet10 child Jacobian or global budget admission is inferred.')
    w.save(args.output, report)
    assert all(w.digest(Path(p)) == h for p, h in inputs.items())
    print(f'{report["no_single_but_legal_red"]}/{report["no_legal_single_edge"]} '
        'tets blocked for every single edge admit a local full-edge template', flush=True)


if __name__ == '__main__':
    main()
