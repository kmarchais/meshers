"""Compare admitted lattice mesh densities before selecting a count-matched trial."""

import argparse
import json
import math
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--actual', type=Path, required=True)
    parser.add_argument('--actual-assessment', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    import numpy as np
    import scorecard

    proof = json.loads(args.actual_assessment.read_text())
    assert proof['fresh_actual_generation'] and not proof['saved_seed_replay']
    assert proof['passes_target_all_models'] and proof['retained_model_count'] == 24
    assert Path(proof['source']).resolve() == args.actual.resolve()
    paths = [Path(__file__), Path(w.__file__), Path(scorecard.__file__), args.actual_assessment]
    paths.extend(p/n for p in (args.native, args.actual) for n in ('mesh.npz', 'result.json'))
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    shutil.copyfile(__file__, args.run/'analysis-driver.py')
    result = {}
    fractions = [0, .01, .1, .5, .9, .99, 1]
    for label, source in (('native', args.native), ('actual', args.actual)):
        row = json.loads((source/'result.json').read_text())
        assert row['case'] == 'smooth-lattice' and row['qualified'] and row['order'] == 2
        assert row['geometry']['geometry_admitted'] and max(row['nodes'], row['elements']) <= w.CAP
        mesh = np.load(source/'mesh.npz')
        points, cells, surface = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
        assert cells.shape == (row['elements'], 10) and surface.shape[1] == 6
        corners = points[cells[:, :4]]
        lengths = np.array([np.linalg.norm(corners[:, i]-corners[:, j], axis=1)
            for i in range(4) for j in range(i+1, 4)]).T
        rms = np.sqrt((lengths**2).mean(axis=1))
        q, _, volume = scorecard.tet_metrics(points, cells[:, :4])
        assert volume.min() > 0 and q.min() >= .1
        face_points = points[surface[:, :3]]
        caps = np.any(np.all(np.isclose(face_points, -.5, atol=1e-12, rtol=0), axis=1)
            | np.all(np.isclose(face_points, .5, atol=1e-12, rtol=0), axis=1), axis=1)
        wall_nodes = np.unique(surface[~caps, :3])
        on_cap = np.any(np.isclose(abs(points[wall_nodes]), .5, atol=1e-12, rtol=0), axis=1)
        curve_nodes = wall_nodes[on_cap]
        curve_cells = np.isin(cells[:, :4], curve_nodes).any(axis=1)
        result[label] = dict(source=str(source), nodes=row['nodes'], tets=row['elements'],
            geometry=row['geometry'], corner_edge_length_quantiles=np.quantile(lengths, fractions).tolist(),
            cell_rms_edge_length_quantiles=np.quantile(rms, fractions).tolist(),
            corner_volume_quantiles=np.quantile(volume, fractions).tolist(),
            corner_quality_quantiles=np.quantile(q, fractions).tolist(),
            volume_weighted_rms_edge=float(volume @ rms / volume.sum()),
            curve_corner_nodes=len(curve_nodes), curve_incident_cells=int(curve_cells.sum()),
            curve_incident_corner_volume_fraction=float(volume[curve_cells].sum()/volume.sum()))
        print(f'{label}: {row["nodes"]} nodes, {row["elements"]} tets, median RMS edge {np.median(rms):.6g}', flush=True)
    native, actual = result['native'], result['actual']
    result.update(quantile_fractions=fractions,
        actual_to_native_element_ratio=actual['tets']/native['tets'],
        target_tets=actual['tets'], target_nodes=actual['nodes'],
        equal_regular_tet_edge_at_actual_count=(6*math.sqrt(2)*native['geometry']['volume']/actual['tets'])**(1/3),
        count_only_uniform_length_scale=(native['tets']/actual['tets'])**(1/3),
        inputs=hashes, mesh_changed=False, stiffness_scores_used_for_size_selection=False,
        new_solve=False, comparative_timing_claim=False,
        scope='Read-only density and count scales on admitted geometry; no selected new mesh or accuracy claim')
    w.save(args.run/'analysis-inputs.json', hashes)
    w.save(args.run/'result.json', result)
    assert all(w.digest(Path(p)) == value for p, value in hashes.items())


if __name__ == '__main__':
    main()
