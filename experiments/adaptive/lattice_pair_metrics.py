"""Compare admitted geometry and corner-quality distributions, without timing claims."""

import argparse
import json
from pathlib import Path
import numpy as np
import microgen_workflow as w
from mmg_quality_localization import quality


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--actual', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    paths = [Path(__file__), w.HERE/'mmg_quality_localization.py']
    results = []
    for run in (args.source, args.actual):
        row = json.loads((run/'result.json').read_text())
        assert row['case'] == 'smooth-lattice' and row['qualified'] and row['geometry']['geometry_admitted']
        mesh = np.load(run/'mesh.npz')
        assert len(mesh['points']) == row['nodes'] <= w.CAP
        assert len(mesh['tetrahedra']) == row['elements'] <= w.CAP
        q = quality(mesh['points'], mesh['tetrahedra'][:, :4])
        assert abs(q.min()-row['geometry']['q_min_corner']) < 1e-12
        results.append(dict(source=str(run), nodes=row['nodes'], elements=row['elements'],
            quality_quantiles=dict(zip(('min','p01','p10','median','p90','p99','max'),
                map(float,np.quantile(q,[0,.01,.1,.5,.9,.99,1])))),
            distance=row['geometry']['sampled_two_sided_max'],
            absolute_relative_volume_error=abs(row['geometry']['relative_volume_error'])))
        paths.extend(run/name for name in ('mesh.npz', 'result.json',
            'smooth-lattice-samples.npz', 'smooth-lattice-specimen.json'))
    for name in ('smooth-lattice-samples.npz', 'smooth-lattice-specimen.json'):
        assert w.digest(args.source/name) == w.digest(args.actual/name)
    w.save(args.output, dict(routes=results, identical_specimen_and_samples=True,
        all_quality_quantiles_no_worse=all(results[0]['quality_quantiles'][k] >= v
            for k,v in results[1]['quality_quantiles'].items()),
        comparative_timing_claim=False, objective_complete=False,
        inputs={str(p):w.digest(p) for p in paths}))
    for row in results:
        print(f'{row["source"]}: qmin {row["quality_quantiles"]["min"]:.6f}, '
            f'qmedian {row["quality_quantiles"]["median"]:.6f}, distance {row["distance"]:.7f}', flush=True)


if __name__ == '__main__':
    main()
