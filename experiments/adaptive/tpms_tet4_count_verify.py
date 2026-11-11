"""Verify returned Tet4 artifacts, frozen specimens and selected repeated pools."""
import argparse
import hashlib
import json
from pathlib import Path

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selection", type=Path, required=True)
    p.add_argument("--studies", type=Path, nargs="+", required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--arms", default="actual,native")
    a = p.parse_args()
    np, *_ = w.imports()
    import scorecard
    selection = json.loads(a.selection.read_text())
    entries = {e["case"]: e for e in selection["cases"]}
    assert entries and len(entries) == len(selection["cases"])
    report = json.loads(a.report.read_text())
    assert {e["case"] for e in report["cases"]} == set(entries)
    assert report["all_orders"] == 1 and report["admission"] == dict(
        distance=.002, absolute_relative_volume=.01, minimum_quality=.1, max_nodes=w.CAP, max_tets=w.CAP)
    artifacts, uniques, fallback_sources = {}, {}, []
    attempts = 0
    measured = {}
    for study in a.studies:
        manifest = json.loads((study/"manifest.json").read_text())
        summary = json.loads((study/"summary.json").read_text())
        assert summary["complete"]
        for source, expected in manifest["inputs"].items():
            path = Path(source)
            if path.exists() and w.digest(path) == expected:
                continue
            archived = study/"driver.py"
            assert path.name == "tpms_tet4_count_benchmark.py" and w.digest(archived) == expected, source
            fallback_sources.append(dict(original=source, archived=str(archived), sha256=expected))
        repeated = manifest["repetitions"] == 5
        for row in summary["rows"]:
            attempts += 1
            directory = Path(row["directory"])
            original = json.loads((directory/"result.json").read_text())
            for key, value in original.items():
                if key in ("status", "qualified") and row["resources"]["returncode"] != 0:
                    continue
                assert row[key] == value, (directory, key)
            entry = entries[row["case"]]
            for suffix in ("-samples.npz", "-specimen.json"):
                path = directory/(row["case"]+suffix)
                expected = Path(entry["fixtures"])/row["case"]/path.name
                assert w.digest(path) == w.digest(expected)
            assert row["config"]["order"] == 1
            if repeated and not row["warmup"]:
                arm = row["arm"]
                assert row["config"] == entry[arm]
                measured.setdefault((row["case"], arm), []).append(row)
            if not row["returned_mesh_verified"]:
                continue
            path = directory/"mesh.npz"
            assert w.digest(path) == row["mesh_sha256"]
            with np.load(path) as saved:
                arrays = {key:saved[key] for key in ("points", "tetrahedra", "surface")}
            points, cells, surface = [arrays[k] for k in ("points", "tetrahedra", "surface")]
            assert points.shape == (row["nodes"], 3) and cells.shape == (row["elements"], 4)
            assert surface.ndim == 2 and surface.shape[1] == 3
            assert max(len(points), len(cells)) <= w.CAP
            assert np.isfinite(points).all() and cells.min() >= 0 and cells.max() < len(points)
            assert surface.min() >= 0 and surface.max() < len(points)
            h = hashlib.sha256()
            for key, array in arrays.items():
                h.update(key.encode())
                h.update(str((array.shape, array.dtype.str)).encode())
                h.update(array.tobytes())
            array_digest = h.hexdigest()
            if array_digest not in uniques:
                q, _, volumes = scorecard.tet_metrics(points, cells)
                nonpositive = int((volumes <= 0).sum())
                if nonpositive:
                    assert not row["qualified"] and "Jacobian" in row["geometry"]["failed_gates"]
                np.testing.assert_allclose(np.quantile(q, [0, .01, .1, .5, .9, .99, 1]),
                    [row["quality_quantiles"][k] for k in ("min", "p01", "p10", "median", "p90", "p99", "max")],
                    atol=2e-14, rtol=0)
                assert abs(float(volumes.sum())-row["geometry"]["volume"]) < 2e-12
                uniques[array_digest] = dict(case=row["case"], arm=row["arm"], nodes=len(points), tets=len(cells),
                    nonpositive_tets=nonpositive)
            artifacts[str(path)] = dict(sha256=row["mesh_sha256"], array_sha256=array_digest)
    repeated_cases = [e for e in entries.values() if e["repeat"]]
    repeated_array_sets = {}
    for entry in repeated_cases:
        for arm in a.arms.split(","):
            pool = measured[(entry["case"], arm)]
            assert len(pool) == 5 and all(r["returned_mesh_verified"] for r in pool)
            repeated_array_sets[entry["case"]+"/"+arm] = len({
                artifacts[str(Path(r["directory"])/"mesh.npz")]["array_sha256"] for r in pool})
    w.save(a.output, dict(passed=True, attempts=attempts, returned_artifacts=len(artifacts),
        distinct_array_sets=len(uniques), repeated_cases=len(repeated_cases), measured_pairs=5*len(repeated_cases),
        repeated_array_set_counts=repeated_array_sets,
        verified_arms=a.arms.split(","),
        source_snapshots_resolved=fallback_sources, artifacts=artifacts, distinct_meshes=uniques,
        report_sha256=w.digest(a.report), selection_sha256=w.digest(a.selection),
        scope="Actual Tet4 arrays, hashes, frozen input copies, positive determinants, counts, volume and quality replay. Existing complete periodic/geometry audits remain in the per-worker evidence. No new FEA or admitted-speed claim."))
    print(f"Verified {attempts} attempts and {len(artifacts)} returned Tet4 artifacts")


if __name__ == "__main__":
    main()
