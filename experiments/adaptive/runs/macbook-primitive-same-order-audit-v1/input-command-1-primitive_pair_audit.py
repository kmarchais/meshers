"""Independent replay of repeated, admitted same-order basic-solid comparisons."""

import argparse
import json
from pathlib import Path
import statistics

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    np, *_ = w.imports()
    from comparison_mesh_audit import observe
    import primitive_family_fixtures as fixtures

    manifest = json.loads((a.source / "manifest.json").read_text())
    summary = json.loads((a.source / "summary.json").read_text())
    assert summary["complete"] and len(summary["rows"]) == 36
    assert manifest["repetitions"] == 5 and manifest["warmups"] == 1
    inputs = {}
    for name, digest in {**manifest["inputs"], **manifest["fixture_hashes"]}.items():
        path = Path(name)
        if not path.exists() or w.digest(path) != digest:
            matches = [p for p in a.source.glob("input-*-" + path.name) if w.digest(p) == digest]
            assert len(matches) == 1, name
            path = matches[0]
        assert w.digest(path) == digest
        inputs[str(path)] = digest
    for f in (Path(__file__), a.source / "summary.json", a.source / "manifest.json",
              w.HERE / "comparison_mesh_audit.py", w.HERE / "quadratic_edge_removal_probe.py"):
        inputs[str(f)] = w.digest(f)
    grouped = {}
    replayed = {}
    for row in summary["rows"]:
        assert row["qualified"] and row["worker"]["returncode"] == 0
        assert row["config"]["order"] == 2 and max(row["nodes"], row["elements"]) <= w.CAP
        path = Path(row["directory"]) / "mesh.npz"
        assert w.digest(path) == row["mesh_sha256"]
        inputs[str(path)] = row["mesh_sha256"]
        key = (row["case"], row["mesh_sha256"])
        if key not in replayed:
            with np.load(path) as mesh:
                points, cells, surface = (mesh[n] for n in ("points", "tetrahedra", "surface"))
                assert cells.shape == (row["elements"], 10) and len(points) == row["nodes"]
                geometry, periodic, quality = observe(
                    path.parent, row["case"], fixtures.case(row["case"], a.fixtures),
                    points, cells, surface, False)
            assert geometry["geometry_admitted"]
            for name in ("sampled_two_sided_max", "volume", "q_min_corner", "minimum_sampled_jacobian"):
                np.testing.assert_allclose(geometry[name], row["geometry"][name], atol=1e-12, rtol=1e-12)
            np.testing.assert_allclose(list(quality.values()), list(row["quality_quantiles"].values()),
                                       atol=1e-12, rtol=1e-12)
            periodic["reason"] = "Complete bounded basic solid; no periodic identification"
            replayed[key] = dict(geometry=geometry, quality=quality, periodic=periodic)
        grouped.setdefault(row["case"], {}).setdefault(row["config"]["arm"], []).append(row)
    results = []
    for name, groups in grouped.items():
        arms = {}
        for arm, rows in groups.items():
            measured = [r for r in rows if not r["warmup"]]
            assert len(rows) == 6 and len(measured) == 5
            assert {r["repeat"] for r in measured} == set(range(5))
            assert len({r["mesh_sha256"] for r in rows}) == 1
            times = [r["generation_s"] for r in measured]
            first = rows[0]
            arms[arm] = dict(
                median_s=statistics.median(times), range_s=[min(times), max(times)], samples=times,
                elements=first["elements"], nodes=first["nodes"], config=first["config"],
                mesh_sha256=first["mesh_sha256"],
                **replayed[(name, first["mesh_sha256"])])
        native, actual = arms["meshers"], arms["microgen"]
        results.append(dict(
            case=name, arms=arms, measured_pairs=5,
            speed_actual_over_native=actual["median_s"] / native["median_s"],
            tets_actual_over_native=actual["elements"] / native["elements"],
            nodes_actual_over_native=actual["nodes"] / native["nodes"],
            distance_actual_over_native=actual["geometry"]["sampled_two_sided_max"] / native["geometry"]["sampled_two_sided_max"],
            quality_native_minus_actual={k: native["quality"][k] - actual["quality"][k] for k in native["quality"]},
        ))
    assert len(results) == 3 and len(replayed) == 6
    assert all(w.digest(Path(f)) == h for f, h in inputs.items())
    w.save(a.run / "result.json", dict(
        rows=results, inputs=inputs, checks_passed=True, unique_mesh_replays=len(replayed),
        attempts_verified=36, FEA_performed=False, same_host=True, same_order=2,
        limitations="Three frozen basic solids. Actual resolution24 first-pass overrides, second-pass omission, explicit repair and straight elevation are charged. Different counts and curvature. No general optimum, FEA convergence or stress claim."))


if __name__ == "__main__":
    main()
