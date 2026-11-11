"""Capture rejected spinodoid stages without weakening the production quality gate.

Run inside macbook_stage.py with an explicit timeout and 8 GiB watchdog.
The existing isolated observer package must reproduce the baseline rejection.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import sys
from unittest.mock import patch

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--package", type=Path, required=True)
    p.add_argument("--selection", type=Path, required=True)
    p.add_argument("--baseline", type=Path, required=True)
    a = p.parse_args()
    a.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(a.package.resolve()))
    np, _, meshers, *_ = w.imports()
    import comparison_benchmark as benchmark
    import scorecard
    import comparison_fixtures
    from cad_bcc_gap_localization import finish_samples, distances, KINDS
    assert meshers._meshers._LINEAR_OBSERVATIONS
    assert Path(meshers.__file__).resolve().is_relative_to(a.package.resolve())
    receipt = json.loads((a.package.parent/"receipt.json").read_text())
    assert w.digest(Path(meshers._meshers.__file__)) == receipt["native_binary_sha256"]
    for name, expected in receipt["inputs"].items():
        if name.endswith((".rs", "Cargo.toml", "Cargo.lock")):
            assert w.digest(w.ROOT/name) == expected, name
    baseline = json.loads(a.baseline.read_text())
    entries = json.loads(a.selection.read_text())["cases"]
    result = dict(complete=False, native_binary_sha256=receipt["native_binary_sha256"],
        baseline_sha256=w.digest(a.baseline), selection_sha256=w.digest(a.selection),
        snapshots_unadmitted=True, diagnostic_timings_only=True, cases=[])
    w.save(a.run/"result.json", result)
    for entry in entries:
        name = entry["case"]
        directory = a.run/name
        directory.mkdir()
        for source in (Path(entry["fixtures"])/name).iterdir():
            if source.is_file():
                shutil.copyfile(source, directory/source.name)
        observations = []
        original = meshers._meshers.generate

        def capture(*args, **kwargs):
            return original(*args, **kwargs, _observations=observations)

        with patch.object(meshers._meshers, "generate", capture):
            benchmark.worker(directory, name, entry["native"], Path(entry["fixtures"]))
        row = json.loads((directory/"result.json").read_text())
        old = next(r for r in baseline["rows"] if r["case"] == name and r["arm"] == "native")
        def rejection(r):
            match = re.search(r"(minimum MMG quality|sampled geometry error) ([^ ]+)", r["error"])
            return match[1], float(match[2])
        current, previous = rejection(row), rejection(old)
        assert current[0] == previous[0] and abs(current[1]-previous[1]) < 1e-11
        if current[0] == "minimum MMG quality":
            assert re.search(r"\(\d+ elements below", row["error"])[0] == re.search(r"\(\d+ elements below", old["error"])[0]
        compiled = meshers.compile_field(comparison_fixtures.case(name, Path(entry["fixtures"])).function)
        stages = []
        for i, state in enumerate(observations):
            path = directory/f"{i:02d}-{state['stage']}.npz"
            arrays = {k:state[k] for k in ("points", "tetrahedra", "surface")}
            np.savez_compressed(path, **arrays)
            points, cells, surface = [arrays[k] for k in ("points", "tetrahedra", "surface")]
            quality, _, volume = scorecard.tet_metrics(points, cells)
            bad = np.flatnonzero(quality < .1)
            locks = np.isclose(abs(points), .5, atol=1e-12, rtol=0)
            wall_nodes = np.unique(surface)
            bad_nodes = cells[bad].ravel()
            item = dict(stage=state["stage"], path=str(path.resolve()), sha256=w.digest(path),
                nodes=len(points), tets=len(cells), minimum_quality=float(quality.min()),
                below_floor=len(bad), nonpositive=int((volume <= 0).sum()),
                bad_tet_ids=bad.tolist(), bad_tet_quality=quality[bad].tolist(),
                bad_corner_on_surface=int(np.isin(bad_nodes, wall_nodes).sum()),
                bad_corner_locks=np.bincount(locks[bad_nodes].sum(axis=1), minlength=4).tolist())
            if entry["native"]["adapt"]:
                queries, locked, face_ids, kinds = finish_samples(np, points, surface, [[-.5]*3, [.5]*3])
                errors, *_ = distances(np, compiled, queries, locked, [[-.5]*3, [.5]*3], [True]*3, 1.)
                worst = int(np.argmax(errors))
                item["finish_error"] = dict(maximum=float(errors[worst]), kind=KINDS[kinds[worst]],
                    point=queries[worst].tolist(), locks=locked[worst].tolist(), face=int(face_ids[worst]),
                    maxima_by_kind={kind:float(errors[kinds == i].max()) for i,kind in enumerate(KINDS)},
                    violations=int((errors > .002).sum()))
            stages.append(item)
        expected = (["clipped", "adapted", "optimized", "geometry-closed", "final-linear"] if entry["native"]["adapt"] else
            ["clipped", "optimized", "geometry-refined", "final-linear"])
        assert [s["stage"] for s in stages] == expected
        if current[0] == "sampled geometry error":
            assert abs(stages[-1]["finish_error"]["maximum"]-current[1]) < 1e-11
        result["cases"].append(dict(case=name, rejection_reproduced=True, stages=stages,
            original_error=old["error"], observed_error=row["error"]))
        w.save(a.run/"result.json", result)
        print(name+": "+", ".join(f"{s['stage']} q={s['minimum_quality']:.5g} bad={s['below_floor']}" for s in stages), flush=True)
    result["complete"] = True
    w.save(a.run/"result.json", result)


if __name__ == "__main__":
    main()
