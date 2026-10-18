"""Measure the quality cost of exact roots on rejected saved wall vertices."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, *_, scorecard = workflow.imports()
    from periodic import periodic_pairs

    source = np.load(args.source / "mesh.npz")
    points, tets, surface = (source[key] for key in ("points", "tetrahedra", "surface"))
    assert tets.shape[1] == 4 and len(tets) <= workflow.CAP
    _, classes, _ = periodic_pairs(points, surface)
    probes = json.loads((args.probe / "result.json").read_text())["classified_probes"]
    rows, seen = [], set()
    for probe in probes:
        barycentric = scorecard.BARY[probe["probe_index"]]
        support = np.flatnonzero(abs(barycentric) > 1e-14)
        if len(support) != 1 or probe["failed"] or probe["distance"] <= .002:
            continue
        node = probe["node_ids"][int(support[0])]
        group = int(classes[node])
        if group in seen:
            continue
        seen.add(group)
        np.testing.assert_allclose(points[node], probe["probe"], atol=1e-14, rtol=0)
        members = np.flatnonzero(classes == group)
        incident = tets[np.isin(tets, members).any(axis=1)]
        moved = points.copy()
        moved[members] += np.array(probe["foot"]) - points[node]
        for axis, locked in enumerate(probe["locked"]):
            if locked:
                np.testing.assert_array_equal(moved[members, axis], points[members, axis])
        before = scorecard.tet_metrics(points, incident)[0]
        after = scorecard.tet_metrics(moved, incident)[0]
        determinants = np.linalg.det(moved[incident[:, 1:]] - moved[incident[:, :1]])
        rows.append({"node": node, "periodic_members": members.tolist(),
            "origin": points[node].tolist(), "target": probe["foot"],
            "locked": probe["locked"], "distance": probe["distance"],
            "incident_elements": len(incident), "minimum_quality_before": float(before.min()),
            "minimum_quality_after": float(after.min()),
            "inverted_elements_after": int((determinants <= 0).sum()),
            "elements_below_quality_after": int((after < .1).sum()),
            "direct_projection_meets_quality_gate": bool(after.min() >= .1 and determinants.min() > 0)})
    assert rows
    workflow.save(args.run / "result.json", {"source": str(args.source), "vertices": rows,
        "scope": "Separate single-periodic-class root moves, not a repaired or admitted mesh",
        "directly_projectable_classes": sum(row["direct_projection_meets_quality_gate"] for row in rows)})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [args.source / "mesh.npz", args.source / "result.json", args.probe / "result.json",
        Path(__file__), Path(workflow.__file__), workflow.HERE / "scorecard.py"]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes,
        "memory_bytes": workflow.MEMORY, "threads": 1, "timeout_s": 120})
    with Job(parts=["projection-stars"]) as job:
        with job.part("projection-stars") as part:
            command = [sys.executable, __file__, "--run", str(args.run), "--source",
                str(args.source), "--probe", str(args.probe), "--worker"]
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=120)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("projection-stars", seconds=record["wall_s"],
                directly_projectable_classes=row["directly_projectable_classes"])
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
