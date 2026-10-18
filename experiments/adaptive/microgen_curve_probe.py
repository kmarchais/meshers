"""Locate sampled boundary failures on required VTK cap vertices and edges."""

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


def on_required_cap(support, cap_nodes, cap_edges):
    return ((len(support) == 1 and int(support[0]) in cap_nodes)
        or (len(support) == 2 and tuple(sorted(map(int, support))) in cap_edges))


def self_check():
    assert on_required_cap([1], {1, 2}, {(1, 2)})
    assert on_required_cap([2, 1], {1, 2}, {(1, 2)})
    assert not on_required_cap([2, 3], {1, 2}, {(1, 2)})
    assert not on_required_cap([1, 2, 3], {1, 2}, {(1, 2)})


def worker(source, run):
    np, *_, scorecard = workflow.imports()
    from cap_quality_bound import caps

    case = json.loads((source / "result.json").read_text())["case"]
    seed = np.load(source / "vtk-seed.npz")
    geometry = workflow.geometry(source, case, seed["points"], seed["tetrahedra"], seed["surface"], detail_count=32)
    cap_faces, _ = caps(seed)
    cap_nodes = set(map(int, cap_faces.ravel()))
    cap_edges = {tuple(sorted((int(face[first]), int(face[second]))))
        for face in cap_faces for first, second in ((0, 1), (1, 2), (2, 0))}
    probes = []
    for probe in geometry["worst_boundary_samples"]:
        face = seed["surface"][probe["face"]]
        coordinates = seed["points"][face]
        np.testing.assert_array_equal(coordinates, probe["coordinates"])
        barycentric = scorecard.BARY[probe["probe_index"]]
        np.testing.assert_allclose(barycentric @ coordinates, probe["probe"], atol=1e-14, rtol=0)
        support = face[np.flatnonzero(abs(barycentric) > 1e-14)]
        fixed = on_required_cap(support, cap_nodes, cap_edges)
        probes.append({**probe, "required_cap_probe": fixed,
            "fixed_probe_exceeds_sampled_gate": fixed and not probe["failed"] and probe["distance"] > 0.002})
    previous = json.loads((source / "vtk-seed.json").read_text())["geometry"]
    for name, value in previous.items():
        assert geometry[name] == value, name
    workflow.save(run / "result.json", {"case": case, "geometry": geometry, "classified_probes": probes,
        "fixed_cap_gate_failures_among_top_probes": sum(probe["fixed_probe_exceeds_sampled_gate"] for probe in probes),
        "scope": "Finite forward projections under the unchanged admission metric; no exact closest-distance bound or exhaustive cap-curve audit"})


def main():
    self_check()
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.source, args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    case = json.loads((args.source / "result.json").read_text())["case"]
    inputs = {str(args.source / name): workflow.digest(args.source / name) for name in
        ("result.json", "vtk-seed.json", "vtk-seed.npz", case + "-samples.npz", case + "-specimen.json")}
    sources = {str(path): workflow.digest(path) for path in (Path(__file__), Path(workflow.__file__))}
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "source_sha256": sources,
        "memory_bytes": workflow.MEMORY, "threads": 1, "worker_timeout_s": 180})
    with Job(parts=["fixed-cap-probes"]) as job:
        with job.part("fixed-cap-probes") as part:
            record = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                "--source", str(args.source), "--worker"], args.run / "worker.log", part, timeout=180)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("fixed-cap-probes", seconds=record["wall_s"],
                sampled_fixed_failures=row["fixed_cap_gate_failures_among_top_probes"])
    assert all(workflow.digest(Path(path)) == value for path, value in {**inputs, **sources}.items())


if __name__ == "__main__":
    main()
