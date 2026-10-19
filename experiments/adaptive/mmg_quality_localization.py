"""Locate failed MMG tetrahedra relative to fixed caps and wall triangles."""

import argparse
import itertools
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def quality(points, cells):
    import numpy as np

    coordinates = points[cells]
    volume = abs(np.linalg.det(coordinates[:, 1:] - coordinates[:, :1])) / 6
    squared = sum(((coordinates[:, first] - coordinates[:, second])**2).sum(axis=1)
        for first, second in itertools.combinations(range(4), 2))
    return (12 * (3 * volume)**(2 / 3) / squared)**1.5


def worker(source, output):
    import numpy as np
    from cap_quality_bound import bounds, self_check

    self_check()
    equilateral = np.array([[0, 0, 0], [1, 0, 0], [.5, np.sqrt(3) / 2, 0],
        [.5, np.sqrt(3) / 6, np.sqrt(2 / 3)]])
    np.testing.assert_allclose(quality(equilateral, np.arange(4)[None]), 1, atol=1e-14, rtol=0)
    row = json.loads((source / "result.json").read_text())
    mesh = np.load(source / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 4 and surface.shape[1] == 3 and len(cells) <= workflow.CAP
    measured = quality(points, cells)
    assert abs(float(measured.min()) - row["geometry"]["q_min_corner"]) < 1e-12
    triangles = points[surface]
    cap = np.zeros(len(surface), dtype=bool)
    for axis in range(3):
        for side in (-.5, .5):
            cap |= np.all(abs(triangles[:, :, axis] - side) < 1e-10, axis=1)
    ceiling, _ = bounds(triangles)
    faces = {tuple(sorted(face)): index for index, face in enumerate(surface)}
    cap_vertices = np.unique(surface[cap])
    wall_vertices = np.unique(surface[~cap])
    bad = np.flatnonzero(measured < .1)
    counts = {"bad_with_cap_face": 0, "bad_with_wall_face": 0,
        "bad_with_obstructing_wall_face": 0, "bad_with_obstructing_cap_face": 0}
    examples = []
    for cell in sorted(bad, key=lambda index: measured[index]):
        incident = [faces[key] for face in itertools.combinations(cells[cell], 3)
            if (key := tuple(sorted(face))) in faces]
        cap_faces = [index for index in incident if cap[index]]
        wall_faces = [index for index in incident if not cap[index]]
        counts["bad_with_cap_face"] += bool(cap_faces)
        counts["bad_with_wall_face"] += bool(wall_faces)
        counts["bad_with_obstructing_cap_face"] += any(ceiling[index] < .1 for index in cap_faces)
        counts["bad_with_obstructing_wall_face"] += any(ceiling[index] < .1 for index in wall_faces)
        if len(examples) < 20:
            examples.append({"element": int(cell), "quality": float(measured[cell]),
                "coordinates": points[cells[cell]].tolist(),
                "cap_face_quality_ceilings": [float(ceiling[index]) for index in cap_faces],
                "wall_face_quality_ceilings": [float(ceiling[index]) for index in wall_faces],
                "nearest_box_plane_distance": float((.5 - abs(points[cells[cell]])).min())})
    workflow.save(output, {"source": str(source), "nodes": len(points), "elements": len(cells),
        "bad_elements": len(bad), "quality_minimum": float(measured.min()),
        "bad_incident_on_cap_vertices": int(np.isin(cells[bad], cap_vertices).any(axis=1).sum()),
        "bad_incident_on_wall_vertices": int(np.isin(cells[bad], wall_vertices).any(axis=1).sum()),
        "cap_faces_below_quality_ceiling": int((cap & (ceiling < .1)).sum()),
        "wall_faces_below_quality_ceiling": int((~cap & (ceiling < .1)).sum()),
        **counts, "worst_examples": examples,
        "scope": "Saved actual MMG output; fixed-face obstructions and adjacency only, no guarantee that grading or interior optimization repairs it"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--sources", nargs="+", type=Path)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.run / (args.worker.name + ".json"))
        return
    assert args.sources
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "cap_quality_bound.py"]
    paths.extend(source / name for source in args.sources for name in ("result.json", "mesh.npz"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 120, "threads": 1})
    with Job(parts=[source.name for source in args.sources]) as job:
        for source in args.sources:
            with job.part(source.name) as part:
                record = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                    "--worker", str(source)], args.run / (source.name + ".log"), part, timeout=120)
                workflow.save(args.run / (source.name + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {source.name}.log"
                result = json.loads((args.run / (source.name + ".json")).read_text())
                job.result(source.name, seconds=record["wall_s"], bad_elements=result["bad_elements"],
                    obstructing_wall_faces=result["wall_faces_below_quality_ceiling"])
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
