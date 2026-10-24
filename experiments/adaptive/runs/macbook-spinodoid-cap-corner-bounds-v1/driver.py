"""Necessary cap quality bounds at fixed convex polygon corners; no edits."""

import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow
from job_parts import worker_part

import numpy as np

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def wedge_bound(previous, center, following):
    first, second = previous - center, following - center
    product = float(np.linalg.norm(first) * np.linalg.norm(second))
    assert product > 0
    dot = float(first @ second)
    cross = float(first[0] * second[1] - first[1] * second[0])
    convex = cross < 0  # The material is on the left of the oriented polygon.
    angle = float(np.degrees(np.arctan2(abs(cross), dot)))
    # For angle a <= 60 degrees, max triangle quality = sqrt(3)*sin(a)/(2-cos(a)).
    # Equal adjacent lengths attain the maximum, by a^2+b^2 >= 2ab.
    ceiling = np.sqrt(3) * abs(cross) / (2 * product - dot) if convex and dot > .5 * product else 1.
    return float(ceiling), convex, angle if convex else 360 - angle


def decimal_wedge_bound(previous, center, following):
    with localcontext() as context:
        context.prec = 80
        first = [Decimal.from_float(float(a)) - Decimal.from_float(float(b))
            for a, b in zip(previous, center, strict=True)]
        second = [Decimal.from_float(float(a)) - Decimal.from_float(float(b))
            for a, b in zip(following, center, strict=True)]
        dot = sum(a * b for a, b in zip(first, second, strict=True))
        cross = first[0] * second[1] - first[1] * second[0]
        product = (sum(a * a for a in first) * sum(a * a for a in second)).sqrt()
        if cross >= 0 or dot <= product / 2:
            return Decimal(1)
        return Decimal(3).sqrt() * abs(cross) / (2 * product - dot)


def polygon_corners(points, triangles):
    triangles = triangles.copy()
    edges = points[triangles[:, 1:]] - points[triangles[:, :1]]
    signed = edges[:, 0, 0] * edges[:, 1, 1] - edges[:, 0, 1] * edges[:, 1, 0]
    assert (signed != 0).all(), "Degenerate cap triangle"
    negative = signed < 0
    triangles[negative, 1], triangles[negative, 2] = triangles[negative, 2].copy(), triangles[negative, 1].copy()
    incidence = defaultdict(list)
    for face in triangles:
        for first, second in ((0, 1), (1, 2), (2, 0)):
            edge = tuple(map(int, face[[first, second]]))
            incidence[tuple(sorted(edge))].append(edge)
    incoming, outgoing = defaultdict(list), defaultdict(list)
    for copies in incidence.values():
        assert len(copies) <= 2
        if len(copies) == 2:
            assert copies[0] == copies[1][::-1], "Incoherent planar cap incidence"
            continue
        first, second = copies[0]
        outgoing[first].append(second)
        incoming[second].append(first)
    rows, ambiguous = [], []
    for node in sorted(set(incoming) | set(outgoing)):
        if len(incoming[node]) != 1 or len(outgoing[node]) != 1:
            ambiguous.append(node)
            continue
        previous, following = incoming[node][0], outgoing[node][0]
        ceiling, convex, angle = wedge_bound(points[previous], points[node], points[following])
        rows.append({"node": node, "previous": previous, "following": following,
            "convex": convex, "interior_angle_deg": angle, "maximum_incident_cap_quality": ceiling})
    return rows, ambiguous


def worker(args):
    workflow.imports()
    from facet_incidence_audit import audit
    from periodic import periodic_pairs

    mesh = np.load(args.source / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    consistency = audit(points, cells, surface)
    assert consistency["oriented_facet_consistency_passed"] and not consistency["nonpositive_corner_determinants"]
    _, _, periodic = periodic_pairs(points, surface)
    planes, witnesses = [], []
    for axis in range(3):
        transverse = [other for other in range(3) if other != axis]
        for side in (-.5, .5):
            selected = np.all(workflow.on_box_plane(points[surface[:, :3], axis], side), axis=1)
            if not selected.any():
                continue
            planar = points[:, transverse]
            corners, ambiguous = polygon_corners(planar, surface[selected, :3])
            blocked = [row for row in corners if row["maximum_incident_cap_quality"] < .1 - 1e-12]
            assert len(witnesses) + len(blocked) <= 512, "Bound extended-precision witness count"
            for row in blocked:
                nodes = [row[key] for key in ("previous", "node", "following")]
                extended = decimal_wedge_bound(*planar[nodes])
                assert extended < Decimal(".1")
                assert abs(float(extended) - row["maximum_incident_cap_quality"]) < 1e-12
                witnesses.append(dict(row, axis=axis, side=side,
                    coordinates=points[nodes].tolist(), decimal80_cap_quality_bound=str(extended)))
            planes.append({"axis": axis, "side": side, "boundary_corners": len(corners),
                "ambiguous_boundary_nodes": ambiguous, "obstructing_corners": len(blocked)})
    case = json.loads((args.source / "result.json").read_text())["case"]
    workflow.save(args.run / "result.json", {"source": str(args.source), "case": case,
        "planes": planes, "witnesses": witnesses, "obstructing_corner_count": len(witnesses),
        "fixed_polygon_corner_obstruction": bool(witnesses), "periodic_source": periodic,
        "qualified": False, "source_mesh_unchanged": True, "mmg_calls": [], "stiffness_solve_performed": False,
        "scope": "Necessary quality bound for triangles inside each fixed convex polygon wedge, combined with the retained-face tet bound; 80-digit replay of floating-point coordinates, not interval certification or impossibility after boundary changes"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.source, args.run = args.source.resolve(), args.run.resolve()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), args.source / "mesh.npz", args.source / "result.json"]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "job_parts.py",
        "test_cap_corner_probe.py", "test_job_parts.py", "facet_incidence_audit.py", "cap_quality_bound.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY, "platform": sys.platform,
        "worker_timeout_s": 120, "maximum_corner_witnesses": 512, "decimal_precision": 80,
        "quality_target": .1, "no_mesh_changes_mmg_or_fea": True})
    with Job(parts=["controls", "corner-bounds"]) as job:
        for stage in ("controls", "corner-bounds"):
            with worker_part(job, stage) as part:
                command = ([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_cap_corner_probe.py",
                    "experiments/adaptive/test_job_parts.py"] if stage == "controls" else
                    [sys.executable, __file__, "--source", str(args.source), "--run", str(args.run), "--worker"])
                record = workflow.bounded(command, args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage == "corner-bounds":
                    values["obstructing_corners"] = json.loads((args.run / "result.json").read_text())["obstructing_corner_count"]
                job.result(stage, **values)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
