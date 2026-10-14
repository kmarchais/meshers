"""Bounded geometry-only surface coarsening trial from immutable inputs."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report
from run import digest, write_json

OUTPUT = Path(__file__).resolve().parents[1]


def owned_process_rss(process):
    """Count the worker and its recursive children, including Windows launchers."""
    import psutil
    total = 0
    try:
        family = [process, *process.children(recursive=True)]
    except psutil.NoSuchProcess:
        return 0
    for member in family:
        try:
            total += member.memory_info().rss
        except psutil.NoSuchProcess:
            pass
    return total


def surface_signature(surface):
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    vertices = np.unique(surface)
    edges = np.unique(np.sort(np.concatenate([surface[:, [0, 1]], surface[:, [1, 2]],
                                             surface[:, [2, 0]]]), axis=1), axis=0)
    ids = np.searchsorted(vertices, edges)
    graph = coo_matrix((np.ones(len(ids)), (ids[:, 0], ids[:, 1])), shape=(len(vertices), len(vertices)))
    components, _ = connected_components(graph, directed=False)
    return {"components": int(components), "euler_characteristic": len(vertices) - len(edges) + len(surface)}


def caps(points, surface, original_ids):
    import numpy as np
    xyz = points[surface]
    mask = np.any(np.all(np.isclose(xyz, -.5, atol=1e-10, rtol=0), axis=1), axis=1)
    mask |= np.any(np.all(np.isclose(xyz, .5, atol=1e-10, rtol=0), axis=1), axis=1)
    triangles = original_ids[surface[mask]]
    rotation = (triangles.argmin(axis=1)[:, None] + np.arange(3)) % 3
    triangles = np.take_along_axis(triangles, rotation, axis=1)
    return triangles[np.lexsort(triangles.T[::-1])]


def wall_guard(case, tolerance):
    import numpy as np
    weights = np.array([[1/3]*3, [1, 0, 0], [0, 1, 0], [0, 0, 1],
                        [.5, .5, 0], [0, .5, .5], [.5, 0, .5]])
    offsets = np.eye(3) * 1e-5

    def check(triangles):
        if not len(triangles):
            return False
        probes = np.einsum("si,tij->tsj", weights, triangles).reshape(-1, 3)
        values = case.field(*probes.T)
        gradient = np.column_stack([(case.field(*(probes + d).T) - case.field(*(probes - d).T)) / 2e-5
                                    for d in offsets])
        size = np.linalg.norm(gradient, axis=1)
        distance = np.abs(np.abs(values) - 1) / np.maximum(size, 1e-14)
        normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        direction = gradient.reshape(-1, len(weights), 3)[:, 0] * np.sign(values.reshape(-1, len(weights))[:, :1])
        return bool(np.isfinite(distance).all() and size.min() > 1e-12 and distance.max() <= tolerance
                    and np.all(np.sum(normals * direction, axis=1) > 0))
    return check


def worker(manifest):
    import numpy as np
    from fixtures import make_case
    from microgen_benchmark.geometry import canonical_mesh, validate
    from periodic import periodic_pairs
    from surface_coarsen import coarsen
    progress_path = OUTPUT / "worker-progress.json"

    def progress(fraction, message, **values):
        write_json(progress_path, {"fraction": fraction, "message": message, "values": values})

    progress(.01, "Load admitted native68 seed")
    seed = json.loads(Path(manifest["seed_result"]).read_text())
    with np.load(manifest["seed_mesh"]) as archive:
        points, tets, surface = [archive[k].copy() for k in ["points", "tetrahedra", "surface"]]
    assert digest(manifest["seed_mesh"]) == seed["mesh_sha256"]
    assert seed["status"] == "ok" and all(seed[k] for k in ["geometry_eligible", "quality_eligible", "periodic_eligible"])
    case = make_case(seed["case"], manifest["specimen"])
    assert case.scale == 1 and not case.implicit
    values = case.field(*points.T)
    levels = np.where(np.abs(np.abs(values) - 1) <= 1e-8, np.sign(values), 0)
    baseline_signature = surface_signature(surface)
    original_caps = caps(points, surface, np.arange(len(points)))
    config = manifest["config"]
    before = time.monotonic()
    new_points, new_tets, retained, metadata = coarsen(
        points, tets, surface, config["remove_nodes"], config["max_edge_length"],
        minimum_quality=max(seed["q_min"], config["minimum_quality"]), levels=levels,
        geometry_check=wall_guard(case, config["local_geometry_tolerance"]),
        checkpoint=lambda fraction, removed, tested, **_: progress(
            .05 + .60 * fraction, "Contract disjoint short wall edges", removed=removed, tested=tested))
    metadata["adaptation_seconds"] = time.monotonic() - before
    np.testing.assert_array_equal(new_points, points[retained])
    np.savez_compressed(OUTPUT / "raw-mesh.npz", points=new_points, tetrahedra=new_tets,
                        retained_vertex_ids=retained)
    progress(.70, "Validate complete geometry and cell overlaps", removed=metadata["removed_nodes"], nodes=len(new_points))
    new_points, new_tets, flipped, unused = canonical_mesh(new_points, new_tets)
    assert flipped == unused == 0, "Coarsening required orientation or unused-node repair"
    with np.load(manifest["geometry_reference"]) as archive:
        metrics, new_surface = validate(new_points, new_tets, case, config["max_validation_pairs"], archive["boundary_points"])
    np.savez_compressed(OUTPUT / "mesh.npz", points=new_points, tetrahedra=new_tets,
                        surface=new_surface, retained_vertex_ids=retained)
    identical_prior_mesh = None
    if manifest.get("expected_mesh"):
        with np.load(manifest["expected_mesh"]) as prior:
            for name, value in [("points", new_points), ("tetrahedra", new_tets),
                                ("surface", new_surface), ("retained_vertex_ids", retained)]:
                np.testing.assert_array_equal(value, prior[name])
        identical_prior_mesh = True
    assert np.array_equal(caps(new_points, new_surface, retained), original_caps), "Coarsening changed oriented cap triangles"
    signature = surface_signature(new_surface)
    assert signature == baseline_signature, "Coarsening changed boundary components or Euler characteristic"
    progress(.90, "Verify unchanged periodic caps and matching", removed=metadata["removed_nodes"], quality=metrics["q_min"])
    pairs, classes, periodic = periodic_pairs(new_points, new_surface, config["periodic_pair_tolerance"])
    np.savez_compressed(OUTPUT / "periodic.npz", classes=classes,
                        **{f"pairs_{i}": p for i, p in enumerate(pairs)})
    volume = manifest["volume_reference"]
    volume_error = abs(metrics["solid_fraction"] - volume["fraction"]) / volume["fraction"]
    volume_with_uncertainty = volume_error + volume["uncertainty"] / volume["fraction"]
    relative_volume_change = abs(metrics["solid_fraction"] - seed["solid_fraction"]) / seed["solid_fraction"]
    geometry_ok = (metrics["sampled_geometry_distance"] is not None
                   and metrics["sampled_geometry_distance"] <= config["geometry_tolerance"]
                   and metrics["sampled_reverse_distance"] <= config["geometry_tolerance"]
                   and metrics["sampled_cap_overshoot"] <= config["geometry_tolerance"]
                   and metrics["sampled_outward_orientation"] and volume_with_uncertainty <= config["volume_tolerance"])
    quality_ok = metrics["q_min"] >= config["minimum_quality"] and metrics["q_min"] >= seed["q_min"] - 1e-12
    admitted = bool(metadata["removed_nodes"] > 0 and geometry_ok and quality_ok
                    and relative_volume_change <= config["maximum_relative_volume_change"])
    result = {"case": seed["case"], "stage": "geometry-only", "fea_performed": False,
              "admitted_geometry": admitted, **metadata, **metrics, **periodic,
              "geometry_eligible": bool(geometry_ok), "quality_eligible": bool(quality_ok), "periodic_eligible": True,
              "oriented_caps_unchanged": True, "retained_coordinates_exact": True,
              "surface_signature": signature, "relative_volume_error_with_uncertainty": volume_with_uncertainty,
              "same_mesh_as_prior": identical_prior_mesh,
              "relative_volume_change": relative_volume_change, "mesh_sha256": digest(OUTPUT / "mesh.npz"),
              "scope": "Geometry-only Python coarsening of an archive; no FEA, adoption or generation speed claim"}
    write_json(OUTPUT / "geometry-result.json", result)
    progress(1., "Geometry-only coarsening evaluated", removed=metadata["removed_nodes"], quality=metrics["q_min"])
    if not admitted:
        raise SystemExit(2)


def main():
    manifest = json.loads((OUTPUT / "manifest.json").read_text())
    if "--worker" in sys.argv:
        try:
            worker(manifest)
        except Exception as exc:
            write_json(OUTPUT / "geometry-error.json", {"error": f"{type(exc).__name__}: {exc}", "fea_performed": False})
            traceback.print_exc()
            raise
        return
    report(0, "Verify frozen coarsening inputs", eta_s=180)
    for group, root in [("source_sha256", OUTPUT / "source"), ("input_sha256", Path())]:
        for name, expected in manifest[group].items():
            assert digest(root / name) == expected, name
    tests = subprocess.run([sys.executable, "-m", "pytest", "-q", *map(str, sorted((OUTPUT / "source").glob("test_*.py")))],
                           capture_output=True, text=True)
    (OUTPUT / "tests.log").write_text(tests.stdout + tests.stderr)
    if tests.returncode:
        raise RuntimeError("Coarsening regressions failed")
    import psutil
    started, peak, stopped = time.monotonic(), 0., None
    with (OUTPUT / "worker.log").open("w") as log:
        process = subprocess.Popen([sys.executable, "-u", __file__, "--worker"], stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            elapsed = time.monotonic() - started
            try:
                peak = max(peak, owned_process_rss(psutil.Process(process.pid)) / 2**20)
            except psutil.NoSuchProcess:
                pass
            if peak > manifest["config"]["memory_limit_mib"] or elapsed > manifest["config"]["worker_timeout_s"]:
                stopped = "memory limit" if peak > manifest["config"]["memory_limit_mib"] else "timeout"
                try:
                    children = psutil.Process(process.pid).children(recursive=True)
                except psutil.NoSuchProcess:
                    children = []
                for child in children:
                    try:
                        child.kill()
                    except psutil.NoSuchProcess:
                        pass
                process.kill()
                process.wait()
                break
            progress_path = OUTPUT / "worker-progress.json"
            if progress_path.exists():
                progress = json.loads(progress_path.read_text())
                f = progress["fraction"]
                report(.04 + .94 * f, progress["message"], eta_s=elapsed * (1-f) / max(f, .02),
                       peak_rss={"value": peak / 1024, "unit": "GiB"}, **progress["values"])
            else:
                report(.04, "Start bounded geometry worker", eta_s=180)
            time.sleep(2)
    for group, root in [("source_sha256", OUTPUT / "source"), ("input_sha256", Path())]:
        for name, expected in manifest[group].items():
            assert digest(root / name) == expected, name
    state = "complete" if process.returncode == 0 else "rejected" if process.returncode == 2 and stopped is None else "failed"
    write_json(OUTPUT / "completion.json", {"state": state, "worker_exit_code": process.returncode,
                "stopped": stopped, "peak_rss_mib": peak, "worker_seconds": time.monotonic()-started,
                "rss_scope": "Sum of owned worker and recursive child processes, sampled every 2 seconds",
                "input_hashes_unchanged": True, "source_hashes_unchanged": True, "fea_performed": False,
                "manifest_sha256": digest(OUTPUT / "manifest.json")})
    report(1, f"Geometry coarsening trial {state}", eta_s=0, peak_rss={"value": peak/1024, "unit": "GiB"})
    if process.returncode:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
