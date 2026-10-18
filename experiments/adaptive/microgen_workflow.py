"""Bounded G3 experiment using unmodified Microgen VTK and MMG defaults."""

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "experiments/adaptive"
CASES = ("gyroid", "schwarz-solid", "spinodoid", "smooth-lattice")
BOUNDS = (-0.5, 0.5) * 3
MEMORY = 8 * 1024**3
CAP = 250_000
TIMEOUT = 1200


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def imports():
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str(ROOT.parent / "periodic-microgen/experiments/periodic-fea"))
    import numpy as np
    import fedoo as fd
    import meshers
    import microgen
    import pyvista as pv
    import scorecard
    return np, fd, meshers, microgen, pv, scorecard


def specimen(name):
    np, _, _, microgen, _, scorecard = imports()
    from microgen.shape import from_field, surface_functions
    from microgen.shape.implicit_ops import smooth_union

    if name == "gyroid":
        shape = microgen.Tpms(surface_functions.gyroid, offset=1.1)
        return scorecard.Case(name, scorecard.gyroid, (-1.0, 1.0)), shape, "sheet"
    if name == "schwarz-solid":
        def field(x, y, z):
            return np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y) + np.cos(2 * np.pi * z)
        shape = microgen.Tpms(surface_functions.schwarz_p, offset=0.0)
        return scorecard.Case(name, field, (-np.inf, 0.0), 0.5), shape, "lower skeletal"
    if name == "spinodoid":
        shape = microgen.Spinodoid(
            k0=4 * np.pi, bandwidth=3, offset=0, cell_size=1,
            resolution=16, seed=20261010, mode_amplitude_threshold=0.05,
        )
        def field(x, y, z):
            arrays = np.broadcast_arrays(x, y, z)
            coordinates = [array.ravel() + 0.5 for array in arrays]
            values = np.empty(coordinates[0].size)
            for start in range(0, len(values), 4096):
                subset = slice(start, start + 4096)
                values[subset] = shape.evaluate(*(axis[subset] for axis in coordinates))
            return values.reshape(arrays[0].shape)
        return scorecard.Case(name, field, (-np.inf, 0.0)), shape, "solid"
    if name != "smooth-lattice":
        raise ValueError(name)
    shape = None
    for corner in scorecard.CORNERS:
        def capsule(x, y, z, corner=corner):
            parameter = np.clip((x * corner[0] + y * corner[1] + z * corner[2]) / 0.75, 0, 1)
            return np.sqrt(sum((coordinate - parameter * end)**2 for coordinate, end in
                               zip((x, y, z), corner, strict=True))) - 0.12
        strut = from_field(capsule, bounds=BOUNDS)
        shape = strut if shape is None else smooth_union(shape, strut, k=0.08)
    def field(x, y, z):
        distance = None
        for corner in scorecard.CORNERS:
            parameter = np.clip((x * corner[0] + y * corner[1] + z * corner[2]) / 0.75, 0, 1)
            radius = np.sqrt((x - parameter * corner[0])**2 +
                             (y - parameter * corner[1])**2 + (z - parameter * corner[2])**2)
            if distance is None:
                distance = radius
            else:
                blend = np.maximum(0.08 - np.abs(distance - radius), 0.0) / 0.08
                distance = np.minimum(distance, radius) - blend**3 * 0.08 / 6
        return distance - 0.12
    return scorecard.Case(name, field, (-np.inf, 0.0)), shape, "solid"


def manifest(run):
    np, fd, meshers, microgen, _, _ = imports()
    from microgen.remesh import remesh_keeping_boundaries_for_fem
    from microgen.shape.shape import Shape
    modules = [microgen, microgen.remesh, microgen.external,
               sys.modules["microgen.shape.tpms"], sys.modules["microgen.shape.spinodoid"],
               sys.modules["microgen.shape.shape"], sys.modules["microgen.shape.implicit_ops"]]
    source = run / "microgen-source"
    source.mkdir(exist_ok=True)
    hashes = {}
    for module in modules:
        path = Path(module.__file__)
        name = module.__name__ + ".py"
        shutil.copyfile(path, source / name)
        hashes[name] = digest(path)
    binary = next(Path(meshers.__file__).parent.glob("_meshers*.so"))
    mmg = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
    mmg.chmod(mmg.stat().st_mode | 0o111)
    versions = {dist.metadata["Name"]: dist.version for dist in importlib.metadata.distributions()}
    data = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "host": os.uname().nodename, "python": sys.version, "packages": versions,
        "native_sha256": digest(binary), "mmg_cli_sha256": digest(mmg),
        "source_sha256": hashes, "script_sha256": digest(Path(__file__)),
        "defaults": {"Tpms": str(inspect.signature(microgen.Tpms)),
                     "Spinodoid": str(inspect.signature(microgen.Spinodoid)),
                     "Shape.generate_volume_mesh": str(inspect.signature(Shape.generate_volume_mesh)),
                     "remesh_keeping_boundaries_for_fem": str(inspect.signature(remesh_keeping_boundaries_for_fem))},
        "material": {"E": 2.5, "nu": 0.25, "lambda": 1, "mu": 1},
        "limits": {"worker_rss_bytes": MEMORY, "worker_timeout_s": TIMEOUT,
                   "tetrahedra": CAP, "solve_nodes": 100_000, "threads": 1},
        "admission": {"sampled_distance": 0.002, "relative_volume": 0.01,
                      "corner_quality": 0.1, "periodic_tolerance": 1e-10},
        "scope": "Pilot. Finite sampled maxima, not Hausdorff bounds. No local-stress or continuum certification.",
    }
    save(run / "manifest.json", data)


def preparation(run, name):
    np, _, _, _, _, scorecard = imports()
    from scipy.stats import qmc
    case, shape, _ = specimen(name)
    probes = qmc.Sobol(3, scramble=True, seed=971).random_base2(10) - 0.5
    if name == "smooth-lattice":
        assert np.allclose(case.function(*probes.T), shape.evaluate(*probes.T), atol=1e-13, rtol=1e-13)
    periodic_errors, material_mismatches, near_boundary_errors = [], [], []
    for axis in range(3):
        low, high = probes.copy(), probes.copy()
        low[:, axis], high[:, axis] = -0.5, 0.5
        low_value, high_value = case.function(*low.T), case.function(*high.T)
        periodic_errors.append(float(np.max(np.abs(low_value - high_value))))
        low_material = (low_value >= case.levels[0]) & (low_value <= case.levels[1])
        high_material = (high_value >= case.levels[0]) & (high_value <= case.levels[1])
        material_mismatches.append(int((low_material != high_material).sum()))
        levels = np.array([level for level in case.levels if np.isfinite(level)])
        near = np.minimum(np.min(abs(low_value[:, None] - levels), axis=1),
                          np.min(abs(high_value[:, None] - levels), axis=1)) < 0.03
        near_boundary_errors.append(float(np.max(abs(low_value[near] - high_value[near]), initial=0)))
    if max(near_boundary_errors) > 1e-10 or any(material_mismatches):
        raise ValueError(f"Periodic material sample mismatch: {material_mismatches}, {near_boundary_errors}")
    estimates = []
    for seed in range(4):
        probes = qmc.Sobol(3, scramble=True, seed=seed).random_base2(18) - 0.5
        values = case.function(*probes.T)
        estimates.append(float(np.mean((values >= case.levels[0]) & (values <= case.levels[1]))))
    walls = scorecard.wall_samples(case, count=2**15)
    caps = []
    probes = qmc.Sobol(2, scramble=True, seed=119).random_base2(12) - 0.5
    for axis in range(3):
        for side in (-0.5, 0.5):
            points = np.empty((len(probes), 3))
            points[:, axis] = side
            points[:, [index for index in range(3) if index != axis]] = probes
            values = case.function(*points.T)
            caps.append(points[(values >= case.levels[0]) & (values <= case.levels[1])])
    np.savez_compressed(run / f"{name}-samples.npz", walls=walls, caps=np.concatenate(caps))
    save(run / f"{name}-specimen.json", {
        "volume": case.volume if case.volume is not None else float(np.mean(estimates)),
        "volume_sobol_estimates": estimates, "volume_scramble_spread": float(np.ptp(estimates)),
        "volume_samples_per_scramble": 2**18, "wall_samples": len(walls),
        "periodic_field_errors": periodic_errors,
        "periodic_material_mismatches": material_mismatches,
        "periodic_near_boundary_field_errors": near_boundary_errors,
        "definition": {"gyroid": "abs(sin X cos Y + sin Y cos Z + sin Z cos X) <= 0.55",
                       "schwarz-solid": "cos X + cos Y + cos Z <= 0",
                       "spinodoid": "Fixed Microgen spectral realization, resolution=16, seed=20261010, k0=4*pi, bandwidth=3, threshold=0.05, offset=0; native [0,1] translated by -0.5",
                       "smooth-lattice": "Eight capsules from center to corners, radius=0.12; sequential Microgen cubic smooth_union k=0.08 in scorecard.CORNERS order"}[name],
    })


def tri_values(barycentric):
    np, *_ = imports()
    first, second, third = np.asarray(barycentric).T
    return np.array([first * (2 * first - 1), second * (2 * second - 1),
                     third * (2 * third - 1), 4 * first * second,
                     4 * second * third, 4 * third * first]).T


def tessellate(points, surface, subdivisions):
    np, *_ = imports()
    if surface.shape[1] == 3:
        return points, surface
    coordinates = [(first, second) for first in range(subdivisions + 1)
                   for second in range(subdivisions + 1 - first)]
    lookup = {pair: index for index, pair in enumerate(coordinates)}
    barycentric = np.array([[first / subdivisions, second / subdivisions,
                            1 - (first + second) / subdivisions] for first, second in coordinates])
    triangles = []
    for first, second in coordinates:
        if first + second < subdivisions:
            triangles.append([lookup[first, second], lookup[first + 1, second], lookup[first, second + 1]])
        if first + second < subdivisions - 1:
            triangles.append([lookup[first + 1, second], lookup[first + 1, second + 1], lookup[first, second + 1]])
    vertices = np.einsum("sn,fnd->fsd", tri_values(barycentric), points[surface])
    faces = (np.array(triangles)[None] + np.arange(len(surface))[:, None, None] * len(barycentric))
    return vertices.reshape(-1, 3), faces.reshape(-1, 3)


def reorient_linear(points, tetrahedra):
    np, *_ = imports()
    corners = points[tetrahedra]
    negative = np.linalg.det(corners[:, 1:] - corners[:, :1]) < 0
    tetrahedra = tetrahedra.copy()
    tetrahedra[negative, 1], tetrahedra[negative, 2] = tetrahedra[negative, 2].copy(), tetrahedra[negative, 1].copy()
    return tetrahedra, int(negative.sum())


def on_box_plane(coordinates, value):
    np, *_ = imports()
    return np.isclose(coordinates, value, atol=1e-12, rtol=0)


def geometry(run, name, points, tets, surface):
    np, fd, _, _, _, scorecard = imports()
    case, _, _ = specimen(name)
    samples = np.load(run / f"{name}-samples.npz")
    reference = json.loads((run / f"{name}-specimen.json").read_text())
    corners = points[surface[:, :3]]
    caps = np.any(np.all(on_box_plane(corners, -0.5), axis=1) |
                  np.all(on_box_plane(corners, 0.5), axis=1), axis=1)
    wall = surface[~caps]
    weights = tri_values(scorecard.BARY) if surface.shape[1] == 6 else scorecard.BARY
    probes = np.einsum("sn,fnd->fsd", weights, points[wall]).reshape(-1, 3)
    locked = on_box_plane(np.abs(probes), 0.5)
    foot = probes.copy()
    finite_levels = np.array([level for level in case.levels if np.isfinite(level)])
    values = case.function(*probes.T)
    levels = finite_levels[np.argmin(np.abs(values[:, None] - finite_levels), axis=1)]
    for _ in range(40):
        gradient = scorecard.gradient(case.function, foot)
        gradient[locked] = 0
        residual = case.function(*foot.T) - levels
        step = residual[:, None] * gradient / np.maximum((gradient**2).sum(1), 1e-30)[:, None]
        foot -= step
        if np.max(np.linalg.norm(step, axis=1)) < 1e-12:
            break
    residual = np.abs(case.function(*foot.T) - levels) / np.maximum(np.linalg.norm(scorecard.gradient(case.function, foot), axis=1), 1e-30)
    failed = (residual > 1e-9) | np.any(np.abs(foot) > 0.5 + 1e-9, axis=1)
    forward = np.linalg.norm(foot - probes, axis=1)
    outside = (values < case.levels[0]) | (values > case.levels[1])
    queries = np.concatenate([samples["walls"], samples["caps"]])
    reverse_values = []
    for divisions in ((4, 8) if surface.shape[1] == 6 else (1,)):
        vertices, faces = tessellate(points, surface, divisions)
        reverse_values.append(scorecard.distance_to_mesh(vertices, faces, queries))
    quality, angles, corner_volumes = scorecard.tet_metrics(points, tets[:, :4])
    if tets.shape[1] == 10:
        element = fd.lib_elements.element_list.get_element("tet10")(15)
        coordinates = points[tets]
        determinants = np.array([np.linalg.det(np.einsum("an,tnb->tab", derivative, coordinates))
                                 for derivative in element.shape_function_derivative(element.xi_pg)])
        volume = float(np.einsum("gt,g->", determinants, element.w_pg))
        extra = element.shape_function_derivative(np.concatenate([element.xi_nd, [[0.25]*3]]))
        minimum = min(float(determinants.min()), *(float(np.linalg.det(np.einsum("an,tnb->tab", derivative, coordinates)).min()) for derivative in extra))
    else:
        volume = float(corner_volumes.sum())
        minimum = float(corner_volumes.min() * 6)
    reverse = float(reverse_values[-1].max())
    row = {"sampled_forward_max": float(forward.max()), "sampled_reverse_max": reverse,
           "sampled_two_sided_max": float(max(forward.max(), reverse)),
           "outward_forward_max": float(forward[outside].max(initial=0)),
           "inward_forward_max": float(forward[~outside].max(initial=0)),
           "projection_failed": int(failed.sum()), "forward_samples": len(probes),
           "reverse_samples": len(queries), "volume": volume,
           "relative_volume_error": volume / reference["volume"] - 1,
           "q_min_corner": float(quality.min()), "q_p01_corner": float(np.quantile(quality, 0.01)),
           "minimum_dihedral_deg_corner": float(angles[np.isfinite(angles)].min()) if np.isfinite(angles).any() else None,
           "minimum_sampled_jacobian": minimum,
           "reverse_tessellation_delta_max": float(np.max(np.abs(reverse_values[-1] - reverse_values[0]))),
           "quadratic_reverse_subdivisions": [4, 8] if tets.shape[1] == 10 else None}
    row["geometry_admitted"] = bool(row["sampled_two_sided_max"] <= 0.002 and
                                    abs(row["relative_volume_error"]) <= 0.01 and
                                    row["q_min_corner"] >= 0.1 and not failed.any() and minimum > 0)
    return row


def worker(run, name, route, cells):
    np, fd, meshers, microgen, pv, scorecard = imports()
    from homogenization import AuditedCG
    from periodic import periodic_pairs
    row = {"case": name, "route": route, "cells": cells, "status": "started"}
    output = run / f"{name}-{route}-{cells}.json"
    save(output, row)
    try:
        started = time.perf_counter()
        case, shape, part = specimen(name)
        row["specimen_setup_s"] = time.perf_counter() - started
        if route.startswith("control-"):
            archive = run / f"{name}-meshers-curved-12.npz"
            data = np.load(archive)
            points, tets, surface = (data[key].copy() for key in ("points", "tetrahedra", "surface"))
            row["replay_source_sha256"] = digest(archive)
            row["workflow_scope"] = "Saved-mesh same-geometry FEA control; no live generation timing"
            if route == "control-linear":
                used = np.unique(tets[:, :4])
                mapping = np.full(len(points), -1, dtype=int)
                mapping[used] = np.arange(len(used))
                points, tets, surface = points[used], mapping[tets[:, :4]], mapping[surface[:, :3]]
            elif route == "control-straight":
                for index, (first, second) in enumerate(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))):
                    points[tets[:, 4 + index]] = (points[tets[:, first]] + points[tets[:, second]]) / 2
            else:
                raise ValueError(route)
        elif route in {"vtk-default", "mmg-default"}:
            if name in {"gyroid", "schwarz-solid"}:
                grid = shape.generate_volume_mesh(type_part=part).triangulate()
            else:
                grid = shape.generate_volume_mesh().triangulate()
                if name == "spinodoid":
                    grid.translate((-0.5,) * 3, inplace=True)
            row["vtk_extraction_s"] = time.perf_counter() - started
            row["seed_nodes"], row["seed_elements"] = grid.n_points, grid.n_cells
            save(output, row)
            if grid.n_cells > CAP:
                raise ValueError("VTK seed exceeds staged element cap before MMG")
            if route == "mmg-default":
                from microgen.remesh import remesh_keeping_boundaries_for_fem
                from unittest.mock import patch
                binary = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
                binary.chmod(binary.stat().st_mode | 0o111)
                os.environ["PATH"] = str(binary.parent) + os.pathsep + os.environ["PATH"]
                row["mmg_arguments"] = {"periodic": True, "mesh_version": 2, "dimension": 3,
                                        "tol": 1e-8, "hausd": None, "hgrad": None,
                                        "hmax": None, "hmin": None, "hsiz": None}
                save(output, row)
                original = microgen.Mmg.mmg3d
                row["mmg_calls"] = []
                def record_call(**parameters):
                    row["mmg_calls"].append({key: value for key, value in parameters.items()
                                              if key not in {"input", "output"}})
                    save(output, row)
                    return original(**parameters)
                with patch.object(microgen.Mmg, "mmg3d", record_call):
                    grid = remesh_keeping_boundaries_for_fem(grid)
            points = np.asarray(grid.points).copy()
            if set(grid.cells_dict) != {int(pv.CellType.TETRA)}:
                raise ValueError("Microgen output is not exclusively Tet4")
            tets = np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()
            poly = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
            surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[poly.faces.reshape(-1, 4)[:, 1:]]
        else:
            order = 1 if route == "meshers-linear" else 2
            mesh = meshers.generate(case.function, bounds=BOUNDS, cells=cells, band=case.band,
                                    periodic=(True, True, True), order=order, threads=1,
                                    optimize_passes=1, minimum_quality=0.1,
                                    geometry_tolerance=1.0, max_tetrahedra=CAP)
            points, tets, surface = np.asarray(mesh.points), np.asarray(mesh.tetrahedra), np.asarray(mesh.surface)
            row["native_diagnostics"] = mesh.diagnostics
            row["native_arguments"] = {"geometry_tolerance": 1.0, "minimum_quality": 0.1,
                                       "optimize_passes": 1, "order": order, "periodic": [True]*3}
        if route.startswith("control-"):
            row["replay_load_s"] = time.perf_counter() - started
            row["generation_s"] = None
        else:
            row["generation_s"] = time.perf_counter() - started
        row["nodes"], row["elements"], row["order"] = len(points), len(tets), 2 if tets.shape[1] == 10 else 1
        save(output, row)
        if len(tets) > CAP:
            raise ValueError("Generated mesh exceeds staged element cap")
        corner = points[tets[:, :4]]
        determinants = np.linalg.det(corner[:, 1:] - corner[:, :1])
        row["input_signed_volume_ordering_diagnostic"] = float(determinants.sum() / 6)
        row["reoriented_tetrahedra"] = 0
        if tets.shape[1] == 4:
            tets, row["reoriented_tetrahedra"] = reorient_linear(points, tets)
        row["nonpositive_or_tiny_tetrahedra"] = int((abs(determinants) <= 1e-18).sum())
        audit_start = time.perf_counter()
        row["geometry"] = geometry(run, name, points, tets, surface)
        row["audit_s"] = time.perf_counter() - audit_start
        np.savez_compressed(output.with_suffix(".npz"), points=points, tetrahedra=tets, surface=surface)
        save(output, row)
        if row["nonpositive_or_tiny_tetrahedra"]:
            raise ValueError("Tiny absolute-volume tetrahedra; no cleanup or cell deletion")
        try:
            _, _, row["periodic"] = periodic_pairs(points + 0.5, surface[:, :3])
        except Exception as error:
            row["periodic_error"] = str(error)
            raise ValueError("Periodic topology unsupported; no nonmatching-BC substitution") from error
        if len(points) > 100_000:
            raise ValueError("Solve node cap exceeded")
        if row["geometry"]["minimum_sampled_jacobian"] <= 0:
            raise ValueError("Invalid sampled Jacobian; solve withheld")
        row["qualified"] = row["geometry"]["geometry_admitted"]
        solve_start = time.perf_counter()
        fd.ModelingSpace("3D")
        fem = fd.Mesh(points + 0.5, tets.copy(), "tet10" if row["order"] == 2 else "tet4")
        assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), fem)
        solver = AuditedCG()
        stiffness = np.asarray(fd.homogen.get_homogenized_stiffness(
            assembly, meshperio=True, solver=solver, rigid_body_constraint="pin"))
        row["solve_s"] = time.perf_counter() - solve_start
        row["stiffness_skew_relative"] = float(np.linalg.norm(stiffness - stiffness.T) / np.linalg.norm(stiffness))
        if row["stiffness_skew_relative"] > 1e-8 or len(solver.records) != 6:
            raise ValueError("Six-load tensor audit failed")
        row["stiffness"] = ((stiffness + stiffness.T) / 2).tolist()
        row["solver_records"] = solver.records
        row["status"] = "solved" if row["qualified"] else "diagnostic-unqualified"
        row["replay_workflow_s" if route.startswith("control-") else "total_workflow_s"] = time.perf_counter() - started
    except Exception as error:
        row.update(status="failed", error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc())
        row["failed_workflow_s"] = time.perf_counter() - started
    save(output, row)


def summarize(run):
    np, *_ = imports()
    from periodic import tensor_difference
    import fea
    import fea_solid
    rows = [json.loads(path.read_text()) for path in sorted(run.glob("*-*-*.json"))
            if path.name != "manifest.json" and "-specimen" not in path.name]
    workers = json.loads((run / "workers.json").read_text())
    for row in rows:
        key = f"{row['case']}-{row['route']}-{row['cells']}"
        row["worker"] = workers.get(key, {})
        if row["worker"].get("returncode"):
            row["status"] = "resource-failure" if row["worker"].get("limit_failure") else "process-failure"
    references = {}
    for name in CASES:
        candidates = {row["cells"]: row for row in rows if row.get("case") == name and
                      row.get("route") == "meshers-curved" and "stiffness" in row}
        if 16 in candidates and 20 in candidates:
            coarse, fine = (np.array(candidates[cells]["stiffness"]) for cells in (16, 20))
            extrapolated = (20**4 * fine - 16**4 * coarse) / (20**4 - 16**4)
            reference = {"stiffness": extrapolated.tolist(), "assumed_order": 4,
                         "fine_to_extrapolated": tensor_difference(fine, extrapolated),
                         "fine_coarse": tensor_difference(coarse, fine),
                         "qualified_inputs": candidates[16].get("qualified", False) and candidates[20].get("qualified", False),
                         "justification": "Two-level p=4 extrapolation; assumed order, not a certified error bound"}
            if 12 in candidates:
                previous = (16**4 * coarse - 12**4 * np.array(candidates[12]["stiffness"])) / (16**4 - 12**4)
                reference["extrapolation_window_difference"] = tensor_difference(previous, extrapolated)
            if name in {"gyroid", "schwarz-solid"}:
                historical = fea.reference() if name == "gyroid" else fea_solid.reference()
                reference["historical_crosscheck"] = tensor_difference(extrapolated, historical)
            references[name] = reference
    for row in rows:
        if "stiffness" in row and row["case"] in references:
            row["mixed_load_reference_disagreement"] = tensor_difference(np.array(row["stiffness"]), np.array(references[row["case"]]["stiffness"]))
    decomposition = {}
    from scipy.linalg import eigvalsh
    for name in CASES:
        selected = {row["route"]: np.array(row["stiffness"]) for row in rows
                    if row.get("case") == name and row.get("cells") == 12 and "stiffness" in row}
        if all(route in selected for route in ("control-linear", "control-straight", "meshers-curved")):
            discretization = selected["control-straight"] - selected["control-linear"]
            geometry_change = selected["meshers-curved"] - selected["control-straight"]
            assert np.allclose(discretization + geometry_change, selected["meshers-curved"] - selected["control-linear"])
            result = {"discretization_delta_C_same_straight_domain": discretization.tolist(),
                      "geometry_delta_C_same_quadratic_order": geometry_change.tolist()}
            if name in references:
                reference = np.array(references[name]["stiffness"])
                result["discretization_relative_eigenvalues"] = eigvalsh(discretization, reference).tolist()
                result["geometry_relative_eigenvalues"] = eigvalsh(geometry_change, reference).tolist()
            decomposition[name] = result
    save(run / "summary.json", {"references": references, "rows": rows,
                                "same_skeleton_decomposition": decomposition,
                                "limits": "Diagnostic tensors are not admitted results. Same-skeleton deltas separate order and mapping changes, not certified continuum errors. Local stress was not measured."})
    return rows


def self_check():
    np, fd, meshers, _, _, _ = imports()
    from homogenization import AuditedCG
    from periodic import tensor_difference
    nodes = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                      [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
    assert np.allclose(tri_values(nodes), np.eye(6))
    assert np.allclose(tri_values([[0.2, 0.3, 0.5]]).sum(), 1)
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0]])
    vertices, faces = tessellate(points, np.arange(6)[None], 4)
    assert len(faces) == 16
    area = np.linalg.norm(np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                                   vertices[faces[:, 2]] - vertices[faces[:, 0]]), axis=1).sum() / 2
    assert abs(area - 0.5) < 1e-14
    element = fd.lib_elements.element_list.get_element("tet10")(15)
    assert abs(element.w_pg.sum() - 1 / 6) < 1e-14
    corners = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    oriented, changed = reorient_linear(corners, np.array([[0, 2, 1, 3], [0, 1, 2, 3]]))
    assert changed == 1 and np.array_equal(oriented, np.array([[0, 1, 2, 3], [0, 1, 2, 3]]))
    stiffness = np.eye(6)
    candidate = stiffness.copy()
    candidate[0, 3] = candidate[3, 0] = 0.1
    assert abs(tensor_difference(candidate, stiffness)["worst_relative_strain_energy_difference"] - 0.1) < 1e-14
    cube = meshers.generate(lambda x, y, z: x * 0 - 1, bounds=BOUNDS, cells=4,
                            order=2, periodic=(True, True, True), threads=1)
    fd.ModelingSpace("3D")
    fem = fd.Mesh(cube.points + 0.5, cube.tetrahedra.copy(), "tet10")
    assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), fem)
    solver = AuditedCG()
    actual = np.asarray(fd.homogen.get_homogenized_stiffness(assembly, meshperio=True, solver=solver, rigid_body_constraint="pin"))
    expected = np.eye(6)
    expected[:3, :3] = 1
    np.fill_diagonal(expected[:3, :3], 3)
    assert np.allclose(actual, expected, atol=1e-10, rtol=1e-10)
    assert len(solver.records) == 6
    print("Interpolation, curved tessellation, volume weights, mixed-load score and native Tet10 fedoo affine patch checks pass")


def bounded(command, logfile, part, timeout=TIMEOUT):
    import psutil
    start = time.monotonic()
    peak = 0
    def limit_address_space():
        resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    with logfile.open("w") as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True, preexec_fn=limit_address_space)
        reason = None
        while process.poll() is None:
            elapsed = time.monotonic() - start
            try:
                family = [psutil.Process(process.pid)] + psutil.Process(process.pid).children(recursive=True)
                rss = sum(member.memory_info().rss for member in family if member.is_running())
                peak = max(peak, rss)
            except psutil.Error:
                rss = 0
            try:
                part.update(min(elapsed / timeout, 0.99), f"{logfile.stem}: {elapsed:.0f}s / {timeout}s cap",
                            rss_gib=rss / 1024**3, peak_gib=peak / 1024**3,
                            elapsed_s=elapsed, timeout_remaining_s=max(timeout - elapsed, 0))
            except Exception:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise
            if peak > MEMORY or elapsed > timeout:
                reason = "memory-cap" if peak > MEMORY else "timeout"
                os.killpg(process.pid, signal.SIGKILL)
                break
            time.sleep(0.25)
        code = process.wait()
    return {"returncode": code, "limit_failure": reason, "peak_rss_bytes": peak,
            "wall_s": time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", nargs=3, metavar=("CASE", "ROUTE", "CELLS"))
    parser.add_argument("--prepare")
    parser.add_argument("--manifest", action="store_true")
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    if args.worker:
        worker(run, args.worker[0], args.worker[1], int(args.worker[2]))
        return
    if args.prepare:
        preparation(run, args.prepare)
        return
    if args.manifest:
        manifest(run)
        return
    if args.summarize:
        summarize(run)
        return
    if args.self_check:
        self_check()
        return
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job
    if run.exists():
        raise SystemExit("Run exists; keep it frozen and select a new directory")
    run.mkdir(parents=True)
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                      RAYON_NUM_THREADS="1", CARGO_BUILD_JOBS="2", VIRTUAL_ENV=str(ROOT / ".venv"),
                      LD_LIBRARY_PATH=str(ROOT / ".local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu"),
                      PATH=str(Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"],
                      RUSTUP_TOOLCHAIN="1.95.0")
    shutil.copyfile(__file__, run / "driver.py")
    requirements = subprocess.check_output(["uv", "pip", "freeze", "--python", sys.executable], text=True)
    (run / "requirements.txt").write_text(requirements)
    tasks = [("build", [str(ROOT / ".venv/bin/maturin"), "develop", "--release", "--locked",
                         "--manifest-path", str(ROOT / "crates/meshers-python/Cargo.toml"),
                         "--features", "experimental-surfaces"], 1200),
             ("manifest", [sys.executable, __file__, "--run", str(run), "--manifest"], 120)]
    tasks.append(("self-check", [sys.executable, __file__, "--run", str(run), "--self-check"], 120))
    for name in CASES:
        tasks.append((name + "-prepare", [sys.executable, __file__, "--run", str(run), "--prepare", name], 600))
        for route, cells in [("vtk-default", 0), ("mmg-default", 0), ("meshers-linear", 24),
                             ("meshers-curved", 12), ("meshers-curved", 16), ("meshers-curved", 20)]:
            tasks.append((f"{name}-{route}-{cells}", [sys.executable, __file__, "--run", str(run),
                                                     "--worker", name, route, str(cells)], TIMEOUT))
    tasks.append(("summarize", [sys.executable, __file__, "--run", str(run), "--summarize"], 120))
    records = {}
    if args.baseline:
        baseline = args.baseline.resolve()
        old_workers = json.loads((baseline / "workers.json").read_text())
        provenance = {}
        for path in baseline.glob("*-*-*.json"):
            row = json.loads(path.read_text())
            if "case" not in row or row.get("route") == "vtk-default":
                continue
            name = path.stem
            for artifact in (path, path.with_suffix(".npz")):
                if artifact.exists():
                    shutil.copyfile(artifact, run / artifact.name)
                    provenance[artifact.name] = {"source": str(artifact), "sha256": digest(artifact)}
            records[name] = old_workers[name] | {"inherited_from": str(baseline)}
        for name in CASES:
            for suffix in ("-specimen.json", "-samples.npz"):
                path = baseline / (name + suffix)
                if path.exists():
                    shutil.copyfile(path, run / path.name)
                    provenance[path.name] = {"source": str(path), "sha256": digest(path)}
            if (run / (name + "-specimen.json")).exists():
                records[name + "-prepare"] = old_workers[name + "-prepare"] | {"inherited_from": str(baseline)}
        save(run / "inherited.json", provenance)
        save(run / "workers.json", records)
        tasks = [(name, command, timeout) for name, command, timeout in tasks
                 if name in {"build", "manifest", "self-check", "summarize"} or name not in records]
        controls = []
        for name in CASES:
            for route in ("control-linear", "control-straight"):
                controls.append((f"{name}-{route}-12", [sys.executable, __file__, "--run", str(run),
                                                        "--worker", name, route, "12"], TIMEOUT))
        tasks[-1:-1] = controls
    with Job(parts={name: 1 for name, _, _ in tasks}) as job:
        for name, command, timeout in tasks:
            with job.part(name) as part:
                part.update(0, "Starting " + name, worker_cap_gib=8, threads=1)
                if name not in {"build", "manifest", "self-check", "summarize"}:
                    case = next(case for case in CASES if name.startswith(case))
                    if "prepare" not in name and records.get(case + "-prepare", {}).get("returncode") != 0:
                        records[name] = {"skipped": "specimen preparation failed"}
                        part.update(1, "Skipped after specimen preparation failure", skipped=1)
                        continue
                    if "-control-" in name and not (run / f"{case}-meshers-curved-12.npz").exists():
                        records[name] = {"skipped": "No native quadratic mesh survived the original gate"}
                        part.update(1, "Skipped; native quality gate rejected source mesh", skipped=1)
                        continue
                result = bounded(command, run / f"{name}.log", part, timeout)
                records[name] = result
                save(run / "workers.json", records)
                if name == "manifest" and args.baseline:
                    current = json.loads((run / "manifest.json").read_text())
                    previous = json.loads((args.baseline / "manifest.json").read_text())
                    assert current["native_sha256"] == previous["native_sha256"]
                    assert current["packages"] == previous["packages"]
                job.result(name, returncode=result["returncode"], seconds=result["wall_s"],
                           peak_gib=result["peak_rss_bytes"] / 1024**3, jsonl=run / "queue-results.jsonl")
                if name in {"build", "manifest", "self-check"} and result["returncode"]:
                    raise RuntimeError(f"{name} failed; inspect {run / (name + '.log')}")
    assert digest(Path(__file__)) == digest(run / "driver.py"), "Driver changed during run"


if __name__ == "__main__":
    main()
