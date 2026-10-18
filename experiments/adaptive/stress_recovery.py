"""Check periodic stress recovery against fixed-domain polynomial enrichment."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job
import microgen_workflow as workflow


def solve(points, tetrahedra, checkpoint=None, *, n_elm_gp=None, element_type=None, solver=None):
    import fedoo as fd
    import numpy as np
    from homogenization import AuditedCG, COMPONENTS

    fd.ModelingSpace("3D")
    if element_type is None:
        element_type = {4: "tet4", 10: "tet10"}[tetrahedra.shape[1]]
    assert fd.lib_elements.element_list.get_element(element_type).n_nodes == tetrahedra.shape[1]
    mesh = fd.Mesh(points.copy(), tetrahedra.copy(), element_type)
    assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(
        fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), mesh, n_elm_gp=n_elm_gp)
    solver = AuditedCG() if solver is None else solver
    fields, tensor = [], np.zeros((6, 6))
    for load, engineering in enumerate(np.eye(6)):
        problem = fd.problem.Linear(assembly)
        if hasattr(solver, "problem"):
            solver.problem = problem
        problem.bc.add(fd.constraint.PeriodicBC("small_strain", meshperio=True))
        problem.bc.add("Dirichlet", [int(np.linalg.norm(points, axis=1).argmin())], "Disp", 0)
        for component, value in zip(COMPONENTS, engineering, strict=True):
            problem.bc.add("Dirichlet", component, float(value))
        problem.set_solver(solver)
        problem.solve()
        fields.append(problem.get_disp().T.copy())
        tensor[:, load] = problem.get_ext_forces("MeanStrain").ravel()
        if checkpoint is not None:
            checkpoint(load + 1, tensor[:, :load + 1], solver.records)
    return np.asarray(fields), (tensor + tensor.T) / 2, solver.records


def stresses(points, tetrahedra, fields, derivative):
    import numpy as np

    jacobian = np.einsum("an,tnb->tab", derivative, points[tetrahedra])
    gradients = np.einsum("tba,an->tbn", np.linalg.inv(jacobian), derivative)
    gradient = np.einsum("ltni,tjn->tlij", fields[:, tetrahedra], gradients)
    strain = (gradient + gradient.swapaxes(-1, -2)) / 2
    return 2 * strain + np.trace(strain, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)


def gram(stress):
    import numpy as np

    trace = np.trace(stress, axis1=-2, axis2=-1)
    return 0.5 * np.einsum("tlij,tmij->tlm", stress, stress) - 0.1 * np.einsum("tl,tm->tlm", trace, trace)


def recover(stress, volumes, tetrahedra, classes):
    import numpy as np

    totals = np.zeros((len(classes),) + stress.shape[1:])
    weights = np.zeros(len(classes))
    np.add.at(totals, classes[tetrahedra.ravel()], np.repeat(stress * volumes[:, None, None, None], 4, axis=0))
    np.add.at(weights, classes[tetrahedra.ravel()], np.repeat(volumes, 4))
    assert np.all(weights[classes] > 0)
    return totals[classes] / weights[classes, None, None, None]


def periodic_displacement_jump_error(points, surface, fields):
    np, *_ = workflow.imports()
    from homogenization import strain_tensor
    from periodic import periodic_pairs

    assert fields.shape == (6, len(points), 3) and np.isfinite(fields).all()
    matching, _, _ = periodic_pairs(points, surface)
    maximum = 0.0
    for load, engineering in enumerate(np.eye(6)):
        for paired in matching:
            shift = points[paired[:, 1]] - points[paired[:, 0]]
            jump = fields[load, paired[:, 1]] - fields[load, paired[:, 0]]
            maximum = max(maximum, float(abs(jump - shift @ strain_tensor(engineering).T).max(initial=0)))
    return maximum


def worker(args):
    np, fd, _, _, _, _ = workflow.imports()
    from periodic import periodic_pairs, tensor_difference
    from scipy.linalg import eigvalsh
    from scipy.stats import spearmanr
    from tet10_integration import integrated_tensors

    started = time.perf_counter()
    meshes, fields, tensors, records, consistency = {}, {}, {}, {}, {}
    source_row = None
    if args.skeleton:
        source_row = json.loads((args.skeleton / "result.json").read_text())
        assert source_row["case"] == args.case and source_row["order"] == 2
        assert source_row["qualified"] and source_row["geometry"]["geometry_admitted"]
        original = np.load(args.skeleton / "mesh.npz")
        corner_tets = original["tetrahedra"][:, :4]
        corner_count = int(corner_tets.max()) + 1
        np.testing.assert_array_equal(np.unique(corner_tets), np.arange(corner_count))
        meshes["linear"] = {"points": original["points"][:corner_count].copy(),
            "tetrahedra": corner_tets.copy(), "surface": original["surface"][:, :3].copy()}
        straight_points = original["points"].copy()
        edge_order = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))
        straight_points[original["tetrahedra"][:, 4:].ravel()] = original["points"][
            corner_tets[:, edge_order]].mean(axis=2).reshape(-1, 3)
        meshes["straight"] = {"points": straight_points,
            "tetrahedra": original["tetrahedra"].copy(), "surface": original["surface"].copy()}
        np.testing.assert_array_equal(straight_points[original["surface"][:, 3:]],
            original["points"][original["surface"][:, :3][:, edge_order[:3]]].mean(axis=2))
    for route in ("linear", "straight"):
        if not args.skeleton:
            label = f"{args.case}-control-{route}-12"
            archive = np.load(args.source / (label + ".npz"))
            meshes[route] = {key: archive[key] for key in ("points", "tetrahedra", "surface")}
        data = meshes[route]
        periodic_pairs(data["points"], data["surface"])
        fields[route], tensors[route], records[route] = solve(data["points"], data["tetrahedra"],
            n_elm_gp=4 if args.skeleton else None)
        if not args.skeleton:
            saved = json.loads((args.source / (label + ".json")).read_text())
            assert tensor_difference(tensors[route], np.array(saved["stiffness"]))["worst_relative_strain_energy_difference"] < 1e-7
        stress, energy, volume = integrated_tensors(data["points"], data["tetrahedra"], fields[route])
        reaction_error = float(np.linalg.norm(stress - tensors[route]) / np.linalg.norm(tensors[route]))
        energy_error = float(np.linalg.norm(energy - tensors[route]) / np.linalg.norm(tensors[route]))
        assert max(reaction_error, energy_error) < 1e-7
        jump_error = periodic_displacement_jump_error(data["points"], data["surface"], fields[route])
        assert jump_error < 1e-9
        assert len(records[route]) == 6 and max(record["relative_residual"] for record in records[route]) < 2e-9
        consistency[route] = {"reaction_relative_error": reaction_error,
            "hill_mandel_relative_error": energy_error, "integrated_volume": volume,
            "periodic_jump_max": jump_error}
        if args.skeleton:
            np.savez_compressed(args.run / (args.case + "-control-" + route + ".npz"), **data)
    linear, quadratic = meshes["linear"], meshes["straight"]
    points, tetrahedra = linear["points"], linear["tetrahedra"]
    assert np.array_equal(points[tetrahedra], quadratic["points"][quadratic["tetrahedra"][:, :4]])
    edges = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))
    assert np.allclose(quadratic["points"][quadratic["tetrahedra"][:, 4:]], points[tetrahedra[:, edges]].mean(axis=2), atol=1e-14)
    _, classes, pairing = periodic_pairs(points, linear["surface"])
    element = fd.lib_elements.element_list.get_element("tet4")(4)
    enriched = fd.lib_elements.element_list.get_element("tet10")(4)
    volumes = np.linalg.det(points[tetrahedra[:, 1:]] - points[tetrahedra[:, :1]]) / 6
    assert np.all(volumes > 0)
    constant_stress = stresses(points, tetrahedra, fields["linear"], element.shape_function_derivative([[0.25] * 3])[0])
    recovered = recover(constant_stress, volumes, tetrahedra, classes)
    constant_check = np.broadcast_to(np.eye(3), constant_stress.shape).copy()
    assert np.allclose(recover(constant_check, volumes, tetrahedra, classes),
                       np.broadcast_to(np.eye(3), recovered.shape), atol=1e-14)
    estimated = np.zeros((len(tetrahedra), 6, 6))
    observed = np.zeros_like(estimated)
    for coordinate, weight in zip(element.xi_pg, element.w_pg, strict=True):
        interpolation = element.shape_function(coordinate[None])[0]
        recovered_stress = np.einsum("n,tnlij->tlij", interpolation, recovered[tetrahedra])
        higher_stress = stresses(quadratic["points"], quadratic["tetrahedra"], fields["straight"],
                                enriched.shape_function_derivative(coordinate[None])[0])
        estimated += 6 * weight * volumes[:, None, None] * gram(recovered_stress - constant_stress)
        observed += 6 * weight * volumes[:, None, None] * gram(higher_stress - constant_stress)
    expected_change = tensors["linear"] - tensors["straight"]
    energy_identity_error = float(np.linalg.norm(observed.sum(axis=0) - expected_change) / np.linalg.norm(expected_change))
    assert energy_identity_error < 1e-7
    inverse = np.linalg.inv(tensors["linear"])
    estimated_mass = np.einsum("ij,tji->t", inverse, estimated)
    observed_mass = np.einsum("ij,tji->t", inverse, observed)
    assert np.all(estimated_mass >= -1e-14) and np.all(observed_mass >= -1e-14)
    magnitude = volumes * np.einsum("ij,tji->t", inverse, gram(constant_stress))
    rankings = {}
    for name, indicator in (("recovery", estimated_mass), ("stress-magnitude", magnitude),
                            ("volume", volumes), ("enrichment-oracle", observed_mass)):
        selected = np.argsort(indicator)[-int(np.ceil(0.2 * len(indicator))):]
        rankings[name] = float(observed_mass[selected].sum() / observed_mass.sum())
    density_correlation = float(spearmanr(estimated_mass / volumes, observed_mass / volumes).statistic)
    worth_testing = (rankings["recovery"] > max(rankings["volume"], rankings["stress-magnitude"])
        and density_correlation > 0.5)
    control_geometry = None
    if args.skeleton:
        control_geometry = workflow.geometry(args.run, args.case, points, tetrahedra, linear["surface"])
        assert control_geometry["q_min_corner"] >= 0.1 and control_geometry["minimum_sampled_jacobian"] > 0
    np.testing.assert_allclose([value["integrated_volume"] for value in consistency.values()],
        volumes.sum(), atol=1e-12, rtol=1e-10)
    length = np.cbrt(6 * volumes)
    coefficient = np.maximum(estimated_mass / (volumes * length**2), 1e-30)
    target = coefficient**(-1 / 5)
    target *= (np.sum(volumes / target**3) / np.sum(volumes / length**3))**(1 / 3)
    nodal_target = np.full(len(points), np.inf)
    np.minimum.at(nodal_target, classes[tetrahedra.ravel()], np.repeat(target, 4))
    nodal_target = nodal_target[classes]
    assert np.isfinite(nodal_target).all() and np.all(nodal_target > 0)
    np.savez_compressed(args.run / (args.case + "-recovery.npz"), points=points, tetrahedra=tetrahedra,
        estimated_mass=estimated_mass, observed_mass=observed_mass, target_length=nodal_target,
        linear_displacement=fields["linear"], straight_displacement=fields["straight"])
    workflow.save(args.run / (args.case + ".json"), {
        "case": args.case, "elements": len(tetrahedra), "linear_nodes": len(points),
        "straight_nodes": len(quadratic["points"]), "periodic": pairing,
        "energy_identity_relative_error": energy_identity_error,
        "top_20_percent_elements_enrichment_energy_capture": rankings,
        "density_spearman": density_correlation,
        "live_sizing_pilot_worth_testing": bool(worth_testing),
        "selector_scope": "Recovery capture beats volume and stress magnitude; density correlation >0.5; not benchmark admission",
        "control_geometry": control_geometry, "field_consistency": consistency,
        "source_skeleton": str(args.skeleton) if args.skeleton else None,
        "source_cells": source_row["cells"] if source_row else 12,
        "recovery_worst_relative_energy_norm": float(np.sqrt(eigvalsh(estimated.sum(axis=0), tensors["linear"]).max())),
        "enrichment_worst_relative_energy_norm": float(np.sqrt(eigvalsh(observed.sum(axis=0), tensors["linear"]).max())),
        "target_length_quantiles": np.quantile(nodal_target, [0, 0.1, 0.5, 0.9, 1]).tolist(),
        "stiffness": {route: tensor.tolist() for route, tensor in tensors.items()},
        "solver_records": records, "workflow_s": time.perf_counter() - started,
        "scope": "Fixed-domain polynomial enrichment only; separately reported control geometry is not a continuum stress reference or admitted workflow ranking."})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-microgen-defaults-v4")
    parser.add_argument("--case")
    parser.add_argument("--skeleton", type=Path)
    args = parser.parse_args()
    if args.case:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    cases = (json.loads((args.skeleton / "result.json").read_text())["case"],) if args.skeleton else ("gyroid", "schwarz-solid")
    if args.skeleton:
        inputs = [args.skeleton / name for name in ("result.json", "mesh.npz", "manifest.json",
            cases[0] + "-samples.npz", cases[0] + "-specimen.json")]
        for path in inputs[-2:]:
            shutil.copyfile(path, args.run / path.name)
    else:
        inputs = list(args.source.glob("*-control-*-12.*"))
    inputs += [Path(workflow.__file__), workflow.HERE / "tet10_integration.py"]
    hashes = {str(path): workflow.digest(path) for path in inputs}
    baseline = (args.skeleton or args.source) / "manifest.json"
    workflow.save(args.run / "manifest.json", {"source_sha256": workflow.digest(Path(__file__)),
        "inputs": hashes, "memory_bytes": workflow.MEMORY, "timeout_s": 600, "threads": 1,
        "baseline_manifest_sha256": workflow.digest(baseline),
        "source_skeleton": str(args.skeleton) if args.skeleton else None,
        "selector": {"recovery_capture_beats": ["volume", "stress-magnitude"], "density_spearman_minimum": 0.5}})
    workers = {}
    with Job(parts=list(cases)) as job:
        for case in cases:
            with job.part(case) as part:
                command = [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source), "--case", case]
                if args.skeleton:
                    command.extend(["--skeleton", str(args.skeleton)])
                workers[case] = workflow.bounded(command, args.run / (case + ".log"), part, timeout=600)
                workflow.save(args.run / "workers.json", workers)
                assert workers[case]["returncode"] == 0, f"Read {case}.log"
                row = json.loads((args.run / (case + ".json")).read_text())
                job.result(case, recovery_capture=row["top_20_percent_elements_enrichment_energy_capture"]["recovery"],
                           density_correlation=row["density_spearman"], seconds=row["workflow_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    workflow.save(args.run / "summary.json", {"inputs_unchanged": True,
        "cases": [json.loads((args.run / (case + ".json")).read_text()) for case in cases]})


if __name__ == "__main__":
    main()
