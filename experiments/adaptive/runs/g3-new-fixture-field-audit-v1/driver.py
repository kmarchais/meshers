"""Bounded validation of the public direct-fedoo conversion on G3."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))

from pueue_progress import Job
import microgen_workflow as workflow


def integrated_tensors(points, tetrahedra, displacements):
    np, fd, *_ = workflow.imports()
    order = 2 if tetrahedra.shape[1] == 10 else 1
    element = fd.lib_elements.element_list.get_element("tet10" if order == 2 else "tet4")(15 if order == 2 else 4)
    coordinates = points[tetrahedra]
    stiffness, energy = np.zeros((6, 6)), np.zeros((6, 6))
    volume = 0.0
    for derivative, weight in zip(element.shape_function_derivative(element.xi_pg), element.w_pg, strict=True):
        jacobian = np.einsum("an,tnb->tab", derivative, coordinates)
        determinants = np.linalg.det(jacobian)
        assert determinants.min() > 0
        weights = weight * determinants
        gradients = np.einsum("tba,an->tbn", np.linalg.inv(jacobian), derivative)
        displacement_gradient = np.einsum("ltni,tjn->ltij", displacements[:, tetrahedra], gradients)
        strain = (displacement_gradient + displacement_gradient.swapaxes(-1, -2)) / 2
        stress = 2 * strain + np.trace(strain, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)
        mean = np.einsum("t,ltij->lij", weights, stress)
        stiffness += mean[:, [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]].T
        energy += np.einsum("t,itkl,jtkl->ij", weights, strain, stress)
        volume += weights.sum()
    return stiffness, energy, float(volume)


def audit_saved(source, output):
    np, *_ = workflow.imports()
    from homogenization import strain_tensor
    from periodic import periodic_pairs, tensor_difference

    row = json.loads((source / "result.json").read_text())
    assert row["qualified"] and len(row["solver_records"]) == 6
    mesh = np.load(source / "mesh.npz")
    fields = np.load(source / "fields.npz")["displacement"]
    points, tetrahedra = mesh["points"], mesh["tetrahedra"]
    assert fields.shape == (6, len(points), 3) and np.isfinite(fields).all()
    stiffness, energy, volume = integrated_tensors(points, tetrahedra, fields)
    expected = np.asarray(row["stiffness"])
    reaction_error = float(np.linalg.norm(stiffness - expected) / np.linalg.norm(stiffness))
    energy_error = float(np.linalg.norm(stiffness - energy) / np.linalg.norm(stiffness))
    skew = float(np.linalg.norm(stiffness - stiffness.T) / np.linalg.norm(stiffness))
    assert max(reaction_error, energy_error, skew) < 1e-7
    disagreement = tensor_difference((stiffness + stiffness.T) / 2, expected)
    assert disagreement["worst_relative_strain_energy_difference"] < 1e-7
    assert abs(volume / row["geometry"]["volume"] - 1) < 1e-10
    pairs, _, _ = periodic_pairs(points, mesh["surface"])
    mids = np.ones(len(points), dtype=bool)
    mids[np.unique(tetrahedra[:, :4])] = False
    jump_max, mid_jump_max = 0.0, 0.0
    for load, engineering in enumerate(np.eye(6)):
        for axis, matching in enumerate(pairs):
            jump = fields[load, matching[:, 1]] - fields[load, matching[:, 0]]
            error = abs(jump - strain_tensor(engineering)[:, axis])
            jump_max = max(jump_max, float(error.max(initial=0)))
            mid_jump_max = max(mid_jump_max, float(error[mids[matching[:, 0]]].max(initial=0)))
    assert max(jump_max, mid_jump_max) < 1e-9
    workflow.save(output, {"source": str(source), "stress_integrated_stiffness": stiffness.tolist(),
        "energy_integrated_stiffness": energy.tolist(), "reaction_relative_error": reaction_error,
        "hill_mandel_relative_error": energy_error, "integrated_stiffness_skew_relative": skew,
        "stress_reaction_mixed_load_difference": disagreement, "volume": volume,
        "periodic_jump_max": jump_max, "periodic_mid_edge_jump_max": mid_jump_max,
        "scope": "Saved-field FE consistency; no continuum or local stress accuracy certification"})


def audit_saved_runs(run, sources):
    run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, run / "driver.py")
    inputs = {str(source / name): workflow.digest(source / name) for source in sources
              for name in ("result.json", "mesh.npz", "fields.npz")}
    workflow.save(run / "manifest.json", {"source_sha256": workflow.digest(Path(__file__)),
        "inputs": inputs, "memory_bytes": workflow.MEMORY, "threads": 1, "worker_timeout_s": 300})
    with Job(parts=[source.name for source in sources]) as job:
        for source in sources:
            with job.part(source.name) as part:
                result = workflow.bounded([sys.executable, __file__, "--saved-worker", str(source),
                    "--run", str(run)], run / (source.name + ".log"), part, timeout=300)
                workflow.save(run / (source.name + "-worker.json"), result)
                assert result["returncode"] == 0, f"Read {source.name}.log"
                row = json.loads((run / (source.name + ".json")).read_text())
                job.result(source.name, reaction_error=row["reaction_relative_error"],
                    jump=row["periodic_jump_max"], seconds=result["wall_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in inputs.items())


def check_case(name, source, output):
    import numpy as np
    import fedoo as fd
    import meshers
    import microgen
    import pyvista as pv
    from microgen.single_mesh import NotOnlyLinearTetrahedraError

    workflow.imports()
    from homogenization import AuditedCG, COMPONENTS, MATERIAL, strain_tensor
    from periodic import periodic_pairs, tensor_difference

    started = time.perf_counter()
    if name.startswith("cube-"):
        mesh = meshers.generate(lambda x, y, z: x * 0 - 1, bounds=workflow.BOUNDS,
                                cells=4, order=int(name[-1]), threads=1,
                                periodic=(True, True, True))
        expected = MATERIAL
    else:
        label = name + "-meshers-curved-12"
        saved = json.loads((source / (label + ".json")).read_text())
        assert saved["qualified"]
        data = np.load(source / (label + ".npz"))
        mesh = meshers.Mesh(points=data["points"], tetrahedra=data["tetrahedra"],
                            surface=data["surface"], boundary_tags=np.empty(0, dtype=np.uint8),
                            periodic_pairs=(), diagnostics={})
        expected = np.asarray(saved["stiffness"])
    assert len(mesh.points) <= 100_000 and len(mesh.tetrahedra) <= 250_000
    original_points = mesh.points.copy()
    original_tets = mesh.tetrahedra.copy()
    pairs, classes, pairing = periodic_pairs(mesh.points, mesh.surface)
    if mesh.periodic_pairs:
        for native, independent in zip(mesh.periodic_pairs, pairs, strict=True):
            assert set(map(tuple, native)) == set(map(tuple, independent))
    fem = mesh.to_fedoo()
    assert np.array_equal(fem.nodes, mesh.points)
    assert np.array_equal(fem.elements, mesh.tetrahedra)
    assert not np.shares_memory(fem.nodes, mesh.points)
    assert not np.shares_memory(fem.elements, mesh.tetrahedra)
    order = 2 if mesh.tetrahedra.shape[1] == 10 else 1
    element = fd.lib_elements.element_list.get_element(fem.elm_type)(15 if order == 2 else 4)
    assert np.allclose(element.shape_function(element.xi_nd), np.eye(mesh.tetrahedra.shape[1]))
    cell_type = pv.CellType.QUADRATIC_TETRA if order == 2 else pv.CellType.TETRA
    grid = pv.UnstructuredGrid({cell_type: mesh.tetrahedra}, mesh.points)
    cell = grid.GetCell(0)
    interpolation_error = 0.0
    for coordinate in element.xi_pg:
        values = [0.0] * mesh.tetrahedra.shape[1]
        vtk_coordinate = (coordinate[2], 1 - coordinate.sum(), coordinate[0])
        cell.InterpolateFunctions(vtk_coordinate, values)
        interpolation_error = max(interpolation_error, float(np.max(abs(
            np.asarray(values) - element.shape_function(coordinate[None])[0]))))
    assert interpolation_error < 1e-14
    conversion = {}
    if order == 2:
        for label, converter in (("SingleMesh", microgen.SingleMesh), ("BoxMesh", microgen.BoxMesh)):
            try:
                converter.from_pyvista(grid)
            except NotOnlyLinearTetrahedraError as error:
                conversion[label] = {"supported": False, "error": str(error)}
            else:
                raise AssertionError("Microgen quadratic conversion behavior changed")
    corner_nodes = np.unique(mesh.tetrahedra[:, :4])
    mids = np.setdiff1d(np.unique(mesh.tetrahedra), corner_nodes)
    if order == 2:
        edges = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))
        assert np.allclose(element.xi_nd[4:], element.xi_nd[edges].mean(axis=1))
        edge_nodes = np.column_stack((np.sort(mesh.tetrahedra[:, edges].reshape(-1, 2), axis=1),
                                      mesh.tetrahedra[:, 4:].ravel()))
        unique = np.unique(edge_nodes, axis=0)
        assert len(unique) == len(np.unique(unique[:, :2], axis=0)) == len(np.unique(unique[:, 2]))
        assert np.array_equal(np.sort(unique[:, 2]), mids)
    coordinates = mesh.points[mesh.tetrahedra]
    minimum_jacobian = min(float(np.linalg.det(np.einsum("an,tnb->tab", derivative, coordinates)).min())
        for derivative in element.shape_function_derivative(
            np.concatenate((element.xi_pg, element.xi_nd, [[0.25] * 3]))))
    assert minimum_jacobian > 0
    fd.ModelingSpace("3D")
    assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(
        fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), fem, n_elm_gp=15 if order == 2 else 4)
    solver = AuditedCG()
    displacements, checks = [], []
    reactions = np.zeros((6, 6))
    center = int(np.linalg.norm(mesh.points, axis=1).argmin())
    for load, engineering in enumerate(np.eye(6)):
        problem = fd.problem.Linear(assembly)
        problem.bc.add(fd.constraint.PeriodicBC("small_strain", meshperio=True))
        problem.bc.add("Dirichlet", [center], "Disp", 0)
        for component, value in zip(COMPONENTS, engineering, strict=True):
            problem.bc.add("Dirichlet", component, float(value))
        problem.set_solver(solver)
        problem.solve()
        displacement = problem.get_disp().T
        assert len(displacement) == len(mesh.points)
        assert np.all(np.isfinite(displacement))
        macro = strain_tensor(engineering)
        jumps, mid_jumps = [], []
        for axis_pairs in pairs:
            error = abs(displacement[axis_pairs[:, 1]] - displacement[axis_pairs[:, 0]]
                        - (mesh.points[axis_pairs[:, 1]] - mesh.points[axis_pairs[:, 0]]) @ macro.T)
            jumps.append(float(error.max()))
            mid_jumps.append(float(error[np.isin(axis_pairs[:, 0], mids)].max(initial=0)))
        assert max(jumps) < 1e-10
        forces = problem.get_ext_forces("Disp", include_mpc=False).T
        balanced = np.zeros_like(forces)
        np.add.at(balanced, classes, forces)
        balance = float(np.linalg.norm(balanced) / np.linalg.norm(forces))
        assert balance < 1e-7
        reactions[:, load] = problem.get_ext_forces("MeanStrain").ravel()
        displacements.append(displacement)
        checks.append({"load": COMPONENTS[load], "jump_max": max(jumps),
                       "mid_edge_jump_max": max(mid_jumps), "force_balance_relative": balance,
                       **solver.records[-1]})
    displacements = np.asarray(displacements)
    stiffness, energy, volume = integrated_tensors(mesh.points, mesh.tetrahedra, displacements)
    reaction_error = float(np.linalg.norm(stiffness - reactions) / np.linalg.norm(stiffness))
    energy_error = float(np.linalg.norm(stiffness - energy) / np.linalg.norm(stiffness))
    assert max(reaction_error, energy_error) < 1e-7
    disagreement = tensor_difference(stiffness, expected)
    assert disagreement["worst_relative_strain_energy_difference"] < 1e-7
    assert np.array_equal(mesh.points, original_points) and np.array_equal(mesh.tetrahedra, original_tets)
    np.savez_compressed(output / (name + "-fields.npz"), displacement=displacements,
                        points=mesh.points, tetrahedra=mesh.tetrahedra)
    workflow.save(output / (name + ".json"), {
        "case": name, "nodes": len(mesh.points), "elements": len(mesh.tetrahedra),
        "order": order, "volume": float(volume), "minimum_sampled_jacobian": minimum_jacobian,
        "vtk_fedoo_interpolation_error": interpolation_error, "periodic": pairing,
        "periodic_mid_edge_pair_counts": [int(np.isin(pair[:, 0], mids).sum()) for pair in pairs],
        "microgen_from_pyvista": conversion, "solver_checks": checks,
        "stiffness": stiffness.tolist(), "reference_disagreement": disagreement,
        "reaction_relative_error": reaction_error, "hill_mandel_relative_error": energy_error,
        "source_mesh_unchanged": True, "workflow_s": time.perf_counter() - started,
        "scope": "Conversion and periodic FE integration; no continuum or local stress accuracy claim."})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-microgen-defaults-v4")
    parser.add_argument("--case")
    parser.add_argument("--saved-fields", nargs="+", type=Path)
    parser.add_argument("--saved-worker", type=Path)
    args = parser.parse_args()
    if args.saved_worker:
        audit_saved(args.saved_worker, args.run / (args.saved_worker.name + ".json"))
        return
    if args.saved_fields:
        audit_saved_runs(args.run, args.saved_fields)
        return
    if args.case:
        check_case(args.case, args.source, args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    sources = [Path(__file__), ROOT / "crates/meshers-python/python/meshers/__init__.py",
               ROOT / "crates/meshers-python/tests/test_fedoo.py",
               ROOT.parent / "periodic-microgen/experiments/periodic-fea/homogenization.py",
               ROOT.parent / "periodic-microgen/experiments/periodic-fea/periodic.py"]
    frozen = {str(path): workflow.digest(path) for path in args.source.iterdir() if path.is_file()}
    workflow.save(args.run / "manifest.json", {"source_sha256": {str(path): workflow.digest(path) for path in sources},
        "input_sha256": frozen, "packages": {name: importlib.metadata.version(name) for name in
        ("fedoo", "microgen", "numpy", "scipy", "pyvista", "vtk", "meshers")},
        "worker_memory_bytes": workflow.MEMORY, "worker_timeout_s": 600, "threads": 1})
    cases = ["cube-1", "cube-2", "gyroid", "schwarz-solid"]
    workers = {}
    with Job(parts=["conversion-tests"] + cases) as job:
        for name in ["conversion-tests"] + cases:
            with job.part(name) as part:
                command = ([sys.executable, "-m", "pytest", "-q", "crates/meshers-python/tests/test_fedoo.py"]
                    if name == "conversion-tests" else [sys.executable, __file__, "--case", name,
                    "--source", str(args.source), "--run", str(args.run)])
                workers[name] = workflow.bounded(command, args.run / (name + ".log"), part, timeout=600)
                workflow.save(args.run / "workers.json", workers)
                assert workers[name]["returncode"] == 0, f"Read {args.run / (name + '.log')}"
                job.result(name, peak_gib=workers[name]["peak_rss_bytes"] / 1024**3,
                           seconds=workers[name]["wall_s"])
    assert frozen == {str(path): workflow.digest(path) for path in args.source.iterdir() if path.is_file()}
    workflow.save(args.run / "summary.json", {"frozen_inputs_unchanged": True,
        "cases": [json.loads((args.run / (name + ".json")).read_text()) for name in cases]})


if __name__ == "__main__":
    main()
