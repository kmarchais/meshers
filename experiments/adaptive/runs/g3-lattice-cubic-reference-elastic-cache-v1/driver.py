"""Reference-only cubic enrichment feasibility; not a Meshers Tet20 API."""

import argparse
import itertools
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def register_cubic():
    import numpy as np
    from fedoo.lib_elements.tetrahedron import ElementTetrahedron, Tet4
    from fedoo.lib_elements.element_list import add_elements

    indices = [tuple(3 if axis == node else 0 for axis in range(4)) for node in range(4)]
    for first, second in itertools.combinations(range(4), 2):
        for weights in ((2, 1), (1, 2)):
            indices.append(tuple(weights[0] if axis == first else weights[1] if axis == second else 0
                for axis in range(4)))
    indices.extend(tuple(0 if axis == missing else 1 for axis in range(4)) for missing in range(4))

    class CubicTetrahedron(ElementTetrahedron):
        name = "tet20_reference"
        default_n_gp = 15
        n_nodes = 20
        multi_indices = np.array(indices)

        def __init__(self, n_elm_gp=15, **kwargs):
            self.xi_nd = self.multi_indices @ Tet4().xi_nd / 3
            self.n_elm_gp = n_elm_gp
            super().__init__(n_elm_gp)

        def factors(self, coordinates):
            barycentric = Tet4().shape_function(np.asarray(coordinates))
            factors, derivatives = [], []
            for index in self.multi_indices:
                values, gradients = [], []
                for axis, degree in enumerate(index):
                    terms = [3 * barycentric[:, axis] - shift for shift in range(degree)]
                    values.append(np.prod(terms, axis=0) / math.factorial(degree) if terms else np.ones(len(coordinates)))
                    gradients.append(sum((3 * np.prod([term for position, term in enumerate(terms)
                        if position != omitted], axis=0) for omitted in range(degree)),
                        np.zeros(len(coordinates))) / math.factorial(degree))
                factors.append(np.array(values).T)
                derivatives.append(np.array(gradients).T)
            return np.array(factors), np.array(derivatives)

        def shape_function(self, coordinates):
            factors, _ = self.factors(coordinates)
            return np.prod(factors, axis=2).T

        def shape_function_derivative(self, coordinates):
            factors, derivatives = self.factors(coordinates)
            barycentric_derivative = Tet4().shape_function_derivative([[.25] * 3])[0]
            gradient = np.zeros((len(coordinates), 3, 20))
            for axis in range(4):
                partial = derivatives[:, :, axis] * np.prod(np.delete(factors, axis, axis=2), axis=2)
                gradient += np.einsum("np,d->pdn", partial, barycentric_derivative[:, axis])
            return gradient

    add_elements(CubicTetrahedron)
    return CubicTetrahedron


def enrich(points, tetrahedra, surface, element):
    import numpy as np
    from fedoo.lib_elements.tetrahedron import Tet10

    assert tetrahedra.shape[1] == 10 and surface.shape[1] == 6
    values = Tet10().shape_function(element.xi_nd)
    ownership, nodes = {}, []
    enriched = np.empty((len(tetrahedra), 20), dtype=np.int64)
    for node, index in enumerate(element.multi_indices):
        locations = np.einsum("n,tnd->td", values[node], points[tetrahedra])
        for cell, corners in enumerate(tetrahedra[:, :4]):
            key = tuple(sorted((int(corner), int(weight)) for corner, weight in zip(corners, index, strict=True) if weight))
            if key not in ownership:
                ownership[key] = len(nodes)
                nodes.append(locations[cell])
            else:
                np.testing.assert_allclose(nodes[ownership[key]], locations[cell], atol=1e-13, rtol=0)
            enriched[cell, node] = ownership[key]
    facets = []
    face_indices = np.array([index for index in element.multi_indices if index[3] == 0])[:, :3]
    assert face_indices.shape == (10, 3)
    for face in surface[:, :3]:
        facets.append([ownership[tuple(sorted((int(corner), int(weight))
            for corner, weight in zip(face, index, strict=True) if weight))] for index in face_indices])
    return np.array(nodes), enriched, np.array(facets)


def basis_checks(element):
    import numpy as np
    from fedoo.lib_elements.tetrahedron import Tet10

    np.testing.assert_allclose(element.shape_function(element.xi_nd), np.eye(20), atol=2e-14, rtol=0)
    coordinates = np.random.default_rng(20261018).dirichlet(np.ones(4), 64)[:, :3]
    values = element.shape_function(coordinates)
    derivatives = element.shape_function_derivative(coordinates)
    np.testing.assert_allclose(values.sum(1), 1, atol=2e-14, rtol=0)
    np.testing.assert_allclose(derivatives.sum(2), 0, atol=2e-14, rtol=0)
    for degrees in itertools.product(range(4), repeat=3):
        if sum(degrees) > 3:
            continue
        nodal = np.prod(element.xi_nd ** np.array(degrees), axis=1)
        expected = np.prod(coordinates ** np.array(degrees), axis=1)
        np.testing.assert_allclose(values @ nodal, expected, atol=2e-14, rtol=0)
        for axis in range(3):
            shifted = list(degrees)
            shifted[axis] = max(0, shifted[axis] - 1)
            expected_derivative = degrees[axis] * np.prod(coordinates ** np.array(shifted), axis=1)
            np.testing.assert_allclose(derivatives[:, axis] @ nodal, expected_derivative, atol=5e-14, rtol=0)
    for degrees in itertools.product(range(5), repeat=3):
        if sum(degrees) <= 4:
            actual = element.w_pg @ np.prod(element.xi_pg ** np.array(degrees), axis=1)
            expected = math.prod(math.factorial(degree) for degree in degrees) / math.factorial(sum(degrees) + 3)
            assert abs(actual - expected) < 1e-14
    parent = Tet10()
    curved = parent.xi_nd + np.random.default_rng(17).normal(size=(10, 3)) * .02
    cubic_nodes = parent.shape_function(element.xi_nd) @ curved
    error = abs(values @ cubic_nodes - parent.shape_function(coordinates) @ curved).max()
    derivative_error = abs(np.einsum("pan,nd->pad", derivatives, cubic_nodes) -
        np.einsum("pan,nd->pad", parent.shape_function_derivative(coordinates), curved)).max()
    assert max(error, derivative_error) < 1e-12
    return {"nodal_and_polynomial_identities": True, "degree_four_quadrature_exact": True,
        "quadratic_map_embedding_error": float(error), "map_derivative_embedding_error": float(derivative_error)}


def worker(run):
    np, _, meshers, *_ = workflow.imports()
    from stress_recovery import solve, periodic_displacement_jump_error
    from tet10_integration import integrated_tensors
    from quadrature_audit import duffy_rule
    from homogenization import MATERIAL, strain_tensor
    from periodic import periodic_pairs, tensor_difference

    element_class = register_cubic()
    element = element_class()
    checks = basis_checks(element)
    workflow.save(run / "basis.json", checks)
    parent = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
        bounds=(-.5, .5) * 3, cells=4, order=2, geometry_tolerance=1,
        minimum_quality=.1, periodic=(True, True, True), threads=1)
    points, tets, surface = enrich(parent.points, parent.tetrahedra, parent.surface, element)
    assert len(points) <= 10_000 and len(tets) <= 2_000
    _, _, periodic = periodic_pairs(points, surface)
    rows = []
    for name, locations, cells, facets, element_type in (
        ("tet10-control", parent.points, parent.tetrahedra, parent.surface, None),
        ("tet20-reference", points, tets, surface, element_class.name)):
        fields, tensor, records = solve(locations, cells, n_elm_gp=15, element_type=element_type)
        stress, energy, volume = integrated_tensors(locations, cells, fields,
            duffy_rule(5), element_type=element_type)
        jump = periodic_displacement_jump_error(locations, facets, fields)
        affine_error = max(float(abs(fields[load] - locations @ strain_tensor(engineering).T).max())
            for load, engineering in enumerate(np.eye(6)))
        difference = tensor_difference(tensor, MATERIAL)
        assert max(difference["worst_relative_strain_energy_difference"],
            np.linalg.norm(stress - tensor), np.linalg.norm(energy - tensor), abs(volume - 1),
            affine_error, jump) < 1e-8
        rows.append({"route": name, "nodes": len(locations), "elements": len(cells),
            "difference": difference, "higher_rule_stress_error": float(np.linalg.norm(stress - tensor)),
            "higher_rule_energy_error": float(np.linalg.norm(energy - tensor)), "volume": volume,
            "periodic_jump_max": jump, "affine_field_error": affine_error, "solver_records": records})
    workflow.save(run / "result.json", {"basis_checks": checks, "patches": rows,
        "cubic_periodic": periodic, "analytic_patch_passed": True,
        "scope": "Reference element feasibility only; no lattice tensor or Meshers Tet20 support inferred"})


def assembly_worker(run, source):
    np, fd, *_ = workflow.imports()
    from periodic import periodic_pairs

    parent = json.loads((source / "result.json").read_text())
    assert parent["case"] == "smooth-lattice" and parent["qualified"]
    assert parent["order"] == 2 and parent["geometry"]["geometry_admitted"]
    archive = np.load(source / "mesh.npz")
    element_class = register_cubic()
    element = element_class()
    basis_checks(element)
    points, cells, surface = enrich(archive["points"], archive["tetrahedra"], archive["surface"], element)
    assert len(points) <= 250_000 and len(cells) <= workflow.CAP
    _, _, periodic = periodic_pairs(points, surface)
    np.savez_compressed(run / "mesh.npz", points=points, tetrahedra=cells, surface=surface)
    workflow.save(run / "prepared.json", {"nodes": len(points), "elements": len(cells),
        "periodic": periodic, "geometry_source": str(source),
        "geometry_mapping": "Cubic basis exactly embeds source quadratic mapping"})
    fd.ModelingSpace("3D")
    mesh = fd.Mesh(points, cells, element_class.name)
    assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(
        fd.constitutivelaw.ElasticIsotrop(2.5, .25)), mesh, n_elm_gp=15)
    problem = fd.problem.Linear(assembly)
    matrix = assembly.get_global_matrix()
    assert problem.mesh is mesh and matrix.shape == (3 * len(points),) * 2
    assert np.isfinite(matrix.data).all() and (matrix.diagonal() > 0).all()
    workflow.save(run / "result.json", {"assembly_only": True, "nodes": len(points),
        "elements": len(cells), "matrix_shape": matrix.shape, "matrix_nonzeros": matrix.nnz,
        "matrix_csr_bytes": matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes,
        "no_fea_solve": True,
        "scope": "Resource feasibility only; no cubic lattice tensor or accuracy claim"})


def reference_worker(run, source, prepared, use_amg=False, elastic_blocks=False):
    np, *_ = workflow.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    from stress_recovery import solve, periodic_displacement_jump_error
    from tet10_integration import integrated_tensors
    from quadrature_audit import duffy_rule
    from periodic import tensor_difference, periodic_pairs

    parent = json.loads((source / "result.json").read_text())
    assert parent["qualified"] and parent["geometry"]["geometry_admitted"]
    assert parent["case"] == "smooth-lattice" and len(parent["solver_records"]) == 6
    mesh = np.load(prepared / "mesh.npz")
    original = np.load(source / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 20 and len(points) <= 250_000 and len(cells) <= workflow.CAP
    element_class = register_cubic()
    element = element_class()
    basis_checks(element)
    _, _, periodic = periodic_pairs(points, surface)
    interpolation = Tet10().shape_function(element.xi_nd)
    original_fields = np.load(source / "fields.npz")["displacement"][:, original["tetrahedra"]]
    embedded = np.zeros((6, len(points), 3))
    assigned = np.zeros(len(points), dtype=bool)
    embedding_error = 0.0
    for node, weights in enumerate(interpolation):
        destinations = cells[:, node]
        expected = np.einsum("n,tnd->td", weights, original["points"][original["tetrahedra"]])
        embedding_error = max(embedding_error, float(abs(points[destinations] - expected).max()))
        values = np.einsum("n,ltnc->ltc", weights, original_fields)
        known = assigned[destinations]
        np.testing.assert_allclose(embedded[:, destinations[known]], values[:, known], atol=1e-12, rtol=0)
        embedded[:, destinations] = values
        assigned[destinations] = True
    assert assigned.all() and embedding_error < 1e-12
    parent_tensor = np.asarray(parent["stiffness"])
    embedded_stress, embedded_energy, volume = integrated_tensors(points, cells, embedded,
        element_type=element_class.name)
    embedding_tensor_error = max(float(np.linalg.norm(tensor - parent_tensor) / np.linalg.norm(parent_tensor))
        for tensor in (embedded_stress, embedded_energy))
    assert embedding_tensor_error < 1e-7
    assert abs(volume / parent["geometry"]["volume"] - 1) < 1e-10
    assert periodic_displacement_jump_error(points, surface, embedded) < 1e-9
    del original_fields, embedded
    row = {"case": parent["case"], "order": 3, "element_type": element_class.name,
        "nodes": len(points), "elements": len(cells), "geometry": parent["geometry"],
        "geometry_source": str(source), "exact_quadratic_mapping_embedding_error": embedding_error,
        "embedded_parent_tensor_error": embedding_tensor_error, "periodic": periodic,
        "solver_quadrature_points": 15, "qualified": False, "status": "embedded-parent-verified",
        "scope": "Saved-mesh reference enrichment, not a live workflow timing or continuum certificate",
        "first_load_wall_cap_s": 240,
        "solver": ("cached-elastic-block-amg-cg" if elastic_blocks else
            "cached-smoothed-aggregation-cg" if use_amg else "jacobi-cg")}
    workflow.save(run / "result.json", row)
    phase = time.perf_counter()

    def first_load_timeout(signum, frame):
        row["status"] = "first-load-timeout"
        workflow.save(run / "result.json", row)
        raise TimeoutError("First cubic reference load exceeds 240-second stage")

    def checkpoint(completed, partial, records):
        row.update(status="solving", completed_loads=completed,
            partial_stiffness_columns=partial.tolist(), solver_records=records)
        if completed == 1:
            signal.setitimer(signal.ITIMER_REAL, 0)
        workflow.save(run / "result.json", row)

    def save_completed_field(completed, displacement):
        np.savez_compressed(run / f"load-{completed}.npz", displacement=displacement)

    solver = None
    if use_amg:
        from amg_reference import AuditedAMG
        solver = AuditedAMG(elastic_blocks=elastic_blocks)
    previous_alarm = signal.signal(signal.SIGALRM, first_load_timeout)
    signal.setitimer(signal.ITIMER_REAL, row["first_load_wall_cap_s"])
    try:
        fields, tensor, records = solve(points, cells, checkpoint=checkpoint, n_elm_gp=15,
            element_type=element_class.name, solver=solver, field_checkpoint=save_completed_field)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_alarm)
    row["solve_s"] = time.perf_counter() - phase
    np.savez_compressed(run / "fields.npz", displacement=fields)
    row.pop("partial_stiffness_columns", None)
    row.update(status="auditing", stiffness=tensor.tolist(), solver_records=records)
    workflow.save(run / "result.json", row)
    stress, energy, volume = integrated_tensors(points, cells, fields, element_type=element_class.name)
    own_error = max(float(np.linalg.norm(value - tensor) / np.linalg.norm(tensor)) for value in (stress, energy))
    jump = periodic_displacement_jump_error(points, surface, fields)
    assert own_error < 1e-7 and jump < 1e-9
    _, higher_energy, _ = integrated_tensors(points, cells, fields, duffy_rule(5),
        element_type=element_class.name)
    sensitivity = tensor_difference(higher_energy, tensor)
    minimum_enrichment = float(np.linalg.eigvalsh(parent_tensor - tensor).min())
    assert minimum_enrichment >= -1e-8 * np.linalg.norm(parent_tensor)
    row.pop("partial_stiffness_columns", None)
    passed = sensitivity["worst_relative_strain_energy_difference"] <= 1e-4
    row.update(status="solved-needs-reference" if passed else "quadrature-sensitivity-rejected",
        qualified=passed, stiffness=tensor.tolist(), solver_records=records,
        own_rule_relative_error=own_error, periodic_jump_max=jump,
        higher_rule_fixed_field_difference=sensitivity, quadrature_sensitivity_limit=1e-4,
        quadrature_sensitivity_passed=passed, minimum_parent_minus_cubic_eigenvalue=minimum_enrichment,
        integrated_volume=volume)
    workflow.save(run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--amg", action="store_true")
    parser.add_argument("--elastic-blocks", action="store_true")
    args = parser.parse_args()
    assert args.prepared is None or args.source is not None
    assert not args.amg or args.prepared is not None
    assert not args.elastic_blocks or args.amg
    if args.worker:
        if args.source is None:
            worker(args.run)
        elif args.prepared is None:
            assembly_worker(args.run, args.source)
        else:
            reference_worker(args.run, args.source, args.prepared, args.amg, args.elastic_blocks)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
        workflow.HERE / "tet10_integration.py", workflow.HERE / "quadrature_audit.py"]
    if args.source is not None:
        paths.extend(args.source / name for name in ("result.json", "mesh.npz"))
    if args.prepared is not None:
        paths.extend((args.prepared / "mesh.npz", args.prepared / "result.json", args.source / "fields.npz"))
    if args.amg:
        paths.append(workflow.HERE / "amg_reference.py")
    hashes = {str(path): workflow.digest(path) for path in paths}
    timeout = 1800 if args.prepared is not None else 600 if args.source is not None else 300
    workflow.save(args.run / "manifest.json", {"sources": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": timeout, "threads": 1,
        "node_cap": 250000 if args.source is not None else 10000,
        "tet_cap": workflow.CAP if args.source is not None else 2000,
        "packages": {name: importlib.metadata.version(name) for name in ("fedoo", "numpy", "scipy", "meshers")},
        "assembly_only": args.source is not None and args.prepared is None,
        "solver": ("cached-elastic-block-amg-cg" if args.elastic_blocks else
            "cached-smoothed-aggregation-cg" if args.amg else "jacobi-cg"),
        "elastic_blocks": args.elastic_blocks,
        "pyamg": importlib.metadata.version("pyamg") if args.amg else None,
        "first_load_wall_cap_s": 240 if args.prepared is not None else None})
    label = "cubic-reference-solve" if args.prepared is not None else "cubic-assembly-memory" if args.source is not None else "basis-and-solid-patch"
    with Job(parts=[label]) as job:
        with job.part(label) as part:
            command = [sys.executable, __file__, "--run", str(args.run), "--worker"]
            if args.source is not None:
                command.extend(["--source", str(args.source)])
            if args.prepared is not None:
                command.extend(["--prepared", str(args.prepared)])
            if args.amg:
                command.append("--amg")
            if args.elastic_blocks:
                command.append("--elastic-blocks")
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=timeout)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result(label, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3,
                nodes=row.get("nodes"), analytic_patch_passed=int(row.get("analytic_patch_passed", False)))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
