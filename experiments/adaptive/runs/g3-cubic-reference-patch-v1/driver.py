"""Reference-only cubic enrichment feasibility; not a Meshers Tet20 API."""

import argparse
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sys

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
        workflow.HERE / "tet10_integration.py", workflow.HERE / "quadrature_audit.py"]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"sources": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 300, "threads": 1, "solve_node_cap": 10000, "tet_cap": 2000})
    with Job(parts=["basis-and-solid-patch"]) as job:
        with job.part("basis-and-solid-patch") as part:
            command = [sys.executable, __file__, "--run", str(args.run), "--worker"]
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("basis-and-solid-patch", seconds=record["wall_s"],
                analytic_patch_passed=int(row["analytic_patch_passed"]))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
