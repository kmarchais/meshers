"""Coupled periodic edge-mode trial energies on an unchanged cubic FE domain."""

import argparse
import importlib.metadata
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

EDGES = tuple(itertools.combinations(range(4), 2))
SOURCE = workflow.HERE / "runs/g3-lattice-cubic-reference-cache-release-v1"
PREPARED = workflow.HERE / "runs/g3-lattice-cubic-assembly-memory-v1"
EDGE_BASELINE = workflow.HERE / "runs/g3-lattice-quartic-edge-witness-higher-v2"


def edge_basis(barycentric):
    import numpy as np

    barycentric = np.asarray(barycentric)
    values = np.zeros((len(barycentric), 6))
    derivatives = np.zeros((len(barycentric), 6, 4))
    for edge, (first, second) in enumerate(EDGES):
        left, right = barycentric[:, first], barycentric[:, second]
        difference = left - right
        values[:, edge] = 16 * left * right * difference**2
        derivatives[:, edge, first] = 16 * right * (difference**2 + 2 * left * difference)
        derivatives[:, edge, second] = 16 * left * (difference**2 - 2 * right * difference)
    return values, derivatives


def mode_basis(barycentric, relax_corners=False):
    import numpy as np

    values, derivatives = edge_basis(barycentric)
    if not relax_corners:
        return values, derivatives
    barycentric = np.asarray(barycentric)
    corners = .5 * barycentric * (3 * barycentric - 1) * (3 * barycentric - 2)
    corner_derivatives = np.zeros((len(barycentric), 4, 4))
    for corner in range(4):
        corner_derivatives[:, corner, corner] = 13.5 * barycentric[:, corner]**2 - 9 * barycentric[:, corner] + 1
    return np.concatenate((values, corners), axis=1), np.concatenate((derivatives, corner_derivatives), axis=1)


def select_modes(cells, classes, physical_edges, relax_corners=False):
    import numpy as np

    physical_edges = {tuple(sorted(map(int, edge))) for edge in physical_edges}
    keys = sorted({tuple(sorted(map(int, classes[list(edge)]))) for edge in physical_edges})
    corner_keys = sorted({int(classes[node]) for edge in physical_edges for node in edge}) if relax_corners else []
    assert 0 < len(keys) + len(corner_keys) <= 512
    lookup = {key: index for index, key in enumerate(keys)}
    corner_lookup = {key: len(keys) + index for index, key in enumerate(corner_keys)}
    identifiers = np.full((len(cells), 10 if relax_corners else 6), -1, dtype=int)
    for cell, corners in enumerate(cells[:, :4]):
        for edge, (first, second) in enumerate(EDGES):
            pair = tuple(sorted((int(corners[first]), int(corners[second]))))
            key = tuple(sorted(map(int, classes[list(pair)])))
            if key in lookup:
                assert pair in physical_edges, "Unexpected winding edge in a periodic orbit"
                identifiers[cell, edge] = lookup[key]
        for corner, node in enumerate(corners):
            if int(classes[node]) in corner_lookup:
                identifiers[cell, 6 + corner] = corner_lookup[int(classes[node])]
    active = np.flatnonzero((identifiers >= 0).any(axis=1))
    assert 0 < len(active) <= 5000
    assert set(identifiers[identifiers >= 0]) == set(range(len(keys) + len(corner_keys)))
    return {"active": active, "identifiers": identifiers[active], "orbits": keys,
            "mode_count": len(keys) + len(corner_keys), "physical_edges": sorted(physical_edges),
            "edge_mode_count": len(keys), "corner_mode_count": len(corner_keys),
            "corner_orbits": corner_keys, "relax_corners": relax_corners}


def strain_basis(gradients):
    import numpy as np

    result = np.zeros(gradients.shape[:2] + (6, 3))
    result[:, :, 0, 0] = gradients[:, :, 0]
    result[:, :, 1, 1] = gradients[:, :, 1]
    result[:, :, 2, 2] = gradients[:, :, 2]
    result[:, :, 3, 0], result[:, :, 3, 1] = gradients[:, :, 1], gradients[:, :, 0]
    result[:, :, 4, 0], result[:, :, 4, 2] = gradients[:, :, 2], gradients[:, :, 0]
    result[:, :, 5, 1], result[:, :, 5, 2] = gradients[:, :, 2], gradients[:, :, 1]
    return result


def assemble(points, cells, fields, selection, element, quadrature):
    import numpy as np
    from fedoo.lib_elements.tetrahedron import Tet4
    from homogenization import MATERIAL
    from stress_recovery import stresses

    selected = cells[selection["active"]]
    identifiers = selection["identifiers"]
    size = 3 * selection["mode_count"]
    matrix, residual, cached = np.zeros((size, size)), np.zeros((size, 6)), []
    local_size = 3 * identifiers.shape[1]
    dofs = (3 * identifiers[:, :, None] + np.arange(3)).reshape(-1, local_size)
    valid = dofs >= 0
    pairs = valid[:, :, None] & valid[:, None, :]
    matrix_indices = (dofs[:, :, None] * size + dofs[:, None, :])[pairs]
    locations, weights = quadrature
    linear = Tet4()
    _, derivatives = mode_basis(linear.shape_function(locations), selection["relax_corners"])
    references = np.einsum("pea,da->ped", derivatives, linear.shape_function_derivative([[.25] * 3])[0])
    for derivative, reference, weight in zip(element.shape_function_derivative(locations), references, weights, strict=True):
        jacobian = np.einsum("an,tnb->tab", derivative, points[selected])
        determinants = np.linalg.det(jacobian)
        assert weight > 0 and determinants.min() > 0
        integration = weight * determinants
        gradients = np.einsum("tba,ea->teb", np.linalg.inv(jacobian), reference)
        gradients *= (identifiers >= 0)[:, :, None]
        strain = strain_basis(gradients)
        stress = stresses(points, selected, fields, derivative)
        voigt = np.stack([stress[:, :, first, second] for first, second in
                          ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))], axis=-1)
        local_matrix = np.einsum("teai,ab,tfbj,t->teifj", strain, MATERIAL, strain, integration).reshape(-1, local_size, local_size)
        local_residual = np.einsum("teai,tla,t->teil", strain, voigt, integration).reshape(-1, local_size, 6)
        np.add.at(matrix.ravel(), matrix_indices, local_matrix[pairs])
        np.add.at(residual, dofs[valid], local_residual[valid])
        cached.append((integration, gradients, stress))
    np.testing.assert_allclose(matrix, matrix.T, atol=1e-12, rtol=1e-12)
    np.linalg.cholesky(matrix)
    return matrix, residual, cached


def coupled_gains(matrix, residual, edge_modes):
    import numpy as np

    split = 3 * edge_modes
    assert 0 < split < len(matrix)
    edges, corners = matrix[:split, :split], matrix[split:, split:]
    cross = matrix[:split, split:]
    edge_residual, corner_residual = residual[:split], residual[split:]
    edge_gain = edge_residual.T @ np.linalg.solve(edges, edge_residual)
    corner_gain = corner_residual.T @ np.linalg.solve(corners, corner_residual)
    coefficients = np.linalg.solve(matrix, residual)
    joint_gain = residual.T @ coefficients
    incremental = joint_gain - corner_gain
    schur = edges - cross @ np.linalg.solve(corners, cross.T)
    np.testing.assert_allclose(schur, schur.T, atol=1e-12, rtol=1e-12)
    np.linalg.cholesky(schur)
    reduced = edge_residual - cross @ np.linalg.solve(corners, corner_residual)
    schur_gain = reduced.T @ np.linalg.solve(schur, reduced)
    np.testing.assert_allclose(incremental, schur_gain, atol=1e-12, rtol=1e-8)
    assert np.linalg.eigvalsh(incremental).min() >= -1e-12
    assert np.linalg.eigvalsh(joint_gain - edge_gain).min() >= -1e-12
    return {"coefficients": coefficients, "edge_gain": edge_gain, "corner_gain": corner_gain,
        "joint_gain": joint_gain, "incremental_gain": incremental,
        "schur_identity_relative_error": float(np.linalg.norm(incremental - schur_gain) / max(np.linalg.norm(incremental), 1e-30))}


def trial_energy(cached, selection, coefficients):
    import numpy as np
    from stress_recovery import gram

    identifiers = selection["identifiers"]
    local = coefficients.reshape(-1, 3, 6)[np.maximum(identifiers, 0)].copy()
    local *= (identifiers >= 0)[:, :, None, None]
    gain, mode_energy = np.zeros((6, 6)), np.zeros((6, 6))
    for weights, gradients, stress in cached:
        gradient = np.einsum("teil,tej->tlij", local, gradients)
        strain = (gradient + gradient.swapaxes(-1, -2)) / 2
        correction = 2 * strain + np.trace(strain, axis1=-2, axis2=-1)[:, :, None, None] * np.eye(3)
        gain += np.einsum("t,tlm->lm", weights, gram(stress) - gram(stress - correction))
        mode_energy += np.einsum("t,tlm->lm", weights, gram(correction))
    return gain, mode_energy


def gauge_at_node(cells, selection, element, coefficients, pin):
    import numpy as np
    from fedoo.lib_elements.tetrahedron import Tet4

    identifiers = selection["identifiers"]
    local = coefficients.reshape(-1, 3, 6)[np.maximum(identifiers, 0)].copy()
    local *= (identifiers >= 0)[:, :, None, None]
    values = mode_basis(Tet4().shape_function(element.xi_nd), selection["relax_corners"])[0]
    result = None
    for cell, node in np.argwhere(cells[selection["active"]] == pin):
        current = np.einsum("e,eil->li", values[node], local[cell])
        if result is not None:
            np.testing.assert_allclose(current, result, atol=1e-11, rtol=0)
        result = current
    return np.zeros((6, 3)) if result is None else result


def trace_audit(points, cells, surface, classes, selection, coefficients, element):
    import numpy as np

    triangles = np.array([[first / 4, second / 4, (4 - first - second) / 4]
        for first in range(5) for second in range(5 - first)])
    assert len(triangles) == 15
    identifiers = selection["identifiers"]
    local = coefficients.reshape(-1, 3, 6)[np.maximum(identifiers, 0)].copy()
    local *= (identifiers >= 0)[:, :, None, None]
    boundary = {tuple(sorted(map(int, face))) for face in surface[:, :3]}
    faces, caps = {}, {}
    continuity = 0.
    for index, cell in enumerate(cells[selection["active"], :4]):
        for face in itertools.combinations(range(4), 3):
            ordered = sorted(face, key=lambda position: int(cell[position]))
            barycentric = np.zeros((15, 4))
            barycentric[:, ordered] = triangles
            values = np.einsum("pe,eil->pli", mode_basis(barycentric, selection["relax_corners"])[0], local[index])
            key = tuple(sorted(map(int, cell[list(face)])))
            if key in faces:
                previous, count = faces[key]
                continuity = max(continuity, float(abs(values - previous).max()))
                faces[key] = previous, count + 1
            else:
                faces[key] = values, 1
            for axis in range(3):
                for side in (-.5, .5):
                    if np.all(workflow.on_box_plane(points[cell[list(face)], axis], side)):
                        ordered = sorted(face, key=lambda position: int(classes[cell[position]]))
                        canonical = tuple(map(int, classes[cell[ordered]]))
                        assert len(set(canonical)) == 3
                        barycentric[:] = 0
                        barycentric[:, ordered] = triangles
                        caps[axis, side, canonical] = np.einsum("pe,eil->pli", mode_basis(barycentric, selection["relax_corners"])[0], local[index])
    for key, (values, count) in faces.items():
        assert count <= 2
        if count == 1 and key not in boundary:
            continuity = max(continuity, float(abs(values).max()))
    jump, nonzero = 0., 0
    for axis, _, canonical in caps:
        low = caps.get((axis, -.5, canonical), np.zeros((15, 6, 3)))
        high = caps.get((axis, .5, canonical), np.zeros((15, 6, 3)))
        jump = max(jump, float(abs(low - high).max()))
        nonzero += int(max(abs(low).max(), abs(high).max()) > 1e-12)
    assert continuity < 1e-10 and jump < 1e-10
    np.testing.assert_array_equal(edge_basis(np.eye(4))[0], np.zeros((4, 6)))
    pin = int(np.linalg.norm(points, axis=1).argmin())
    translation = gauge_at_node(cells, selection, element, coefficients, pin)
    return {"continuity_max": continuity, "periodic_quartic_trace_max": jump,
            "quartic_face_samples": 15, "nonzero_cap_trace_entries": nonzero,
            "pin_node": pin, "uniform_translation_to_add_to_trial": translation.tolist(),
            "pin_gauge_max": float(abs(translation).max()),
            "pin_scope": "A uniform translation restores the source pin without changing strains, energies, continuity or periodic jumps"}


def nodal_quadratic_field(points, cells, selection, element, coefficients):
    import numpy as np
    from fedoo.lib_elements.tetrahedron import Tet4

    barycentric = Tet4().shape_function(element.xi_nd)
    values = np.array([4 * barycentric[:, first] * barycentric[:, second] for first, second in EDGES]).T
    if selection["relax_corners"]:
        values = np.concatenate((values, mode_basis(barycentric, True)[0][:, 6:]), axis=1)
    identifiers = selection["identifiers"]
    local = coefficients.reshape(-1, 3, 6)[np.maximum(identifiers, 0)].copy()
    local *= (identifiers >= 0)[:, :, None, None]
    sampled = np.einsum("pe,teil->tlpi", values, local)
    fields = np.zeros((6, len(points), 3))
    assigned = np.zeros(len(points), dtype=bool)
    for cell, field in zip(cells[selection["active"]], sampled, strict=True):
        np.testing.assert_allclose(fields[:, cell[assigned[cell]]], field[:, assigned[cell]], atol=1e-12, rtol=0)
        fields[:, cell] = field
        assigned[cell] = True
    return fields


def controls(run, relax_corners=False):
    np, _, meshers, *_ = workflow.imports()
    from cubic_reference import register_cubic, enrich
    from homogenization import strain_tensor
    from periodic import periodic_pairs
    from quadrature_audit import duffy_rule

    parent = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
        bounds=workflow.BOUNDS, cells=4, order=2, periodic=(True, True, True), minimum_quality=.1, threads=1)
    element = register_cubic()()
    points, cells, surface = enrich(parent.points, parent.tetrahedra, parent.surface, element)
    pairs, classes, _ = periodic_pairs(points, surface)
    candidates = cells[:, :4]
    interior = np.flatnonzero((abs(points[candidates]) < .5 - 1e-10).all(axis=(1, 2)))
    assert len(interior)
    chosen = candidates[interior[0]]
    all_edges = {tuple(sorted((int(cell[first]), int(cell[second])))) for cell in candidates for first, second in EDGES}
    selected_keys = {tuple(sorted(map(int, classes[chosen[list(edge)]]))) for edge in EDGES}
    for planes in (1, 2):
        matches = [edge for edge in sorted(all_edges) if
            np.all(workflow.on_box_plane(points[list(edge), :planes], -.5)) and
            (planes == 2 or not np.any(np.all(workflow.on_box_plane(abs(points[list(edge), 1:]), .5), axis=0)))]
        assert matches
        selected_keys.add(tuple(sorted(map(int, classes[list(matches[0])]))))
    selected_edges = [edge for edge in all_edges if tuple(sorted(map(int, classes[list(edge)]))) in selected_keys]
    selection = select_modes(cells, classes, selected_edges, relax_corners)
    affine = np.array([points @ strain_tensor(load).T for load in np.eye(6)])
    matrix, residual, _ = assemble(points, cells, affine, selection, element, duffy_rule(5))
    affine_residual_norm = float(np.linalg.norm(residual))
    assert np.linalg.norm(residual) < 1e-11 * max(np.linalg.norm(matrix), 1)
    random = np.random.default_rng(20261020).normal(size=(len(matrix), 6))
    trace = trace_audit(points, cells, surface, classes, selection, random, element)
    assert trace["nonzero_cap_trace_entries"] > 0
    nonzero_gauges = [float(abs(gauge_at_node(cells, selection, element, random, int(node))).max())
        for node in cells[selection["active"][0]]]
    assert max(nonzero_gauges) > 1e-6
    fields = nodal_quadratic_field(points, cells, selection, element, random)
    for matching in pairs:
        np.testing.assert_allclose(fields[:, matching[:, 1]], fields[:, matching[:, 0]], atol=1e-12, rtol=0)
    matrix, residual, cached = assemble(points, cells, fields, selection, element, duffy_rule(5))
    coefficients = np.linalg.solve(matrix, residual)
    expected = residual.T @ coefficients
    direct, mode = trial_energy(cached, selection, coefficients)
    np.testing.assert_allclose(direct, expected, atol=1e-12, rtol=1e-8)
    np.testing.assert_allclose(mode, coefficients.T @ matrix @ coefficients, atol=1e-12, rtol=1e-10)
    assert np.linalg.norm(expected) > 1e-8
    coupled = coupled_gains(matrix, residual, selection["edge_mode_count"]) if relax_corners else None
    higher, _, _ = assemble(points, cells, fields, selection, element, duffy_rule(6))
    np.testing.assert_allclose(matrix, higher, atol=1e-12, rtol=1e-11)
    overlap = matrix.copy()
    for orbit in range(selection["mode_count"]):
        overlap[3 * orbit:3 * orbit + 3, 3 * orbit:3 * orbit + 3] = 0
    assert np.linalg.norm(overlap) > 1e-6
    workflow.save(run / "controls.json", {"passed": True, "nodes": len(points), "selected_tets": len(selection["active"]),
        "edge_orbits": selection["edge_mode_count"], "corner_orbits": selection["corner_mode_count"],
        "mode_dofs": len(matrix), "affine_residual_norm": affine_residual_norm,
        "trial_identity_relative_error": float(np.linalg.norm(direct - expected) / np.linalg.norm(expected)),
        "nonzero_pin_gauge_control": max(nonzero_gauges), "off_orbit_coupling_norm": float(np.linalg.norm(overlap)),
        "trace": trace, "affine125_vs216_matrix_relative_error": float(np.linalg.norm(matrix - higher) / np.linalg.norm(matrix)),
        "schur_identity_relative_error": None if coupled is None else coupled["schur_identity_relative_error"]})


def lattice(run, higher_rules=False, relax_corners=False):
    np, *_ = workflow.imports()
    from cubic_reference import register_cubic
    from periodic import periodic_pairs
    from periodic_crease_audit import features
    from quadrature_audit import duffy_rule
    from scipy.linalg import eigvalsh

    row = json.loads((SOURCE / "result.json").read_text())
    assert row["qualified"] and row["geometry"]["geometry_admitted"]
    assert row["completed_loads"] == 6 and row["quadrature_sensitivity_passed"]
    parent_path = Path(row["geometry_source"])
    parent_path = parent_path if parent_path.is_absolute() else workflow.ROOT / parent_path
    parent, archive = np.load(parent_path / "mesh.npz"), np.load(PREPARED / "mesh.npz")
    points, cells, surface = (archive[name] for name in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    before, after = parent["tetrahedra"][:, :4].ravel(), cells[:, :4].ravel()
    mapping = dict(zip(map(int, before), map(int, after), strict=True))
    np.testing.assert_array_equal(np.array([mapping[int(node)] for node in before]), after)
    np.testing.assert_allclose(parent["points"][list(mapping)], points[list(mapping.values())], atol=1e-12, rtol=0)
    selected_edges = [tuple(mapping[node] for node in edge) for edge in features(parent)["curve_edges"]]
    _, classes, periodic = periodic_pairs(points, surface)
    selection = select_modes(cells, classes, selected_edges, relax_corners)
    fields = np.load(SOURCE / "fields.npz")["displacement"]
    assert fields.shape == (6, len(points), 3)
    assert abs(fields[:, np.linalg.norm(points, axis=1).argmin()]).max() < 1e-10
    reference = np.asarray(row["stiffness"])
    element = register_cubic()()
    rules, previous = [], None
    orders = (5, 6) if higher_rules else (4, 5)
    for order in orders:
        matrix, residual, cached = assemble(points, cells, fields, selection, element, duffy_rule(order))
        coupled = coupled_gains(matrix, residual, selection["edge_mode_count"]) if relax_corners else None
        coefficients = np.linalg.solve(matrix, residual) if coupled is None else coupled["coefficients"]
        solve_residual = float(np.linalg.norm(matrix @ coefficients - residual) / max(np.linalg.norm(residual), 1e-30))
        assert solve_residual < 1e-10
        gain = residual.T @ coefficients
        direct, mode = trial_energy(cached, selection, coefficients)
        np.testing.assert_allclose(direct, gain, atol=1e-12, rtol=1e-7)
        np.testing.assert_allclose(mode, coefficients.T @ matrix @ coefficients, atol=1e-12, rtol=1e-10)
        spectrum = eigvalsh(gain, reference)
        assert spectrum.min() >= -1e-12
        signal_gain = gain if coupled is None else coupled["incremental_gain"]
        signal_spectrum = eigvalsh(signal_gain, reference)
        extra = {}
        if coupled is not None:
            baseline = np.load(EDGE_BASELINE / f"positive{order**3}-coefficients.npz")
            split = 3 * selection["edge_mode_count"]
            np.testing.assert_allclose(matrix[:split, :split], baseline["mode_stiffness"], atol=1e-12, rtol=1e-10)
            np.testing.assert_allclose(residual[:split], baseline["mode_residual"], atol=1e-12, rtol=1e-8)
            saved = json.loads((EDGE_BASELINE / "result.json").read_text())
            saved_rule = next(rule for rule in saved["rules"] if rule["points"] == order**3)
            saved_gain = np.asarray(saved_rule["energy_reduction_gram"])
            regression = float(np.linalg.norm(coupled["edge_gain"] - saved_gain) / np.linalg.norm(saved_gain))
            assert regression < 1e-7
            extra = {"edge_only_regression_relative_error": regression,
                "schur_identity_relative_error": coupled["schur_identity_relative_error"],
                "edge_only_energy_reduction_gram": coupled["edge_gain"].tolist(),
                "corner_only_energy_reduction_gram": coupled["corner_gain"].tolist(),
                "edge_only_relative_reduction_spectrum": eigvalsh(coupled["edge_gain"], reference).tolist(),
                "corner_only_relative_reduction_spectrum": eigvalsh(coupled["corner_gain"], reference).tolist()}
        trace = trace_audit(points, cells, surface, classes, selection, coefficients, element)
        change = None if previous is None else float(abs(eigvalsh(signal_gain - previous, reference)).max() / max(abs(signal_spectrum).max(), 1e-30))
        rules.append({"points": order**3, "energy_reduction_gram": gain.tolist(),
            "signal_reduction_gram": signal_gain.tolist(), "signal_relative_reduction_spectrum": signal_spectrum.tolist(),
            "mode_solve_relative_residual": solve_residual,
            "mixed_load_relative_reduction_spectrum": spectrum.tolist(), "quadrature_relative_change": change,
            "trial_identity_relative_error": float(np.linalg.norm(direct - gain) / max(np.linalg.norm(gain), 1e-30)), "trace": trace, **extra})
        np.savez_compressed(run / f"positive{order**3}-coefficients.npz", coefficients=coefficients,
                            mode_stiffness=matrix, mode_residual=residual)
        previous = signal_gain
    passed = rules[-1]["quadrature_relative_change"] <= .01
    signal = max(rules[-1]["signal_relative_reduction_spectrum"]) >= 5e-4
    workflow.save(run / "result.json", {"source": str(SOURCE), "prepared": str(PREPARED),
        "quadrature_orders": list(orders),
        "edge_orbits": selection["edge_mode_count"], "corner_orbits": selection["corner_mode_count"],
        "relax_curve_corners": relax_corners, "selected_tets": len(selection["active"]), "mode_dofs": len(matrix),
        "periodic": periodic, "rules": rules, "quadrature_sensitivity_passed": passed,
        "predeclared_curve_signal_passed": bool(passed and signal), "geometry_unchanged": True,
        "retained_models_unchanged": True,
        "scope": "Coupled quartic edge-mode trial energy on a frozen quadratic domain, not global equilibrium, homogenized stiffness, continuum or local-stress certification"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--lattice", action="store_true")
    parser.add_argument("--higher-rules", action="store_true")
    parser.add_argument("--relax-curve-corners", action="store_true")
    parser.add_argument("--worker", choices=("controls", "lattice"))
    args = parser.parse_args()
    if args.relax_curve_corners and not args.higher_rules:
        parser.error("Corner coupling requires the validated higher-rule pair")
    if args.worker:
        if args.worker == "controls":
            controls(args.run, args.relax_curve_corners)
        else:
            lattice(args.run, args.higher_rules, args.relax_curve_corners)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), workflow.HERE / "test_quartic_edge_probe.py"]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "cubic_reference.py", "stress_recovery.py",
                 "periodic_crease_audit.py", "quadrature_audit.py"))
    import meshers._meshers
    paths.extend((Path(meshers._meshers.__file__),
        workflow.ROOT / "crates/meshers-python/python/meshers/__init__.py",
        workflow.ROOT / "crates/meshers-python/python/meshers/expression.py"))
    paths.extend(workflow.ROOT.parent / "periodic-microgen/experiments/periodic-fea" / name
                 for name in ("homogenization.py", "periodic.py"))
    if args.lattice:
        row = json.loads((SOURCE / "result.json").read_text())
        parent = Path(row["geometry_source"])
        parent = parent if parent.is_absolute() else workflow.ROOT / parent
        paths.extend((SOURCE / "result.json", SOURCE / "fields.npz", SOURCE / "manifest.json", PREPARED / "mesh.npz",
                      parent / "mesh.npz", parent / "result.json"))
        frozen = {Path(path).resolve(): value for path, value in json.loads((SOURCE / "manifest.json").read_text())["sources"].items()}
        for path in (PREPARED / "mesh.npz", parent / "mesh.npz", parent / "result.json"):
            assert workflow.digest(path) == frozen[path.resolve()]
        audit = workflow.HERE / "runs/g3-lattice-cubic-assessment-v1"
        checked = json.loads((audit / "result.json").read_text())
        assert checked["fresh_process_audit_passed"] and Path(checked["source"]).resolve() == SOURCE.resolve()
        audited = {Path(path).resolve(): value for path, value in json.loads((audit / "manifest.json").read_text())["inputs"].items()}
        for path in (SOURCE / "result.json", SOURCE / "fields.npz", PREPARED / "mesh.npz"):
            assert workflow.digest(path) == audited[path.resolve()]
        paths.extend((audit / "result.json", audit / "manifest.json"))
        paths.extend(workflow.HERE / "runs" / name for name in (
            "g3-lattice-curve-grading-reference-comparison-v1.json", "g3-lattice-reference-sensitivity-v2.json",
            "g3-lattice-microgen-reference-comparison-v1.json"))
        if args.relax_curve_corners:
            paths.extend(EDGE_BASELINE / name for name in ("result.json", "manifest.json",
                "positive125-coefficients.npz", "positive216-coefficients.npz"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY, "threads": 1,
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "fedoo", "meshers")},
        "worker_timeout_s": 300, "maximum_scalar_modes": 512, "maximum_support_tets": 5000,
        "relax_curve_corners": args.relax_curve_corners,
        "quadrature_points": [125, 216] if args.higher_rules else [64, 125],
        "maximum_relative_gain_quadrature_change": .01,
        "minimum_mixed_load_signal": 5e-4, "equilibrium_solve": False, "lattice_selected": args.lattice})
    stages = ["unit-tests", "controls"] + (["lattice"] if args.lattice else [])
    with Job(parts=stages) as job:
        for stage in stages:
            with job.part(stage) as part:
                command = ([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_quartic_edge_probe.py"]
                    if stage == "unit-tests" else [sys.executable, __file__, "--run", str(args.run), "--worker", stage])
                if stage != "unit-tests" and args.higher_rules:
                    command.append("--higher-rules")
                if stage != "unit-tests" and args.relax_curve_corners:
                    command.append("--relax-curve-corners")
                record = workflow.bounded(command, args.run / (stage + ".log"), part, timeout=300)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage == "lattice":
                    row = json.loads((args.run / "result.json").read_text())
                    values.update(signal_passed=int(row["predeclared_curve_signal_passed"]),
                        quadrature_passed=int(row["quadrature_sensitivity_passed"]),
                        energy_reduction_ppm=max(row["rules"][-1]["mixed_load_relative_reduction_spectrum"]) * 1e6,
                        incremental_signal_ppm=max(row["rules"][-1]["signal_relative_reduction_spectrum"]) * 1e6,
                        quadrature_change_percent=row["rules"][-1]["quadrature_relative_change"] * 100)
                job.result(stage, **values)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
