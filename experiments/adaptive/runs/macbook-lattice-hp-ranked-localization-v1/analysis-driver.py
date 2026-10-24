"""Locate the audited ranked increment and compare frozen wall residuals."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy.linalg import eigvalsh

import microgen_workflow as workflow
from hp_coupled_preflight import inherited_inputs
from hp_shared_domain import ROOT, SECOND_SOLVE, CUBIC, CUBIC_SOLVE


def increment_bases(old, new):
    """Basis order is old H, original P, lifted old base, new base."""
    previous = np.vstack((-old, old, np.eye(6), np.zeros((6, 6))))
    current = np.vstack((np.zeros((6, 6)), new, -new, np.eye(6)))
    return previous, current, previous - current


def wall_residuals(points, faces, function, order):
    """Positive surface quadrature and sampled normal projection, no error bound."""
    from periodic_crease_audit import project_box_curve

    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    area, squared, signed, maximum = (np.zeros(len(faces)) for _ in range(4))
    for r, wr in zip(nodes, weights, strict=True):
        for s, ws in zip(nodes, weights, strict=True):
            a, b, c = (1 - r) * (1 - s), r, (1 - r) * s
            values = workflow.tri_values([[a, b, c]])[0]
            derivative = np.array([[-(4*a-1), 4*b-1, 0, 4*(a-b), 4*c, -4*c],
                [-(4*a-1), 0, 4*c-1, -4*b, 4*b, 4*(a-c)]])
            probe = np.einsum("n,fnd->fd", values, points[faces])
            tangent = np.einsum("an,fnd->fad", derivative, points[faces])
            normal = np.cross(tangent[:, 0], tangent[:, 1])
            density = np.linalg.norm(normal, axis=1)
            assert density.min() > 0
            foot, _ = project_box_curve(function, probe)
            distance = np.linalg.norm(foot - probe, axis=1)
            measure = wr * ws * (1 - r) * density
            area += measure
            squared += measure * distance**2
            signed += measure * np.einsum("fd,fd->f", foot - probe, normal / density[:, None])
            maximum = np.maximum(maximum, distance)
    return dict(area=area, squared_distance_integral=squared, signed_normal_integral=signed,
        maximum_sampled_projection=maximum)


def quadratic_lift(points, cells, old_cells, fields, parents, maps):
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    quadratic, linear = Tet10(), Tet4()
    nodes, first = np.unique(cells.ravel(), return_index=True)
    np.testing.assert_array_equal(nodes, np.arange(len(points)))
    child, local = np.divmod(first, 10)
    barycentric = linear.shape_function(quadratic.xi_nd)
    result = np.empty((6, len(points), 3))
    for start in range(0, len(points), 2048):
        selection = slice(start, start + 2048)
        locations = np.einsum("na,nab->nb", barycentric[local[selection]], maps[child[selection]]) @ linear.xi_nd
        values = quadratic.shape_function(locations)
        result[:, selection] = np.einsum("tn,ltnd->ltd", values, fields[:, old_cells[parents[child[selection]]]])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_cross_energy import integrate
    from hp_coupled_audit import cubic_nodal_values
    from hp_shared_domain import periodic_jumps
    from nested_refinement import interpolation_check, FACES
    from nested_preflight import fixture_inputs
    from periodic_crease_audit import features
    from quadrature_audit import duffy_rule
    from fedoo.lib_elements.tetrahedron import Tet10

    audit = json.loads((args.audit / "result.json").read_text())
    assert audit["passed"] and audit["retained_model_count"] == 24 and not audit["reference_promoted"]
    solved_path = Path(audit["source"])
    solved = json.loads((solved_path / "result.json").read_text())
    prepared, coupling = Path(solved["prepared"]), Path(solved["coupling"])
    old_audit_path = Path(audit["previous_audit"])
    old_audit = json.loads((old_audit_path / "result.json").read_text())
    old_source = Path(old_audit["source"])
    paths = [Path(__file__), args.audit / "result.json", solved_path / "hp-fields.npz",
        old_source / "hp-fields.npz", prepared / "prolongation.npz", coupling / "root-maps.npz"]
    paths.extend(workflow.HERE / name for name in ("hp_cross_energy.py", "hp_coupled_audit.py",
        "periodic_crease_audit.py", "quadrature_audit.py", "nested_refinement.py", "nested_preflight.py",
        "test_hp_ranked_localization.py"))
    paths.extend(fixture_inputs(ROOT))
    inputs = {**inherited_inputs(args.audit), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_hp_ranked_localization.py", args.output.parent / "test_hp_ranked_localization.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_ranked_localization.py"], check=True, timeout=45)
    mesh, previous_mesh, cubic = np.load(prepared / "mesh.npz"), np.load(SECOND_SOLVE / "mesh.npz"), np.load(CUBIC / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    transfer = np.load(prepared / "prolongation.npz")
    old_h = np.load(SECOND_SOLVE / "fields.npz")["displacement"]
    h = quadratic_lift(points, cells, previous_mesh["tetrahedra"], old_h, transfer["parents"], transfer["parent_barycentric_maps"])
    reproduction = interpolation_check(previous_mesh["points"], previous_mesh["tetrahedra"], old_h,
        dict(points=points, tetrahedra=cells, displacement=h, parents=transfer["parents"], parent_barycentric_maps=transfer["parent_barycentric_maps"]))
    old_fields, new_fields = np.load(old_source / "hp-fields.npz"), np.load(solved_path / "hp-fields.npz")
    previous, current, delta = increment_bases(old_fields["mode_coefficients"], new_fields["mode_coefficients"])
    extras = np.concatenate((transfer["displacement"], new_fields["base_displacement"]))
    p = np.load(CUBIC_SOLVE / "fields.npz")["displacement"]
    maps = np.load(coupling / "root-maps.npz")
    nodal_p = cubic_nodal_values(points, cells, cubic["tetrahedra"], p, maps["parents"], maps["parent_barycentric_maps"])
    nodal_basis = np.concatenate((h, nodal_p, extras))
    periodic_jumps(points, surface, np.einsum("ij,ind->jnd", previous, nodal_basis))
    periodic_jumps(points, surface, np.einsum("ij,ind->jnd", current, nodal_basis))
    baseline = np.asarray(old_audit["rules"][-1]["energy_tensor"])
    groups = features(mesh)["groups"]
    rows, scalars = [], []
    for count in (15, 125, 216):
        element = Tet10(15)
        quadrature = (element.xi_pg, element.w_pg) if count == 15 else duffy_rule(5 if count == 125 else 6)
        energy, volume, minimum, local, volumes = integrate(points, cells, h, cubic["tetrahedra"], p,
            maps["parents"], maps["parent_barycentric_maps"], quadrature, additional_fields=extras, local_coefficients=delta)
        old_energy, new_energy = previous.T @ energy @ previous, current.T @ energy @ current
        np.testing.assert_allclose(local.sum(axis=0), delta.T @ energy @ delta, atol=1e-12)
        identity = float(abs(eigvalsh(local.sum(axis=0) - old_energy + new_energy, baseline)).max())
        assert identity < (1e-7 if count == 15 else 1e-4)
        assert abs(eigvalsh(new_energy - np.asarray(audit["rules"][(15, 125, 216).index(count)]["energy_tensor"]), baseline)).max() < 1e-7
        assert abs(eigvalsh(old_energy - baseline, baseline)).max() < 1e-4
        scalar = np.einsum("ij,tji->t", np.linalg.inv(baseline), local)
        assert scalar.min() >= -1e-12 * scalar.sum() and scalar.sum() > 0
        order = np.argsort(scalar)[::-1]
        np.savez_compressed(args.output.parent / f"local-{count}.npz", energy_tensors=local, volumes=volumes,
            six_load_normalized_trace=scalar)
        row = dict(points=count, volume=volume, minimum_sampled_jacobian=minimum,
            galerkin_identity_error=identity, correction_energy_tensor=local.sum(axis=0).tolist(),
            trace=float(scalar.sum()), local_archive_sha256=workflow.digest(args.output.parent / f"local-{count}.npz"),
            groups={name: dict(element_fraction=float(mask.mean()), correction_fraction=float(scalar[mask].sum()/scalar.sum())) for name, mask in groups.items()},
            top_fractions={str(f): float(scalar[order[:max(1, int(np.ceil(f*len(cells))))]].sum()/scalar.sum()) for f in (.01, .05, .1, .2)})
        rows.append(row); scalars.append(scalar)
        workflow.save(args.output.parent / f"rule-{count}.json", row)
        print(f"Ranked increment localized at {count} points; identity {identity:.9g}", flush=True)
    precision = float(np.linalg.norm(scalars[-1]-scalars[-2])/np.linalg.norm(scalars[-1]))
    assert precision < .01
    root = np.load(ROOT / "mesh.npz")
    owners = {tuple(sorted(face)): index for index, cell in enumerate(root["tetrahedra"][:, :4]) for face in cell[FACES]}
    triangle = root["points"][root["surface"][:, :3]]
    caps = np.any(np.all(workflow.on_box_plane(triangle, -.5), axis=1) | np.all(workflow.on_box_plane(triangle, .5), axis=1), axis=1)
    wall = root["surface"][~caps]
    wall_owners = np.array([owners[tuple(sorted(face))] for face in wall[:, :3]])
    case, _, _ = workflow.specimen("smooth-lattice")
    residuals = [wall_residuals(root["points"], wall, case.function, order) for order in (5, 6)]
    geometry_precision = float(np.linalg.norm(residuals[-1]["squared_distance_integral"]-residuals[-2]["squared_distance_integral"])/np.linalg.norm(residuals[-1]["squared_distance_integral"]))
    root_score = np.bincount(maps["parents"], weights=scalars[-1], minlength=len(root["tetrahedra"]))
    wall_cells = np.unique(wall_owners)
    square, area = (np.bincount(wall_owners, weights=residuals[-1][key], minlength=len(root_score)) for key in ("squared_distance_integral", "area"))
    rms = np.sqrt(square[wall_cells]/area[wall_cells])
    worst_geometry = wall_cells[np.argsort(rms)[-max(1, int(np.ceil(.1*len(wall_cells)))):]]
    np.savez_compressed(args.output.parent / "geometry-overlay.npz", root_score=root_score, wall_owners=wall_owners,
        root_parents=maps["parents"], **residuals[-1])
    geometry = dict(positive_surface_rule_points=[25, 36], squared_residual_vector_change=geometry_precision,
        precision_gate=.01, precision_passed=geometry_precision < .01,
        wall_area=float(area.sum()), area_weighted_projection_rms=float(np.sqrt(square.sum()/area.sum())),
        maximum_sampled_projection=float(residuals[-1]["maximum_sampled_projection"].max()),
        wall_owner_correction_fraction=float(root_score[wall_cells].sum()/root_score.sum()),
        worst_geometry_decile_wall_correction_fraction=float(root_score[worst_geometry].sum()/root_score[wall_cells].sum()),
        root_wall_cells=len(wall_cells), geometry_overlay_sha256=workflow.digest(args.output.parent / "geometry-overlay.npz"))
    workflow.save(args.output, dict(passed=geometry["precision_passed"], source_audit=str(args.audit),
        ranked_increment_localization=True, reproduction=reproduction, rules=rows, geometry=geometry,
        positive_rule_local_trace_change=precision, geometry_or_fields_modified=False, reference_promoted=False,
        scope="Resolved ranked discrete correction and sampled wall-projection residual co-location on the unchanged parent23 domain; no remaining FE/geometry error estimate, stress certificate, native improvement or workflow claim."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
