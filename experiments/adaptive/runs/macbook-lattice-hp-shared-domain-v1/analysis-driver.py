"""Verify saved h- and p-enriched fields share one unchanged quadratic domain."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np

RUNS = workflow.HERE / "runs"
ROOT = RUNS / "g3-lattice-cubic-parent23-fields-v1"
FIRST = RUNS / "g3-lattice-nested23-preflight-v1"
FIRST_SOLVE = RUNS / "g3-lattice-nested23-level1-v1"
SECOND = RUNS / "macbook-lattice-nested23-level2-preflight-v1"
SECOND_SOLVE = RUNS / "macbook-lattice-nested23-level2-v1"
CUBIC = RUNS / "g3-lattice-cubic23-assembly-memory-v1"
CUBIC_SOLVE = RUNS / "g3-lattice-cubic23-reference-v1"


def compose(first_parents, first_maps, second_parents, second_maps):
    return first_parents[second_parents], second_maps @ first_maps[second_parents]


def periodic_jumps(points, surface, fields):
    from homogenization import strain_tensor
    from periodic import periodic_pairs

    pairs, _, report = periodic_pairs(points, surface)
    maximum = 0.
    for load, engineering in enumerate(np.eye(6)):
        for axis, matching in enumerate(pairs):
            jump = fields[load, matching[:, 1]] - fields[load, matching[:, 0]]
            maximum = max(maximum, float(abs(jump - strain_tensor(engineering)[:, axis]).max(initial=0)))
    assert maximum < 1e-9
    return {"maximum_six_load_jump": maximum, **report}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from cubic_reference import register_cubic, basis_checks
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4
    from nested_refinement import interpolation_check

    directories = (ROOT, FIRST, FIRST_SOLVE, SECOND, SECOND_SOLVE, CUBIC, CUBIC_SOLVE)
    paths = [Path(__file__), workflow.HERE / "test_hp_shared_domain.py"]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "cubic_reference.py", "nested_refinement.py"))
    for directory in directories:
        paths.extend(directory / name for name in ("manifest.json", "result.json", "mesh.npz", "fields.npz", "prolongation.npz", "audit-inputs.json")
            if (directory / name).exists())
    transfer = workflow.ROOT.parent / "meshers-handoff-20261024/data-manifest.json"
    paths.append(transfer)
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_hp_shared_domain.py", args.output.parent / "test_hp_shared_domain.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_shared_domain.py"], check=True, timeout=30)
    transferred = json.loads(transfer.read_text())
    frozen = {str((workflow.ROOT.parent / repository["name"] / item["path"]).resolve()): item["sha256"]
        for repository in transferred["repositories"] for item in repository["files"]}
    frozen_checked = 0
    for path in paths:
        if path.parent.name.startswith("g3-") and path.suffix == ".npz":
            assert inputs[str(path)] == frozen[str(path.resolve())]
            frozen_checked += 1
    for solved, prepared in ((CUBIC_SOLVE, CUBIC), (FIRST_SOLVE, FIRST), (SECOND_SOLVE, SECOND)):
        manifest = json.loads((solved / "manifest.json").read_text())
        recorded = manifest.get("inputs", manifest.get("sources"))
        linked = []
        for name, expected in recorded.items():
            path = Path(name.replace("/home/kmarchais/dev/", str(workflow.ROOT.parent) + "/", 1))
            path = path if path.is_absolute() else workflow.ROOT / path
            if path.resolve() == (prepared / "mesh.npz").resolve():
                linked.append(expected)
        assert linked == [inputs[str(prepared / "mesh.npz")]], "Solve/prepared mesh provenance differs"
    meshes = {directory: dict(np.load(directory / "mesh.npz")) for directory in directories if (directory / "mesh.npz").exists()}
    root, first, second, cubic = (meshes[directory] for directory in (ROOT, FIRST, SECOND, CUBIC))
    for prepared, solved in ((FIRST, FIRST_SOLVE), (SECOND, SECOND_SOLVE)):
        for key in ("points", "tetrahedra", "surface"):
            np.testing.assert_array_equal(meshes[prepared][key], meshes[solved][key])
    fields = {directory: np.load(directory / "fields.npz")["displacement"] for directory in (ROOT, FIRST_SOLVE, SECOND_SOLVE, CUBIC_SOLVE)}
    for directory, field in fields.items():
        mesh = cubic if directory == CUBIC_SOLVE else meshes[directory]
        assert field.shape == (6, len(mesh["points"]), 3) and np.isfinite(field).all()
        assert len(mesh["points"]) <= 250000 and len(mesh["tetrahedra"]) <= workflow.CAP
        row = json.loads((directory / "result.json").read_text())
        assert row["qualified"] and len(row["solver_records"]) == 6
    first_map, second_map = (dict(np.load(directory / "prolongation.npz")) for directory in (FIRST, SECOND))
    first_check = interpolation_check(root["points"], root["tetrahedra"], fields[ROOT], {**first, **first_map})
    second_check = interpolation_check(first["points"], first["tetrahedra"], fields[FIRST_SOLVE], {**second, **second_map})
    parents, maps = compose(first_map["parents"], first_map["parent_barycentric_maps"],
        second_map["parents"], second_map["parent_barycentric_maps"])
    assert maps.min() >= 0 and maps.max() <= 1
    np.testing.assert_allclose(maps.sum(axis=2), 1., atol=1e-15)
    element = Tet10()
    vertices = maps @ element.xi_nd[:4]
    ratios = np.linalg.det(vertices[:, 1:] - vertices[:, :1]) / np.linalg.det(element.xi_nd[1:4] - element.xi_nd[:1])
    assert ratios.min() > 0
    np.testing.assert_allclose(np.bincount(parents, weights=ratios, minlength=len(root["tetrahedra"])), 1., atol=1e-14)
    barycentric = Tet4().shape_function(element.xi_nd)
    maximum = 0.
    for start in range(0, len(parents), 512):
        selection = slice(start, start + 512)
        mapped = np.einsum("pa,cab->cpb", barycentric, maps[selection])
        values = element.shape_function(mapped.reshape(-1, 4) @ element.xi_nd[:4]).reshape(-1, 10, 10)
        expected = np.einsum("cpn,cnd->cpd", values, root["points"][root["tetrahedra"][parents[selection]]])
        actual = second["points"][second["tetrahedra"][selection]]
        maximum = max(maximum, float(abs(expected - actual).max()))
    assert maximum < 1e-12
    cubic_element = register_cubic()()
    basis_checks(cubic_element)
    assert len(cubic["tetrahedra"]) == len(root["tetrahedra"])
    np.testing.assert_allclose(cubic["points"][cubic["tetrahedra"][:, :4]], root["points"][root["tetrahedra"][:, :4]], atol=1e-13)
    embedded = np.einsum("pn,cnd->cpd", element.shape_function(cubic_element.xi_nd), root["points"][root["tetrahedra"]])
    cubic_error = float(abs(embedded - cubic["points"][cubic["tetrahedra"]]).max())
    assert cubic_error < 1e-12
    periodic = {"nested": periodic_jumps(second["points"], second["surface"], fields[SECOND_SOLVE]),
        "cubic": periodic_jumps(cubic["points"], cubic["surface"], fields[CUBIC_SOLVE])}
    np.savez_compressed(args.output.parent / "composed-maps.npz", parents=parents, parent_barycentric_maps=maps)
    workflow.save(args.output, {"passed": True, "root": str(ROOT), "nested": str(SECOND_SOLVE), "cubic": str(CUBIC_SOLVE),
        "cubic_mesh": str(CUBIC / "mesh.npz"), "frozen_npz_inputs_verified": frozen_checked,
        "nested_nodes": len(second["points"]), "nested_tets": len(parents), "cubic_nodes": len(cubic["points"]),
        "cubic_tets": len(cubic["tetrahedra"]), "first_level_field_reproduction": first_check,
        "second_level_field_reproduction": second_check, "maximum_composed_geometry_error": maximum,
        "maximum_cubic_geometry_error": cubic_error, "periodic": periodic,
        "composed_maps_sha256": workflow.digest(args.output.parent / "composed-maps.npz"),
        "unisolvent_geometry_checks": "All ten Tet10 child nodes and all twenty cubic parent nodes; affine child maps preserve polynomial degrees",
        "source_meshes_and_fields_unchanged": True, "new_equilibrium_solve": False, "reference_promoted": False,
        "scope": "Shared physical-domain and field-provenance preflight; different equilibrium fields are not asserted equal or combined into a new reference."})
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())
    print(f"Shared-domain preflight passes: nested {maximum:.3g}, cubic {cubic_error:.3g}")


if __name__ == "__main__":
    main()
