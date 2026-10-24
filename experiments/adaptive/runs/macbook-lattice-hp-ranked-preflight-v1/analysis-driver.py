"""One measured-score nested preflight; no field or reference promotion."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np

from hp_coupled_preflight import inherited_inputs
from hp_shared_domain import SECOND_SOLVE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--localization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from nested_refinement import refine, interpolation_check
    from nested_preflight import fixture_inputs
    from periodic_crease_audit import features
    from hp_shared_domain import periodic_jumps

    localization = json.loads((args.localization / "result.json").read_text())
    assert localization["passed"] and max(localization["positive_rule_local_energy_relative_changes"].values()) < .01
    audit_path = Path(localization["source_audit"])
    audit = json.loads((audit_path / "result.json").read_text())
    source = Path(audit["source"])
    solved = json.loads((source / "result.json").read_text())
    consistent = Path(solved["consistent_source"])
    paths = [Path(__file__), args.localization / "result.json", args.localization / "analysis-inputs.json",
        source / "hp-fields.npz", source / "result.json"]
    paths.extend(workflow.HERE / name for name in ("nested_refinement.py", "test_nested_refinement.py",
        "nested_preflight.py", "periodic_crease_audit.py", "hp_coupled_preflight.py", "microgen_workflow.py"))
    paths.extend(args.localization / f"local-{count}.npz" for count in (125, 216))
    fixtures = fixture_inputs(SECOND_SOLVE)
    paths.extend(fixtures)
    inputs = {**inherited_inputs(args.localization), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    for name in ("nested_refinement.py", "test_nested_refinement.py", "hp_coupled_preflight.py"):
        shutil.copyfile(workflow.HERE / name, args.output.parent / name)
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_nested_refinement.py"], check=True, timeout=60)
    mesh = dict(np.load(SECOND_SOLVE / "mesh.npz"))
    groups = features(mesh)["groups"]
    baseline = np.load(consistent / "coupling.npz")["energy"][:6, :6]
    previous = None
    scalar_precision = {}
    for count in (125, 216):
        row = next(row for row in localization["rules"] if row["points"] == count)
        path = args.localization / f"local-{count}.npz"
        assert workflow.digest(path) == row["local_archive_sha256"]
        local = np.load(path)
        tensors = local["energy_tensors"]
        np.testing.assert_allclose(tensors.sum(axis=0), row["correction_energy_tensor"], atol=1e-12)
        scalars = {"six_load_normalized_trace": np.einsum("ij,tji->t", np.linalg.inv(baseline), tensors)}
        scalars.update({name: np.einsum("i,tij,j->t", np.asarray(localization[key]), tensors, np.asarray(localization[key]))
            for name, key in (("largest_correction", "gain_witness"), ("worst_retained_model", "worst_model_witness"))})
        for name, values in scalars.items():
            np.testing.assert_allclose(values, local[name], atol=1e-15)
            for label, mask in groups.items():
                np.testing.assert_allclose(values[mask].sum() / values.sum(), row["summaries"][name]["groups"][label]["energy_fraction"], atol=1e-12)
            if previous is not None:
                scalar_precision[name] = float(np.linalg.norm(values - previous[name]) / np.linalg.norm(values))
                assert scalar_precision[name] < .01
        previous = scalars
    scores = previous["six_load_normalized_trace"]
    assert scores.min() >= 0
    original = np.load(SECOND_SOLVE / "fields.npz")["displacement"]
    modal = np.load(source / "hp-fields.npz")
    base = original - np.einsum("ij,ind->jnd", modal["mode_coefficients"], modal["response_displacement"])
    result = refine(mesh["points"], mesh["tetrahedra"], mesh["surface"], base, element_scores=scores)
    reproduction = interpolation_check(mesh["points"], mesh["tetrahedra"], base, result)
    for path in fixtures:
        shutil.copyfile(path, args.output.parent / path.name)
    np.savez_compressed(args.output.parent / "mesh.npz", **{key: result[key] for key in ("points", "tetrahedra", "surface")})
    np.savez_compressed(args.output.parent / "prolongation.npz", **{key: result[key]
        for key in ("displacement", "parents", "parent_barycentric_maps")})
    geometry = workflow.geometry(args.output.parent, "smooth-lattice", result["points"], result["tetrahedra"], result["surface"])
    old = json.loads((SECOND_SOLVE / "result.json").read_text())["geometry"]
    volume_change = geometry["volume"] / old["volume"] - 1
    assert abs(volume_change) < 1e-10 and geometry["geometry_admitted"]
    periodic = periodic_jumps(result["points"], result["surface"], result["displacement"])
    multiplicity = np.bincount(result["parents"], minlength=len(mesh["tetrahedra"]))
    selected = multiplicity == 2
    assert np.all((multiplicity == 1) | selected)
    pin = int(np.linalg.norm(mesh["points"], axis=1).argmin())
    np.testing.assert_array_equal(result["points"][pin], mesh["points"][pin])
    workflow.save(args.output, {"passed": True, "geometry": geometry, "periodic_base_field": periodic,
        "physical_nodes": len(result["points"]), "physical_tets": len(result["tetrahedra"]),
        "selection": result["diagnostics"], "reproduction": reproduction, "relative_volume_change": volume_change,
        "selected_parent_fraction": float(selected.mean()), "selected_parent_score_fraction": float(scores[selected].sum() / scores.sum()),
        "localization_scalar_replays": scalar_precision, "source_coupled_solve": str(source),
        "preserved_source_pin_node": pin, "new_nearest_node": int(np.linalg.norm(result["points"], axis=1).argmin()),
        "fixed_cubic_function_span_preserved": True, "source_augmented_space_included": True,
        "extra_scalar_mode_unknowns": 6, "new_equilibrium_solve": False, "reference_promoted": False,
        "scope": "Measured six-load score per added node on complete periodic stars; exact physical-domain and FE-base prolongation, with the unchanged cubic function span; no accuracy or workflow claim."})
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
    print(f"Ranked preflight passes: {len(result['points'])} nodes, {len(result['tetrahedra'])} tets")


if __name__ == "__main__":
    main()
