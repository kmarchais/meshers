"""Localize resolved full modal correction energy; never estimate continuum error."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np
from scipy.linalg import eigh

from hp_coupled_preflight import inherited_inputs
from hp_shared_domain import SECOND_SOLVE, CUBIC_SOLVE, CUBIC


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_cross_energy import integrate
    from quadrature_audit import duffy_rule
    from periodic_crease_audit import features
    from fedoo.lib_elements.tetrahedron import Tet10

    audit = json.loads((args.audit / "result.json").read_text())
    assert audit["passed"] and not audit["reference_promoted"]
    solved_path = Path(audit["source"])
    solved = json.loads((solved_path / "result.json").read_text())
    consistent_path = Path(solved["consistent_source"])
    consistent = json.loads((consistent_path / "result.json").read_text())
    cross = json.loads((Path(consistent["source"]) / "result.json").read_text())
    shared = Path(cross["preflight"])
    paths = [Path(__file__), args.audit / "result.json", args.audit / "analysis-inputs.json", solved_path / "hp-fields.npz"]
    paths.extend(workflow.HERE / name for name in ("hp_cross_energy.py", "test_hp_cross_energy.py",
        "periodic_crease_audit.py", "hp_coupled_preflight.py", "quadrature_audit.py"))
    inputs = {**inherited_inputs(args.audit), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    for name in ("hp_cross_energy.py", "test_hp_cross_energy.py", "hp_coupled_preflight.py"):
        shutil.copyfile(workflow.HERE / name, args.output.parent / name)
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_cross_energy.py"], check=True, timeout=45)
    mesh, cubic = np.load(SECOND_SOLVE / "mesh.npz"), np.load(CUBIC / "mesh.npz")
    assert len(mesh["points"]) <= 250000 and len(mesh["tetrahedra"]) <= workflow.CAP
    h = np.load(SECOND_SOLVE / "fields.npz")["displacement"]
    p = np.load(CUBIC_SOLVE / "fields.npz")["displacement"]
    fields = np.load(solved_path / "hp-fields.npz")
    responses, coefficients = fields["response_displacement"], fields["mode_coefficients"]
    maps = np.load(shared / "composed-maps.npz")
    baseline = np.load(consistent_path / "coupling.npz")["energy"][:6, :6]
    expected_gain = baseline - np.asarray(audit["rules"][0]["energy_tensor"])
    values, vectors = eigh(expected_gain, baseline)
    gain_witness = vectors[:, -1]
    worst = max(audit["models"], key=lambda row: row["worst_relative_difference"])
    model = np.asarray(json.loads(Path(worst["source"]).read_text())["models"][worst["index"]]["stiffness"])
    errors, witnesses = eigh(np.asarray(audit["rules"][-1]["energy_tensor"]) - model, model)
    model_witness = witnesses[:, np.argmax(abs(errors))]
    directions = {"largest_correction": gain_witness, "worst_retained_model": model_witness}
    groups = features(mesh)["groups"]
    correction_basis = np.vstack((coefficients, -coefficients, coefficients))
    records, scalars_by_rule = [], []
    for count in (15, 125, 216):
        element = Tet10(15)
        quadrature = (element.xi_pg, element.w_pg) if count == 15 else duffy_rule(5 if count == 125 else 6)
        energy, volume, minimum, local, volumes = integrate(mesh["points"], mesh["tetrahedra"], h,
            cubic["tetrahedra"], p, maps["parents"], maps["parent_barycentric_maps"], quadrature,
            additional_fields=responses, local_coefficients=correction_basis)
        total = local.sum(axis=0)
        np.testing.assert_allclose(total, correction_basis.T @ energy @ correction_basis, atol=1e-12)
        identity = float(abs(eigh(total - expected_gain, baseline, eigvals_only=True)).max())
        assert identity < (1e-7 if count == 15 else 1e-4)
        scalars = {"six_load_normalized_trace": np.einsum("ij,tji->t", np.linalg.inv(baseline), local)}
        scalars.update({name: np.einsum("i,tij,j->t", direction, local, direction) for name, direction in directions.items()})
        summaries = {}
        for name, scalar in scalars.items():
            assert scalar.sum() > 0
            if count != 15:
                assert scalar.min() > -1e-12 * scalar.sum()
            order = np.argsort(scalar)[::-1]
            summaries[name] = {"total": float(scalar.sum()),
                "top_element_energy_fractions": {str(fraction): float(scalar[order[:max(1, int(np.ceil(fraction * len(scalar))))]].sum() / scalar.sum())
                    for fraction in (.001, .01, .05, .1, .2)},
                "groups": {label: {"element_fraction": float(mask.mean()),
                    "energy_fraction": float(scalar[mask].sum() / scalar.sum())} for label, mask in groups.items()},
                "top_elements": [{"element": int(index), "energy": float(scalar[index]),
                    "parent": int(maps["parents"][index]), "corner_nodes": mesh["tetrahedra"][index, :4].tolist()}
                    for index in order[:20]]}
        np.savez_compressed(args.output.parent / f"local-{count}.npz", energy_tensors=local, volumes=volumes, **scalars)
        records.append({"points": count, "minimum_sampled_jacobian": minimum, "volume": volume,
            "local_archive_sha256": workflow.digest(args.output.parent / f"local-{count}.npz"),
            "galerkin_identity_error_normalized_to_baseline": identity, "correction_energy_tensor": total.tolist(),
            "summaries": summaries})
        scalars_by_rule.append(scalars)
        workflow.save(args.output.parent / f"rule-{count}.json", records[-1])
        print(f"Localized complete correction at {count} points; identity {identity:.9g}", flush=True)
    precision = {name: float(np.linalg.norm(scalars_by_rule[-1][name] - scalars_by_rule[-2][name])
        / np.linalg.norm(scalars_by_rule[-1][name])) for name in scalars_by_rule[-1]}
    assert max(precision.values()) <= .01
    workflow.save(args.output, {"passed": True, "source_audit": str(args.audit), "rules": records,
        "positive_rule_local_energy_relative_changes": precision, "precision_gate": .01,
        "gain_witness": gain_witness.tolist(), "worst_model_witness": model_witness.tolist(),
        "worst_model": worst, "field_or_geometry_modified": False, "reference_promoted": False,
        "scope": "Resolved discrete correction-energy localization in all six loads and two declared mixed directions, not remaining continuum error or local stress certification."})
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
