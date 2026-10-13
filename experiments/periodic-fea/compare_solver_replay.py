"""Check solver agreement on the identical mesh, including volume-weighted stresses."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from accuracy_metrics import field_score
from homogenization import element_geometry
from periodic import tensor_difference
from residual_analysis import stress_and_strain
from run import digest, write_json
from verify import verify_run


def compare(output):
    document = json.loads((output / "results.json").read_text())
    assert document["state"] == "complete"
    rows = document["rows"]
    assert [row["workflow"] for row in rows] == ["solver-cg-replay", "solver-pardiso-replay"]
    assert all(row["status"] == "ok" and all(row[key] for key in
               ["geometry_eligible", "quality_eligible", "periodic_eligible", "homogenization_audited"]) for row in rows)
    assert rows[0]["mesh_sha256"] == rows[1]["mesh_sha256"] == rows[0]["input_mesh_sha256"]
    audit_hash = digest(output / "artifact-verification.json")
    verify_run(output, write_report=False)
    assert digest(output / "artifact-verification.json") == audit_hash
    mesh = np.load(output / rows[0]["artifact"] / "mesh.npz")
    volume, gradients = element_geometry(SimpleNamespace(points=mesh["points"], tetrahedra=mesh["tetrahedra"]))
    stresses = []
    for row in rows:
        loads = []
        for field in row["homogenization"]["displacement_fields"]:
            path = output / row["artifact"] / field["artifact"]
            assert digest(path) == field["sha256"]
            with np.load(path) as data:
                stress, _ = stress_and_strain(data["displacement"], mesh["tetrahedra"], gradients)
            loads.append(stress[:, [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]])
        stresses.append(np.stack(loads, axis=-1) * np.sqrt(volume)[:, None, None])
    local = field_score(stresses[1], stresses[0])
    tensor = tensor_difference(rows[1]["homogenization"]["stiffness_engineering_voigt"],
                               rows[0]["homogenization"]["stiffness_engineering_voigt"])
    assert tensor["worst_relative_strain_energy_difference"] < 1e-7
    assert local["worst_relative_stress_energy_norm_difference"] < 1e-7
    times = []
    for row in rows:
        h = row["homogenization"]
        blocks = [h["strain_controlled_loads"], h["automatic_stress_controlled_solves"]]
        assert [len(block) for block in blocks] == [6, 6]
        if row["workflow"] == "solver-pardiso-replay":
            for block in blocks:
                assert len({r["matrix_sha256"] for r in block}) == 1
                assert [r["factorizations"] for r in block] == [1] * 6
        times.append({"workflow": row["workflow"], "homogenization_s": h["seconds"],
                      "reduced_solves_s": sum(r["seconds"] for block in blocks for r in block),
                      "peak_rss_mib": row["peak_rss_mib"]})
    result = {"state": "complete", "source_sha256": digest(Path(__file__)),
              "results_sha256": digest(output / "results.json"), "unchanged_audit_sha256": audit_hash,
              "stiffness": tensor, "volume_weighted_local_stress": local, "timings": times,
              "scope": "One-thread single calibration samples; same saved mesh, twelve audited systems per solver; no meshing speed claim"}
    target = output / "solver-comparison.json"
    if target.exists():
        raise FileExistsError(target)
    write_json(target, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.run), indent=2))
