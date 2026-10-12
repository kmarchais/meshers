"""Reconstruct periodic stiffness and equilibrium from saved displacement fields."""

import hashlib
import json
from pathlib import Path

import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_run(output):
    document = json.loads((output / "results.json").read_text(encoding="utf-8"))
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert sha(output / "manifest.json") == document["manifest_sha256"]
    for name, expected in manifest["source_sha256"].items():
        assert sha(output / "source" / name) == expected, name
    for reference in document["references"].values():
        assert sha(reference["path"]) == reference["sha256"]
    verified = []
    for row in document["rows"]:
        if not row.get("homogenization_audited"):
            continue
        directory = output / row["artifact"]
        assert sha(directory / "mesh.npz") == row["mesh_sha256"]
        mesh, pairing = np.load(directory / "mesh.npz"), np.load(directory / "periodic.npz")
        points, tets = mesh["points"], mesh["tetrahedra"]
        vertices = points[tets]
        jacobian = vertices[:, 1:] - vertices[:, :1]
        volumes = np.linalg.det(jacobian) / 6
        assert volumes.min() > 0
        gradients = np.empty((len(tets), 4, 3))
        gradients[:, 1:] = np.linalg.inv(jacobian).transpose(0, 2, 1)
        gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
        columns, strains, stresses, force_errors, jump_errors = [], [], [], [], []
        fields = row["homogenization"]["displacement_fields"]
        assert len(fields) == 6
        for index, field_record in enumerate(fields):
            field_path = directory / field_record["artifact"]
            assert sha(field_path) == field_record["sha256"]
            with np.load(field_path) as field:
                displacement, engineering = field["displacement"], field["engineering"]
            assert np.array_equal(engineering, np.eye(6)[index])
            gradient = np.einsum("tia,tib->tab", displacement[tets], gradients)
            strain = (gradient + gradient.transpose(0, 2, 1)) / 2
            stress = 2 * strain + np.trace(strain, axis1=1, axis2=2)[:, None, None] * np.eye(3)
            average = np.einsum("t,tij->ij", volumes, stress)
            columns.append(average[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]])
            forces = np.zeros_like(points)
            element_forces = np.einsum("tab,tib,t->tia", stress, gradients, volumes)
            np.add.at(forces, tets.ravel(), element_forces.reshape(-1, 3))
            balanced = np.zeros_like(forces)
            np.add.at(balanced, pairing["classes"], forces)
            force_errors.append(float(np.linalg.norm(balanced) / np.linalg.norm(forces)))
            xx, yy, zz, xy, xz, yz = engineering
            macro = np.array([[xx, xy / 2, xz / 2], [xy / 2, yy, yz / 2], [xz / 2, yz / 2, zz]])
            jump_errors.append(
                max(
                    float(
                        np.max(
                            np.abs(
                                displacement[pairing[f"pairs_{axis}"][:, 1]]
                                - displacement[pairing[f"pairs_{axis}"][:, 0]]
                                - (
                                    points[pairing[f"pairs_{axis}"][:, 1]]
                                    - points[pairing[f"pairs_{axis}"][:, 0]]
                                )
                                @ macro.T
                            )
                        )
                    )
                    for axis in range(3)
                )
            )
            strains.append(strain)
            stresses.append(stress)
        stiffness = np.asarray(columns).T
        reported = np.asarray(row["homogenization"]["stiffness_engineering_voigt"])
        difference = float(np.linalg.norm(stiffness - reported) / np.linalg.norm(reported))
        energy = np.einsum("t,itkl,jtkl->ij", volumes, strains, stresses)
        hill_mandel = float(np.linalg.norm(energy - stiffness) / np.linalg.norm(stiffness))
        assert difference < 1e-10, difference
        assert hill_mandel < 1e-7, hill_mandel
        assert max(force_errors) < 1e-7, force_errors
        assert max(jump_errors) < 1e-10, jump_errors
        verified.append(
            {
                "artifact": row["artifact"],
                "fields": 6,
                "stiffness_reconstruction_relative_difference": difference,
                "hill_mandel_relative_difference": hill_mandel,
                "maximum_periodic_force_balance_error": max(force_errors),
                "maximum_periodic_displacement_jump_error": max(jump_errors),
            }
        )
    result = {
        "verified_rows": len(verified),
        "rows": verified,
        "scope": "Saved fields, source/reference/mesh hashes, stiffness, energy and periodic equilibrium",
    }
    (output / "artifact-verification.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result
