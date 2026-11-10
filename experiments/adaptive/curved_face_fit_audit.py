"""Replay frozen face-fit evidence and probe a denser independent surface grid."""

import argparse
import json
from pathlib import Path
import shutil

import numpy as np

import comparison_fixtures as fixtures
from comparison_mesh_audit import observe
import microgen_workflow as w
from periodic_trace import periodic_pairs
import scorecard


def audit(source, output):
    w.imports()
    selection = json.loads((source / "selection.json").read_text())
    checks = []
    hashes = {
        str(p): w.digest(p)
        for p in (
            Path(__file__),
            source / "selection.json",
            w.HERE / "comparison_mesh_audit.py",
            w.HERE / "microgen_workflow.py",
        )
    }
    for entry in selection["cases"]:
        name = entry["case"]
        run = source / name
        saved = json.loads((run / "result.json").read_text())
        for filename, digest in saved["inputs"].items():
            path = Path(filename)
            if w.digest(path) != digest:
                path = run / ("input-" + path.name)
            assert path.exists() and w.digest(path) == digest, filename
            hashes[str(path)] = digest
        hashes[str(run / "result.json")] = w.digest(run / "result.json")
        assert w.digest(run / "mesh.npz") == saved["mesh_sha256"]
        hashes[str(run / "mesh.npz")] = saved["mesh_sha256"]
        with np.load(run / "mesh.npz") as m:
            points, cells, surface = (m[k].copy() for k in ("points", "tetrahedra", "surface"))
        with np.load(Path(entry["source"]) / "mesh.npz") as m:
            old = m["points"].copy()
            np.testing.assert_array_equal(cells, m["tetrahedra"])
            np.testing.assert_array_equal(surface, m["surface"])
        np.testing.assert_array_equal(old[np.unique(cells[:, :4])], points[np.unique(cells[:, :4])])
        target = fixtures.case(name, Path(entry["fixtures"]))
        replay = observe(
            run,
            name,
            target,
            points,
            cells,
            surface,
            fixtures.periodic(name),
            periodic_auditor=periodic_pairs,
        )
        assert replay[0]["geometry_admitted"] == saved["admitted"]
        for key in ("sampled_two_sided_max", "relative_volume_error", "minimum_sampled_jacobian"):
            np.testing.assert_allclose(
                replay[0][key], saved["after"]["geometry"][key], atol=1e-12, rtol=1e-12
            )
        # This grid is independent of degree6 fitting and degree11 proxy checks.
        # Original benchmark samples/gates are unchanged; dense values are supplementary.
        previous_bary = scorecard.BARY
        previous_specimen = w.specimen
        degree = 13
        dense = (
            np.array(
                [(i, j, degree - i - j) for i in range(degree + 1) for j in range(degree + 1 - i)],
                dtype=float,
            )
            / degree
        )
        scorecard.BARY = np.unique(np.vstack((previous_bary, dense)), axis=0)
        w.specimen = lambda requested: (
            (target, None, None) if requested == name else previous_specimen(requested)
        )
        try:
            before = w.geometry(
                run,
                name,
                old,
                cells,
                surface,
                reverse_subdivisions=(8, 16),
                box_projection_fallback=True,
            )
            after = w.geometry(
                run,
                name,
                points,
                cells,
                surface,
                reverse_subdivisions=(8, 16),
                box_projection_fallback=True,
            )
        finally:
            scorecard.BARY = previous_bary
            w.specimen = previous_specimen
        row = dict(
            case=name,
            original_admission_reproduced=True,
            counts_and_corner_coordinates_unchanged=True,
            dense_before=before,
            dense_after=after,
            dense_boundary_ratio=after["sampled_two_sided_max"] / before["sampled_two_sided_max"],
            dense_projections_passed=before["projection_failed"] == after["projection_failed"] == 0,
        )
        w.save(output / (name + ".json"), row)
        checks.append(row)
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    shutil.copyfile(__file__, output / "input-curved_face_fit_audit.py")
    w.save(
        output / "result.json",
        dict(
            checks=checks,
            inputs=hashes,
            source_unchanged=True,
            passed=all(
                r["dense_projections_passed"] and r["original_admission_reproduced"] for r in checks
            ),
            original_geometry_gates_unchanged=True,
            dense_sampling_not_global_certificate=True,
            new_FE_solve=False,
            reference_promoted=False,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    assert not (args.run / "result.json").exists()
    audit(args.source, args.run)
