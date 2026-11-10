"""Independent dense replay, fixed-space checks and evidence for fresh native fitting."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    a = p.parse_args()
    np, *_ = w.imports()
    import comparison_fixtures as fixtures
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs
    import scorecard

    manifest = json.loads((a.source / "manifest.json").read_text())
    resolved = {}
    for name, digest in manifest["inputs"].items():
        path = Path(name)
        if not path.exists() or w.digest(path) != digest:
            candidates = [p for p in a.source.glob("input-*-" + path.name) if w.digest(p) == digest]
            assert len(candidates) == 1, name
            path = candidates[0]
        assert w.digest(path) == digest
        resolved[str(path)] = digest
    rows = json.loads((a.source / "summary.json").read_text())["rows"]
    report = []
    for entry in manifest["selection"]:
        name = entry["case"]
        pair = [
            next(r for r in rows if r["case"] == name and r["config"]["curved_fit"] == flag)
            for flag in (False, True)
        ]
        assert pair[0]["frozen_baseline_arrays_equal"]
        meshes = []
        for row in pair:
            run = Path(row["directory"])
            for filename in (
                "mesh.npz",
                "result.json",
                name + "-samples.npz",
                name + "-specimen.json",
            ):
                path = run / filename
                resolved[str(path)] = w.digest(path)
            with np.load(run / "mesh.npz") as mesh:
                meshes.append([mesh[k].copy() for k in ("points", "tetrahedra", "surface")])
        before, after = meshes
        np.testing.assert_array_equal(before[1], after[1])
        np.testing.assert_array_equal(before[2], after[2])
        corners = np.unique(before[1][:, :4])
        np.testing.assert_array_equal(before[0][corners], after[0][corners])
        assert before[0].shape == after[0].shape
        locked = np.abs(before[0]) == 0.5
        np.testing.assert_array_equal(before[0][locked], after[0][locked])
        wall_corners = before[0][before[2][:, :3]]
        cap = np.any(
            np.all(w.on_box_plane(wall_corners, -0.5), axis=1)
            | np.all(w.on_box_plane(wall_corners, 0.5), axis=1),
            axis=1,
        )
        wall = before[2][~cap]
        mids = wall[:, 3:].ravel()
        length = np.linalg.norm(
            before[0][wall[:, [1, 2, 0]]] - before[0][wall[:, :3]], axis=2
        ).ravel()
        limits = np.full(len(before[0]), np.inf)
        np.minimum.at(limits, mids, 0.02 * length)
        displacement = np.linalg.norm(after[0] - before[0], axis=1)
        assert np.all(displacement <= limits + 1e-12)
        fixed = np.ones(len(before[0]), dtype=bool)
        fixed[mids] = False
        np.testing.assert_array_equal(before[0][fixed], after[0][fixed])
        target = fixtures.case(name, Path(entry["fixtures"]))
        run = Path(pair[1]["directory"])
        replay = observe(
            run, name, target, *after, fixtures.periodic(name), periodic_auditor=periodic_pairs
        )
        assert replay[0]["geometry_admitted"] == pair[1]["qualified"]
        for key in ("sampled_two_sided_max", "relative_volume_error", "minimum_sampled_jacobian"):
            np.testing.assert_allclose(
                replay[0][key], pair[1]["geometry"][key], atol=1e-12, rtol=1e-12
            )
        old_bary, old_specimen = scorecard.BARY, w.specimen
        dense = (
            np.array([(i, j, 13 - i - j) for i in range(14) for j in range(14 - i)], dtype=float)
            / 13
        )
        scorecard.BARY = np.unique(np.vstack((old_bary, dense)), axis=0)
        w.specimen = lambda _: (target, None, None)
        try:
            measured = [
                w.geometry(
                    run, name, *mesh, reverse_subdivisions=(8, 16), box_projection_fallback=True
                )
                for mesh in meshes
            ]
        finally:
            scorecard.BARY, w.specimen = old_bary, old_specimen
        assert all(r["projection_failed"] == 0 for r in measured)
        result = dict(
            case=name,
            counts_and_corners_fixed=True,
            cap_coordinates_fixed=True,
            only_wall_midpoints_move=True,
            maximum_move_bound_ratio=float(np.max(displacement / limits)),
            frozen_baseline_reproduced=True,
            standard_reproduced=True,
            dense_before=measured[0],
            dense_after=measured[1],
            dense_ratio=measured[1]["sampled_two_sided_max"] / measured[0]["sampled_two_sided_max"],
            before=pair[0],
            after=pair[1],
        )
        w.save(a.run / (name + ".json"), result)
        report.append(result)
    assert all(w.digest(Path(f)) == h for f, h in resolved.items())
    w.save(
        a.run / "result.json",
        dict(
            rows=report,
            inputs=resolved,
            checks_passed=True,
            FEA_performed=False,
            reference_promoted=False,
            global_geometry_certificate=False,
        ),
    )


if __name__ == "__main__":
    main()
