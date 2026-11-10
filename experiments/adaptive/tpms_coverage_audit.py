"""Explicit periodic-cell controls and localization of the first coverage failures."""

import argparse
from functools import partial
import json
from pathlib import Path

import microgen_workflow as w


def controls():
    import numpy as np
    from periodic_trace import periodic_pairs

    check = partial(periodic_pairs, declared_bounds=([-0.5] * 3, [0.5] * 3))
    # A closed interior tetrahedron has no material trace on any cell face.
    island = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]])
    surface = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    pairs, _, _ = check(island, surface)
    assert all(len(p) == 0 for p in pairs)
    try:
        periodic_pairs(island, surface)
    except ValueError:
        pass
    else:
        raise AssertionError("Historical full-span rule changed")
    # A one-sided material cap cannot be excused as an empty opposite trace.
    one_side = island.copy()
    one_side[[0, 1, 2], 2] = -0.5
    try:
        check(one_side, surface)
    except ValueError as error:
        assert "Unmatched periodic cap node" in str(error)
    else:
        raise AssertionError("One-sided cap accepted")
    # Extruded triangle: matching x caps and empty y/z traces.
    points = np.array(
        [
            [-0.5, 0.0, 0.0],
            [-0.5, 0.1, 0.0],
            [-0.5, 0.0, 0.1],
            [0.5, 0.0, 0.0],
            [0.5, 0.1, 0.0],
            [0.5, 0.0, 0.1],
        ]
    )
    caps = np.array([[0, 1, 2], [3, 4, 5]])
    assert check(points, caps)[2]["periodic_pair_counts"] == [3, 0, 0]
    bad = points.copy()
    bad[4, 1] += 1e-5
    try:
        check(bad, caps)
    except ValueError:
        pass
    else:
        raise AssertionError("Mismatched cap accepted")
    outside = island.copy()
    outside[0, 0] = -0.6
    try:
        check(outside, surface)
    except ValueError:
        pass
    else:
        raise AssertionError("Outside-cell node accepted")
    return dict(
        empty_opposite_traces=True,
        old_full_span_unchanged=True,
        one_sided_cap_rejected=True,
        nonempty_cap_bijection=True,
        displaced_cap_rejected=True,
        outside_cell_rejected=True,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument(
        "--cases",
        default="tpms-neovius-upper,tpms-honeycomb-upper,tpms-neovius-sheet,tpms-honeycomb-sheet",
    )
    a = p.parse_args()
    np, *_ = w.imports()
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs
    import tpms_family_fixtures as fixtures

    w.save(a.run / "controls.json", controls())
    hashes = {
        str(f): w.digest(f)
        for f in (Path(__file__), a.source / "summary.json", w.HERE / "periodic_trace.py")
    }
    rows = []
    for row in json.loads((a.source / "summary.json").read_text())["rows"]:
        if row["case"] not in a.cases.split(","):
            continue
        directory = Path(row["directory"])
        path = directory / "mesh.npz"
        if not path.exists():
            continue
        hashes[str(path)] = w.digest(path)
        with np.load(path) as mesh:
            points, cells, surface = (mesh[k] for k in ("points", "tetrahedra", "surface"))
        target = fixtures.case(row["case"], a.fixtures)
        geometry, periodic, quality = observe(
            directory,
            row["case"],
            target,
            points,
            cells,
            surface,
            True,
            periodic_auditor=partial(periodic_pairs, declared_bounds=([-0.5] * 3, [0.5] * 3)),
        )
        details = None
        if row["config"]["arm"] == "meshers":
            old = w.specimen
            w.specimen = lambda _: (target, None, None)
            try:
                details = w.geometry(
                    directory,
                    row["case"],
                    points,
                    cells,
                    surface,
                    detail_count=4,
                    reverse_subdivisions=(8, 16),
                    box_projection_fallback=True,
                )
            finally:
                w.specimen = old
        rows.append(
            dict(
                case=row["case"],
                arm=row["config"]["arm"],
                source=str(directory),
                geometry=geometry,
                periodic=periodic,
                quality=quality,
                details=details,
                source_diagnostics=row.get("diagnostics"),
                old_result_unchanged=True,
            )
        )
    assert all(w.digest(Path(f)) == h for f, h in hashes.items())
    w.save(
        a.run / "result.json",
        dict(
            rows=rows,
            inputs=hashes,
            controls_passed=True,
            scope="New coverage diagnostic with explicitly declared physical cell; historical gate and results unchanged",
        ),
    )


if __name__ == "__main__":
    main()
