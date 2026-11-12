"""Audit saved public-backend results and render the existing comparison table."""

import argparse
import hashlib
import importlib.util
import json
import html
from pathlib import Path
import sys

import numpy as np

import microgen_workflow as w
from tpms_tet4_count_report import ARMS, page, stats, median


QUALITY = ("min", "max", "median", "mean", "p10", "p90")


def collect(run, part):
    manifest = json.loads((run / "manifest.json").read_text())
    result = json.loads((run / "result.json").read_text())
    selection = json.loads((run / "selection.json").read_text())
    assert result["complete"]
    assert manifest["script_sha256"] == w.digest(run / "driver.py")
    assert manifest["selection_sha256"] == w.digest(run / "selection.json")
    for name, digest in manifest["audit_source_sha256"].items():
        assert w.digest(run / "audit-source" / name) == digest
    source = run / "audit-source/scorecard.py"
    spec = importlib.util.spec_from_file_location("_saved_scorecard_" + w.digest(source), source)
    saved = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = saved
    spec.loader.exec_module(saved)
    tet_metrics = saved.tet_metrics

    verified = []
    for index, row in enumerate(result["rows"]):
        part.update(
            index / len(result["rows"]),
            "Replaying saved arrays and metrics",
            meshes_verified=len(verified),
        )
        directory = Path(row["directory"])
        worker = json.loads((directory / "result.json").read_text())
        for key in ("case", "arm", "qualified", "workflow_s", "arrays_sha256"):
            assert row.get(key) == worker.get(key), (directory, key)
        assert row["timing_valid"] == (
            row["resource"]["returncode"] == 0
            and row.get("workflow_s") is not None
            and abs(row["resource"]["suspension_or_clock_adjustment_s"]) <= 0.1
        )
        if row["qualified"]:
            assert row["status"] == "admitted" and row["resource"]["returncode"] == 0
            assert row["geometry"]["geometry_admitted"] and not row["geometry"]["failed_gates"]
            assert row["geometry"]["distance_audit_complete"]
            assert row["geometry"]["sampled_two_sided_max"] <= 0.002
            assert abs(row["geometry"]["relative_volume_error"]) <= 0.01
            assert row["geometry"]["q_min_corner"] >= 0.1
            assert row["geometry"]["minimum_sampled_jacobian"] > 0
            assert row["geometry"]["oriented_facet_consistency_passed"]
            assert max(row["elements"], row["nodes"]) <= 250000
            if next(e["periodic"] for e in selection["cases"] if e["case"] == row["case"]):
                assert row["periodic"]["complete_existing_node_classes"]
                assert row["periodic"]["periodic_triangles_match"]
                assert (
                    row["native_pairs_match_independent_audit"] if row["arm"] == "native" else True
                )
        if row.get("mesh_sha256"):
            mesh = directory / "mesh.npz"
            assert w.digest(mesh) == row["mesh_sha256"]
            with np.load(mesh) as data:
                points, cells, surface = (data[k] for k in ("points", "tetrahedra", "surface"))
            digest = hashlib.sha256()
            for array in (points, cells, surface):
                digest.update(str((array.shape, array.dtype.str)).encode())
                digest.update(array.tobytes())
            assert digest.hexdigest() == row["arrays_sha256"]
            assert (len(points), len(cells)) == (row["nodes"], row["elements"])
            q, _, volumes = tet_metrics(points, cells)
            np.testing.assert_allclose(
                volumes.sum(), row["geometry"]["volume"], rtol=1e-11, atol=1e-14
            )
            if row["qualified"]:
                assert volumes.min() > 0
            expected = dict(
                zip(
                    ("min", "p01", "p10", "median", "p90", "p99", "max"),
                    np.quantile(q, [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1]),
                    strict=True,
                )
            )
            expected["mean"] = q.mean()
            for key, value in expected.items():
                np.testing.assert_allclose(value, row["quality_stats"][key], rtol=1e-12, atol=1e-14)
            if row.get("audit_reused_from"):
                previous = Path(row["audit_reused_from"])
                assert w.digest(previous) == row["audit_reused_sha256"]
                old = json.loads(previous.read_text())
                for key in (
                    "arrays_sha256",
                    "fixture_hashes",
                    "audit_protocol_sha256",
                    "case",
                    "arm",
                    "periodic_required",
                    "geometry",
                    "periodic",
                    "quality_stats",
                    "qualified",
                ):
                    assert old[key] == row[key], (directory, key)
            verified.append(
                dict(
                    directory=str(directory),
                    arrays_sha256=row["arrays_sha256"],
                    quality_recomputed=True,
                    tet_volume_recomputed=True,
                )
            )

    cases = []
    for entry in selection["cases"]:
        name = entry["case"]
        specimen = json.loads(
            (Path(entry["fixtures"]) / name / (name + "-specimen.json")).read_text()
        )
        exact = specimen.get("volume_reference") == "analytic" or name in (
            "schwarz-solid",
            "schwarz-d-solid",
        )
        reference = dict(
            value=specimen["volume"],
            interval=[specimen["volume"]] * 2
            if exact
            else [min(specimen["volume_sobol_estimates"]), max(specimen["volume_sobol_estimates"])],
            spread=0 if exact else specimen["volume_scramble_spread"],
            kind="analytic or exact half-cell"
            if exact
            else "frozen empirical Sobol, not a certified interval",
        )
        arms = {}
        for arm in ARMS:
            pool = [
                r
                for r in result["rows"]
                if r["case"] == name and r["arm"] == arm and not r.get("warmup", False)
            ]
            returned = [
                r
                for r in pool
                if r.get("mesh_sha256") and r["status"] in ("admitted", "admission-rejected")
            ]
            geometry = [r["geometry"] for r in returned]
            facets = bool(returned) and all(
                g["oriented_facet_consistency_passed"] for g in geometry
            )
            arms[arm] = dict(
                source_kind="not run"
                if not pool
                else "five measured repeats"
                if len(pool) == 5
                else "single bounded screen",
                attempts=len(pool),
                returned=len(returned),
                admitted=sum(r["qualified"] for r in pool),
                all_admitted=bool(pool) and all(r["qualified"] for r in pool),
                seconds=stats([r["workflow_s"] for r in returned if r["timing_valid"]]),
                tets=stats([r["elements"] for r in returned]),
                nodes=stats([r["nodes"] for r in returned]),
                quality={k: stats([r["quality_stats"][k] for r in returned]) for k in QUALITY},
                reverse=stats([g.get("sampled_reverse_max") for g in geometry]),
                two_sided=stats([g.get("sampled_two_sided_max") for g in geometry]),
                volume=stats([g["volume"] for g in geometry]),
                volume_error_percent=stats([
                    100 * abs(g["relative_volume_error"]) for g in geometry
                ]),
                volume_diagnostic_only=bool(returned) and not facets,
                facets_pass=facets,
                periodic_required=entry["periodic"],
                periodic_pass=bool(returned)
                and all(
                    not entry["periodic"]
                    or (
                        r["periodic"].get("complete_existing_node_classes")
                        and r["periodic"].get("periodic_triangles_match")
                        and not r["periodic"].get("error")
                    )
                    for r in returned
                ),
                failed_gates=sorted({
                    g for r in pool for g in r.get("geometry", {}).get("failed_gates", [])
                }),
                execution_failures=[
                    dict(
                        status=r["status"],
                        error=r.get("error"),
                        limit_failure=r["resource"].get("limit_failure"),
                    )
                    for r in pool
                    if r["status"] not in ("admitted", "admission-rejected")
                ],
                config={
                    "name": "native-fixed-count" if arm == "native" else "actual-default-vtk-mmg",
                    **entry[arm],
                },
            )
        actual, native = (arms[a] for a in ARMS)
        verdict = "unavailable"
        if native["volume"] and actual["volume"]:
            differences = [
                abs(median(native["volume"]) - v) - abs(median(actual["volume"]) - v)
                for v in reference["interval"]
            ]
            verdict = (
                "native smaller"
                if max(differences) < -1e-15
                else "native larger"
                if min(differences) > 1e-15
                else "reference-sensitive"
            )
            if not native["facets_pass"] or not actual["facets_pass"]:
                verdict = "reference-sensitive"
        cases.append(
            dict(
                case=name,
                arms=arms,
                reference=reference,
                volume_ranking_over_reference_interval=verdict,
                admitted_pair=native["all_admitted"] and actual["all_admitted"],
                count_ratios_native_over_actual={
                    k: median(native[k]) / median(actual[k]) if native[k] and actual[k] else None
                    for k in ("tets", "nodes")
                },
                attempt_history=[
                    dict(
                        arm=r["arm"],
                        config=r["config"],
                        status=r["status"],
                        qualified=r["qualified"],
                        failed_gates=r.get("geometry", {}).get("failed_gates", []),
                    )
                    for r in result["rows"]
                    if r["case"] == name
                ],
            )
        )
    paired = [c for c in cases if all(c["arms"][a]["returned"] for a in ARMS)]
    counts = {
        key: sum(
            median(c["arms"]["native"][key]) < median(c["arms"]["actual"][key]) for c in paired
        )
        for key in ("seconds", "tets", "nodes", "reverse")
    }
    counts.update({
        f"quality_{key}": sum(
            median(c["arms"]["native"]["quality"][key])
            > median(c["arms"]["actual"]["quality"][key])
            for c in paired
        )
        for key in QUALITY
    })
    return dict(
        cases=cases,
        source_run=str(run),
        source_sha256={
            n: w.digest(run / n)
            for n in ("result.json", "manifest.json", "selection.json", "driver.py")
        },
        verified_meshes=verified,
        paired_returned=len(paired),
        native_raw_better_counts=counts,
        quality_rows=[("Quality " + k, k) for k in QUALITY],
        time_label="Constructor + backend s",
        measurement_note="Fresh 12 November M4 measurements, one numerical thread. This section retains "
        + (
            "five interleaved measured pairs after one excluded warmup pair."
            if manifest["repetitions"] == 5
            else "one bounded attempt per route per specimen. These are screening timings, not a repeated performance claim."
        ),
        timing_note="Charged cost is construction plus the selected public backend. Imports, fixed-fixture preparation and prior density calibration, independent identity/audit work and recorded instrumentation I/O are excluded.",
        new_fea=False,
        reference_promoted=False,
        all_four_original_cases_retained_in_study=True,
        all_24_lattice_models_retained=True,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen", type=Path, required=True)
    p.add_argument("--repeated", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--controls", type=Path, nargs="*", default=[])
    p.add_argument("--history", type=Path, nargs="*", default=[])
    p.add_argument("--stiffness-pair", type=Path)
    p.add_argument("--reference-diagnostics", type=Path)
    a = p.parse_args()
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    runs = [("screen", a.screen), ("repeated", a.repeated)] + [
        (f"native-control-{i}", run) for i, run in enumerate(a.controls)
    ]
    with Job(parts=[name for name, _ in runs]) as job:
        reports = []
        history = []
        for directory in a.history:
            result = json.loads((directory / "result.json").read_text())
            assert result["complete"]
            history.extend(
                {**r, "history_source_sha256": w.digest(directory / "result.json")}
                for r in result["rows"]
            )
        for name, run in runs:
            with job.part(name) as part:
                part.update(0, "Replaying saved arrays and metrics", meshes_verified=0, force=True)
                d = collect(run.resolve(), part)
                d["title"] = "Public Microgen backend comparison, " + (
                    "30 frozen coarse Tet4 recipes"
                    if name == "screen"
                    else "repeated Gyroid and adaptive spinodoids"
                )
                d["scope"] = (
                    "The screen retains original Gyroid, Schwarz P, spinodoid and ordered smooth-union lattice, 16 catalog TPMS sheets, Schwarz D solid, anisotropic spinodoid, three graded Gyroids and five primitive controls. Coarse recipes are historical fixed-count selections, not the latest best accurate Meshers recipes."
                    if name == "screen"
                    else "This follow-up changes only the two native spinodoid recipes to previously selected adaptive settings. Actual spinodoids explicitly use sampling32 with MMG defaults. The flagship coarse Gyroid remains a rejected accuracy control. All gates and caps stay unchanged."
                )
                if name.startswith("native-control"):
                    admitted = sum(c["arms"]["native"]["all_admitted"] for c in d["cases"])
                    d["title"] = f"Final native feasibility: {admitted}/{len(d['cases'])} admitted"
                    for case in d["cases"]:
                        case["attempt_history"] += [
                            dict(
                                arm=r["arm"],
                                config=r["config"],
                                status=r["status"],
                                qualified=r["qualified"],
                                failed_gates=r.get("geometry", {}).get("failed_gates", []),
                                directory=r["directory"],
                                error=r.get("error"),
                                source_sha256=r["history_source_sha256"],
                            )
                            for r in history
                            if r["case"] == case["case"]
                        ]
                    d["measurement_note"] = (
                        "Fresh 12 November M4 native-only verification, one numerical thread. "
                        + (
                            "Five measured generations after one excluded warmup per specimen. "
                            if all(c["arms"]["native"]["attempts"] == 5 for c in d["cases"])
                            else "One measured generation per specimen. "
                        )
                        + "Actual Microgen+MMG is not rerun in this section. These are feasibility and repeatability checks, not comparative performance claims."
                    )
                    d["scope"] = (
                        "Frozen final recipes cover the original four specimens, the TPMS catalog, anisotropic spinodoid, linear thickness/density grading, periodic grading and bounded primitives. The lattice explicitly uses historical min(f, 0.02); its zero sublevel set and ordered union are unchanged. Existing public Microgen primitive SDFs are used directly. Adaptive Tet4 recognition retains unit-distance box/cylinder branches; scalar defaults, bands and explicit-gradient policies stay unchanged. Earlier failures are retained in each case history. Strict study gates and caps stay unchanged; internal locked-corner roundoff handling is documented in project tracking. This is not universal geometry support: mapped/clipped TPMS and sampled OffsetGrading remain explicitly unsupported by the adapter."
                    )
                if name.startswith("native-control") and len(d["cases"]) == 1:
                    d["title"] = "Original Gyroid adaptive count control, native only"
                    d["scope"] = "One frozen original Gyroid sheet specimen using the predeclared catalog adaptive recipe. This is a separate native count tradeoff, not a replacement for the 30-case selection or the frozen film metrics. Both meshes pass mesh admission; the separate six-load stiffness comparison does not establish FEA equivalence."
                reports.append(d)
                job.result(
                    name,
                    verified_meshes=len(d["verified_meshes"]),
                    paired_returned=d["paired_returned"],
                    admitted_pairs=sum(c["admitted_pair"] for c in d["cases"]),
                )
        stiffness = None
        stiffness_html = ""
        if a.stiffness_pair:
            run = a.stiffness_pair
            assessment = json.loads((run / "assessment.json").read_text())
            assert assessment["complete"] and assessment["new_fea"]
            assert not assessment["reference_promoted"] and not assessment["reference_accuracy_assessed"]
            for path, digest in assessment["source_inputs"].items():
                assert w.digest(Path(path)) == digest
            records = []
            for name in assessment["cases"]:
                row = json.loads((run / name / "result.json").read_text())
                audit = json.loads((run / name / "field-audit.json").read_text())
                positive = json.loads((run / name / "positive-quadrature-audit.json").read_text())
                assert w.digest(run / name / "fields.npz") == row["fields_sha256"]
                assert row["completed_loads"] == 6 and row["field_audit_passed"] and row["positive_quadrature_passed"]
                assert max(audit[k] for k in ("reaction_relative_error", "hill_mandel_relative_error", "periodic_jump_max")) < 1e-7
                records.append(dict(name=name, elements=row["elements"], nodes=row["nodes"],
                    saved_mesh_solve_s=row["saved_mesh_solve_s"], stiffness=row["stiffness"],
                    field_audit=audit, positive_quadrature=positive))
            sys.path.insert(0, str(w.ROOT.parent / "periodic-microgen/experiments/periodic-fea"))
            from periodic import tensor_difference
            difference = tensor_difference(np.array(records[1]["stiffness"]), np.array(records[0]["stiffness"]))
            np.testing.assert_allclose(difference["relative_energy_difference_eigenvalues"],
                assessment["pairwise_mixed_load_difference"]["relative_energy_difference_eigenvalues"], rtol=1e-12)
            stiffness = dict(assessment=assessment, records=records,
                assessment_sha256=w.digest(run / "assessment.json"))
            worst = 100 * difference["worst_relative_strain_energy_difference"]
            table_rows = "".join(
                f"<tr><td>{html.escape(r['name'])}</td><td>{r['elements']:,}</td><td>{r['nodes']:,}</td>"
                f"<td>{r['saved_mesh_solve_s']:.3f}</td><td>{r['field_audit']['periodic_jump_max']:.3g}</td></tr>"
                for r in records)
            stiffness_html = (
                '<section id="stiffness"><h2>Gyroid Tet4 stiffness tradeoff</h2>'
                f'<p>The admitted adaptive mesh uses 52.1% fewer tetrahedra, but its full six-load stiffness differs by <strong>{worst:.3f}%</strong> from the uniform mesh across mixed strains. '
                'This exceeds the study 1% target. Neither mesh is a reference for the other. The domains are different admitted polyhedral approximations and the spaces are not nested.</p>'
                '<table><thead><tr><th>Saved native mesh</th><th>Tets</th><th>Nodes</th><th>Six-load solve seconds, single run</th><th>Periodic field jump</th></tr></thead><tbody>'
                + table_rows + '</tbody></table><p>Both pass independent reaction/energy, positive-quadrature and periodic-field checks. '
                'Solve times exclude generation and audits and are not repeated workflow timings. No continuum reference, local-stress claim or reference promotion follows.</p></section>')
        reference_diagnostics = None
        reference_html = ""
        if a.reference_diagnostics:
            path = a.reference_diagnostics / "assessment.json"
            reference_diagnostics = json.loads(path.read_text())
            assert reference_diagnostics["complete"] and not reference_diagnostics["reference_promoted"]
            assert not reference_diagnostics["spinodoid_reference_assessed"]
            for source, digest in reference_diagnostics["inputs"].items():
                assert w.digest(Path(source)) == digest
            rows_html = ""
            for name, record in reference_diagnostics["records"].items():
                low, high = (100 * x for x in record["declared_model_range"])
                color = "#62e3a6" if record["target_passes_all"] else "#ff8792"
                rows_html += f'<tr><td>{html.escape(name)}</td><td style="color:{color}">{low:.3f}–{high:.3f}%</td><td style="color:{color}">{"PASS" if record["target_passes_all"] else "FAIL"}</td></tr>'
            reference_html = (
                '<section id="reference"><h2>Retained stiffness reference checks, native Tet4</h2>'
                '<p>Full mixed-load tensor differences. Gyroid and Schwarz P use the previously qualified curved Tet10 16/20 reference evidence, with orders 2–6 and two length proxies declared before assessment. '
                'Lattice retains all 24 historical model tensors. These empirical reference assumptions are not continuum bounds. No reference is promoted.</p>'
                '<table><thead><tr><th>Saved native mesh</th><th>Disagreement across retained models</th><th>1% target under every model</th></tr></thead><tbody>'
                + rows_html + '<tr><td>spinodoid</td><td>Unassessed: no justified reference selected</td><td>UNKNOWN</td></tr></tbody></table>'
                '<p>Passing meshing gates and FE consistency does not imply stiffness accuracy. Actual Microgen+MMG has no new admitted FEA comparator here. Generation counts, timing and quality tables remain mesh-only evidence.</p></section>')
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.with_suffix(".json").write_text(
            json.dumps(
                {"sections": reports, "stiffness_tradeoff": stiffness, "reference_diagnostics": reference_diagnostics, "new_fea": stiffness is not None, "reference_promoted": False},
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        sections = [
            f'<section id="section-{i}">' + page(d).split("<body>", 1)[1].split("</body>", 1)[0] + "</section>"
            for i, d in enumerate(reports)
        ]
        template = page(reports[0]).split("<body>", 1)[0]
        palette = ":root{--foreground:#e4edf7;--card:#142438;--border:#30445b;--radius:10px;--muted-foreground:#b7c8d9;--success:#62e3a6;--destructive:#ff8792;--warning:#ffd16b;--secondary:#263950;--chart-1:#52d5b8;--chart-2:#8eafff}body{background:#091421!important;padding:24px;max-width:1700px;margin:auto!important}h1{font-size:25px}a{color:#8eafff}"
        a.output.with_suffix(".html").write_text(
            template.replace(
                "</head>",
                '<title>Meshers through Microgen, verified 12 November</title><meta name="viewport" content="width=device-width,initial-scale=1"><style>'
                + palette
                + "</style></head>",
            )
            + "<body><h1>Meshers through Microgen</h1><p>Fresh public API measurements. Raw metric colors do not imply a qualified winner. Expand a case for all six quality statistics, periodicity, failures and reference limits. Meshing sections assess mesh admission only. The separate Gyroid stiffness check compares two saved native meshes; no stiffness reference is promoted.</p>"
            + '<nav><a href="#section-0">Coarse comparison</a> · <a href="#section-1">Repeated comparisons</a> · <a href="#section-2">Final 30-case native results</a> · <a href="#stiffness">Gyroid stiffness tradeoff</a> · <a href="#reference">Stiffness reference gaps</a></nav>'
            + '<p style="padding:16px;background:#142438;border:1px solid #30445b;border-radius:10px"><strong style="color:#62e3a6">Final native wheel: 30/30 specimens admitted, 150/150 measured generations pass.</strong> <span style="color:#ffd16b">Coarse paired screen: 0 admitted matched-accuracy pairs. Broad comparative speedup remains unproven.</span></p>'
            + reference_html
            + stiffness_html
            + "".join(sections)
            + "</body></html>"
        )


if __name__ == "__main__":
    main()
