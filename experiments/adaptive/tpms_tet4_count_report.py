"""Report frozen Tet4 count selections without choosing candidates by outcomes."""
import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import statistics


ROOT = Path(__file__).resolve().parents[2]
ARMS = ("actual", "native")


def finite(value):
    return isinstance(value, (float, int)) and math.isfinite(value)


def stats(values):
    values = [float(v) for v in values if finite(v)]
    return dict(median=statistics.median(values), minimum=min(values), maximum=max(values), n=len(values)) if values else None


def identity(config):
    return json.dumps(config, sort_keys=True, separators=(",", ":"))


def collect(args):
    inputs = {}

    def read(path):
        path = Path(path).resolve()
        content = path.read_bytes()
        inputs[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = hashlib.sha256(content).hexdigest()
        return json.loads(content)

    def study(path):
        path = Path(path).resolve()
        if not (path/"summary.json").exists():
            path = path/"study"
        manifest = read(path/"manifest.json")
        summary = read(path/"summary.json")
        assert summary["complete"], f"Study is incomplete: {path}"
        for row in summary["rows"]:
            directory = Path(row["directory"])
            for name in ("result.json", "worker.json"):
                if (directory/name).exists():
                    read(directory/name)
        return manifest, summary["rows"]

    manifest, pilot = study(args.pilot)
    base = manifest["selection"]
    selected = read(args.selection) if args.selection else base
    entries = {e["case"]: e for e in selected["cases"]}
    assert set(entries) == {e["case"] for e in base["cases"]}, "Report selection must retain every pilot case"
    correction = study(args.correction)[1] if args.correction else []
    repeated = study(args.repeated)[1] if args.repeated else []
    all_rows = pilot + correction + repeated

    def arm_rows(entry, arm):
        expected = identity(entry[arm])
        pool = [r for r in repeated if r["case"] == entry["case"] and r["arm"] == arm
                and not r.get("warmup", False) and r.get("repeat", 0) >= 0
                and identity(r["config"]) == expected]
        if pool:
            return pool, "repeated"
        chosen = entry.get("chosen_native_directory" if arm == "native" else "actual_directory")
        if chosen:
            pool = [r for r in pilot + correction if Path(r["directory"]).resolve() == Path(chosen).resolve()]
            assert len(pool) == 1 and identity(pool[0]["config"]) == expected, "Selected directory/config mismatch"
            return pool, "selected pilot/correction"
        pool = [r for r in correction + pilot if r["case"] == entry["case"] and r["arm"] == arm
                and identity(r["config"]) == expected]
        assert pool, f"Missing selected attempt: {entry['case']} {arm}"
        return pool[:1], "pilot/correction failure"

    def aggregate(rows, source_kind, reference):
        # A killed/failed worker can leave seed counts or partial timers. They
        # are never reported as measurements of a returned mesh.
        returned = [r for r in rows if r.get("returned_mesh_verified", False)]
        geometry = [r.get("geometry", {}) for r in returned]
        facets_ok = bool(returned) and all(g.get("oriented_facet_consistency_passed", False) for g in geometry)
        periodic_ok = bool(returned) and all(not r.get("periodic_required", True) or (
            r.get("periodic", {}).get("complete_existing_node_classes", False)
            and r.get("periodic", {}).get("periodic_triangles_match", False) and not r.get("periodic", {}).get("error"))
            for r in returned)
        failures = sorted({str(gate) for r in rows for gate in r.get("geometry", {}).get("failed_gates", [])})
        errors = []
        for r in rows:
            if (not r.get("returned_mesh_verified", False) or r.get("resources", {}).get("returncode", 0) != 0
                    or r.get("status") not in ("admitted", "admission-rejected")):
                errors.append(dict(status=r.get("status"), error=r.get("error"),
                    limit_failure=r.get("resources", {}).get("limit_failure"), directory=r["directory"]))
        volumes = [g.get("volume") for g in geometry]
        volume_errors = [100*abs(v/reference["value"]-1) for v in volumes if finite(v)]
        return dict(source_kind=source_kind, attempts=len(rows), returned=len(returned),
            admitted=sum(bool(r.get("qualified")) for r in rows), all_admitted=bool(rows) and all(r.get("qualified", False) for r in rows),
            seconds=stats([r.get("generation_s") for r in returned]), tets=stats([r.get("elements") for r in returned]),
            nodes=stats([r.get("nodes") for r in returned]), quality={key:stats([r.get("quality_quantiles", {}).get(key) for r in returned])
                for key in ("min", "p01", "median", "max")},
            reverse=stats([g.get("sampled_reverse_max") for g in geometry]),
            two_sided=stats([g.get("sampled_two_sided_max") for g in geometry]),
            complete_distance_audits=sum(bool(g.get("distance_audit_complete")) for g in geometry),
            volume=stats(volumes), volume_error_percent=stats(volume_errors),
            volume_diagnostic_only=bool(returned) and not facets_ok, facets_pass=facets_ok, periodic_pass=periodic_ok,
            periodic_required=rows[0].get("periodic_required", True),
            failed_gates=failures, execution_failures=errors,
            source_directories=[r["directory"] for r in rows], config=rows[0]["config"])

    cases = []
    for original in base["cases"]:
        entry = entries[original["case"]]
        name = entry["case"]
        specimen_path = Path(entry["fixtures"])/name/(name+"-specimen.json")
        specimen = read(specimen_path)
        if specimen.get("volume_reference") == "analytic":
            reference = dict(value=specimen["volume"], interval=[specimen["volume"]]*2,
                spread=0., kind="analytic solid volume")
        elif name in ("schwarz-solid", "schwarz-d-solid"):
            assert specimen["volume"] == .5
            reference = dict(value=.5, interval=[.5, .5], spread=0., kind="exact half-cell",
                reason="The frozen zero-level Schwarz P/D solid has a measure-preserving sign-reversal symmetry.")
        else:
            estimates = specimen["volume_sobol_estimates"]
            reference = dict(value=specimen["volume"], interval=[min(estimates), max(estimates)],
                spread=specimen["volume_scramble_spread"], kind="frozen empirical Sobol",
                scope="Scramble range is empirical sensitivity, not a certified reference interval.")
        result = dict(case=name, fixture_kind=entry["fixture_kind"], reference=reference,
            selection_reason=entry.get("selection"), arms={})
        for arm in ARMS + (("before",) if "before" in entry else ()):
            rows, source_kind = arm_rows(entry, arm)
            assert all(r["config"]["order"] == 1 for r in rows)
            result["arms"][arm] = aggregate(rows, source_kind, reference)
        actual, native = [result["arms"][arm] for arm in ARMS]
        result["count_ratios_native_over_actual"] = {key:(native[key]["median"]/actual[key]["median"]
            if native[key] and actual[key] and actual[key]["median"] else None) for key in ("tets", "nodes")}
        verdict = "unavailable"
        if actual["volume"] and native["volume"]:
            av, nv = actual["volume"]["median"], native["volume"]["median"]
            differences = [abs(nv-r)-abs(av-r) for r in reference["interval"]]
            verdict = ("native smaller" if max(differences) < -1e-15 else "native larger"
                if min(differences) > 1e-15 else "equal" if max(map(abs, differences)) <= 1e-15
                else "reference-sensitive")
        result["volume_ranking_over_reference_interval"] = verdict
        result["admitted_pair"] = actual["all_admitted"] and native["all_admitted"]
        result["attempt_history"] = [dict(arm=r["arm"], config=r["config"], status=r.get("status"),
            returned_mesh=r.get("returned_mesh_verified", False), qualified=r.get("qualified", False),
            directory=r["directory"], warmup=r.get("warmup", False), error=r.get("error"),
            failed_gates=r.get("geometry", {}).get("failed_gates", []),
            limit_failure=r.get("resources", {}).get("limit_failure")) for r in all_rows if r["case"] == name]
        cases.append(result)
    source = Path(__file__).resolve()
    inputs[str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
    runtime_path = ROOT/selected.get("runtime_file", "handoff/macbook/tpms-tet4-runtime-20261111.json")
    runtime = read(runtime_path) if runtime_path.exists() else None
    return dict(cases=cases, inputs=inputs, runtime=runtime, case_count=len(cases), all_orders=1,
        title=selected.get("title", "Tet4 comparison across 18 fixed TPMS specimens"),
        scope=selected.get("scope", "16 frozen approximately30%-volume sheet specimens and the original Schwarz P/D half solids. The sheet Gyroid differs from the original35.6%-volume Gyroid."),
        admission=dict(distance=.002, absolute_relative_volume=.01, minimum_quality=.1, max_nodes=250000, max_tets=250000),
        all_four_original_cases_retained=True, all_24_lattice_models_retained=True, new_fea=False,
        production_defaults_changed=False, successful_speed_ratios_from_rejected_meshes=False,
        claim_limits="Colors compare raw observed metrics. A pair with either route rejected is not an admitted matched-accuracy performance win. Counts and ratios remain explicit.")


def number(value, digits=5):
    if value is None:
        return "—"
    if isinstance(value, int) or float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:.{digits}g}"


def median(stat):
    return stat["median"] if stat else None


def comparison_style(native, actual, higher=False, verdict=None):
    if native is None or actual is None or native == actual:
        return ""
    if verdict == "reference-sensitive":
        return "color:var(--warning)"
    if verdict == "equal":
        return ""
    better = native > actual if higher else native < actual
    color = "success" if better else "destructive"
    strength = min(abs(math.log(max(native, 1e-300)/max(actual, 1e-300))), 2.)/2
    return f"color:var(--{color});background:color-mix(in srgb,var(--{color}) {8+18*strength:.1f}%,transparent)"


def page(data):
    overview, details, chart = [], [], []
    for case in data["cases"]:
        actual, native = [case["arms"][arm] for arm in ARMS]
        label = case["case"].removeprefix("tpms-").replace("-", " ")
        if case["case"] == "schwarz-solid":
            label = "Schwarz P original solid"
        elif case["case"] == "schwarz-d-solid":
            label = "Schwarz D original solid"
        if case.get("comparison_group"):
            label += " · "+case["comparison_group"]
        cells = []
        specs = [("tets", False), ("nodes", False), ("seconds", False), ("quality.min", True),
                 ("quality.median", True), ("reverse", False), ("two_sided", False), ("volume_error_percent", False)]
        for key, higher in specs:
            def value(arm):
                result = arm
                for part in key.split("."):
                    result = result[part]
                return median(result)
            av, nv = value(actual), value(native)
            verdict = case["volume_ranking_over_reference_interval"] if key == "volume_error_percent" else None
            style = comparison_style(nv, av, higher, verdict)
            labels = ("ACT ", "MESH ") if key == "tets" else ("", "")
            samples = (f'<small>n={actual["seconds"]["n"] if actual["seconds"] else 0}/{actual["attempts"]}</small>',
                       f'<small>n={native["seconds"]["n"] if native["seconds"] else 0}/{native["attempts"]}</small>') if key == "seconds" else ("", "")
            cells.append(f'<td><span>{labels[0]}{number(av)}{samples[0]}</span><span style="{style}">{labels[1]}{number(nv)}{samples[1]}</span></td>')
        statuses = []
        for arm_name, arm in (("ACT", actual), ("MESH", native)):
            status = "PASS" if arm["all_admitted"] else "REJECT" if arm["returned"] else "NO MESH"
            color = "success" if arm["all_admitted"] else "destructive"
            statuses.append(f'<span style="color:var(--{color})">{arm_name} {status}<small>{arm["returned"]}/{arm["attempts"]} returned</small></span>')
        ratios = case["count_ratios_native_over_actual"]
        ratio = f'Tet ratio {number(ratios["tets"],3)}; node ratio {number(ratios["nodes"],3)}'
        source_kind = f'ACT {actual["source_kind"]}; MESH {native["source_kind"]}'
        source_kind += '; '+actual['config']['name']
        overview.append(f'<tr><th>{html.escape(label)}<small>{html.escape(ratio)}</small><small>{html.escape(source_kind)}</small></th>'+"".join(cells)+f'<td>{"".join(statuses)}</td></tr>')
        av, nv = median(actual["seconds"]), median(native["seconds"])
        if av is not None and nv is not None:
            scale = max(av, nv, 1e-12)
            tracks = "".join(f'<div class="track"><i style="width:{100*v/scale:.4f}%;background:var(--{color})"></i></div><b>{v:.4g}s</b>'
                for v, color in ((av, "chart-2"), (nv, "chart-1")))
            chart.append(f'<div class="timing"><span>{html.escape(label)}</span><div class="tracks">{tracks}</div></div>')
        else:
            chart.append(f'<div class="timing"><span>{html.escape(label)}</span><small>No two returned-mesh timings</small></div>')
        rows = []
        for title, key in (("Quality minimum", "min"), ("Quality p01", "p01"), ("Quality median", "median"), ("Quality maximum", "max")):
            rows.append(f'<tr><th>{title}</th><td>{number(median(actual["quality"][key]))}</td><td>{number(median(native["quality"][key]))}</td></tr>')
        for title, key in (("Common reverse distance", "reverse"), ("Two-sided sampled maximum", "two_sided"), ("Volume discrepancy %", "volume_error_percent")):
            rows.append(f'<tr><th>{title}</th><td>{number(median(actual[key]))}</td><td>{number(median(native[key]))}</td></tr>')
        prose = []
        for arm_name, arm in (("Actual "+actual["config"]["name"], actual), ("Meshers", native)):
            timing = arm["seconds"]
            time = (f'{number(timing["median"])}s median [{number(timing["minimum"])}, {number(timing["maximum"])}], n={timing["n"]}' if timing else "no returned-mesh generation time")
            errors = "; ".join(str(x.get("limit_failure") or x.get("error") or x.get("status")) for x in arm["execution_failures"])
            diagnostic = "; volume diagnostic because oriented facets fail" if arm["volume_diagnostic_only"] else ""
            periodic = "not required" if not arm["periodic_required"] else "passes" if arm["periodic_pass"] else "fails or unobserved"
            prose.append(f'<p><strong>{arm_name}</strong>: {html.escape(time)}. Returned {arm["returned"]}/{arm["attempts"]}; admitted {arm["admitted"]}/{arm["attempts"]}. Periodic {periodic}. '
                f'Failed gates: {html.escape(", ".join(arm["failed_gates"]) or "none recorded")}{html.escape(diagnostic)}. '
                f'{html.escape(errors[:400])}</p>')
        ref = case["reference"]
        ref_text = (f'Volume reference {number(ref["value"],10)}; {ref["kind"]}; empirical scramble spread {number(ref["spread"])}. '
                    f'Native volume ranking over the reference interval: {case["volume_ranking_over_reference_interval"]}.')
        history = sorted({f'{attempt["arm"]} cells{attempt["config"].get("cells", "default")}: {attempt["status"]}'
            + (f' ({", ".join(attempt["failed_gates"])})' if attempt["failed_gates"] else "")
            for attempt in case["attempt_history"] if not attempt["qualified"]})
        history_text = "Attempt history, including superseded pilots and warmups: " + "; ".join(history)
        details.append(f'<details><summary>{html.escape(label)} · detailed metrics and failures</summary><table><tr><th>Metric</th><th>Actual</th><th>Meshers</th></tr>{"".join(rows)}</table>{"".join(prose)}<p>{html.escape(ref_text)}</p><p>{html.escape(history_text)}</p></details>')
    admitted = sum(c["admitted_pair"] for c in data["cases"])
    return f'''<!doctype html><html><head><meta charset="utf-8"><style>
body{{margin:0;font-family:var(--font-sans,system-ui);color:var(--foreground)}}section{{background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:16px;margin-bottom:12px}}
h3{{font-size:16px;margin:0 0 10px}}p,small{{font-size:12px;line-height:1.5;color:var(--muted-foreground)}}table{{width:100%;border-collapse:collapse;font-size:11px;font-variant-numeric:tabular-nums}}td,th{{text-align:right;padding:7px 3px;border-bottom:1px solid var(--border)}}th:first-child{{text-align:left}}td span{{display:block;padding:2px}}th small,.overview td small{{display:block;font-size:9px;font-weight:normal}}.scroll{{overflow-x:auto}}.overview{{min-width:660px}}.timing{{display:grid;grid-template-columns:180px 1fr;gap:10px;align-items:center;margin:12px 0;font-size:12px}}.tracks{{display:grid;grid-template-columns:1fr 65px;gap:4px 8px}}.track,i{{height:9px;display:block}}.track{{background:var(--secondary)}}b{{font-size:11px;text-align:right;font-variant-numeric:tabular-nums}}details{{border-top:1px solid var(--border);padding:9px 0}}summary{{cursor:pointer;font-size:12px}}@media(max-width:450px){{.timing{{grid-template-columns:118px 1fr}}}}
</style></head><body><section><h3>{html.escape(data['title'])}</h3>
<p>Every numeric pair shows ACT, actual Microgen+MMG, above MESH, Meshers. The actual route is named per row; modified routes remain separate from defaults. Generation n counts returned-mesh timings over measured attempts, after excluding warmups. A partial return pool is not a complete five-pair result. Green means a better raw native value; red means worse. Amber volume values depend on the empirical reference interval. Counts were selected before repeated timings; the ratios remain visible. Strictly admitted pairs: {admitted}/{len(data['cases'])}.</p>
<p>Fresh 11 November measurements on the MacBook M4 with one numerical thread. Returned pairs use five interleaved measured pairs after one excluded warmup pair. Cases without two returned meshes retain their bounded pilot failures.</p>
<p>{html.escape(data['scope'])} All meshes here are Tet4, with four nodes per tetrahedron.</p>
<div class="scroll"><table class="overview"><tr><th>Case / native-to-actual count ratios</th><th>Tets</th><th>Nodes</th><th>Generation s</th><th>q min</th><th>q median</th><th>Reverse distance</th><th>Full distance</th><th>Volume %</th><th>Strict checks</th></tr>{''.join(overview)}</table></div>
<p>Reverse distance uses the same frozen wall and cap queries in both arms; smaller is better. It is one side of the geometry audit. The strict 0.002 gate uses the full two-sided audit when earlier checks permit it. Rejected or missing meshes stay rejected. Invalid-facet volume scores are diagnostic, not qualified physical-volume accuracy.</p></section>
<section><h3>Recorded generation times</h3><p><span style="color:var(--chart-2)">Actual route</span> appears above <span style="color:var(--chart-1)">Meshers</span>. Each pair has its own bar scale. Values are medians after excluding warmups; detailed ranges and sample counts are below. Timeout worker durations are never shown as mesh generation times.</p>{''.join(chart)}</section>
<section><h3>Quality, periodicity and failure details</h3>{''.join(details)}</section>
<p>Colors compare observed raw metrics. If either route rejects, the pair does not establish an admitted matched-accuracy speedup. These attempts do not prove universal geometry support or failure. No FEA score transfers, no default is promoted, and all four original study cases and all 24 lattice models remain in scope.</p></body></html>'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilot", type=Path, required=True)
    p.add_argument("--correction", type=Path)
    p.add_argument("--repeated", type=Path)
    p.add_argument("--selection", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    data = collect(args)
    args.output.with_suffix(".json").write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    args.output.with_suffix(".html").write_text(page(data))


if __name__ == "__main__":
    main()
