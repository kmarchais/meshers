"""Build a compact, hash-linked report of completed overnight controls."""

import argparse
import hashlib
import html
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "experiments/adaptive/runs"


def collect():
    inputs = {}

    def read(relative):
        path = ROOT / relative
        content = path.read_bytes()
        inputs[relative] = hashlib.sha256(content).hexdigest()
        return json.loads(content)

    def run(name, file="result.json"):
        return read(f"experiments/adaptive/runs/{name}/{file}")

    old = run("macbook-pmy-native-background-v4", "one-round/result.json")
    new = run("macbook-pmy-rollback-native-v1", "edge/result.json")
    preflight = run("macbook-pmy-second-edge-preflight-v1")
    cad = run("macbook-cad-ridge-general-v1")["policies"]
    restoration = run("macbook-cad-ridge-quality-adaptive-v1")
    ranking = run("macbook-gyroid-volume-rankings-v1")
    trace = run("macbook-pmy-saved-trace-replay-v1")
    linear = run("macbook-pmy-linear-failure-v1")
    diagonal = run("macbook-pmy-boundary-diagonal-v1")
    bcc = run("macbook-cad-bcc-verified-v1")
    repair = run("macbook-pmy-interior-quality-repair-v3")
    coupled = run("macbook-pmy-coupled-interior-v1")
    bcc_gap = run("macbook-cad-bcc-gap-localization-v1")
    payoff = run("macbook-native-volume-probe-nine-v4")
    assert payoff["complete"] and len(payoff["rows"]) == 9
    assert all(row["complete"] for row in (old, new, preflight, restoration, ranking, trace))
    assert trace["trace_contract_verified"]

    def metrics(row):
        geometry = row["geometry"]
        return dict(nodes=row["nodes"], tets=row["elements"],
                    distance=geometry["sampled_two_sided_max"],
                    volume_error=abs(geometry["relative_volume_error"]),
                    q_min=geometry["q_min_corner"],
                    admitted=geometry["geometry_admitted"])

    volume = []
    for row in ranking["cases"]:
        scores = row["saved_mesh_scores"]
        native = next(item for item in payoff["rows"] if item["case"] == row["case"])
        volume.append(dict(case=row["case"], reference=row["estimates"]["gauss128"],
                           empirical_reference_interval=row["empirical_reference_interval"],
                           errors={"baseline": scores["baseline"]["gauss128_absolute_volume_error"],
                                   "native_fit": native["volume_evidence"]["before_absolute_error"],
                                   "signed_envelope": native["volume_evidence"]["after_absolute_error"]},
                           native_correction_verdict=native["volume_evidence"]["empirical_ranking"]["verdict"]))
    return dict(inputs=inputs, pmy=dict(before=metrics(old), after=metrics(new),
                trace_contract_verified=True, second_preflight=preflight["background"]["after"],
                witness_summary=preflight["witness_summary"],
                linear_stages=[{key: stage[key] for key in ("stage", "minimum_quality", "below_quality_gate")}
                               for stage in linear["stages"]],
                locally_feasible_diagonals=diagonal["locally_feasible_count"],
                interior_repair={key:repair[key] for key in ("initial_bad_elements", "final_bad_elements",
                    "accepted_moves", "final_quality", "invariant_guards")},
                coupled_accepted_components=coupled["accepted_components"]),
                cad=dict(before=metrics(cad["sorted"]), after=metrics(cad["minimum"]),
                         quality_before=cad["minimum"]["original_quality"],
                         quality_after=cad["minimum"]["candidate_quality"],
                         restoration_accepted=restoration["restoration_accepted"],
                         restoration_stop=restoration["stop_reason"],
                         bcc_gap={key:bcc_gap[key] for key in ("complete", "final_native_error",
                             "error_reproduction_absolute_difference", "stages")},
                         bcc={key:bcc[key] for key in ("fixture_verified", "native_admitted", "error", "maximum_field_at_cad_wall")}),
                volume=volume, native_payoff=[dict(case=row["case"],
                    seconds=row["timing"]["total_s"]["median"], admitted=row["admitted"],
                    dense_admitted=row["dense_admitted"], boundary_nonworse=row["standard_boundary_nonworse"] and row["dense_boundary_nonworse"],
                    volume_evidence=row["volume_evidence"], moved_nodes=row["native"]["repeats"][0]["moved_nodes"])
                    for row in payoff["rows"]], new_fea=False, new_workflow_timing=False,
                all_four_original_cases_retained=True, all_24_lattice_models_retained=True,
                report_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def page(data):
    def number(value):
        return f"{value:.7g}"

    def bar(label, value, scale, color):
        return (f'<div class="bar"><span>{html.escape(label)}</span><div class="track">'
                f'<i style="width:{100*value/scale:.5f}%;background:var(--{color})"></i>'
                f'</div><b>{number(value)}</b></div>')

    def cell(new, old, *, larger=False):
        better = new > old if larger else new < old
        color = "success" if better else "destructive" if new != old else "foreground"
        return f'<td style="color:var(--{color})">{number(new)}</td>'

    pmy, cad = data["pmy"], data["cad"]
    scale = max(pmy["before"]["distance"], pmy["after"]["distance"])*1.05
    pmy_rows = []
    for label, key in (("Tetrahedra", "tets"), ("Nodes", "nodes"), ("Minimum corner quality", "q_min"),
                       ("Frozen-reference volume error", "volume_error")):
        before, after = pmy["before"][key], pmy["after"][key]
        pmy_rows.append(f'<tr><td>{label}</td><td>{number(before)}</td>{cell(after,before,larger=key=="q_min")}</tr>')
    cad_rows = []
    for label, key in (("Boundary maximum", "distance"), ("Absolute relative volume error", "volume_error"),
                       ("Minimum corner quality", "q_min"), ("Nodes", "nodes"), ("Tetrahedra", "tets")):
        before, after = cad["before"][key], cad["after"][key]
        cad_rows.append(f'<tr><td>{label}</td><td>{number(before)}</td>{cell(after,before,larger=key=="q_min")}</tr>')
    for key in ("p10", "median", "p90"):
        before, after = cad["quality_before"][key], cad["quality_after"][key]
        cad_rows.append(f'<tr><td>Corner quality {key}</td><td>{number(before)}</td>{cell(after,before,larger=True)}</tr>')
    volume_rows = []
    for row in data["volume"]:
        scores = row["errors"]
        volume_rows.append(f'<tr><td>{html.escape(row["case"])}</td><td>{number(scores["baseline"])}</td>'
                           +cell(scores["native_fit"], scores["baseline"])
                           +cell(scores["signed_envelope"], scores["native_fit"])+"</tr>")
    payoff_rows = []
    for row in data["native_payoff"]:
        evidence = row["volume_evidence"]
        gain = "Reference unresolved"
        if "before_absolute_error" in evidence:
            before, after = evidence["before_absolute_error"], evidence["after_absolute_error"]
            gain = (f"{100*abs(before-after)/before:.2f}% {'smaller' if after < before else 'larger'}"
                    if before and after != before else "Unchanged")
        status = "Pass" if row["admitted"] and row["dense_admitted"] and row["boundary_nonworse"] else "Check failed"
        payoff_rows.append(f'<tr><td>{html.escape(row["case"])}</td><td>{row["seconds"]:.4f}</td><td>{gain}</td><td>{status}</td></tr>')
    return f'''<!doctype html><html><head><meta charset="utf-8"><style>
body{{margin:0;color:var(--foreground);font-family:var(--font-sans,system-ui)}}
section{{background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:16px;margin-bottom:12px}}
h3{{font-size:16px;margin:0 0 9px}}p{{font-size:12px;line-height:1.5;color:var(--muted-foreground)}}
table{{width:100%;border-collapse:collapse;font-size:12px}}td,th{{text-align:left;padding:6px 3px;border-bottom:1px solid var(--border);font-variant-numeric:tabular-nums}}
.bar{{display:grid;grid-template-columns:150px minmax(20px,1fr) 82px;gap:8px;align-items:center;margin:10px 0;font-size:12px}}
.track,.track i{{height:17px;display:block}}.track{{background:var(--secondary)}}b{{text-align:right;font-variant-numeric:tabular-nums}}
.good{{color:var(--success)}}.bad{{color:var(--destructive)}}
@media(max-width:450px){{.bar{{grid-template-columns:108px minmax(20px,1fr) 72px;font-size:11px}}table{{font-size:11px}}}}
</style></head><body>
<section><h3>PMY: snapping recovered, boundary accuracy still fails</h3>
<p>Restoring invalid snap moves lets the supplied background finish native meshing. Exact Rust input bits verify the restoration. Raw audited distances use unit-cell coordinates, smaller is better.</p>
{bar('Before refinement',pmy['before']['distance'],scale,'destructive')}
{bar('Edge pass + recovery',pmy['after']['distance'],scale,'warning')}
{bar('Admission limit',.002,scale,'success')}
<p>The 57% boundary reduction uses 25.6% more tets, with lower minimum quality. This is not a matched-count win. Other admission checks pass, but <span class="bad">geometry remains rejected.</span></p>
<table><tr><th>PMY tradeoff</th><th>Before</th><th>After</th></tr>{''.join(pmy_rows)}</table>
<p>Severity-first refinement improves the background estimate but leaves 22 final linear tets below q=0.1. A morning interior-only repair removes four failures; 18 remain and the worst quality stays 0.071331. All four subsequent coupled repairs fail and are rejected. The next hypothesis is an interior connectivity pass, which the current nonadaptive path lacks.</p></section>
<section><h3>CAD ridge correction: better geometry, mixed quality</h3>
<p>Sixteen field-selected ridge motions preserve counts and pass ordinary admission. This is a saved-mesh control, with no workflow speed or FEA claim. Green means better than the previous column; red means worse.</p>
<table><tr><th>Crossed-cylinder metric</th><th>Before</th><th>After</th></tr>{''.join(cad_rows)}</table>
<p>Strict sorted-quality restoration remained unsuccessful after bounded cavity growth. That research guard is separate from ordinary admission. The new actual BCC fixture passes independent CAD checks, but native linear meshing rejects boundary error 0.005759 against 0.002. Morning snapshots reproduce that error exactly: adaptation introduces off-wall corners and closing refinement worsens the error. Every failing sample is away from locked cap axes. Junction handling and wall-corner drift are the next targets; fixed-corner Tet10 elevation cannot remove those corner residuals.</p></section>
<section><h3>Volume: sharper references change two rankings</h3>
<p>Absolute volume error against separate Gauss128 estimates. Existing native fitting is compared with its baseline; the signed correction is compared with the fitted mesh. Empirical reference intervals support these ranking directions, without certified coverage. Frozen historical scores are preserved.</p>
<table><tr><th>Gyroid specimen</th><th>Baseline</th><th>Native fit</th><th>Signed correction</th></tr>{''.join(volume_rows)}</table>
<p>The signed correction remains experimental. These volume changes do not establish stiffness or local-stress accuracy.</p></section>
<section><h3>Native volume correction: measured cost and limits</h3>
<p>Added correction-stage seconds, median of three repetitions after one warmup. Direct Rust fixture evaluators, not complete production generation timings. Volume gains use exact half-cell references or the separate empirical Gyroid estimates. Other volume rankings remain unresolved. Pass requires ordinary and denser admission and no worse sampled boundary maximum.</p>
<table><tr><th>Case</th><th>Added s</th><th>Volume error</th><th>Checks</th></tr>{''.join(payoff_rows)}</table></section>
<p>11 November 2026. All original specimens, all 24 lattice tensors and admission gates remain in scope. No new FEA solve or end-to-end comparative timing appears in this report. The full Meshers versus Microgen+MMG objective remains open.</p>
</body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = collect()
    args.output.with_suffix(".json").write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    args.output.with_suffix(".html").write_text(page(data))


if __name__ == "__main__":
    main()
