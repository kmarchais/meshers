"""Collect completed coverage experiments without rewriting their frozen evidence."""

import argparse
import copy
import json
from pathlib import Path

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--screen", type=Path, action="append", required=True)
    p.add_argument("--audit", type=Path, required=True)
    p.add_argument("--native-audit", type=Path, required=True)
    p.add_argument("--primitive-pairs", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    inputs = {}
    attempts = {}
    for root in a.screen:
        summary = root / "summary.json"
        data = json.loads(summary.read_text())
        assert data["complete"]
        inputs[str(summary)] = w.digest(summary)
        for row in data["rows"]:
            row = copy.deepcopy(row)
            if row["worker"]["returncode"] != 0:
                row["status"] = row["worker"].get("limit_failure") or "worker-failed"
                row["qualified"] = False
            attempts.setdefault(row["case"], []).append(row)
    audit_path = a.audit / "result.json"
    data = json.loads(audit_path.read_text())
    assert data["controls_passed"]
    inputs[str(audit_path)] = w.digest(audit_path)
    for check in data["rows"]:
        for row in attempts[check["case"]]:
            if row["directory"] == check["source"]:
                row["historical_qualified"] = row["qualified"]
                row["qualified"] = check["geometry"]["geometry_admitted"]
                row["geometry"], row["periodic"] = check["geometry"], check["periodic"]
                row["status"] = "admitted" if row["qualified"] else "admission-rejected"
                row["declared_cell_reaudit"] = str(audit_path)
    assert sum(name.startswith("tpms-") for name in attempts) == 48
    assert len(attempts) in (48, 53)
    coverage = []
    for name, rows in attempts.items():
        native = [r for r in rows if r["config"]["arm"] == "meshers"]
        actual = [r for r in rows if r["config"]["arm"] == "microgen"]
        selected = {}
        for label, candidates in (("native", native), ("actual", actual)):
            assert candidates
            selected[label] = next(
                (r for r in reversed(candidates) if r["qualified"]), candidates[-1]
            )
        coverage.append(dict(case=name, attempts=rows, **selected))
    native_path = a.native_audit / "result.json"
    native = json.loads(native_path.read_text())
    assert native["checks_passed"]
    inputs[str(native_path)] = w.digest(native_path)
    previous_path = w.ROOT / "handoff/macbook/comparison-findings-20261026.json"
    previous = json.loads(previous_path.read_text())
    inputs[str(previous_path)] = w.digest(previous_path)
    pairs_path = a.primitive_pairs / "result.json"
    pairs = json.loads(pairs_path.read_text())
    assert pairs["checks_passed"]
    inputs[str(pairs_path)] = w.digest(pairs_path)
    report = dict(
        date="2026-11-10",
        inputs=inputs,
        coverage=coverage,
        native_fitting=native["rows"],
        prior_same_order_pairs=previous["admitted_same_order_ratios"],
        prior_lattice_stiffness=previous["lattice_stiffness"],
        primitive_pairs=pairs["rows"],
        admitted={
            label: sum(r[label]["qualified"] for r in coverage) for label in ("native", "actual")
        },
        tested_TPMS_specimens=48,
        tested_new_specimens=len(coverage),
        inventoried_entries=80,
        limits=[
            "Selected native recipes, not a measured automatic fallback workflow.",
            "Actual default and explicitly modified routes remain labelled; rejected attempts are retained.",
            "New coverage contains different element orders. Only the separate audited paired comparisons support same-order ratios.",
            "Single generation observations; no new FE solves, reference promotion or stress certification.",
            "Paired empty traces use declared physical cell bounds. Old observer defaults and results are unchanged.",
            "All original specimens and 24 lattice tensors remain in scope and unchanged.",
        ],
    )
    a.output.parent.mkdir(parents=True, exist_ok=True)
    w.save(a.output.with_suffix(".json"), report)
    chart_data = json.dumps(report).replace("</", "<\\/")
    html = """<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<style>
*{box-sizing:border-box}body{margin:0;font:14px var(--font-sans,system-ui);color:var(--foreground)}
h2{font-size:18px;margin:20px 0 8px}p{line-height:1.5;margin:8px 0;color:var(--muted-foreground)}
button{background:var(--secondary);color:var(--foreground);border:1px solid var(--border);border-radius:6px;padding:7px 12px;cursor:pointer}
.scroll{overflow-x:auto;max-width:100%}table{border-collapse:collapse;font-size:12px;width:100%;white-space:nowrap}td,th{padding:7px 8px;text-align:right;border-bottom:1px solid var(--border)}td:first-child,th:first-child{text-align:left;white-space:normal;min-width:140px}th{font-weight:600;color:var(--muted-foreground)}
.good{color:var(--success)}.bad{color:var(--destructive)}.muted{color:var(--muted-foreground)}
.bars{display:grid;gap:7px}.barrow{display:grid;grid-template-columns:minmax(110px,1fr) minmax(140px,2fr) 75px;gap:9px;align-items:center;font-size:12px}.track{height:16px;background:var(--secondary);border-radius:3px;overflow:hidden}.bar{height:100%;background:var(--success)}
.stats{display:flex;gap:12px;flex-wrap:wrap}.stat{padding:12px 16px;background:var(--card);border-radius:var(--radius,8px)}.stat strong{font-size:25px;display:block}
small{font-size:12px}.active{border-color:var(--primary)}section[hidden]{display:none}
</style>
<div class="stats" id="stats"></div>
<p>Fresh measurements on the MacBook M4, one numerical thread. Admission requires boundary ≤ 0.002, absolute volume error ≤ 1%, corner q ≥ 0.1, positive sampled curved Jacobians, coherent facets and complete periodic matching where required. The 53 new specimens are separate from the earlier nine fixed cases. Counts include explicitly modified actual routes.</p>
<nav><button id="fitTab" onclick="tab('fit')">New native fitting</button> <button id="coverageTab" onclick="tab('coverage')">Expanded coverage</button> <button id="priorTab" class="active" onclick="tab('prior')">Admitted actual comparisons</button></nav>
<section id="coverage" hidden><h2>Coverage and remaining failures</h2><p>Green means admitted, red means rejected. Times belong to the selected recipe and exclude rejected attempts. Element order is shown with status. These single observations do not establish speed superiority. Number tints compare geometry and quality to their gates. Hover over a row for the recipe and failure details.</p><div class="scroll"><table id="coverageTable"></table></div></section>
<section id="fit" hidden><h2>Boundary error reduction at unchanged counts</h2><p>Positive values improve the original sampled metric. Red values regress. The native option remains off by default. Corner quality is unchanged; volume and curved quality still need their own checks.</p><div class="bars" id="fitBars"></div><h2>Fresh generation and dense verification</h2><div class="scroll"><table id="fitTable"></table></div><p>Single timing pairs, including fitting. Dense checks add degree13 surface samples; neither sample set certifies global geometry. No old FEA score transfers to a changed curved domain.</p></section>
<section id="prior"><h2>New admitted pairs, 10 November</h2><p>Five interleaved measured pairs each, after an excluded warmup. Both use Tet10 and pass all applicable gates. The actual route explicitly uses resolution24, first-pass overrides, repair and straight elevation. This compares those recipes, not the best possible Microgen settings.</p><div class="bars" id="newBars"></div><div class="scroll"><table id="newTable"></table></div><p>Native minimum, p01 and median quality improve on all three. The actual maximum q is higher on all three; its sphere q99 is also higher. All quality quantiles are retained in the numeric report.</p><h2>Earlier verified comparisons, 26 October</h2><p>Five interleaved measured pairs each, after an excluded warmup, on this Mac. Both routes pass admission. Actual routes are explicitly modified. These historical measurements remain separate from the new fitting results.</p><div class="bars" id="priorBars"></div><div class="scroll"><table id="priorTable"></table></div><h2>Lattice stiffness still behind</h2><div class="scroll"><table id="stiffnessTable"></table></div><p>The stiffness statistic retains all 24 model tensors. The actual stiffness mesh includes extra refinement and differs from the mesh-generation timing route. No speed ratio above is a matched-stiffness claim.</p></section>
<p><small>The full objective remains open. CAD strut lattices, mapped/clipped TPMS and other untested inventory entries still need equivalent physical fixtures and actual-route comparisons. FEA convergence is not filled in by these geometry tests.</small></p>
<script>const data=DATA;
const el=id=>document.getElementById(id),fmt=(v,n=3)=>v==null?'—':Number(v).toLocaleString(undefined,{maximumFractionDigits:n});
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function tab(name){for(const n of ['coverage','fit','prior']){el(n).hidden=n!==name;el(n+'Tab').classList.toggle('active',n===name)}}
el('stats').innerHTML=`<div class="stat"><strong>${data.admitted.native}/${data.tested_new_specimens}</strong>new native specimens admitted</div><div class="stat"><strong>${data.admitted.actual}/${data.tested_new_specimens}</strong>new actual specimens admitted</div><div class="stat"><strong>80</strong>inventory entries; coverage incomplete</div>`;
function heat(v,g,up=false){if(v==null)return '';let ratio=up?g/Math.max(v,1e-30):v/g;return `background:color-mix(in srgb,var(${ratio<=1?'--success':'--destructive'}) ${Math.min(30,8+22*Math.abs(1-ratio))}%,transparent)`}
function cells(r){let g=r.geometry||{},ok=r.qualified;return `<td class="${ok?'good':'bad'}">${ok?'pass':escape(r.status)} · Tet${r.config.order===2?10:4}</td><td>${fmt(r.elements,0)} / ${fmt(r.nodes,0)}</td><td>${fmt(r.generation_s)}</td><td style="${heat(g.sampled_two_sided_max,.002)}">${fmt(g.sampled_two_sided_max,6)}</td><td style="${heat(g.relative_volume_error==null?null:Math.abs(g.relative_volume_error),.01)}">${fmt(g.relative_volume_error==null?null:Math.abs(g.relative_volume_error)*100,4)}</td><td style="${heat(g.q_min_corner,.1,true)}">${fmt(g.q_min_corner,3)}</td>`}
el('coverageTable').innerHTML='<thead><tr><th rowspan="2">Microgen field / part</th><th colspan="6">Meshers, selected recipe</th><th colspan="6">Actual Microgen + MMG</th></tr><tr>'+['status','tets / nodes','seconds','boundary','volume %','q min','status','tets / nodes','seconds','boundary','volume %','q min'].map(s=>'<th>'+s+'</th>').join('')+'</tr></thead><tbody>'+data.coverage.map(r=>`<tr title="${escape([r.native.config.name,r.native.error,JSON.stringify(r.native.geometry?.failed_gates),r.actual.config.name,r.actual.error,JSON.stringify(r.actual.geometry?.failed_gates)].filter(Boolean).join('; '))}"><td>${escape(r.case.replace('tpms-',''))}</td>${cells(r.native)}${cells(r.actual)}</tr>`).join('')+'</tbody>';
const fitting=data.native_fitting.map(r=>({...r,reduction:100*(1-r.after.geometry.sampled_two_sided_max/r.before.geometry.sampled_two_sided_max)}));
el('fitBars').innerHTML=fitting.map(r=>`<div class="barrow"><span>${escape(r.case)}</span><div class="track"><div class="bar" style="width:${Math.max(1,Math.abs(r.reduction)/45*100)}%;background:var(${r.reduction>=0?'--success':'--destructive'})"></div></div><span class="${r.reduction>=0?'good':'bad'}">${fmt(r.reduction,2)}%</span></div>`).join('');
el('fitTable').innerHTML='<tr><th>Case</th><th>Before s</th><th>Fitted s</th><th>Dense error change</th><th>Fixed tets / nodes</th></tr>'+fitting.map(r=>`<tr><td>${escape(r.case)}</td><td>${fmt(r.before.generation_s)}</td><td>${fmt(r.after.generation_s)}</td><td class="${r.dense_ratio<=1?'good':'bad'}">${fmt((r.dense_ratio-1)*100,2)}%</td><td>${fmt(r.after.elements,0)} / ${fmt(r.after.nodes,0)}</td></tr>`).join('');
const label=r=>(r.case==='schwarz-solid'?'Schwarz P':'Lattice')+' '+(r.native.includes('uniform37')?'Uniform37 ':r.order===2?'curved ':'')+'Tet'+(r.order===2?'10':'4');
const maxSpeed=Math.max(...data.primitive_pairs.map(r=>r.speed_actual_over_native));
el('newBars').innerHTML=data.primitive_pairs.map(r=>`<div class="barrow"><span>${escape(r.case.replace('primitive-',''))}</span><div class="track"><div class="bar" style="width:${r.speed_actual_over_native/maxSpeed*100}%"></div></div><span>${fmt(r.speed_actual_over_native,1)}× faster</span></div>`).join('');
el('newTable').innerHTML='<tr><th>Solid</th><th>Native / actual s</th><th>Native / actual tets</th><th>Boundary error reduction</th><th>Native / actual q min</th></tr>'+data.primitive_pairs.map(r=>`<tr><td>${escape(r.case.replace('primitive-',''))}</td><td>${fmt(r.arms.meshers.median_s)} / ${fmt(r.arms.microgen.median_s)}</td><td>${fmt(r.arms.meshers.elements,0)} / ${fmt(r.arms.microgen.elements,0)}</td><td class="good">${fmt(r.distance_actual_over_native,1)}×</td><td>${fmt(r.arms.meshers.quality.min)} / ${fmt(r.arms.microgen.quality.min)}</td></tr>`).join('');
el('priorBars').innerHTML=data.prior_same_order_pairs.map(r=>`<div class="barrow"><span>${label(r)}</span><div class="track"><div class="bar" style="width:${r.speed_actual_over_native/46*100}%"></div></div><span>${fmt(r.speed_actual_over_native,2)}× faster</span></div>`).join('');
el('priorTable').innerHTML='<tr><th>Route</th><th>Actual / native tets</th><th>Actual / native distance</th><th>Native minus actual q99</th></tr>'+data.prior_same_order_pairs.map(r=>`<tr><td>${label(r)}</td><td>${fmt(r.tets_actual_over_native)}×</td><td>${fmt(r.distance_actual_over_native)}×</td><td class="${r.quality_native_minus_actual.p99>=0?'good':'bad'}">${fmt(r.quality_native_minus_actual.p99,5)}</td></tr>`).join('');
el('stiffnessTable').innerHTML='<tr><th>Saved route</th><th>Tets / nodes</th><th>Worst disagreement</th></tr>'+data.prior_lattice_stiffness.rows.map(r=>`<tr><td>${escape(r.route)}</td><td>${fmt(r.elements,0)} / ${fmt(r.nodes,0)}</td><td class="${r.worst_model_disagreement<=.01?'good':'bad'}">${fmt(r.worst_model_disagreement*100,6)}%</td></tr>`).join('');
</script>""".replace("DATA", chart_data)
    a.output.with_suffix(".html").write_text(html)


if __name__ == "__main__":
    main()
