"""Collect the fixed nine-case experiment without modifying frozen evidence."""

import json
from pathlib import Path

import microgen_workflow as w


def collect():
    root = w.HERE / "runs/macbook-face-fit-nine-v1"
    audit_path = w.HERE / "runs/macbook-face-fit-nine-audit-v1/result.json"
    audit = json.loads(audit_path.read_text())
    assert audit["passed"] and len(audit["checks"]) == 9
    inputs = {str(audit_path): w.digest(audit_path)}
    for filename, digest in audit["inputs"].items():
        path = Path(filename)
        candidates = [path, audit_path.parent / ("input-" + path.name)]
        candidates += [root / c["case"] / ("input-" + path.name) for c in audit["checks"]]
        path = next((p for p in candidates if p.is_file() and w.digest(p) == digest), None)
        assert path is not None, filename
        inputs[str(path)] = digest
    rows = []
    for check in audit["checks"]:
        path = root / check["case"] / "result.json"
        record = json.loads(path.read_text())
        a, b = (record[k]["geometry"] for k in ("before", "after"))
        assert record["admitted"] and record["before"]["quality"] == record["after"]["quality"]
        rows.append(
            dict(
                case=record["case"],
                source=str(path.relative_to(w.ROOT)),
                elements=record["elements"],
                nodes=record["nodes"],
                fit_seconds=record["fit_seconds"],
                boundary_before=a["sampled_two_sided_max"],
                boundary_after=b["sampled_two_sided_max"],
                boundary_reduction_percent=100
                * (1 - b["sampled_two_sided_max"] / a["sampled_two_sided_max"]),
                dense_boundary_before=check["dense_before"]["sampled_two_sided_max"],
                dense_boundary_after=check["dense_after"]["sampled_two_sided_max"],
                dense_boundary_reduction_percent=100 * (1 - check["dense_boundary_ratio"]),
                abs_relative_volume_before=abs(a["relative_volume_error"]),
                abs_relative_volume_after=abs(b["relative_volume_error"]),
                minimum_jacobian_before=a["minimum_sampled_jacobian"],
                minimum_jacobian_after=b["minimum_sampled_jacobian"],
                quality=record["after"]["quality"],
                admitted=True,
                original_admitted=a["geometry_admitted"],
            )
        )
    base_path = w.ROOT / "handoff/macbook/comparison-findings-20261026.json"
    base = json.loads(base_path.read_text())
    inputs[str(base_path)] = w.digest(base_path)
    for filename, digest in base["lattice_stiffness"]["inputs"].items():
        assert w.digest(Path(filename)) == digest
        inputs[filename] = digest
    resources = []
    for name in (
        "macbook-face-fit-gyroid-v1",
        "macbook-face-fit-gyroid-cap-v1",
        "macbook-face-fit-gyroid-p8-v1",
        "macbook-face-fit-nine-v1",
        "macbook-face-fit-nine-audit-v1",
    ):
        path = w.HERE / "runs" / name / "worker.json"
        worker = json.loads(path.read_text())
        assert worker["returncode"] == 0 and worker["limit_failure"] is None
        inputs[str(path)] = w.digest(path)
        resources.append(
            dict(run=name, seconds=worker["wall_s"], peak_gib=worker["peak_rss_bytes"] / 1024**3)
        )
    return dict(
        date="2026-11-10",
        rows=rows,
        inputs=inputs,
        resources=resources,
        lattice_comparison_unchanged=base["lattice_stiffness"],
        original_fixed_cases_included=True,
        all_24_models_retained=True,
        original_geometry_gates_unchanged=True,
        native_production_change=False,
        new_FE_solve=False,
        reference_promoted=False,
        full_objective_complete=False,
        scope="Saved-mesh geometry fitting; added fitting time is a single observation, not a repeated full-workflow benchmark. Corner quality is unchanged; volume and curved Jacobian tradeoffs remain. No Microgen comparison or continuum certificate is inferred.",
    )


def render(data):
    return """<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Curved-face fitting experiment, 10 November 2026</title>
<style>
body{font-family:var(--font-sans,system-ui);color:var(--foreground,CanvasText);margin:0;font-size:14px}
figure{margin:0}figcaption{font-size:18px;font-weight:650;margin-bottom:8px}
p{line-height:1.5;margin:8px 0 16px}.muted{color:var(--muted-foreground,GrayText)}
.controls{display:flex;gap:8px;flex-wrap:wrap;margin:14px 0}
button{font:inherit;padding:7px 10px;border:1px solid var(--border,GrayText);border-radius:var(--radius,6px);background:var(--secondary,ButtonFace);color:inherit;cursor:pointer}
button[aria-pressed=true]{outline:2px solid var(--primary,Highlight)}
.row{display:grid;grid-template-columns:minmax(135px,1fr) minmax(100px,1.4fr) 57px;gap:10px;align-items:center;min-height:48px;border-bottom:1px solid var(--border,GrayText)}
.label{font-size:13px}.count{font-size:11px;color:var(--muted-foreground,GrayText)}
.track{height:15px;border-left:1px solid var(--border,GrayText);position:relative}
.bar{height:100%;background:var(--success,green);border-radius:0 3px 3px 0}
.number{text-align:right;font-variant-numeric:tabular-nums;font-weight:600;color:var(--success,green)}
.axis{display:flex;justify-content:space-between;font-size:11px;color:var(--muted-foreground,GrayText);margin-top:8px}
details{margin-top:16px}summary{cursor:pointer;font-weight:600}table{width:100%;border-collapse:collapse;font-size:12px;margin-top:10px}
th,td{padding:8px 4px;text-align:right;border-bottom:1px solid var(--border,GrayText)}th:first-child,td:first-child{text-align:left}
.worse{color:var(--destructive,red)}.better{color:var(--success,green)}
@media(max-width:400px){.row{gap:6px;grid-template-columns:135px minmax(65px,1fr) 50px}.label{font-size:12px}}
</style>
<figure aria-label="Boundary-error reductions for nine fixed specimens">
<figcaption>Smaller boundary error, unchanged mesh counts</figcaption>
<p class="muted">Meshers seed → experimental curved-face fit. Nine fixed specimens, 10 November 2026. All fitted meshes pass admission; corner quality is unchanged.</p>
<div class="controls"><button id="standard" aria-pressed="true">Standard audit</button><button id="dense" aria-pressed="false">Independent denser audit</button></div>
<div id="chart"></div><div class="axis"><span>Reduction relative to the original mesh · higher is better</span><span>Bars: 0–45%</span></div>
<p id="scope" class="muted"></p>
</figure>
<details><summary>Volume and added-time tradeoffs</summary>
<p class="muted">Absolute volume error percentages use the retained references. Tiny differences can be below their sampling uncertainty. Added fit times are single observations, including diagnostic holdout work.</p>
<table><thead><tr><th>Specimen</th><th>Volume error %<br>before → after</th><th>Added fit s</th></tr></thead><tbody id="tradeoffs"></tbody></table>
</details>
<p class="muted">Saved-mesh geometry experiment. No fresh Microgen+MMG timing, FE solve or production change. Lattice stiffness disagreement remains 1.213247% for the best saved native experiment versus 0.994671% for the modified actual star-refined route.</p>
<script>
const rows=__DATA__;
const names={'gyroid':'Gyroid sheet','schwarz-solid':'Schwarz P solid','schwarz-d-solid':'Schwarz D solid','spinodoid':'Fixed spinodoid','spinodoid-anisotropic':'Anisotropic spinodoid','smooth-lattice':'Smooth-union lattice','graded-gyroid-linear':'Linear thickness Gyroid','graded-gyroid-periodic':'Periodic grade Gyroid','graded-gyroid-density-linear':'Linear density Gyroid'};
const ordered=Object.keys(names).map(k=>rows.find(r=>r.case===k));
function draw(dense){
 document.getElementById('standard').setAttribute('aria-pressed',String(!dense));
 document.getElementById('dense').setAttribute('aria-pressed',String(dense));
 document.getElementById('chart').innerHTML=ordered.map(r=>{
 const v=r[dense?'dense_boundary_reduction_percent':'boundary_reduction_percent'];
 return `<div class="row"><div class="label">${names[r.case]}<div class="count">${r.elements.toLocaleString()} tets${!r.original_admitted?' · admission recovered':''}</div></div><div class="track"><div class="bar" style="width:${v/45*100}%"></div></div><div class="number">${v.toFixed(2)}%</div></div>`;
 }).join('');
 document.getElementById('scope').textContent=dense?'Degree13 supplementary grid confirms improvement in all nine cases. No failed projections; the original gates remain unchanged.':'The fit uses its own degree6 grid. These standard audit samples are not used to choose node movements.';
}
document.getElementById('standard').onclick=()=>draw(false);
document.getElementById('dense').onclick=()=>draw(true);
document.getElementById('tradeoffs').innerHTML=ordered.map(r=>`<tr><td>${names[r.case]}</td><td>${(100*r.abs_relative_volume_before).toFixed(7)} → ${(100*r.abs_relative_volume_after).toFixed(7)}</td><td class="worse">+${r.fit_seconds.toFixed(3)}</td></tr>`).join('');
draw(false);
</script></html>""".replace("__DATA__", json.dumps(data["rows"]))


if __name__ == "__main__":
    data = collect()
    w.save(w.ROOT / "handoff/macbook/face-fit-findings-20261110.json", data)
    (w.ROOT / "handoff/macbook/face-fit-20261110.html").write_text(render(data))
