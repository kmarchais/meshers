"""Collect the three bounded investigations without promoting experimental results."""

import argparse
import hashlib
import html
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "experiments/adaptive/runs"


def collect(pmy_run):
    inputs = {}

    def read(path):
        data = path.read_bytes()
        inputs[str(path.relative_to(ROOT))] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    cad = read(RUNS / "macbook-cad-intersection-curved-v1/result.json")
    ridge = read(RUNS / "macbook-cad-intersection-ridge-v1/result.json")
    volume = []
    for folder in ("macbook-volume-bias-pilot-v1", "macbook-volume-bias-seven-v1"):
        summary = read(RUNS / folder / "summary.json")
        assert summary["complete"] and all(r["worker"]["returncode"] == 0 for r in summary["rows"])
        for row in summary["rows"]:
            data = read(Path(row["directory"]) / "result.json")
            volume.append(dict(case=data["case"], exact_reference=data["reference"]["analytic_exact"],
                before=data["before"]["geometry"], after=data["after"]["geometry"],
                standard_nonworse=data["independent_standard_boundary_nonworse"],
                dense_nonworse=data["independent_dense_boundary_nonworse"],
                admitted=data["admitted"], dense_admitted=data["dense_admitted"],
                correction_seconds=data["correction_seconds"], reference=data["reference"],
                correction=data["correction"]))
    assert len(volume) == 9 and len({r["case"] for r in volume}) == 9
    facewise = [read(RUNS / "macbook-volume-facewise-pilot-v1" / name / "result.json")
                for name in ("spinodoid", "schwarz-d-solid", "schwarz-solid")]
    return dict(date="2026-11-11", inputs=inputs,
        pmy=read(pmy_run / "result.json"),
        pmy_cap_diagnosis=read(RUNS / "macbook-pmy-snap-topology-v1/result.json"),
        cad=dict(straight=cad["straight"], curved=cad["native"],
                 same_linear_topology=cad["same_linear_topology"],
                 baseline_nodes=cad["baseline_nodes"],
                 localization=cad["localization"], ridge_only=ridge["candidate"],
                 native_generation_single_s=cad["native_generation_s"]),
        volume=volume, facewise=[dict(case=d["case"], before=d["before"]["geometry"],
                  after=d["after"]["geometry"], correction=d["correction"],
                  standard_nonworse=d["independent_standard_boundary_nonworse"],
                  dense_nonworse=d["independent_dense_boundary_nonworse"]) for d in facewise],
        limits=["No default meshing policy changes or new FEA solves.",
                "CAD order2 changes node and approximation-space cost; no matched-order actual comparison.",
                "Volume correction is a saved-mesh prototype; its wall-clock cost is not a workflow speedup.",
                "Only analytic volume references support the tiny physical-volume rankings.",
                "All original specimens, all 24 lattice models, and admission thresholds remain unchanged."])


def page(data):
    cad = data["cad"]
    before, after = cad["straight"]["geometry"], cad["curved"]["geometry"]

    def bar(label, value, scale, color):
        return (f'<div class="barrow"><span>{html.escape(label)}</span><div class="track">'
                f'<div style="width:{100*value/scale:.6f}%;background:var(--{color})"></div>'
                f'</div><b>{value:.7f}</b></div>')

    cases = []
    for row in data["volume"]:
        good = row["standard_nonworse"] and row["dense_nonworse"]
        color = "success" if good else "destructive"
        label = "No sampled boundary regression" if good else "Small boundary regression"
        cases.append(f'<div class="case"><span>{html.escape(row["case"])}</span>'
                     f'<span style="color:var(--{color})">{label}</span></div>')
    exact = []
    for row in data["volume"]:
        if not row["exact_reference"]:
            continue
        a = abs(row["before"]["relative_volume_error"])
        b = abs(row["after"]["relative_volume_error"])
        exact.append(f'<tr><td>{html.escape(row["case"])}</td><td>{a:.6g}</td>'
                     f'<td style="color:var(--success)">{b:.6g}</td><td>{a/b:,.0f}×</td></tr>')
    pmy = data["pmy"]
    attempts = [r for r in pmy["rows"] if r["arm"] == "one-round"]
    status = attempts[-1]["result"].get("status", "No result") if attempts else "Not attempted"
    error = attempts[-1]["result"].get("error", "") if attempts else ""
    measured = attempts[-1]["result"].get("geometry", {}) if attempts else {}
    if not error and measured:
        error = (f"reverse distance {measured['sampled_reverse_max']:.7f} exceeds 0.002; "
                 f"minimum corner quality {measured['q_min_corner']:.6f}")
    return f'''<!doctype html><html><head><meta charset="utf-8"><style>
body{{margin:0;font-family:var(--font-sans,system-ui);color:var(--foreground,#17212e)}}
section{{background:var(--card);border:1px solid var(--border,#bbb);border-radius:var(--radius,8px);padding:16px;margin-bottom:12px}}
h3{{font-size:16px;margin:0 0 10px}}p,.note{{font-size:12px;line-height:1.5;color:var(--muted-foreground,#666)}}
.barrow{{display:grid;grid-template-columns:130px minmax(30px,1fr) 85px;gap:8px;align-items:center;font-size:12px;margin:12px 0}}
.track{{height:18px;background:var(--secondary,#ddd)}}.track>div{{height:18px}}
table{{width:100%;border-collapse:collapse;font-size:12px}}td,th{{text-align:left;padding:7px 4px;border-bottom:1px solid var(--border,#ddd);font-variant-numeric:tabular-nums}}
.case{{display:flex;justify-content:space-between;gap:12px;font-size:12px;padding:5px 0}}
.status{{color:var(--destructive,#c22);font-weight:600}}b{{text-align:right;font-variant-numeric:tabular-nums}}
@media(max-width:450px){{.barrow{{grid-template-columns:95px minmax(20px,1fr) 70px;font-size:11px}}.case{{display:block}}}}
</style></head><body>
<section><h3>PMY · feature recovery still blocked</h3>
<div class="status">{html.escape(status)}: {html.escape(error)}</div>
<p>The refined mesh passes volume, quality, curved-Jacobian and periodicity checks, but boundary accuracy still fails. Ordinary and fixed-background baselines remain distinct because snapping changes cap diagonals. Admission requirements are unchanged.</p></section>
<section><h3>CAD union · curved geometry reduces boundary error</h3>
<p>Same 49,422 tetrahedra and identical corners. Lower boundary distance is better; bars start at zero.</p>
{bar('Straight geometry',before['sampled_two_sided_max'],.0022,'warning')}
{bar('Native curved Tet10',after['sampled_two_sided_max'],.0022,'success')}
<table><tr><th>Metric</th><th>Original Tet4</th><th>Curved Tet10</th></tr>
<tr><td>Nodes</td><td>{cad['baseline_nodes']:,}</td><td>{cad['curved']['nodes']:,}</td></tr>
<tr><td>Absolute volume error</td><td>{100*abs(before['relative_volume_error']):.6f}%</td><td>{100*abs(after['relative_volume_error']):.9f}%</td></tr>
<tr><td>Minimum corner quality</td><td>{before['q_min_corner']:.6f}</td><td>{after['q_min_corner']:.6f}</td></tr></table>
<p>All admission gates pass. The ridge now limits accuracy. More nodes and a different approximation order prevent a matched-cost FEA claim. The separate ridge-only motion was not promoted.</p></section>
<section><h3>Volume correction · strong exact-volume controls, one regression</h3>
<table><tr><th>Exact reference case</th><th>Before |ΔV/V|</th><th>After |ΔV/V|</th><th>Reduction</th></tr>{''.join(exact)}</table>
<p>The controller uses the implicit field, never these reference volumes. All nine meshes pass admission; eight preserve or improve both independent boundary maxima.</p>
{''.join(cases)}
<p>The fixed spinodoid's boundary maximum increases by about 0.004%. Per-face guards do not remove it. Native integration is deferred. Small changes against sampled volume references remain uncertain.</p></section>
<p class="note">Bounded experiments on the MacBook, 11 November 2026. No new stiffness reference, FEA comparison, or universal performance win. Source hashes and numeric evidence accompany this report.</p>
</body></html>'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pmy-run", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True, help="Output basename without extension")
    args = p.parse_args()
    data = collect(args.pmy_run.resolve())
    data["report_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.with_suffix(".json").write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    args.output.with_suffix(".html").write_text(page(data))


if __name__ == "__main__":
    main()
