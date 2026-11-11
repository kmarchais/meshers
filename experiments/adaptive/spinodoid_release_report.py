"""Report frozen old/new native and actual spinodoid workflow measurements."""
import argparse
import hashlib
import html
import json
from pathlib import Path

from tpms_tet4_count_report import collect, comparison_style, median, number


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilot", type=Path, required=True)
    p.add_argument("--correction", type=Path)
    p.add_argument("--repeated", type=Path, required=True)
    p.add_argument("--selection", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    data = collect(a)
    rows, chart, findings = [], [], []
    for case in data["cases"]:
        before, native, actual = [case["arms"][key] for key in ("before", "native", "actual")]
        assert before["all_admitted"] and native["all_admitted"]
        assert all(arm["returned"] == arm["attempts"] == 5 for arm in (before,native,actual))
        assert not actual["all_admitted"]
        metric = lambda arm, key: median(arm[key])
        change = {key:100*(metric(native,key)/metric(before,key)-1) for key in ("seconds", "tets", "nodes")}
        findings.append(dict(case=case["case"], change_percent=change,
            admitted_native_before_after=before["all_admitted"] and native["all_admitted"],
            actual_admitted=actual["all_admitted"],
            quality_minimum_change=median(native["quality"]["min"])-median(before["quality"]["min"])))
        rows.append(f'<tr><th colspan="10">{html.escape(case["case"])} · change in tets {change["tets"]:.1f}%, nodes {change["nodes"]:.1f}%, time {change["seconds"]:.1f}%</th></tr>')
        for arm_name, label in (("before", "Previous Meshers"), ("native", "Improved Meshers"), ("actual", "Actual Microgen+MMG, sampling32")):
            arm = case["arms"][arm_name]
            cells = []
            for key, higher in (("tets",False),("nodes",False),("seconds",False),("quality.min",True),
                ("quality.median",True),("reverse",False),("two_sided",False),("volume_error_percent",False)):
                def value(obj):
                    for part in key.split("."):
                        obj = obj[part]
                    return median(obj)
                av, bv = value(arm), value(before)
                verdict = None
                if key == "volume_error_percent" and arm_name == "native":
                    nv, ov = metric(native,"volume"),metric(before,"volume")
                    differences = [abs(nv-r)-abs(ov-r) for r in case["reference"]["interval"]]
                    if min(differences) < 0 < max(differences):
                        verdict = "reference-sensitive"
                style = comparison_style(av,bv,higher,verdict) if arm_name == "native" else ""
                cells.append(f'<td style="{style}">{number(av)}</td>')
            passed = arm["all_admitted"]
            rows.append(f'<tr><th>{label}</th>'+''.join(cells)+f'<td class="{"pass" if passed else "fail"}">{"PASS" if passed else "REJECT"}</td></tr>')
        times = [metric(case["arms"][key],"seconds") for key in ("before","native","actual")]
        for label, seconds, color in zip(("Previous Meshers","Improved Meshers","Actual sampling32"),times,("chart-2","chart-1","chart-3"),strict=True):
            chart.append(f'<div class="bar"><span>{html.escape(case["case"])} · {label}</span><div><i style="width:{100*seconds/max(times):.4f}%;background:var(--{color})"></i></div><b>{seconds:.4g}s</b></div>')
    data["findings"] = findings
    data["inputs"][str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    data["interpretation"] = "Before/after native comparisons use admitted Tet4 meshes. Actual sampling32 workflows remain rejected and support only raw comparisons. Native minimum quality decreases versus the previous finer recipes. No new FEA result, reference, production-default change or universal superiority claim."
    a.output.with_suffix(".json").write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    page = '''<!doctype html><html><head><meta charset="utf-8"><style>
body{margin:0;font-family:var(--font-sans);color:var(--foreground)}section{padding:16px;background:var(--card);border:1px solid var(--border);border-radius:var(--radius);margin-bottom:12px}h3{font-size:16px;margin:0 0 10px}p{font-size:12px;color:var(--muted-foreground);line-height:1.5}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;min-width:680px;font-size:11px;font-variant-numeric:tabular-nums}td,th{padding:8px 3px;border-bottom:1px solid var(--border);text-align:right}th:first-child{text-align:left}th[colspan]{padding-top:16px}.pass{color:var(--success)}.fail{color:var(--destructive)}.bar{display:grid;grid-template-columns:220px 1fr 65px;gap:8px;align-items:center;margin:14px 0;font-size:12px}.bar div{background:var(--secondary);height:10px}i{display:block;height:10px}b{text-align:right;font-size:11px}@media(max-width:450px){.bar{grid-template-columns:140px 1fr 60px}}
</style></head><body><section><h3>Spinodoid Tet4 release candidate</h3>
<p>Fresh same-host measurements, one numerical thread. Each route has five interleaved measured runs after one excluded warmup. Both native recipes preserve the frozen Fourier fields and pass distance0.002, absolute volume1%, quality0.1, positive Jacobians, coherent facets and complete periodic matching.</p>
<p>Colors compare improved Meshers with previous Meshers. Green is better, red is worse; amber volume ranking depends on the frozen empirical reference interval. The actual Microgen route uses sampling32 and default public MMG on the same FFT16 realization. Both actual routes remain rejected, so their timings do not establish admitted matched-accuracy speedups.</p>
<div class="scroll"><table><tr><th>Route</th><th>Tets</th><th>Nodes</th><th>Generation s</th><th>q min</th><th>q median</th><th>Reverse distance</th><th>Full sampled distance</th><th>Volume error %</th><th>Gates</th></tr>'''+''.join(rows)+'''</table></div>
<p>The improved recipes use fewer elements, but their minimum quality is lower than the previous native recipes and stays above0.1. Full sampled distance is absent for actual meshes rejected before that audit; the common reverse distance remains visible. Volume from invalid actual facets is diagnostic.</p></section>
<section><h3>Generation time</h3>'''+''.join(chart)+'''</section><p>The production fix restores unfinished scalar wall projections when remeshing state is rebuilt. The measured recipes also change seed resolution and interior gradation. These improvements must not all be attributed to the source fix alone. Public defaults stay unchanged. No FEA convergence or new stiffness reference is inferred.</p></body></html>'''
    a.output.with_suffix(".html").write_text(page)
    print(json.dumps(findings,indent=2))


if __name__ == "__main__":
    main()
