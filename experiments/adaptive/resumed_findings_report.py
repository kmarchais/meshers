"""Collect the resumed controls without turning local feasibility into a win."""

import argparse
import hashlib
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "experiments/adaptive/runs"


def collect():
    hashes = {}

    def read(relative):
        path = RUNS / relative
        data = path.read_bytes()
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    pmy = read("macbook-pmy-residual-probe-v1/result.json")
    interpolation = read("macbook-pmy-interpolation-control-v1/result.json")
    splits = read("macbook-cad-ridge-cavity-v2/result.json")
    chain = read("macbook-cad-ridge-chain-v2/result.json")
    cad_before = read("macbook-cad-intersection-curved-v1/result.json")
    volume_probe = read("macbook-volume-residual-probe-v1/result.json")
    volume = []
    for folder in ("macbook-volume-signed-envelope-pilot-v1", "macbook-volume-signed-envelope-six-v1"):
        summary = read(folder + "/summary.json")
        assert summary["complete"] and all(row["worker"]["returncode"] == 0 for row in summary["rows"])
        for row in summary["rows"]:
            data = read(str(Path(row["directory"]).relative_to(RUNS)) + "/result.json")
            assert data["inputs_unchanged"]
            volume.append(data)
    assert len(volume) == 9 and len({row["case"] for row in volume}) == 9
    assert all(row["complete"] for row in (pmy, interpolation, splits, chain, volume_probe))
    return dict(inputs=hashes, pmy=pmy, interpolation=interpolation, cad_splits=splits,
                cad_chain=chain, cad_before=cad_before, volume_diagnosis=volume_probe,
                volume=volume, production_change=False, new_fea=False,
                report_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def page(data):
    def number(value):
        return "n/a" if value is None else f"{value:.7g}"

    def bar(label, value, scale, color):
        return (f'<div class="bar"><span>{html.escape(label)}</span><div class="track">'
                f'<i style="width:{min(100, 100*value/scale):.5f}%;background:var(--{color})"></i>'
                f'</div><b>{number(value)}</b></div>')

    pmy, preflight, chain = data["pmy"], data["interpolation"], data["cad_chain"]
    distances = [row["neighborhood"]["backgrounds"]["refined"]["boundary"]["nearest_distance"]
                 for row in pmy["witnesses"]]
    after = preflight["background"]["after"]
    cad = chain["full_admission"].get("geometry", {})
    before = data["cad_before"]["native"]["geometry"]
    quality = chain["local"]["candidate_global_quality"]
    volume_rows, exact_rows = [], []
    for row in data["volume"]:
        good = (row["admitted"] and row["dense_admitted"] and
                row["independent_standard_boundary_nonworse"] and row["independent_dense_boundary_nonworse"])
        label, color = ("Pass; boundary nonworse", "success") if good else ("Inspect regression", "destructive")
        volume_rows.append(f'<tr><td>{html.escape(row["case"])}</td><td style="color:var(--{color})">{label}</td></tr>')
        if row["reference"]["analytic_exact"]:
            old = abs(row["before"]["geometry"]["relative_volume_error"])
            new = abs(row["after"]["geometry"]["relative_volume_error"])
            exact_rows.append(f'<tr><td>{html.escape(row["case"])}</td><td>{number(old)}</td>'
                              f'<td class="good">{number(new)}</td><td>{100*(1-new/old):.1f}%</td></tr>')
    local_before = data["cad_before"]["localization"]["reverse_groups"]["ridges"]["maximum"]
    local_after = chain["local"]["frozen_witness_distance_to_local_surface"]
    status = "Full admission passes" if chain["admitted"] else "Rejected by admission"
    return f'''<!doctype html><html><head><meta charset="utf-8"><style>
body{{margin:0;color:var(--foreground);font-family:var(--font-sans,system-ui)}}
section{{background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:16px;margin-bottom:12px}}
h3{{font-size:16px;margin:0 0 9px}}p{{font-size:12px;line-height:1.5;color:var(--muted-foreground)}}
table{{width:100%;border-collapse:collapse;font-size:12px}}td,th{{text-align:left;padding:6px 3px;border-bottom:1px solid var(--border);font-variant-numeric:tabular-nums}}
.bar{{display:grid;grid-template-columns:145px minmax(20px,1fr) 82px;gap:8px;align-items:center;margin:10px 0;font-size:12px}}
.track,.track i{{height:17px;display:block}}.track{{background:var(--secondary)}}b{{text-align:right;font-variant-numeric:tabular-nums}}
.good{{color:var(--success)}}.bad{{color:var(--destructive)}}.warn{{color:var(--warning)}}
@media(max-width:450px){{.bar{{grid-template-columns:108px minmax(20px,1fr) 72px;font-size:11px}}table{{font-size:11px}}}}
</style></head><body>
<section><h3>PMY: real boundary gaps; subdivision still fails</h3>
<p>Complete curved-surface exclusion confirms all {len(pmy['witnesses'])} selected gaps. The refined background already places its affine wall too far away. Distances use unit-cell coordinates; smaller is better.</p>
{bar('Admission limit',.002,.012,'info')}
{bar('Best of 8 affine gaps',min(distances),.012,'warning')}
{bar('Worst affine gap',max(distances),.012,'destructive')}
<p>The field-only geometry rule detects all eight neighborhoods. Its center-point subdivision fits the count cap at {after['nodes']:,} nodes and {after['tets']:,} tets, but creates <span class="bad">{after['below_quality_gate']:,} tets below q = 0.1</span>. Minimum q is {after['minimum_quality']:.5f}. No native meshing or solve was launched on this rejected background.</p></section>
<section><h3>CAD: reuse an edge chain, avoid thin inserted wedges</h3>
<p>Two insertion orders failed quality and curved Jacobians. The following control reuses three existing corners and checks all 47 affected tets. <span class="{'good' if chain['admitted'] else 'bad'}">{status}.</span></p>
{bar('Saved ridge witness',local_before,.0013,'warning')}
{bar('Local chain distance',local_after,.0013,'success' if chain['admitted'] else 'destructive')}
<table><tr><th>Global metric</th><th>Before</th><th>After</th></tr>
<tr><td>Boundary maximum</td><td>{number(before['sampled_two_sided_max'])}</td><td>{number(cad.get('sampled_two_sided_max'))}</td></tr>
<tr><td>Absolute volume error</td><td>{100*abs(before['relative_volume_error']):.7g}%</td><td>{number(100*abs(cad['relative_volume_error']))+'%' if cad else 'n/a'}</td></tr>
<tr><td>Minimum corner q</td><td>{number(before['q_min_corner'])}</td><td>{number(quality['min'])}</td></tr>
<tr><td>Nodes / tets</td><td>72,793 / 49,422</td><td>{chain['local']['nodes']:,} / {chain['local']['elements']:,}</td></tr></table>
<p>The local witness is an observer, not an optimization target. The p10, median and p90 corner-quality values decrease slightly. This is a saved-mesh feasibility control; no production implementation, generation speedup or matched-cost FEA win is claimed.</p></section>
<section><h3>Volume: protect both signed error lobes</h3>
<p>The old absolute-maximum guard let a positive lobe improve while a negative lobe worsened. Separate signed bounds use the same grid and displacement budget.</p>
<table><tr><th>Exact volume case</th><th>Before |ΔV/V|</th><th>After</th><th>Reduction</th></tr>{''.join(exact_rows)}</table>
<p>Errors above are dimensionless. These gains are smaller than the aggressive control's gains. Counts, corners and corner quality stay fixed.</p>
<table><tr><th>Case</th><th>Independent standard + dense checks</th></tr>{''.join(volume_rows)}</table>
<p>Frozen volume scores worsen for linear-thickness graded Gyroid and anisotropic spinodoid. Those references are sampled, so the physical ranking remains unresolved. Finite grids do not certify continuous extrema. Prototype correction time excludes audits and is not a workflow speedup.</p></section>
<p>MacBook controls, 11 November 2026. All original specimens, all 24 lattice tensors and admission gates remain unchanged. No new FEA reference or universal Meshers versus Microgen+MMG win.</p>
</body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = collect()
    args.output.with_suffix(".json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    args.output.with_suffix(".html").write_text(page(data))


if __name__ == "__main__":
    main()
