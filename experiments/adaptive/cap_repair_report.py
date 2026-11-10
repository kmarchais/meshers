"""Read completed cap-repair and CAD experiments into an auditable report."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / 'experiments/adaptive/runs'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cad-modified', type=Path, action='append', default=[])
    args = parser.parse_args()
    inputs = {}

    def read(path):
        path = path.resolve()
        data = path.read_bytes()
        inputs[str(path)] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    old = read(ROOT / 'handoff/macbook/native-fit-and-coverage-20261110.json')
    before = {r['case']: r['native'] for r in old['coverage']}
    baseline_frd = read(RUNS / 'macbook-schoen-frd-eight-pass-v1/tpms-schoen-frd-upper-native-adaptive8/result.json')
    frd = read(RUNS / 'macbook-frd-native-cap-repair-v1/summary.json')
    regression = read(RUNS / 'macbook-native-cap-repair-regression-v1/summary.json')
    assert frd['complete'] and regression['complete']
    diagnostic = read(RUNS / 'macbook-frd-coupled-probe-v1/probe/result.json')
    pmy = read(RUNS / 'macbook-pmy-failure-probe-v1/result.json')
    neck = read(RUNS / 'macbook-pmy-neck-control-v1/result.json')
    barycenter = read(RUNS / 'macbook-pmy-barycenter-control-v1/result.json')
    global_preflight = read(RUNS / 'macbook-pmy-global-preclip-preflight-v1/result.json')
    cad = read(RUNS / 'macbook-cad-union-pair-v1/summary.json')
    assert cad['complete'] and pmy['complete'] and diagnostic['admitted']
    rows = [dict(case='Schoen FRD upper, 8 passes', before=baseline_frd, after=frd['rows'][0])]
    for row in regression['rows']:
        previous = before[row['case']]
        rows.append(dict(case=row['case'].removeprefix('tpms-') + ', 4 passes', before=previous, after=row))
    report = dict(inputs=inputs, native_change=rows, fixed_count_diagnostic=dict(
        elements=diagnostic['elements'], nodes=diagnostic['nodes'],
        before=diagnostic['baseline_forward']['maximum'], after=diagnostic['geometry']['sampled_two_sided_max'],
        admitted=diagnostic['admitted'], exact_native_reproduction_verified=False),
        pmy=dict(forward=pmy['forward_witnesses'], reverse=pmy['reverse_witnesses'],
                 local_star_trials=pmy['local_star_trials'], adaptive_attempt=pmy['adaptive_attempt']),
        pmy_neck=neck, pmy_barycenter=barycenter, pmy_global_preflight=global_preflight, native_commit="d49a6b4", cad=cad, cad_modified=[read(p / 'summary.json') for p in args.cad_modified],
        limits=['All original specimens and 24 lattice tensors remain unchanged and in scope.',
                'No new FEA solve or reference promotion. Mesh admission does not certify local stress.',
                'Eight-pass FRD is an explicitly selected recipe; the four-pass result remains rejected.',
                'Single fresh generation observations; no stable timing ratio or dominance claim from a rejected actual route.',
                'One CustomLattice fixture does not establish coverage of all strut lattices or mapped TPMS.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    html = '''<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<style>
body{margin:0;font:14px var(--font-sans,system-ui);color:var(--foreground)}h2{font-size:17px;margin:20px 0 10px}p{line-height:1.5;color:var(--muted-foreground)}
.good{color:var(--success)}.bad{color:var(--destructive)}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:12px}td,th{text-align:right;padding:7px;border-bottom:1px solid var(--border);white-space:nowrap}td:first-child,th:first-child{text-align:left;white-space:normal;min-width:150px}button{font:inherit;background:var(--secondary);color:var(--foreground);border:1px solid var(--border);border-radius:5px;padding:6px 10px;cursor:pointer}.bars{height:122px}small{color:var(--muted-foreground)}
</style>
<h2>FRD recovered by native cap repair</h2>
<p>A wall corner and a neighboring cap node move together. Both periodic sides share the move. The native regression uses the measured 20-tet neighborhood and preserves its connectivity and counts.</p>
<div id="bars" class="bars"></div><div class="scroll"><table id="frd"></table></div>
<p>The eight-pass result passes every unchanged admission gate. Four passes still fail at distance 0.002195. Timing values are single observations.</p>
<button onclick="document.getElementById('regression').hidden=!document.getElementById('regression').hidden">Show adaptive regression results</button>
<div id="regression" hidden class="scroll"><table id="reg"></table></div>
<h2>PMY still needs a topology repair</h2>
<p id="pmy"></p>
<h2>Exact crossed-cylinder union</h2>
<p>Same finite-cylinder union, analytic volume, independent ridge and cap samples. Tet4 on both sides. Actual CAD/Gmsh and actual CAD/Gmsh/MMG remain distinct routes.</p>
<div class="scroll"><table id="cad"></table></div><p id="cadNote"></p>
<p><small>≥ marks a reverse-only lower bound after early rejection. Distance ≤ 0.002 · absolute volume error ≤ 1% · q ≥ 0.1 · positive sampled Jacobians · coherent facets · complete periodic caps and node classes · nodes and tets ≤ 250,000. The full Microgen comparison and stiffness objective remain open.</small></p>
<script>const data=DATA;
const el=id=>document.getElementById(id),f=(x,n=6)=>x==null?'—':Number(x).toLocaleString(undefined,{maximumFractionDigits:n}),esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const row=data.native_change[0], b=row.before,a=row.after;
const heat=(v,limit,up=false)=>v==null?'':`background:color-mix(in srgb,var(${(up?v>=limit:Math.abs(v)<=limit)?'--success':'--destructive'}) 15%,transparent)`;
function cells(r){let g=r.geometry||{},dist=g.sampled_two_sided_max??g.sampled_reverse_max,prefix=g.sampled_two_sided_max==null&&dist!=null?'≥ ':'';return `<td class="${r.qualified?'good':'bad'}">${r.qualified?'pass':esc(r.status||'rejected')}</td><td style="${heat(dist,.002)}">${prefix}${f(dist)}</td><td style="${heat(g.relative_volume_error,.01)}">${f(g.relative_volume_error==null?null:Math.abs(g.relative_volume_error)*100,5)}</td><td style="${heat(g.q_min_corner,.1,true)}">${f(g.q_min_corner,4)}</td><td>${f(r.elements,0)} / ${f(r.nodes,0)}</td><td>${f(r.generation_s,4)}</td>`}
const header='<tr><th>Recipe</th><th>Admission</th><th>Boundary</th><th>Volume %</th><th>q min</th><th>Tets / nodes</th><th>Seconds</th></tr>';
el('frd').innerHTML=header+`<tr><td>Before, 8 passes</td>${cells(b)}</tr><tr><td>Native repair, 8 passes</td>${cells(a)}</tr>`;
el('reg').innerHTML=header+data.native_change.slice(1).map(r=>`<tr><td>${esc(r.case)}</td>${cells(r.after)}</tr>`).join('');
let scale=340/.003;el('bars').innerHTML=`<svg viewBox="0 0 580 118" width="100%" height="118" role="img" aria-label="FRD boundary error falls from ${b.geometry.sampled_two_sided_max} to ${a.geometry.sampled_two_sided_max}; gate 0.002"><text x="0" y="27" fill="currentColor">Before</text><text x="0" y="65" fill="currentColor">Native repair</text><rect x="110" y="10" width="${b.geometry.sampled_two_sided_max*scale}" height="24" fill="var(--destructive)"/><rect x="110" y="48" width="${a.geometry.sampled_two_sided_max*scale}" height="24" fill="var(--success)"/><line x1="${110+.002*scale}" x2="${110+.002*scale}" y1="0" y2="80" stroke="var(--warning)" stroke-dasharray="4 3"/><text x="${110+.002*scale-42}" y="103" fill="currentColor">gate 0.002</text></svg>`;
let pm=data.pmy;el('pmy').textContent=`Newton overstates the worst forward witness (${f(pm.forward[0].historic.distance)} versus a feasible foot at ${f(pm.forward[0].closest_feasible.distance)}). Nevertheless, all eight reverse witnesses exclude the entire curved mesh within radius 0.002. All eight tested full midpoint projections invert local elements. The coarse failure is real. A local preclip control adds two background samples and removes 91 false-material path samples while preserving outer facets, but its clipped elements still fail q ≥ 0.1. A general face/tet-center rule also removes all 91 false samples and keeps background q ≥ 0.253, but clipping still needs quality repair. The complete periodic background also passes at 10,881 nodes / 56,292 tets with q ≥ 0.122 and exact matching. Native material clipping and quality repair still need verification; PMY remains unadmitted.`;
let candidates=[];for(const r of data.cad.rows){if(r.direct)candidates.push(['Actual CAD/Gmsh pilot',r.direct]);if(r.final)candidates.push([r.arm==='native'?'Meshers':'Actual CAD/Gmsh/MMG',r.final]);}for(const control of data.cad_modified)for(const r of control.rows){if(r.direct)candidates.push(['Actual CAD/Gmsh, size .03',r.direct]);if(r.final)candidates.push([r.control?.mmg_requested_both_passes?'Actual + MMG, explicit tolerances':'Actual + MMG, default passes',r.final]);}
candidates=candidates.filter(([name,r],i)=>!candidates.slice(0,i).some(([n,q])=>n===name&&q.mesh_sha256===r.mesh_sha256));
el('cad').innerHTML=header+candidates.map(([name,r])=>`<tr><td>${esc(name)}</td>${cells(r)}</tr>`).join('');
el('cadNote').textContent='The default actual pilot fails. The explicit size .03 / declared-cell control admits direct CAD/Gmsh with fewer elements and better geometry than Meshers, but lower minimum quality. Its unchanged public MMG passes then fail volume. Explicit hausd/hmax in both MMG passes produces an admitted 200,057-tet mesh with better boundary and volume accuracy than the 49,422-tet native mesh, but lower minimum quality. These remain separate routes. Timings are single pilots; no matched-FEA claim follows.';
</script>'''.replace('DATA', json.dumps(report).replace('</', '<\\/'))
    args.output.with_suffix('.html').write_text(html)


if __name__ == '__main__':
    main()
