"""Build an offline dashboard from the frozen, independently audited findings.

Run with the standard library only. --check verifies the committed output is current.
This reads evidence; it never launches meshing or FE workers.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = "handoff/macbook/comparison-findings-20261026.json"
FOLLOW = "handoff/macbook/periodic-grading-findings-20261026.json"
LOOKUP = "handoff/macbook/density-lookup-findings-20261026.json"
CASES = [
    ("gyroid", "Gyroid sheet", "TPMS"),
    ("schwarz-solid", "Schwarz P solid", "TPMS"),
    ("schwarz-d-solid", "Schwarz D solid", "TPMS"),
    ("spinodoid", "Fixed spinodoid", "Spinodoids"),
    ("spinodoid-anisotropic", "Anisotropic spinodoid", "Spinodoids"),
    ("smooth-lattice", "Smooth-union lattice", "Implicit lattices"),
    ("graded-gyroid-linear", "Gyroid · linear thickness parameter", "Graded TPMS"),
    ("graded-gyroid-density-linear", "Gyroid · linear density 27–45%", "Graded TPMS"),
    ("graded-gyroid-periodic", "Gyroid · periodic grading", "Graded TPMS"),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_data():
    inputs = {}

    def read(name):
        path = ROOT / name
        inputs[str(path.relative_to(ROOT))] = digest(path)
        return json.loads(path.read_text())

    base, follow, lookup = read(BASE), read(FOLLOW), read(LOOKUP)
    assert not base["audit_receipt"]["full_objective_complete"]
    assert not follow["full_objective_complete"]
    assert not lookup["full_objective_complete"]
    for key in ("finite_difference", "automatic", "curved_fitting"):
        study = lookup[key]
        audit = ROOT / study["audit_source"]
        assert digest(audit) == study["audit_sha256"]
        inputs[study["audit_source"]] = study["audit_sha256"]
        for source in study["resolved_inputs"]:
            path = Path(source["resolved_path"])
            assert digest(path) == source["sha256"], path
            name = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
            inputs[name] = source["sha256"]
    for source in base["sources"]:
        assert digest(ROOT / source["run"] / "summary.json") == source["sha256"]
    for source in follow["checks"]:
        assert digest(ROOT / source["run"] / "summary.json") == source["summary_sha256"]

    rows = []
    for snapshot, source in (("baseline", base), ("followup", follow)):
        for index, r in enumerate(source["rows"]):
            newer = snapshot == "followup"
            geometries = [r.get("geometry") or {}] if newer else []
            directories = r["directories"] if newer else r["source_rows"]
            for directory in directories:
                raw = read(directory + "/result.json")
                if not newer:
                    geometries.append(raw.get("geometry", {}))
            good = r["admitted"] if newer else r["qualified"]
            admitted = good == r["samples"] and r["samples"] >= 5
            geometry = r.get("geometry") or {}
            med = lambda v: None if v is None else v["median"]
            quality = r["quality_quantiles"] if newer else {k: med(v) for k, v in r["quality"].items()}

            def observation(key, complete=False):
                values = [g[key] for g in geometries if g.get(key) is not None
                          and (not complete or g.get("distance_audit_complete"))]
                return statistics.median(values) if values else None

            volume = abs(geometry["relative_volume_error"]) if newer and geometry else med(r.get("abs_relative_volume"))
            superseded = (snapshot == "baseline" and r["case"] == "schwarz-d-solid" and r["config"]["arm"] == "meshers") or (
                newer and r["case"] == "graded-gyroid-density-linear" and admitted)
            rows.append(dict(
                id=f"{snapshot}-{index}", case=r["case"], name=r["name"], config=r["config"],
                arm=r["config"]["arm"], order=r["config"]["order"], snapshot=snapshot,
                historical=superseded, samples=r["samples"], passes=good, admitted=admitted,
                generation=r.get("generation_median_s") if newer else med(r["generation_s"]),
                time_range=r.get("generation_range_s") if newer else
                    ([r["generation_s"]["minimum"], r["generation_s"]["maximum"]] if r["generation_s"] else None),
                nodes=r["nodes"] if newer else med(r["nodes"]),
                elements=r["elements"] if newer else med(r["elements"]),
                distance=geometry.get("sampled_two_sided_max") if newer else med(r["sampled_distance"]),
                forward=observation("sampled_forward_max", True),
                reverse=observation("sampled_reverse_max", True),
                reverse_only=observation("sampled_reverse_max") if not admitted else None,
                volume_pct=None if volume is None else 100 * volume,
                quality=quality or {},
                failed_gates=geometry.get("failed_gates", []) if newer else r["failed_gates"],
                errors=r["errors"], source=FOLLOW if newer else BASE, directories=directories,
                native_sha256=follow["native_sha256"] if newer else follow["baseline_native_sha256"],
            ))

    for pair in lookup['automatic']['pairs']:
        stats=pair['native_lookup'];geometry=stats['geometry'];order=pair['order']
        config=next(c for c in lookup['automatic']['configs'] if c['order']==order and c.get('native_interp_auto'))
        directories=lookup['automatic']['directories'][config['name']]
        assert stats['samples']==5 and len(directories)==5
        rows.append(dict(id=f'lookup-{order}',case='graded-gyroid-density-linear',name=config['name'],config=config,
            arm='meshers',order=order,snapshot='lookup',historical=False,samples=5,passes=5,admitted=True,
            generation=stats['generation_median_s'],time_range=stats['generation_range_s'],nodes=stats['nodes'],
            elements=stats['elements'],distance=geometry['sampled_two_sided_max'],forward=geometry['sampled_forward_max'],
            reverse=geometry['sampled_reverse_max'],reverse_only=None,volume_pct=100*abs(geometry['relative_volume_error']),
            quality=stats['quality'],failed_gates=[],errors=[],source=LOOKUP,directories=directories,
            native_sha256=lookup['native_binary_sha256']))

    cases = []
    for name, label, family in CASES:
        if name == "graded-gyroid-density-linear":
            specimen = follow["density_fixture"]
        else:
            specimen = read(f"experiments/adaptive/runs/macbook-comparison-fixtures-v1/{name}/{name}-specimen.json")
        spread = specimen.get("volume_scramble_spread", 0)
        cases.append(dict(id=name, label=label, family=family,
                          periodic=name not in ("graded-gyroid-linear", "graded-gyroid-density-linear"),
                          volume_spread_pct=100 * spread / specimen["volume"],
                          definition=specimen.get("definition", name)))

    pairs = []
    for p in base["admitted_same_order_ratios"]:
        native = next(r for r in rows if r["case"] == p["case"] and r["name"] == p["native"] and r["snapshot"] == "baseline")
        actual = next(r for r in rows if r["case"] == p["case"] and r["name"] == p["actual"] and r["snapshot"] == "baseline")
        assert native["admitted"] and actual["admitted"] and native["order"] == actual["order"]
        summary = read(p["run"] + "/summary.json")
        times = {}
        for key in ("native", "actual"):
            samples = [r for r in summary["rows"] if not r["warmup"] and r["config"]["name"] == p[key]]
            assert len(samples) == 5 and all(r["qualified"] for r in samples)
            times[key] = statistics.median(r["generation_s"] for r in samples)
        assert abs(times["actual"] / times["native"] - p["speed_actual_over_native"]) < 1e-12
        pairs.append(dict(p, native_id=native["id"], actual_id=actual["id"], times=times,
                          count_gap_pct=100 * abs(native["elements"] / actual["elements"] - 1)))

    fea = base["lattice_stiffness"]["rows"]
    for r in fea:
        assert len(r["all24_model_disagreements"]) == 24
        assert r["worst_model_disagreement"] == max(r["all24_model_disagreements"])
        assert r["passes_one_percent"] == (r["worst_model_disagreement"] <= .01)
    assert len(rows) == 51 and len(pairs) == 5
    assert all(any(r["case"] == c[0] and r["order"] == order and r["admitted"] for r in rows)
               for c in CASES for order in (1, 2))
    return dict(date="2026-10-26", host="Apple M4 · 24 GiB · one numerical thread",
                cases=cases, rows=rows, pairs=pairs, fea=fea, lookup=lookup,
                inputs=inputs, sources=[BASE, FOLLOW, LOOKUP], model_hashes=base["lattice_stiffness"]["retained_model_tensor_hashes"],
                full_objective_complete=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = build_data()
    template = HERE / "comparison_dashboard.template.html"
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    html = template.read_text().replace("__DASHBOARD_DATA__", payload)
    assert "__DASHBOARD_DATA__" not in html
    output = HERE / "comparison-dashboard.html"
    if args.check:
        assert output.read_text() == html, "Dashboard is stale; rerun the builder"
    else:
        output.write_text(html)
    print(json.dumps(dict(output=str(output), rows=len(data["rows"]), pairs=len(data["pairs"]),
                          cases=len(data["cases"]), input_files=len(data["inputs"]), sha256=digest(output))))


if __name__ == "__main__":
    main()
