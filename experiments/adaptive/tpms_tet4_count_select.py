"""Choose one count correction, then freeze repeats by count proximity only."""
import argparse
import json
import math
from pathlib import Path

import microgen_workflow as w


def returned(row):
    return bool(row and row.get("returned_mesh_verified") and row.get("elements") and row.get("nodes"))


def deviation(native, actual):
    return max(abs(math.log(native[key]/actual[key])) for key in ("elements", "nodes"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selection", type=Path, required=True)
    p.add_argument("--pilot", type=Path, required=True)
    p.add_argument("--correction", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    selection = json.loads(a.selection.read_text())
    pilot = json.loads((a.pilot/"summary.json").read_text())
    assert pilot["complete"]
    corrections = []
    if a.correction:
        corrected = json.loads((a.correction/"summary.json").read_text())
        assert corrected["complete"]
        corrections = corrected["rows"]
    output = []
    decisions = []
    for source in selection["cases"]:
        entry = dict(source)
        rows = [r for r in pilot["rows"] if r["case"] == entry["case"]]
        actual = next(r for r in rows if r["arm"] == "actual")
        native = next(r for r in rows if r["arm"] == "native")
        row = dict(case=entry["case"], chosen_native_directory=native["directory"],
                   actual_directory=actual["directory"], repeat=False)
        if not a.correction:
            if returned(native) and returned(actual) and deviation(native, actual) > math.log(1.2):
                old = native["config"]["cells"]
                factor = min((actual[k]/native[k])**(1/3) for k in ("elements", "nodes"))
                cells = max(6, min(24, round(old*factor)))
                if cells != old:
                    entry["native"] = dict(native["config"], cells=cells)
                    entry["count_correction"] = dict(source=native["directory"], actual=actual["directory"],
                        old_cells=old, new_cells=cells, rule="One cubic correction using the smaller element/node budget ratio")
                    output.append(entry)
                    row.update(action="correct-count-once", old_cells=old, new_cells=cells)
                else:
                    row.update(action="no-distinct-rounded-correction")
            else:
                row.update(action="retain-pilot-or-record-no-output")
        else:
            candidates = [native]+[r for r in corrections if r["case"] == entry["case"] and r["arm"] == "native"]
            produced = [r for r in candidates if returned(r)]
            if returned(actual) and produced:
                chosen = min(produced, key=lambda r: (deviation(r, actual), r["elements"]))
                entry["native"] = chosen["config"]
                row.update(chosen_native_directory=chosen["directory"], repeat=True,
                    elements_ratio=chosen["elements"]/actual["elements"], nodes_ratio=chosen["nodes"]/actual["nodes"],
                    count_deviation=deviation(chosen, actual), action="repeat-frozen-count-nearest-produced")
            else:
                row.update(action="no-returned-pair-no-repeated-timing")
            entry.update(row)
            output.append(entry)
        decisions.append(row)
    result = dict(selection, cases=output, decisions=decisions,
        selection_policy="Count proximity only, independent of timing, distance, quality or admission. One correction when joint log count deviation exceeds log1.2; this20percent trigger is a research scheduling rule, not an admission gate.",
        source_selection_sha256=w.digest(a.selection), pilot_sha256=w.digest(a.pilot/"summary.json"),
        correction_sha256=w.digest(a.correction/"summary.json") if a.correction else None)
    w.save(a.output, result)
    if a.correction:
        repeat = dict(result, cases=[e for e in output if e["repeat"]])
        w.save(a.output.with_name(a.output.stem+"-repeated.json"), repeat)
    print(f"Selected {len(output)} cases; actions: "+str({v:sum(d['action']==v for d in decisions) for v in sorted({d['action'] for d in decisions})}))


if __name__ == "__main__":
    main()
