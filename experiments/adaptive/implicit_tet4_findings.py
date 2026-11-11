"""Join separate default and modified-route reports without pooling their timings."""
import argparse
import hashlib
import json
from pathlib import Path

from tpms_tet4_count_report import median, page


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--defaults", type=Path, required=True)
    p.add_argument("--modified", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    reports = {key: json.loads(path.read_text()) for key, path in
        (("defaults", a.defaults), ("explicit modified routes", a.modified))}
    cases, groups = [], {}
    for label, report in reports.items():
        assert report["all_orders"] == 1
        paired = [c for c in report["cases"] if all(c["arms"][arm]["returned"] for arm in ("actual", "native"))]
        metrics = {}
        for key, higher in (("seconds", False), ("quality.min", True), ("quality.median", True),
                            ("reverse", False), ("two_sided", False), ("tets", False), ("nodes", False)):
            def value(c, arm):
                result = c["arms"][arm]
                for part in key.split("."):
                    result = result[part]
                return median(result)
            observed = [c for c in paired if all(value(c, arm) is not None for arm in ("actual", "native"))]
            better = [c["case"] for c in observed if (value(c, "native") > value(c, "actual") if higher else
                value(c, "native") < value(c, "actual"))]
            worse = [c["case"] for c in observed if (value(c, "native") < value(c, "actual") if higher else
                value(c, "native") > value(c, "actual"))]
            metrics[key] = dict(observed=len(observed), native_better=better, native_worse=worse)
        dominant = []
        for case in paired:
            if not case["admitted_pair"]:
                continue
            comparisons = []
            for key, higher in (("seconds", False), ("tets", False), ("nodes", False), ("two_sided", False),
                ("volume_error_percent", False), ("quality.min", True), ("quality.median", True)):
                values = []
                for arm in ("actual", "native"):
                    value = case["arms"][arm]
                    for part in key.split("."):
                        value = value[part]
                    values.append(median(value))
                actual, native = values
                comparisons.append(native is not None and actual is not None and
                    (native >= actual if higher else native <= actual))
            if all(comparisons):
                dominant.append(case["case"])
        groups[label] = dict(route_pairs=len(report["cases"]), returned_pairs=len(paired),
            admitted_pairs=[c["case"] for c in paired if c["admitted_pair"]], raw_metrics=metrics)
        groups[label]["admitted_native_no_worse_on_seven_primary_metrics"] = dominant
        cases.extend(dict(case, comparison_group=label) for case in report["cases"])
    assert len({c["case"] for c in cases}) == 11 and len(cases) == 19
    result = dict(title="Tet4: spinodoids, composed fields and basic solids", cases=cases,
        case_count=len(cases), specimen_count=11, groups=groups, all_orders=1,
        scope="Eleven frozen specimens, nineteen separately labeled route comparisons. Defaults are retained for all eleven. Eight supplemental comparisons use explicitly modified actual routes. The two linear grades and five bounded solids require no periodic identification; both spinodoids, the ordered lattice and periodic grade retain complete periodic gates.",
        new_fea=False, admission=reports["defaults"]["admission"],
        sources={str(path.resolve()):hashlib.sha256(path.read_bytes()).hexdigest() for path in (a.defaults, a.modified)},
        limits="Raw comparisons cannot establish admitted superiority when either route rejects. Count mismatches remain explicit. No universal arbitrary-field, CAD-lattice or mapped-TPMS coverage claim; no production meshing change.")
    a.output.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    a.output.with_suffix(".html").write_text(page(result))
    print(json.dumps(groups, indent=2))


if __name__ == "__main__":
    main()
