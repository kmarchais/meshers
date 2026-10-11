"""Copy compact benchmark evidence into version control without preview payloads."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    document = json.loads(args.source.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or "provenance" not in document:
        parser.error("Use a common benchmark results.json, not a standalone probe")
    result = {
        "evidence_schema": 1,
        "original_results_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "local_artifact_directory": args.source.parent.name,
        "scope": "Measured runs only. Counts and times include failed configurations. Local output retains source snapshots, meshes and fields.",
        **{
            key: document[key]
            for key in ["schema_version", "state", "provenance", "manifest_sha256", "elapsed_s"]
            if key in document
        },
        "cases": [{k: v for k, v in c.items() if k != "target_preview"} for c in document["cases"]],
        "rows": [{k: v for k, v in r.items() if k != "preview"} for r in document["rows"]],
        "qualification": document.get("qualification"),
    }
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    args.destination.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
