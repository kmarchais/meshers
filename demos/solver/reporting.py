"""Write solver measurements as newline-delimited JSON."""

import json
import sys


def write_json_line(record: object) -> None:
    """Emit one complete measurement for readers following the solve."""
    sys.stdout.write(json.dumps(record) + "\n")
    sys.stdout.flush()
