"""A standalone HTML report, also renderable from retained JSON without compute."""

import json
from pathlib import Path


def render(document):
    payload = json.dumps(document, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    template = Path(__file__).with_name("report.html").read_text(encoding="utf-8")
    return template.replace("__BENCHMARK_DATA__", payload)
