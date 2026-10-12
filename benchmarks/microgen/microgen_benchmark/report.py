"""A standalone HTML report, also renderable from retained JSON without compute."""

import hashlib
import json
from pathlib import Path


def render(document):
    template = Path(__file__).with_name("report.html").read_text(encoding="utf-8")
    rendered = {**document, "render_template_sha256": hashlib.sha256(template.encode()).hexdigest()}
    payload = json.dumps(rendered, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    return template.replace("__BENCHMARK_DATA__", payload)
