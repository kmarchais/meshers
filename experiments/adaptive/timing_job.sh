#!/usr/bin/env bash
# Queue job for the pueue gpu group: freeze a release build of HEAD, then run timing.py
# at one thread and at every core. Results: timing-threads<N>.json, timing-results.jsonl.
set -euo pipefail
cd "$(dirname "$0")"
root=$(git rev-parse --show-toplevel)
if [ -n "$(git status --porcelain -- "$root/crates")" ]; then
  echo "crates/ has uncommitted changes; the timed build must be a commit" >&2
  exit 1
fi
sha=$(git rev-parse --short=8 HEAD)
build="builds/$sha"
threads=$(nproc)
pueue-progress 0 "building $sha (release)" || true
if [ ! -d "$build/meshers" ]; then
  wheels=$(mktemp -d)
  uv run --no-sync maturin build --release --locked --features experimental-surfaces -m "$root/crates/meshers-python/Cargo.toml" -o "$wheels"
  uv run --no-sync python -c "import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" \
    "$(ls "$wheels"/*.whl)" "$build"
  rm -rf "$wheels"
fi
export MESHERS_BUILD="$sha"
rm -f timing-results.jsonl
uv run --no-sync python timing.py --threads 1 --part 0/2
uv run --no-sync python timing.py --threads "$threads" --part 1/2
