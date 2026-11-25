---
icon: lucide/git-pull-request
---

# Build and contribute

<ol class="mesh-flow">
  <li><strong>Rust core</strong>Generate and validate the tetrahedra.</li>
  <li><strong>Python bindings</strong>Build the wheel and exercise the public API.</li>
  <li><strong>Documentation</strong>Compile the browser demo and both API references.</li>
</ol>

Install Python 3.11+, Rust stable and uv.

```sh
uv sync --locked --all-packages --group docs
```

## Check the Python bindings

```sh
uv run ruff check
uv run ruff format --check
uv run ty check
uv run pytest crates/meshers-python/tests tools/tests
```

## Check Rust and build the wheel

```sh
cargo fmt --all -- --check
cargo check --all-targets
cargo clippy --all-targets -- -D warnings
cargo test --release
cargo fmt --manifest-path crates/meshers-python/Cargo.toml -- --check
cargo clippy --manifest-path crates/meshers-python/Cargo.toml --all-targets -- -D warnings
uv build --package meshers --wheel
```

## Preview the docs

Install Node.js 22 or newer. Build the browser mesher, viewer and API references together.

```sh
rustup target add wasm32-unknown-unknown
uv run --group docs tools/build_docs.py
uv run python -m http.server --directory site 8766
```

Pull requests build a review artifact. On `main`, the publish job deploys the successful documentation build to GitHub Pages. [CI workflow](https://github.com/kmarchais/meshers/blob/main/.github/workflows/ci.yml)
