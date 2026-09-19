# Contributing

## Development setup

Use uv, Python 3.12+ for development tooling, and a current Rust stable toolchain.
The wheel continues to target Python 3.10+. Maturin builds the extension because
`uv_build` currently supports pure Python packages only.

```sh
uv sync --locked --all-packages --group docs
uv run ruff check
uv run ruff format --check
uv run ty check
uv run pytest crates/meshers-python/tests tools/tests
cargo fmt --all -- --check
cargo check --all-targets
cargo clippy --all-targets -- -D warnings
cargo test --release
cargo fmt --manifest-path crates/meshers-python/Cargo.toml -- --check
cargo clippy --manifest-path crates/meshers-python/Cargo.toml --all-targets -- -D warnings
uv build --package meshers --wheel
```

`uv.lock` pins development tools. Ruff enables `ALL` and preview. Configuration
documents exceptions for conflicting formatting rules, dynamic expression
operators, tests and downstream integration examples. The active
Python package, examples, tests and demo are checked. Historical records outside the release documentation are excluded from the
user guide. ty checks the Python package, including the FFI stubs.

Rust uses rustfmt, Cargo's type checking and Clippy with warnings denied. The core
forbids unsafe code. The extension's executable-memory and ABI operations remain
isolated with explicit safety comments.

## Build both API references

```sh
uv run --group docs python tools/build_docs.py
uv run python -m http.server --directory site 8766
```

The build script runs Zensical/mkdocstrings and rustdoc, then copies the Rust
reference into the site. `uv run --group docs zensical serve` provides guide/Python
live preview; use the full build for the bundled Rust reference. CI uploads the
site as an artifact. It does not publish automatically.
