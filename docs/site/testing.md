# Testing and coverage

<ol class="mesh-flow">
  <li><strong>Cell</strong>Check signed volume and shape.</li>
  <li><strong>Boundary</strong>Check surface labels and periodic correspondence.</li>
  <li><strong>Solve</strong>Measure a physical error against a reference.</li>
</ol>

| Layer | What to verify |
| --- | --- |
| Rust core | Geometry, topology, positive volumes, periodic correspondence and failure cases |
| Python bindings | Fields and callbacks, cancellation, exceptions, array ownership and exported connectivity |
| Integration | Retained features, boundary labels and physics on the mesh a downstream solver actually reads |

Run the public examples and regression tests:

```sh
uv sync --locked --all-packages
uv run pytest crates/meshers-python/tests tools/tests
```

## Collect coverage

```sh
uv run pytest crates/meshers-python/tests --cov=meshers --cov-branch --cov-report=term-missing --cov-report=xml:coverage/python.xml
rustup component add llvm-tools-preview
cargo install cargo-llvm-cov --version 0.9.1 --locked
cargo llvm-cov --workspace --release --lcov --output-path coverage/rust.lcov
```

The Coverage workflow uploads Python branch/statement and Rust workspace reports. The PyO3 extension and native JIT need separate instrumentation. [Mesh acceptance criteria](quality.md) · [Performance gates](performance.md)
