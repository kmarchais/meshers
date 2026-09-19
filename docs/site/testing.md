# Testing and coverage

## Separate the layers

- **Rust core:** geometry, topology, periodicity, numerical invariants and errors.
- **Python API:** field preparation, callbacks, ownership, exceptions and export.
- **Integration:** independent geometry providers through the public field API.
  Record geometry failures rather than count accepted callbacks as support.

Python coverage measures Python statements and branches. It does not instrument
Rust code. Rust core coverage does not include the PyO3 binding or native JIT;
those need separate instrumented integration runs before claiming full coverage.

## Run and measure

```sh
uv sync --locked --all-packages
uv run pytest crates/meshers-python/tests --cov=meshers --cov-branch --cov-report=term-missing --cov-report=xml:coverage/python.xml
rustup component add llvm-tools-preview
cargo install cargo-llvm-cov --version 0.9.1 --locked
cargo llvm-cov --workspace --release --lcov --output-path coverage/rust.lcov
```

The Coverage workflow runs these measurements and uploads reports. Remote execution
must be verified after the workflow is pushed. No percentage is claimed from a
workflow definition, and no arbitrary global coverage threshold is imposed yet.
Establish a baseline, inspect uncovered behavior and add targeted regression tests.

## Geometry acceptance

- Check finite coordinates, valid connectivity, positive volumes and a closed
  manifold boundary for solid output.
- Compare analytic volumes or reference geometry, container containment and
  forward/inverse coordinate mappings where available.
- Check periodic correspondence only when requested and valid for the field.
- Evaluate quality distributions and refinement convergence; one coarse fixture
  cannot qualify every parameter combination in a geometry family.
- Test unsupported inputs, cancellation and callback failures alongside successful cases.

Meshing settings are part of each fixture. Disabling lattice snapping may avoid
topological pinches for composed fields, but can reduce minimum element quality.
Keep that tradeoff visible and do not suppress topology or geometry-error checks.


## Executable examples and regression gates

`uv run pytest crates/meshers-python/tests tools/tests` runs the API suite,
user examples and performance-comparison unit tests. Example exports are checked
for finite points, valid connectivity and positive volumes. The separate
[performance workflow](performance.md) compares real base/head builds and uploads
raw samples. Wheel tests run against installed artifacts outside the source tree.

Inspect each CI run's coverage reports for current results; historical test
counts and percentages are not release guarantees.
