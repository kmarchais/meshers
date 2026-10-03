# Implicit meshing integration build

This branch combines value-only field evaluation, bulk field sampling, cheaper
signal checks, inlined min/max and squared fields, scale-normalized optimization,
and protection against removing sampled interior components.

The default build preserves the eight-constraint API and tetrahedral meshing.
Optional Cargo features add `wide-constraints` (32 named constraints) and
`experimental-surfaces` (surface-only triangles). Build the Python wheel with:

```
uv tool run --from maturin maturin build --release --manifest-path crates/meshers-python/Cargo.toml --features wide-constraints,experimental-surfaces
```

`meshers.generate_surface(field, cells=24, quality=False)` provides fast extraction.
`quality=True` runs the balanced triangle improvement preset. Sampling can miss
small features. Surface checks do not certify self-intersections or FEA accuracy.

The optional triangle implementation builds on the existing local quad/triangle
work from `codex/native-quad-surface`; the overnight work adds the Python facade,
roundoff handling and benchmarked improvement preset.

The matching Microgen integration is on `codex/meshers-implicit-port-20261003`.
Select `backend="meshers"` and choose FEA volume meshing or preview surface meshing.
The benchmark report and retained build artifacts are under the research worktree
`overnight-20261003/experiments/overnight-20261003`.

The measured TPMS gains do not establish a universal implicit FEA mesher.
Ellipsoid accuracy and graded Split-P quality remain open limitations.
