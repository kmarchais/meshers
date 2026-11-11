# Changelog

## 0.1.0

- Rust/Python direct tetrahedral meshing of bounded implicit functions, TPMS
  bands, callable grading and blends, with compiled expressions and gradients.
- Named multi-constraint intersections through `generate_intersection`, retaining
  surface identities and intersection curves during construction.
- Coordinate maps and verified translation/rotation periodic end correspondence,
  with explicit transforms and metadata preserved in VTKHDF.
- Integrated constrained quality optimization and strict minimum-quality gates.
- Stronger boundary vertex-link validation and periodic snapping rollback.
- Consistent analytic gradients for built-in and compiled TPMS fields; fixes the
  macOS native/compiled comparison without relaxing its test tolerances.
- Robust classification/order of coincident cuts; fixes the original Microgen
  cylindrical fixture through separate TPMS/container constraints.
- Initial manufactured linear-elasticity refinement and affine patch verification.
- Wheel/source-distribution validation, executable examples and performance regression CI.

Scope: implicit tetrahedral meshing. CAD/BRep, general remeshing, universal
minimum-quality guarantees and a solver-qualified replacement for every
Microgen VTK/MMG workflow are not claimed. The constrained path currently uses
one worker and geometry acceptance is sampled rather than certified.
