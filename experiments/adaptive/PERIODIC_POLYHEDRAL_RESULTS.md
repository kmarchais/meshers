# Periodic affine feature recovery, 25 October 2026

The rotated square column now passes the unchanged adaptive recipe. The scalar
route previously lost longitudinal corners despite matching periodic end caps.
Compiled affine constraints now remain available when every plane is invariant
along each requested periodic axis. This improves failure coverage; it does not
establish another Microgen+MMG timing or stiffness result.

## Exact failure and repaired result

The physical column has square half-width 0.2, rotation 0.31 rad about z, and
z bounds [-0.5, 0.5]. Its exact volume is 0.16. Only z periodicity is requested.
Both scalar runs use cells 24, four adaptation rounds, one optimization round,
snap 0.2, distance 0.002, quality 0.1 and the 250,000 node/tet caps.

| Measured mesh | Original scalar | Equivalent named seed | Repaired adaptive scalar |
|---|---:|---:|---:|
| Admission | rejected | admitted | admitted |
| Tetrahedra / nodes | 65827 /15383 | 31738 /7103 | 2663 /655 |
| Minimum MMG quality | 0.120698 | 0.166498 | 0.302251 |
| Independent sampled reverse distance | 0.003073819 | 6.21e-17 | 6.21e-17 |
| Independent sampled forward distance | 0.001999665 | 5.72e-17 | 5.82e-17 |
| Relative volume error | -0.0145014% | 0 | 2.22e-14% |
| Periodic node pairs / cap triangles per end | 195 /317 | 144 /227 | 36 /51 |

The original reverse witness is exactly on a longitudinal corner line. Its
native forward estimate alone did not detect the miss. Quality, volume,
positive Jacobians, oriented facets, topology and periodic caps already passed;
none is relaxed by this fix. The named seed uses the identical domain and
shared geometry/quality settings, without claiming adaptive work. It is a
representation control, not an actual Microgen comparator.

## Implementation boundary

A compiled convex affine graph qualifies only when every plane coefficient
on every enabled periodic axis is exactly zero. No tolerance is used to infer
invariance. Unsupported fields, explicit gradients and other unsupported modes
keep their scalar route. Nonconvex affine region unions remain nonperiodic.
The native interface independently enforces invariant planes, and the original
field and size-trace periodic preflights remain active.

Nonadaptive dispatch forwards the requested periodic axes to the existing
constrained cutter. Adaptive dispatch reuses the affine seed, torus adaptation,
feature projection and periodic unfolding already present in the core. The
regressions cover both paths, exact wall geometry, cap-facet pairing, the native
noninvariant-plane guard and rejection of a nonperiodic size callback.

This does not support arbitrary periodic branch permutations or non-affine
normal jumps at the lattice's box joins. The fixed lattice, its smooth-union
order and the original Fourier specimen are unchanged.

## Verification and provenance

Jobs 259–260 preserve the original rejection and named control. Build 261 and
control 262 establish the repaired route. Suite 263 preserves a new test's wrong
expected exception type; the existing ValueError sizing rejection was correct.
Corrected suite 264 passes 205 tests with one optional skip. Core and binding
Clippy, Python lint and typing pass. Native release suite 265 passes 94 tests,
including the native invariant-plane guard. A test-only namespace error was
fixed before the passing checks; no production gate changed.

Kevin paused while 265 was running. It completed successfully, with no successor.
Kevin then explicitly resumed from the saved pause checkpoint. Final current-
source build 266 reproduces the exact measured native hash
7685573407b7276e230077215bacbbc1d1ffb3e17fc3f404e74c65d969892e52.
The dc707190 predecessor and both Python packages are preserved locally.

Fresh controls 267–270 admit Gyroid, Schwarz P, spinodoid and lattice. Independent
audit 271 replays their actual geometry, quality, Jacobians, facets and periodicity,
and proves exact equality of coordinates/connectivity/surface after canonical
relabeling to the previous admitted meshes. Fixed specimen and sample hashes
also match. Tet4 counts remain 72,757, 90,792, 72,528 and 218,239 respectively.

The same audit independently replays the rejected column, equivalent named
seed and repaired candidate from saved arrays, including source/native hashes,
exact target queries, a second forward sampling grid, volume balance, topology,
periodic node bijection and matching cap facets. It confirms the original
failure rather than reclassifying it. See
[the audit](../../handoff/macbook/periodic-polyhedral-change-audit-20261025.json)
and [native receipt](../../handoff/macbook/periodic-polyhedral-native-candidate-20261025.json).

No FE solve, reference promotion, discarded lattice model, new comparative
workflow timing or full-goal completion follows. Actual Gyroid/spinodoid
comparators and lattice convergence across all 24 models remain unresolved.
The next diagnosis should exercise periodic joins whose active affine faces
change across the boundary, preserving the measured normal jump. Another
curve-sizing sweep is not selected.
