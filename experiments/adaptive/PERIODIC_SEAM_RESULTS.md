# Periodic normal-jump investigation, 25 October 2026

The invariant-column fix is verified and committed as6c843b6. The next reduced
input exposes an unresolved native failure. No native algorithm changed during
this investigation, and no stiffness/reference or timing claim follows.

## Fixed physical domain and gates

The tapered prism is the unit-box restriction of

`max(abs(x) + 0.3*abs(z) - 0.35, abs(y) - 0.2) <= 0`.

Only z periodicity is requested. Its exact volume is0.22. Opposite cap traces
are identical, but the active affine wall planes exchange across the join.
The one-sided wall-normal jump is33.398488 degrees. The specimen also has a
crease atz=0; the current failure has not yet been localized to a specific one.
Do not smooth either feature away or infer that all lattice FE error comes
from its periodic seams.

Twelve fixture controls pass. They check the closed exact surface, oriented
volume, boundary topology, paired cap triangulations, scalar trace identity,
finite-difference one-sided gradients and equivalence of the mapped domain.
A mismatched diagonal in the hand-constructed target cap was caught and
corrected before any queued mesh run. No frozen result or sample was changed.

The admission gates remain distance0.002,absolute volume1%,quality0.1,
positive Jacobians,oriented facets,periodic nodes/caps and250000nodes/tets.
All controls use cells24,snap0.2 and the same exact target queries. Failed
native calls return no mesh; no quality, cap or distance metrics are invented
for their discarded internal outputs.

## Preserved results

| Run | Declared route | Result |
|---|---|---|
| 272 | Scalar adaptive, four adaptation rounds, one optimization round | Rejects native sampled geometry error0.0106425143 |
| 273 | Equivalent mapped-column representation, one optimization round | Rejects sampled surface error0.00234375 |
| 274 | Same mapped representation, zero optimization rounds | Rejects minimum quality0.0459686167 |
| 275 | Manufactured exact-domain center-star mesh, without native generation | All mesh gates pass:13nodes,20tets,qmin0.1711701955 |

The mapped parameter column has four independent planes and physical map
`(x*(0.35-0.3*abs(z))/0.2, y, z)`. Its Jacobian determinant is positive and its
periodic end transform is exactly a unit z translation. It preserves the
negative physical domain. However, the map itself has a crease atz=0. The
one-round error is consistent with feature loss during optimization, but the
zero-round quality rejection prevents claiming this route is an admitted
native workaround. Stop this particular mapped recipe; no pass-count or grid
sweep is selected.

The separate manufactured mesh cones the exact convex-cell boundary to its
interior center. Saved-array audit276 recomputes its quality,positive
Jacobians,oriented facets,volume,boundary balance and periodic pairing. It has
four periodic node pairs and two matching cap triangles per end. Independent
sampled distance is below1e-16 and volume error is at roundoff. This proves a
geometry/quality-feasible periodic mesh exists for this input. It is not a
Meshers mesh, live workflow, reference tensor or native failure fix.

See [the independent audit](../../handoff/macbook/periodic-wedge-investigation-audit-20261025.json).
The replay also checks all three actual rejection records, their source/native
hashes, unchanged gates, identical exact queries and absence of resource-limit
failures. Runs272–274 remain rejected. The current native image is still the
verified76855734 build from the invariant-column phase.

## Mechanism boundary and next diagnostic

The compiled maximum has six affine planes. Their combined scalar traces match
on the periodic ends, but individual z-dependent planes do not match there.
The current periodic specialization therefore correctly opts out under its
invariant-plane contract. The named cutter independently checks each field's
periodic traces. Merely enabling the specialization or bypassing that check
would not implement branch correspondence or validate the cap triangulation.

Next inspect one nonperiodic six-plane seed in physical coordinates, preserving
all original periodic admission checks in the external diagnostic. Read actual
counts, quality, end-node/cap correspondence and facet labels before choosing
any native change. This is a seed preflight, not permission to turn periodicity
off in a claimed production result. It should distinguish incompatible seed
triangulations from missing image-dependent feature correspondence. No further
mapped recipe, curve-sizing sweep or FE solve is selected.

All four study geometries and all24 lattice tensors remain in scope. The last
four-case replay267–271 proves exact unchanged meshes on the same current
native image. No production source changed after it. Admitted actual Gyroid
and spinodoid comparators and lattice reference convergence remain unresolved.
