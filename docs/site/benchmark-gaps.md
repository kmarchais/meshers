---
icon: lucide/goal
---

# Open benchmark targets

CI exercises selected passing recipes. The Linux coverage study of 18 November 2026 also records rejected meshes and missed solver targets. The cavity replay uses main at `fc54032b`. [Benchmark index](benchmark-index.md)

## Where to work next

| Case | Observed result | Pass condition |
| --- | --- | --- |
| **Lost cavity** · `fc54032b` | A four-cell periodic cube fills a spherical void. The mesh volume is 1.000000 versus 0.993629. The material still has one connected component. | Preserve the void and meet the 0.01 wall tolerance, or reject the unresolved geometry. Check material at the cavity centre. [Replay below](#reproduce-the-lost-cavity) |
| **Seeded spinodoid** · eight-cell smoke recipe | Tightening adaptive sizing still rejects the mesh. The 20-cell uniform control passes, but the 24-cell control falls to 0.051 minimum shape quality. | Admit the eight-cell input within 250,000 tetrahedra, sampled wall distance ≤ 0.01 in both directions and minimum shape quality ≥ 0.1 and relative volume error plus reference uncertainty ≤ 10%. Start with the failing cap-edge sample. [Witness and replay](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/README.md) · [Issue #4](https://github.com/kmarchais/meshers/issues/4) |
| **Manufactured elasticity** · Tet4 | Uniform gyroid reaches 5.78% worst energy error at 24 cells. Adaptive BCC lattice reaches 10.49% at 12 cells. The finer BCC attempts fail geometry admission before the solve. | Pass geometry admission, then reach ≤ 2% worst energy error on the same 18 displacement fields and pass the six affine patch loads. Record degrees of freedom and solve cost. [Refinement records and replay](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/README.md#refinement-at-the-2-manufactured-energy-target) |

### Reproduce the lost cavity

Run with the wheel built from `fc54032b`. The radius is 0.115 in a unit cube. The analytic material volume is `1 - 4πr³/3`. Both builds in the [periodicity review](https://github.com/kmarchais/meshers/pull/21) returned 768 tetrahedra and filled the centre.

```python
import meshers as mrs
import numpy as np


def cavity(x, y, z):
    dx = (x - 0.125 + 0.5) % 1 - 0.5
    return np.tanh(10000 * (0.115**2 - dx**2 - (y - 0.5) ** 2 - (z - 0.5) ** 2))


mesh = mrs.generate(
    cavity,
    bounds=(0, 1) * 3,
    cells=4,
    periodic=(True,) * 3,
    minimum_quality=0.1,
    expected_components=1,
    geometry_tolerance=0.01,
    threads=1,
    compile_expression=False,
)
tets = mesh.points[mesh.tetrahedra]
volume = np.linalg.det(tets[:, 1:] - tets[:, :1]).sum() / 6
print(len(tets), volume)
```

## Periodicity candidate

On Linux, [PR #21 at `b721a33f`](https://github.com/kmarchais/meshers/pull/21) returned **66/66** periodic meshes. **48/66** met every geometry limit and **47/66** did so on the requested route. Base `fc54032b` returned **37/66** and **5/66** met every limit.

| Remaining gap | What would pass |
| --- | --- |
| All three adaptive requests for the archived 50% spinodoid fall back to uniform. Both routes at 16 and 24 cells miss the dihedral-angle limit. | Keep adaptation and satisfy the full geometry audit, including angles from 10° to 165°. [Suite and baseline](https://github.com/kmarchais/meshers/tree/b721a33f/demos/solver/periodic-bench) |
| The 12-cell uniform thin gyroid improves from 66.01% to 0.53% absolute volume error but misses the 0.5% limit. Minimum dihedral falls from 10.89° to 9.05°. | Meet volume, angle and two-sided wall limits together. Report element cost alongside accuracy (6,252 → 24,176 in this run). |
| The lost cavity reports zero reverse coverage error. | Distinguish a lost wall from a wall-free solid. [Coverage check](https://github.com/kmarchais/meshers/blob/b721a33f/crates/meshers-core/src/implicit.rs#L2907) |
| The benchmark exits successfully with an empty baseline. | Require a matching baseline for every requested case, including its quality floor. [Baseline comparison](https://github.com/kmarchais/meshers/blob/b721a33f/demos/solver/periodic_bench.py#L477-L502) |

## Platform-specific comparison

The repaired Microgen Schwarz P route at `hmax=0.035` left one tetrahedron at shape quality 0.09628 in three macOS/Apple M4 runs. The floor is 0.1. Four x86-64 configurations passed. Replay the same setting and isolate the cause before attributing the difference to the OS or CPU. [Matched-cost trials and replays](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/schwarz-matched-cost-20261118/README.md)

Close a row with its original input, measured result, build and replay command. Keep the failing record for regression testing.
