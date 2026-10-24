# Recovered Tet4 coverage, 24 October 2026

Kevin prioritized failure coverage before optimization. The first resumed check
reproduced the lattice uniform32 rejection at sampled distance 0.004628834.
Quality, volume, positive Jacobians, oriented facets and periodicity passed.

The recent comparison omitted historical admitted adaptive32 recipes for both
lattice and spinodoid. Restoring those recipes gives three admitted generations
each on the Mac. No native code, fixed specimen, threshold or budget changed.
The older failing settings remain rejected.

Lattice returns 218,239 tetrahedra and 48,077 nodes. Sampled distance is
0.001768701, minimum quality 0.171902 and signed volume error -0.429123%.
Spinodoid returns 72,528 tetrahedra and 16,732 nodes. Distance is 0.001850419,
minimum quality 0.232844 and signed volume error +0.044600%.
Every geometry, quality, Jacobian, facet and periodic admission gate passes.

Both use a 32-cell seed, Tet4, adaptation, native geometry tolerance 0.002,
adaptation tolerance 0.002, minimum length 0.002, maximum length 0.1,
gradation 0.3, four adaptation rounds, quality guard 0.3 and one optimization
pass. The unchanged limits are 250,000 nodes/tets, 8 GiB and one numerical
thread. Spinodoid uses its exact frozen compiled 52-mode field. Lattice uses
the original sequential smooth union with the explicit positive-void clamp
0.02. The clamp preserves the material domain and zero set. Independent
auditing still uses the original unclamped field and frozen samples.

These recipes differ from the rejected settings in multiple parameters.
Recovery alone does not isolate a single cause. Code inspection confirms that
native geometry tolerance constrains surface relocation as well as final
acceptance. A loose internal acceptance value is therefore not merely a
reporting choice. No speculative sizing sweep or native patch was needed.

Jobs 116–117 perform the initial reproductions. Jobs 118–119 repeat each twice
with `--require-admitted`. That option saves the complete summary before
failing if any output is inadmissible. Job 120 independently replays all four
acceptance meshes, including hashes, counts, quality quantiles, volume,
distance, Jacobians, facets, periodic nodes/caps and absence of unused nodes.
Five focused recipe/Fourier controls pass.

Job 121 checks the rejection path with the known failing uniform32 lattice.
It exits 1 after saving its summary, as required. Its failure is an intentional
acceptance-command control, not an interrupted computation.

The run directories are `runs/macbook-lattice-tet4-coverage-v1` and
`runs/macbook-spinodoid-tet4-coverage-v1`. The audit is
[tet4-coverage-audit-20261024.json](../../handoff/macbook/tet4-coverage-audit-20261024.json).
Reproduce either check from the repository root with the recorded native
environment and the shared queue, using the manifest command. Select
`meshers-tet4-historical32` with `--strengthened --require-admitted`.

All four fixed specimens now have admitted linear routes, using the already
replayed Gyroid and Schwarz controls. This does not prove arbitrary-input
support or an across-the-board Microgen/MMG advantage. The recovered lattice
uses more tetrahedra than the admitted actual control and has lower median
quality. No additional speedup is claimed from these acceptance checks.
Actual Gyroid/spinodoid controls remain inadmissible in tested routes. Raw
lattice scalar periodicity and true crease handling remain input limitations.
The full four-case objective and lattice stiffness reference remain open.
