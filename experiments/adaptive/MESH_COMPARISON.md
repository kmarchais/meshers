# Mesh-only comparison protocol

Kevin clarified the primary objective on 24 October 2026: generate a mesh faster,
with better element quality and geometry fidelity, using fewer or equivalent
numbers of elements. A smaller mesh that also improves the other metrics is a
stronger result. FE solve time and stiffness-reference investigations do not
replace this comparison.

Keep all four fixed specimens: Gyroid sheet, Schwarz P solid, FFT16 spinodoid
with the same 52 Fourier modes, and the same sequential cubic smooth-union
lattice. The former source data and hashes remain frozen.

Primary comparisons use Tet4 on both sides. A separate Tet10 tier compares
native curved elements against explicitly elevated straight elements from the
actual Microgen/MMG output. Report both tetrahedra and nodes, and state the
geometry representation. Never silently combine element orders.

Use the actual Microgen VTK extraction and public two-pass MMG workflow as the
default baseline. Also report declared stronger routes. The Schwarz strong
route omits the second call and repairs interior vertices. The lattice strong
route includes an interior repair. Spinodoid denser sampling keeps the original
FFT16 realization; it does not construct a new random specimen. Every extra
generation, repair and elevation cost is charged.

Measure specimen/field preparation, extraction, meshing and output conversion
in fresh workers on this Mac. Exclude imports, archive writes, independent
validation, measured instrumentation overhead and FE. Record failed-attempt
time separately from successful mesh generation. Meshers, numerical libraries
and VTK use one thread. All heavy workers use the shared one-worker queue,
8 GiB limit, 300-second per-worker timeout and 250,000 node/tet caps. Read each
case pilot before selecting its successor. Use alternating order, one warmup
and five measured repetitions for the final admitted comparisons.

Report generation time, tetrahedra, nodes, minimum/1st/5th/median quality,
fraction below 0.1, sampled two-sided surface distance, volume error, positive
sampled Jacobians, coherent facets and periodic nodes/caps. Quality uses the
same normalized corner metric on both arms. Curved Jacobian positivity is an
additional gate, not a claim of globally certified validity. Geometry distance
is sampled, not a certified Hausdorff bound. Keep the shared box-constrained
observer and identical thresholds for both arms.

Admission remains distance 0.002, absolute relative volume 1%, quality 0.1,
positive sampled Jacobians, coherent facets and periodic topology. Label all
failed gates. Failed competitors do not establish an admitted workflow speedup.
Raw metrics and failed attempts remain visible. Distinguish clear domination
from trade-offs; do not force equal counts where one result uses fewer elements
and improves every relevant metric. Volume estimates have sampling uncertainty,
so small apparent differences alone do not establish superior volume fidelity.

The first pilots are diagnostics, not final timing statistics. They deliberately
retain incorrect-recipe attempts and adapter failures, with their original
source/configuration snapshots, rather than replacing them with later successes.
