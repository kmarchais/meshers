# Difficult periodic specimens on G3

## Selecting the next change

The saved finer sizing pilots put relocation at 0.69 to 1.23 seconds within
26.22 to 37.77 second total workflows. Eliminating that phase entirely would
save about 3%. The current Gyroid volume bias is +0.20% to +0.26%, and Schwarz
volume error is below 0.004%. All pass the fixed 1% volume gate. These costs and
errors remain optimization opportunities, but neither explains the missing
spinodoid and lattice coverage.

## Adaptive spinodoid

Job 26 tests the existing adaptive route on the fixed spectral spinodoid from
the real Microgen comparison. It starts from 16 cells, uses four adaptation
rounds, geometry tolerance 0.002, minimum quality 0.1, minimum length 0.01 and
maximum length 0.1. It fails after 38.09 seconds while unfolding the torus into
the box, reporting `a tetrahedron crosses a periodic plane`. No mesh is returned
or admitted. Records are in `runs/g3-spinodoid-adaptive16-v1`.

The admission runner keeps the 8 GiB worker ceiling, 600-second meshing timeout
and 250,000-tet output cap. It measures geometry and periodic matching before
admission. An admitted result still needs a justified stiffness reference and
the full six-load comparison.

Job 27 rebuilds with a detailed unfolding error and reproduces the failure.
The offending vertex has wall mask 1 and canonical y=0.4999788240285261. Its
tetrahedron uses the image y=-0.5000211759714739 beside three vertices just
inside the lower face. This is a finite crossing, not roundoff at the face.
The rebuild passes 20 core tests. Records are in
`runs/g3-spinodoid-adaptive16-v2`.

The projection helper currently wraps a projected root into the cell before
checking bounds. A periodic sine-wave regression tests whether this permits
a wall vertex to cross a constrained box plane. It also checks that projection
to a nearby root in the same cell remains available.

Job 28 confirms the regression fails under the original helper while 20 other
core tests pass. The fix checks projected coordinates before any periodic wrap,
and applies the same rule to interior relocation targets. Job 29 passes all 21
core tests and reruns the original spinodoid. Unfolding now completes. The next
failure is sampled geometry error 0.007612 above tolerance 0.002, after 44.87
seconds. Admission still fails. Records are in
`runs/g3-periodic-projection-red-v1` and `runs/g3-spinodoid-adaptive16-v3`.
