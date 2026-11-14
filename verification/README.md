# Numerical verification cases

Put analytical reference solutions, reproducible case drivers, and refinement studies here.
These are planned cases; none is implemented yet. Reuse solver routines from
`src/texas_transport/` rather than copying transport algorithms into case drivers.

Development order:

1. `infinite_medium_square`: uniform absorption-only material/source, with
   Sigma_s = Sigma_f = 0 and exact incoming angular flux Q / (4*pi*Sigma_t).
2. `infinite_medium_quads`: the same homogeneous test on the arbitrary-quadrilateral
   Texas geometry; replace geographic material/source assignments for this test.
3. `attenuation`: nonconstant pure-absorber solution to test directional transport.
4. `manufactured`: prescribed nonconstant angular flux, consistent forcing/inflow,
   and spatial/angular refinement. Expect first-order spatial accuracy for step differencing.
5. `balance_and_coupling`: face/global balance, scattering/fission coupling,
   inflow/outflow face patterns, and sweep-cycle convergence where applicable.

The constant infinite-medium solution tests consistency, not spatial convergence order.
Keep small case inputs/reference values in Git. Save generated fields, error tables,
and figures under `results/<case-name>/`.
Production Texas eigenvalue and population-source runs belong in workflow scripts,
not in the homogeneous infinite-medium verification case.
