# Reusable solver code

This is a package scaffold. No transport solver has been implemented yet.
Keep reusable numerical routines here, separate from scripts that read files or run cases.
Package installation/entry points should be added when the first solver modules are implemented.

Planned modules (not yet files):

| Module | Responsibility |
| --- | --- |
| `mesh.py` | Shared cell/face representation, square grids, mesh import, areas, normals, neighbors, region IDs. |
| `quadrature.py` | Discrete directions and weights with an explicit normalization. |
| `materials.py` | 2 MeV cross-section input, density/composition conversion, units, and consistency checks. |
| `sources.py` | Population IDs mapped to physical source strengths; scattering/fission source construction. |
| `boundaries.py` | Vacuum, prescribed analytical inflow, and eventually incident-beam boundary data. |
| `sweep.py` | Upwind step cell balance, ordinate-dependent dependencies, and cycle handling. |
| `fixed_source.py` | Scattering/fission iteration, convergence, and source-driven flux calculation. |
| `eigenvalue.py` | Source-free k-effective iteration and normalization. |

Start with a square-grid mesh and the absorption-only infinite-medium case.
Use the same cell/face interface for the arbitrary quadrilateral solver afterward.
Add a transient module only after the steady-state solver is verified.
Avoid file downloads, plotting windows, or full calculations as side effects of imports.
