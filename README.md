# Deep in the Flux of Texas

### A 2D Neutron Transport Study Using Elevation-Based Materials and Population-Based Sources

**NUEN 430 — Fall 2026 Final Project**<br>
**Team:** Joseangel Lima-Tapia, Madalyn Oliver, and Luke Wood<br>
**Status:** Geographic preprocessing implemented; transport solver and verification pending<br>
**Last updated:** November 13, 2026

## Overview

This project uses a map of Texas as the basis for a two-dimensional, single-energy-group neutron transport study, with **all cross sections evaluated at 2 MeV**. Terrain elevation determines hypothetical material regions, and population density determines where external neutron sources will be placed.

The repository currently contains the **geographic preprocessing workflow**: download elevation tiles, simplify the Texas boundary, generate an unstructured quadrilateral Gmsh mesh, sample elevation and population density at cell centers, and export material/source IDs for use by future solver scripts. It also includes functions for visualizing the mesh and mapped data.

Following the latest discussion with our professor, the planned core solver uses the **discrete ordinates ($S_N$) method with step differencing (first-order upwind finite volume) on the existing unstructured quadrilateral mesh**. This replaces the earlier Cartesian diamond-difference plan. Transport sweeps, scattering iteration, cross sections, and verification cases have not yet been implemented in this repository.

> **Modeling note:** Texas provides a recognizable geometry for a numerical experiment. The material labels and neutron sources are hypothetical. Geographic distances and transport cross sections will need consistent units or nondimensionalization before solving the transport equation.

## Repository contents

| File | Purpose |
| --- | --- |
| [`Quads_Meshing.py`](Quads_Meshing.py) | Builds the Texas quadrilateral mesh, samples raster data, assigns material/source IDs, exports arrays, and defines plotting functions. |
| [`elevation_download.py`](elevation_download.py) | Downloads the elevation GeoTIFFs listed in `data.txt` into `Elevation/`. |
| [`data.txt`](data.txt) | List of 191 USGS elevation-tile download URLs. |
| [`cb_2023_us_state_500k.zip`](cb_2023_us_state_500k.zip) | 2023 U.S. Census cartographic state boundaries at 1:500,000 scale; read directly from the ZIP. |
| [`usa_pd_2020_1km.tif`](usa_pd_2020_1km.tif) | Population-density raster used by the preprocessing script. |
| [`tx_mesh.msh`](tx_mesh.msh) | Included Gmsh mesh; regenerated when `Quads_Meshing.py` runs. |
| [`.gitignore`](.gitignore) | Excludes downloaded `Elevation/` data and Python bytecode. |

The workflow additionally creates:

```text
Elevation/          # Downloaded elevation GeoTIFFs; ignored by Git
mesh_params.npz     # Generated mesh and material/source arrays
```

There are currently no solver, test, or results directories. `mesh_params.npz` is generated locally and is not included in the repository; it is also not currently excluded by `.gitignore`.

## Setup and run

Use Python 3 with NumPy, GeoPandas, Shapely, Gmsh, meshio, Rasterio, Matplotlib, and Requests. Dependency versions have not yet been pinned in a requirements file.

The following commands are for **PowerShell**, run from the repository root. They use a virtual environment without requiring activation:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install numpy geopandas shapely gmsh meshio rasterio matplotlib requests
```

On other platforms, use the equivalent virtual-environment Python executable (`.venv/bin/python`). All scripts resolve input and output paths relative to the working directory, so run them from the repository root.

### 1. Download elevation data

```powershell
.\.venv\Scripts\python.exe elevation_download.py
```

The downloader reads all nonempty lines in `data.txt`, creates `Elevation/`, and downloads the tiles sequentially. Existing destination files are skipped. Each download is streamed to a temporary `.tif.blah` file and renamed to `.tif` after the transfer completes. This step requires internet access and can take substantial time and disk space.

The boundary ZIP and population raster are already included in the repository. Keep their filenames unchanged unless you also edit the paths in `Quads_Meshing.py`.

### 2. Generate the mesh and region assignments

```powershell
.\.venv\Scripts\python.exe Quads_Meshing.py
```

This regenerates `tx_mesh.msh`, samples the elevation and population rasters, prints cell counts and material/source distributions, and writes `mesh_params.npz`. It also prints the number of cells missing population data. The script rejects meshes with triangular surface cells or no quadrilaterals.

The script executes its preprocessing at the top level, including when imported. Load `mesh_params.npz` from another script to reuse the exported arrays without repeating the meshing and raster sampling.

### 3. Display maps

The plotting functions are defined but are not called automatically. To run preprocessing once and then inspect maps interactively:

```powershell
.\.venv\Scripts\python.exe -i Quads_Meshing.py
```

At the Python prompt, call any of these functions:

```python
mesh_map()    # Quadrilateral mesh edges
elev_map()    # Sampled elevation
pop_map()     # Sampled population density
mat_map()     # Material IDs
source_map()  # Source IDs
```

The functions display Matplotlib windows using `plt.show()`; they do not automatically save figures.

## Geographic preprocessing

`Quads_Meshing.py` performs the following steps:

1. Read the state boundaries, select Texas, and project coordinates to **EPSG:5070** in meters.
2. Keep the largest mainland polygon and rebuild its exterior, removing interior holes and smaller detached pieces.
3. Apply positive/negative buffers to close bays and smooth boundary details, then simplify the outline.
4. Create a first-order Gmsh mesh using recombination and subdivision to produce quadrilateral surface cells. Physical groups are named `interior` and `outer_boundary`.
5. Compute each quadrilateral's area centroid and transform those coordinates into each raster's coordinate system.
6. Sample elevation from `Elevation/*.tif` and population density from `usa_pd_2020_1km.tif`, then assign material and source IDs.
7. Save the mesh geometry and region IDs in `mesh_params.npz`.

The current geometry settings are edited directly in `Quads_Meshing.py`:

| Parameter | Current value | Purpose |
| --- | --- | --- |
| `bay_radius` | 10,000 m | Buffer distance for closing bays. |
| `ridge_radius` | 7,000 m | Buffer distance for smoothing smaller boundary features. |
| `simplify_tolerance` | 10,000 m | Allowed geometric simplification tolerance. |
| `h` | 30,000 m | Requested Gmsh characteristic mesh size. |

The simplified boundary approximates the Texas mainland. The requested mesh size does not imply uniform square cells. Raster values are sampled at cell centroids, rather than averaged over each cell's area.

### Elevation-based material regions

| Material ID | Elevation $h$ (m) | Hypothetical material label |
| --- | --- | --- |
| 0 | $h < 100$ | Water |
| 1 | $100 \le h < 300$ | Graphite |
| 2 | $300 \le h < 700$ | Boron |
| 3 | $700 \le h < 1500$ | Low-enriched uranium (5%) |
| 4 | $h \ge 1500$ | Weapons-grade plutonium |

These IDs currently label regions only; no cross sections are assigned. The planned model includes absorption, isotropic scattering, and fission production, with all cross sections evaluated at 2 MeV. The data source and material densities/compositions remain undecided: either read macroscopic cross sections directly or calculate them from microscopic cross sections and specified number densities ($\Sigma_x=\sum_i N_i\sigma_{x,i}$).

**Current missing-data behavior:** Elevation starts as `NaN`. Any cell that remains unsampled falls through the material-assignment `else` branch to material ID 4. The script does not stop if `Elevation/` is absent or incomplete, so confirm elevation coverage before interpreting the material map. Missing elevation is not explicitly flagged in the exported IDs.

### Population-based source regions

| Source ID | Population density $\rho$ (people/km²) | Intended external source |
| --- | --- | --- |
| 0 | $\rho < 100$, or missing population data | No source: $Q=0$ |
| 1 | $100 \le \rho < 1000$ | Weak: $Q=Q_{\mathrm{weak}}$ |
| 2 | $\rho \ge 1000$ | Strong: $Q=Q_{\mathrm{strong}}$ |

Only the integer source IDs are currently stored. Numerical values for $Q_{\mathrm{weak}}$ and $Q_{\mathrm{strong}}$ remain to be selected. Population density is used to position hypothetical sources; it does not represent actual neutron production.

## Exported mesh parameters

`mesh_params.npz` contains five NumPy arrays, with one shared cell ordering:

| Key | Shape | Contents |
| --- | --- | --- |
| `nodes` | `(N_nodes, 2)` | Node coordinates in EPSG:5070 meters. |
| `quads` | `(N_cells, 4)` | Zero-based node indices for each quadrilateral. |
| `centers` | `(N_cells, 2)` | Quadrilateral area-centroid coordinates in EPSG:5070 meters. |
| `cell_mat` | `(N_cells,)` | Integer material IDs, 0–4. |
| `cell_source` | `(N_cells,)` | Integer source IDs, 0–2. |

Load them in another script with:

```python
import numpy as np

with np.load("mesh_params.npz", allow_pickle=False) as data:
    nodes = data["nodes"]
    quads = data["quads"]
    centers = data["centers"]
    cell_mat = data["cell_mat"]
    cell_source = data["cell_source"]
```

The export does not include the raw sampled elevation/population arrays, cross sections, source strengths, boundary physical groups, or transport connectivity. The Gmsh physical groups remain in `tx_mesh.msh`.

## Planned transport model

The core project first solves a **source-free, one-group fission eigenvalue problem** to determine $k_{\mathrm{eff}}$ for the elevation-based Texas material map. Establish $k_{\mathrm{eff}}<1$ with numerical convergence before running the **steady-state, source-driven** calculation, which includes absorption, isotropic scattering, population-based external sources, and fission-produced neutrons. Subcriticality is a required result to check, not an assumed property of the prescribed materials.

The model uses a **single energy group with all cross sections evaluated at 2 MeV**, rather than averaged over an energy spectrum. These fixed cross sections will be used throughout the transport calculation.

For discrete direction $\boldsymbol\Omega_m=(\mu_m,\eta_m,\ldots)$:

$$
\mu_m\frac{\partial\psi_m}{\partial x}
+\eta_m\frac{\partial\psi_m}{\partial y}
+\Sigma_t(x,y)\psi_m
=\frac{\Sigma_s(x,y)}{4\pi}\phi(x,y)
+\frac{\nu\Sigma_f(x,y)}{4\pi}\phi(x,y)
+\frac{Q(x,y)}{4\pi},
$$

with scalar flux approximated by

$$
\phi(x,y)\approx\sum_m w_m\psi_m(x,y),
\qquad \sum_m w_m=4\pi.
$$

Here $\psi_m$ is angular flux, $\Sigma_t$ and $\Sigma_s$ are macroscopic total and scattering cross sections, and $Q$ is the external source. Fission neutrons remain in the same single energy group. Step-difference directional solves and source iteration couple the angular fluxes through scattering and fission. For the eigenvalue problem, set $Q=0$ and divide the fission production term by $k_{\mathrm{eff}}$.

The planned transport domain is a **convex hull enclosing Texas**. Cells inside Texas retain the elevation-based material assignments; the region inside the hull but outside Texas is void, with zero cross sections and no volumetric neutron sources. Apply vacuum inflow on the outer hull for the eigenvalue and steady-state application calculations. Texas/void interfaces are internal transport interfaces: neutrons may leave Texas, cross the void, and reenter Texas. The existing preprocessing script currently meshes the simplified Texas outline only; hull construction, void meshing, region IDs, face connectivity, and hull boundary tagging remain to be implemented.

### Step differencing on quadrilaterals

For a fixed ordinate, a general quadrilateral can have one incoming edge and three outgoing edges, or the reverse. The usual Cartesian diamond-difference closure assumes paired faces and cannot simply be transferred to these cells. Generalized diamond schemes exist, but step differencing provides a simpler closure: every outgoing face uses the cell-average angular flux, while incoming faces use upwind neighbor values or prescribed boundary inflow. No face pairing is needed. This follows the face-based balance and upwind interpolation in [NIST's finite-volume documentation](https://pages.nist.gov/fipy/en/master/numerical/discret.html) and [numerical schemes documentation](https://pages.nist.gov/fipy/en/master/numerical/scheme.html).

For cell $K$, let $A_K$ be its area, $\ell_f$ an edge length, and $\mathbf n_{Kf}$ the outward unit normal. Define $a_{mf}=(\mu_m,\eta_m)\cdot\mathbf n_{Kf}$: negative is incoming, positive is outgoing, and zero is tangential with no streaming contribution. Integrating our transport equation over a cell and applying the step closure gives the planned update:

$$
\psi_{Km}=\frac{A_K q_{Km}+\sum_{f:a_{mf}<0}\ell_f|a_{mf}|\psi^{\mathrm{in}}_{fm}}
{A_K\Sigma_{t,K}+\sum_{f:a_{mf}>0}\ell_f a_{mf}},
\qquad q_{Km}=\frac{(\Sigma_{s,K}+\nu\Sigma_{f,K})\phi_K+Q_K}{4\pi}.
$$

Set $\psi^{\mathrm{out}}_{fm}=\psi_{Km}$ on every outgoing edge. In the motivating example, one known inflow and the cell source determine the cell flux and all three outflows. Use one consistent upwind flux on each shared edge so that leakage cancels between neighboring cells. Nonnegative sources and inflows yield a nonnegative local update when the denominator is positive.

The tradeoff is **first-order spatial accuracy and numerical diffusion**, which can smear localized sources and material interfaces. Measure this error through refinement; second-order spatial convergence is no longer the target. Ordinary step differencing is distinct from step-characteristics integration. An unstructured neutron-transport application of first-order step differencing is described in [Jareteg's thesis](https://www.skc.kth.se/polopoly_fs/1.1028292.1605188289!/KlasJareteg-PhD-Thesis_compressed.pdf).

Local closure does not guarantee a cycle-free global sweep. Build a dependency graph for each ordinate and topologically order it where possible. If cycles occur, use a converged treatment such as lagged cyclic face values with iteration or a coupled solve for the cyclic cells. [OpenSn's discretization documentation](https://open-sn.github.io/opensn/theory/discretization.html) discusses sweep graphs and unstructured-mesh cycles. Validate positive areas, outward normals, and cell geometry, including convexity for the intended simple sweep implementation.

DG-P1 on an unstructured mesh, white-reflective boundaries, and Crank–Nicolson transient transport remain optional extensions after the core solver is verified.

### Infinite-medium verification and development order

Start with the step-difference solver on a grid of squares, then repeat the infinite-medium test on the arbitrary-quadrilateral Texas mesh before applying the heterogeneous Texas material/source map. For this test, use spatially uniform material coefficients and a uniform isotropic volumetric source throughout the test domain; the geographic mesh supplies geometry only. Prescribe the exact infinite-medium **incoming angular flux** at every boundary face instead of the application's vacuum inflow. Do not retain the heterogeneous materials or void region in this homogeneous test.

With the README convention that $Q$ is the angle-integrated external source, our homogeneous balance gives

$$
\phi_\infty=\frac{Q}{\Sigma_t-\Sigma_s-\nu\Sigma_f},
\qquad \psi_{m,\infty}=\frac{\phi_\infty}{4\pi},
$$

provided the denominator is positive. The initial infinite-medium test explicitly sets $\Sigma_s=\Sigma_f=0$. Its analytical scalar flux is $\phi_\infty=Q/\Sigma_t$, and the prescribed incoming angular flux is $\psi_{m,\infty}=Q/(4\pi\Sigma_t)$ under the README source convention. Test scattering and fission coupling separately afterward. If $q$ denotes the **total angular source**, including scattering and fission, the equivalent angular-flux relation is $\psi_\infty=q/\Sigma_t$. Thus the meaning of the source and whether the quantity is scalar or angular must be explicit. [OpenSn's pure-absorber example](https://open-sn.github.io/opensn/tutorials/foundations/first_1d_fixed_source.html) uses the scalar-flux relation $Q/\Sigma_t$ in the absorption-only case.

This constant solution checks boundary normalization, geometric balance, and source coupling; it does not measure spatial convergence order. Retain attenuation, nonconstant manufactured solutions, spatial/angular refinement, and particle-balance checks. Use homogeneous test coefficients with a finite positive solution even if the eventual finite Texas configuration is subcritical by leakage.

## Progress and verification plan

- [x] Include the Texas boundary dataset and population-density raster.
- [x] Add a USGS elevation URL list and a downloader that skips existing files.
- [x] Implement mainland boundary simplification and quadrilateral meshing.
- [x] Implement centroid-based raster sampling and material/source classification.
- [x] Export mesh geometry and region IDs for use by other scripts.
- [x] Add mesh, elevation, population, material, and source plotting functions.
- [ ] Validate elevation coverage and explicitly handle missing elevation before assigning materials.
- [ ] Extend the mesh to a convex hull around Texas, with source-free void cells outside Texas and vacuum conditions on the hull.
- [ ] Validate quad geometry and build cell areas, edge lengths, outward normals, face neighbors, and direction-dependent sweep graphs.
- [ ] Detect sweep cycles and implement a converged treatment if cycles occur.
- [ ] Obtain all material cross sections at 2 MeV, select source strengths, and establish consistent transport length units.
- [ ] Build the square-grid step-difference solver and pass the infinite-medium test using exact incoming boundary fluxes.
- [ ] Extend to arbitrary quadrilaterals and pass the same homogeneous test on the Texas mesh.
- [ ] Implement scattering/fission iteration and the source-free eigenvalue solve; establish subcriticality before the source-driven Texas run.
- [ ] Verify pure-absorber attenuation against an exponential solution.
- [ ] Verify a manufactured solution with consistent forcing and boundary data.
- [ ] Test one-inflow/three-outflow and three-inflow/one-outflow cells, shared-face conservation, and positivity for nonnegative inputs.
- [ ] Run spatial/angular refinement studies; check for approximately first-order spatial convergence on smooth cases and quantify numerical smearing.
- [ ] Track source-iteration convergence and document its stopping criterion.
- [ ] Generate Texas scalar-flux maps and material/source parameter studies.
- [ ] If time permits, add a time-dependent incident-beam problem, select its time integration, and verify temporal accuracy. DG-P1 remains a possible comparison.

Verification will precede the full Texas transport study. Track particle balance, iteration residuals, and flux positivity for nonnegative physical inputs. Manufactured forcing may be signed, so positivity is not required for every manufactured-solution test.

## Course milestones

| Due date (2026) | Deliverable |
| --- | --- |
| October 30 | Project proposal (submitted) |
| **November 20** | **Verification checkpoint:** error-versus-mesh-size or time-step plot with observed order |
| **December 7, 11:59 p.m.** | Final report, slides, and Git repository link |
| December 10 | Presentation and live Q&A defense |

## Data attribution and reproducibility

- **State boundaries:** U.S. Census Bureau, *2023 Cartographic Boundary File (SHP), State and Equivalent Entities for United States, 1:500,000*, included as `cb_2023_us_state_500k.zip` with its metadata.
- **Elevation:** USGS elevation GeoTIFFs; the exact tile download URLs are recorded in `data.txt`. Downloaded tiles stay local in the Git-ignored `Elevation/` directory.
- **Population density:** `usa_pd_2020_1km.tif`, included in the repository. Its original download source and attribution still need to be recorded.

Dependency versions, transport parameters, verification tolerances, and figure-generation commands will be recorded as the solver is developed. The final report will include the course-required disclosure of AI use.
