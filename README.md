# Deep in the Flux of Texas

### A 2D Neutron Transport Study Using Elevation-Based Materials and Population-Based Sources

**NUEN 430 — Fall 2026 Final Project**<br>
**Team:** Joseangel Lima-Tapia, Madalyn Oliver, and Luke Wood<br>
**Status:** Geographic preprocessing implemented; transport solver and verification pending<br>
**Last updated:** November 13, 2026

## Overview

This project uses a map of Texas as the basis for a two-dimensional, single-energy-group neutron transport study, with **all cross sections evaluated at 2 MeV**. Terrain elevation determines hypothetical material regions, and population density determines where external neutron sources will be placed.

The repository currently contains the **geographic preprocessing workflow**: download elevation tiles, simplify the Texas boundary, generate an unstructured quadrilateral Gmsh mesh, sample elevation and population density at cell centers, and export material/source IDs for use by future solver scripts. It also includes functions for visualizing the mesh and mapped data.

The planned core solver remains the **discrete ordinates ($S_N$) method with diamond-difference spatial discretization on a Cartesian grid**, following the November 5 instructor feedback. Transport sweeps, scattering iteration, cross sections, and verification cases have not yet been implemented in this repository. The current quadrilateral mesh is an unstructured geographic mesh; connecting the preprocessing work to the planned Cartesian solver remains an implementation task.

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

These IDs currently label regions only; no cross sections are assigned. For the planned nonmultiplying baseline, all regions will have absorption/scattering coefficients with $\nu\Sigma_f=0$.

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

The core project is a **steady-state, one-group, fixed-source** problem with absorption, isotropic scattering, and vacuum incoming angular flux on the outer Cartesian boundary. Fission multiplication is excluded from the baseline.

The model uses a **single energy group with all cross sections evaluated at 2 MeV**, rather than averaged over an energy spectrum. These fixed cross sections will be used throughout the transport calculation.

For discrete direction $\boldsymbol\Omega_m=(\mu_m,\eta_m,\ldots)$:

$$
\mu_m\frac{\partial\psi_m}{\partial x}
+\eta_m\frac{\partial\psi_m}{\partial y}
+\Sigma_t(x,y)\psi_m
=\frac{\Sigma_s(x,y)}{4\pi}\phi(x,y)
+\frac{Q(x,y)}{4\pi},
$$

with scalar flux approximated by

$$
\phi(x,y)\approx\sum_m w_m\psi_m(x,y),
\qquad \sum_m w_m=4\pi.
$$

Here $\psi_m$ is angular flux, $\Sigma_t$ and $\Sigma_s$ are macroscopic total and scattering cross sections, and $Q$ is the external source. Diamond-difference directional sweeps and source iteration will couple the angular fluxes through isotropic scattering.

The planned Cartesian domain will contain a rasterized Texas mask, with a specified nonmultiplying background and zero external source outside Texas. The current unstructured quad mesh does not provide the rectangular cells or sweep ordering assumed by this method; the Cartesian mapping and solver connectivity still need to be built.

DG-P1 on an unstructured mesh, white-reflective boundaries, and Crank–Nicolson transient transport remain optional extensions after the core solver is verified.

## Progress and verification plan

- [x] Include the Texas boundary dataset and population-density raster.
- [x] Add a USGS elevation URL list and a downloader that skips existing files.
- [x] Implement mainland boundary simplification and quadrilateral meshing.
- [x] Implement centroid-based raster sampling and material/source classification.
- [x] Export mesh geometry and region IDs for use by other scripts.
- [x] Add mesh, elevation, population, material, and source plotting functions.
- [ ] Validate elevation coverage and explicitly handle missing elevation before assigning materials.
- [ ] Build the Cartesian Texas/material/source grid required by the planned solver.
- [ ] Obtain all material cross sections at 2 MeV, select source strengths, and establish consistent transport length units.
- [ ] Implement steady-state $S_N$ diamond-difference sweeps and isotropic-scattering source iteration.
- [ ] Verify pure-absorber attenuation against an exponential solution.
- [ ] Verify a manufactured solution with consistent forcing and boundary data.
- [ ] Run spatial/angular refinement studies and report observed errors and convergence order.
- [ ] Track source-iteration convergence and document its stopping criterion.
- [ ] Generate Texas scalar-flux maps and material/source parameter studies.
- [ ] If time permits, add and verify transient transport or compare with DG-P1.

Verification will precede the full Texas transport study. Any negative angular fluxes encountered with diamond differencing will be examined through mesh and optical-thickness studies.

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
