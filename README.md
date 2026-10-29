# Deep in the Flux of Texas
### Modeling Neutron Behavior Across Heterogeneous Media

**NUEN 430 — Fall 2026 Final Project**  
**Team:** Joseangel Lima-Tapia, Madalyn Oliver, and Luke Wood

## Overview

This project aims to develop a **two-dimensional, one-group, time-dependent neutron transport solver** on a triangular mesh shaped like Texas. Elevation data will determine the material assigned to each mesh element, while population-density data will determine the strength of an illustrative external neutron source. We will study the resulting steady-state angular/scalar flux and the response to a time-limited incident beam along part of the boundary.

Texas is used as a computational geometry and data-mapping exercise; the materials and sources are hypothetical, **not a physical model of neutron transport across the real state**.

## Planned numerical method

| Component | Planned approach |
| --- | --- |
| Geometry and mesh | 2D Texas outline; triangular mesh in Gmsh format |
| Spatial discretization | **P1 discontinuous Galerkin (DG-P1)** finite elements |
| Interface treatment | Upwind numerical fluxes |
| Angular discretization | Discrete ordinates (**S_N**) |
| Time integration | **Crank–Nicolson** |
| Baseline boundary | **White (diffuse) reflection**; zero net boundary current |
| Transient excitation | Prescribed incident beam on a selected boundary segment, then switched off |
| Implementation | Python, NumPy, SciPy; Gmsh for mesh work |

The model will include streaming, total interactions, isotropic scattering, and an optional one-group fission-production term. The transport unknown is the **angular flux**; scalar flux maps will be computed by angular quadrature. After spatial and angular discretization, the time-dependent system will have the form

$$
M\frac{d\boldsymbol{\Psi}}{dt}+A\boldsymbol{\Psi}=\mathbf{F}(t),
$$

where $\boldsymbol{\Psi}$ collects angular-flux degrees of freedom. A Crank–Nicolson step is

$$
\left(M+\frac{\Delta t}{2}A\right)\boldsymbol{\Psi}^{n+1}
=\left(M-\frac{\Delta t}{2}A\right)\boldsymbol{\Psi}^{n}
+\frac{\Delta t}{2}\left(\mathbf{F}^{n+1}+\mathbf{F}^{n}\right).
$$

## Elevation-based materials

These elevation thresholds are **proposed model inputs** and may be revised when data are processed.

| Material ID | Elevation (m) | Material |
| --- | --- | --- |
| 0 | $h < 100$ | Water |
| 1 | $100 \le h < 300$ | Graphite |
| 2 | $300 \le h < 700$ | Boron |
| 3 | $700 \le h < 1500$ | Low-enriched uranium (5%) |
| 4 | $h \ge 1500$ | Weapons-grade plutonium |

## Population-based external sources

Source strengths are hypothetical and will be assigned based on population density.

| Source ID | Population density (people/km²) | Source |
| --- | --- | --- |
| 0 | $\rho < 100$ | None ($Q=0$) |
| 1 | $100 \le \rho < 1000$ | Weak ($Q_{\mathrm{weak}}$) |
| 2 | $\rho \ge 1000$ | Strong ($Q_{\mathrm{strong}}$) |

Numerical values for $Q_{\mathrm{weak}}$ and $Q_{\mathrm{strong}}$, the material cross sections, and the data sources still need to be finalized.

## Current status

- [x] Defined project scope and selected DG-P1, $S_N$, and Crank–Nicolson methods.
- [x] Created an initial Texas-shaped **2D triangular mesh** (`texas_2d_triangles.msh`) and editable geometry (`texas_2d.geo`).
- [x] Created a separate **toy diffusion prototype** to test mesh loading, localized sources, and animation output. **This is not the transport solver.**
- [ ] Obtain and map elevation and population-density datasets.
- [ ] Assign material cross sections and source strengths by element.
- [ ] Implement and verify steady-state DG-P1 / $S_N$ transport.
- [ ] Implement white-reflective boundaries and the incident beam.
- [ ] Implement and verify Crank–Nicolson transient transport.
- [ ] Produce convergence studies, flux maps, and an animation.

## Repository layout (proposed)

```text
.
├── README.md
├── mesh/
│   ├── texas_2d.geo
│   └── texas_2d_triangles.msh
├── data/                 # Elevation/population data or download instructions
├── src/                  # Planned transport solver and mesh/data utilities
├── tests/                # Planned manufactured and benchmark problems
├── results/              # Generated figures, tables, and animations
└── prototype/            # Optional: earlier diffusion-only demonstration
```

*The folder structure above is a plan; not all folders or files exist yet.*

## Mesh preview and regeneration

Open `texas_2d_triangles.msh` directly in **Gmsh** to inspect the initial triangular mesh. If the editable geometry is in `mesh/`, regenerate a mesh with:

```bash
gmsh mesh/texas_2d.geo -2 -format msh2 -o mesh/texas_refined.msh
```

The supplied mesh uses a projected Texas coordinate system measured in **kilometers**, while nuclear macroscopic cross sections are typically specified in inverse centimeters. Coordinate and coefficient units must therefore be made consistent (or the numerical problem nondimensionalized). The existing outline is simplified for computational use.

## Verification plan

1. **Manufactured solution:** prescribe a known angular flux, derive the corresponding source/boundary data, and evaluate error versus mesh refinement.
2. **Pure absorber benchmark:** compare a beam in a homogeneous absorbing domain with analytical exponential attenuation.
3. **Time-step refinement:** check temporal convergence of Crank–Nicolson for a smooth transient test.
4. **Steady/transient consistency:** after the extra boundary beam is removed, check that a stable system relaxes toward the steady solution with the remaining sources.
5. **Optional external comparison:** compare a matched steady-state case against **OpenSn**.

## Planned deliverables

- Texas mesh with elevation-based material and population-based source maps.
- Steady-state scalar neutron-flux map.
- Animation of a transient boundary beam and its subsequent relaxation.
- Mesh-refinement error plot and accompanying error table.
- Time-step refinement study and steady-state comparison.
- Reproducible code, final report, and presentation.

## Reproducibility

Detailed setup commands and figure-generation instructions will be added as the **transport solver** is implemented. The existing diffusion prototype can be kept separately for illustration, but its output should not be presented as a verified transport result.

## Credits and references

Developed for **NUEN 430 (Fall 2026)**. Data sources, nuclear cross-section references, external software citations, and the course-required AI-use statement will be added as the project progresses.
