# Deep in the Flux of Texas
### A 2D Neutron Transport Study Using Elevation-Based Materials and Population-Based Sources

**NUEN 430 — Fall 2026 Final Project**  
**Team:** Joseangel Lima-Tapia, Madalyn Oliver, and Luke Wood  
**Status:** In development — scope updated after instructor feedback (November 5, 2026)

## Overview

This project uses a map of Texas as the basis for a **two-dimensional, one-group neutron transport simulation**. Terrain elevation will be used to assign hypothetical material regions, and population density will be used to place weak, strong, or zero external neutron sources. The primary goal is to implement, verify, and study a working transport solver before adding advanced features.

The **core solver** will use the **discrete ordinates ($S_N$) method with diamond-difference spatial discretization on a Cartesian grid**. We will first solve steady-state, fixed-source transport with absorption and isotropic scattering. Once the solver is verified, we will apply it to the Texas map and explore how the material layout and source distribution affect the scalar neutron flux.

> **Modeling note:** Texas provides a recognizable geometry for a numerical experiment. The assigned materials and neutron sources are hypothetical, not a physical prediction of neutron behavior over real geographic distances. The computational domain/length scales and macroscopic cross sections will be made consistent or nondimensionalized.

## Revised project scope

Following instructor feedback, the project is divided into a required **core** and **optional extensions**.

| Component | Core project | Optional extension (only after verification) |
| --- | --- | --- |
| Spatial grid | Structured 2D Cartesian grid covering a map of Texas | Unstructured triangular Gmsh mesh |
| Spatial method | **Diamond differencing** | **DG-P1** with upwind numerical fluxes |
| Angular method | **$S_N$ discrete ordinates** | Higher angular orders / comparison with DG-P1 |
| Problem type | **Steady-state, one-group, fixed-source transport** | Time-dependent transport |
| Interactions | Absorption and isotropic scattering; **no fission multiplication** | Carefully controlled fission study only if subcriticality is established |
| Boundaries | **Vacuum incoming flux** on the outer Cartesian boundary | Prescribed incident beam; more advanced reflection models |
| Time integration | Not required for the steady-state core | **Crank–Nicolson** with time-step refinement |
| Outputs | Verified flux solutions, error/convergence plots, Texas scalar-flux map | Transient flux animation; DG-P1 comparisons |

The previously proposed **DG-P1 + $S_N$ + white-reflective boundaries + Crank–Nicolson** formulation is **not the baseline project**. We will consider DG-P1 **only if the diamond-difference solver is working and verified first**.

### Why the change?

The original formulation combined multiple challenging methods: discontinuous finite elements on an unstructured mesh, angular transport, reflective boundary coupling, and implicit transients. A Cartesian diamond-difference $S_N$ solver is more closely aligned with the course methods and provides a tractable path to the required verification and results.

We are also avoiding fission multiplication in the baseline problem. Fissile material regions combined with zero-leakage reflective boundaries could make the system critical or supercritical, in which case a finite steady-state fixed-source solution would not exist. **For the core project, $\nu\Sigma_f=0$ everywhere.** Labels such as LEU and plutonium are retained as hypothetical region categories; modeling their actual fission physics is out of scope for the baseline.

## Mathematical model

For each discrete angular direction $\boldsymbol\Omega_m=(\mu_m,\eta_m,\ldots)$, the steady-state, one-group transport equation is

$$
\mu_m\frac{\partial\psi_m}{\partial x}
+\eta_m\frac{\partial\psi_m}{\partial y}
+\Sigma_t(x,y)\psi_m
=\frac{\Sigma_s(x,y)}{4\pi}\phi(x,y)
+\frac{Q(x,y)}{4\pi},
$$

where the scalar flux is approximated by angular quadrature,

$$
\phi(x,y)\approx\sum_m w_m\psi_m(x,y),
\qquad \sum_m w_m=4\pi.
$$

Here $\psi_m$ is the angular flux, $\Sigma_t$ and $\Sigma_s$ are the total and scattering macroscopic cross sections, and $Q$ is the external source. The physical domain is two-dimensional in space; the angular quadrature retains direction dependence. There is **no fission source term** in the baseline solver.

**Diamond differencing** balances angular flux across each rectangular cell and approximates a cell-centered flux by averages of incoming and outgoing face fluxes. Directional transport sweeps will be coupled through **source iteration** for isotropic scattering. Spatial and angular refinement, together with source-iteration convergence criteria, will be documented.

### Boundaries and Texas map representation

The baseline implementation will use a **rectangular Cartesian computational grid** containing a rasterized Texas outline. Inside Texas, elevation and population data determine the materials and sources. Cells outside the Texas mask will be assigned a specified nonmultiplying background material and no external sources. The **outer rectangular boundary** will use vacuum conditions (zero incoming angular flux).

This is deliberately simpler than a reflective boundary that follows every irregular segment of the state outline. A staircase representation of the Texas border or white-reflective boundaries may be investigated only after the core solver is verified. A time-limited incident boundary beam is likewise an optional extension, not a prerequisite for the steady-state study.

## Elevation-based material regions

The following thresholds come from our original proposal and may be adjusted when we choose and process the elevation dataset.

| Material ID | Elevation (m) | Hypothetical material label |
| --- | --- | --- |
| 0 | $h < 100$ | Water |
| 1 | $100 \le h < 300$ | Graphite |
| 2 | $300 \le h < 700$ | Boron |
| 3 | $700 \le h < 1500$ | Low-enriched uranium (5%) |
| 4 | $h \ge 1500$ | Weapons-grade plutonium |

For the **nonmultiplying baseline**, these are material-region identifiers with assigned absorption/scattering coefficients; they are **not complete physical models of fissile nuclear fuel**. Cross-section datasets and unit conventions remain to be finalized.

## Population-based source regions

| Source ID | Population density (people/km²) | External source |
| --- | --- | --- |
| 0 | $\rho < 100$ | No source: $Q=0$ |
| 1 | $100 \le \rho < 1000$ | Weak: $Q=Q_{\mathrm{weak}}$ |
| 2 | $\rho \ge 1000$ | Strong: $Q=Q_{\mathrm{strong}}$ |

The numerical values of $Q_{\mathrm{weak}}$ and $Q_{\mathrm{strong}}$ are still to be selected. The population density is a **data-driven way to position hypothetical neutron sources**, not a claim that populated areas produce these neutron source strengths.

## Verification and numerical study

Verification comes **before** the full Texas application.

1. **Pure-absorber attenuation:** For a homogeneous non-scattering case, compare the computed directional flux with the known exponential attenuation solution.
2. **Manufactured solution:** Prescribe a smooth angular-flux field, derive consistent forcing and boundary data, and measure numerical error.
3. **Mesh refinement:** Plot error versus Cartesian mesh size and estimate the observed spatial convergence order. Examine angular-order dependence as appropriate.
4. **Source-iteration convergence:** Track residual/change per iteration and explain the stopping criterion for scattering problems.
5. **Texas study:** Map materials and population-based sources and interpret how heterogeneity and source contrast change the steady-state scalar flux.
6. **Only if time permits:** Introduce a Crank–Nicolson transient, verify temporal order by time-step refinement, and animate an incident beam. An external-code comparison or DG-P1 implementation is also a stretch goal.

Diamond differencing can produce nonphysical negative angular fluxes in difficult, optically thick cells; such behavior, if encountered, will be reported and addressed through mesh/parameter studies rather than hidden.

## Current progress

- [x] Defined the Texas geography, hypothetical materials, and population-based source categories.
- [x] Created an initial **triangular Gmsh Texas mesh** for exploratory work; it is not the baseline diamond-difference grid.
- [x] Built a separate **toy diffusion prototype** for visualizing localized sources; it is not a neutron transport solver.
- [x] Updated the numerical approach in response to instructor feedback.
- [ ] Create a Cartesian grid with a rasterized Texas/material/source map.
- [ ] Implement steady-state one-group $S_N$ diamond-difference sweeps.
- [ ] Add isotropic scattering and source iteration.
- [ ] Verify pure-absorber attenuation and a manufactured solution.
- [ ] Run spatial/angular refinement and source-iteration convergence tests.
- [ ] Generate steady-state Texas flux maps and parameter-study results.
- [ ] **Optional:** Implement Crank–Nicolson transient transport and an animation.
- [ ] **Optional:** Revisit DG-P1 on the triangular mesh after the diamond-difference solver is verified.

## Proposed repository organization

```text
.
├── README.md
├── src/                 # Cartesian S_N / diamond-difference solver (planned)
├── data/                # Texas raster/material/source maps or retrieval instructions
├── tests/               # Pure absorber and manufactured-solution tests (planned)
├── results/             # Flux maps, convergence plots, tables (planned)
├── legacy_mesh/         # Earlier triangular Gmsh mesh, if retained
└── prototypes/          # Earlier diffusion visualization, if retained
```

*This layout is a suggestion, not a claim that these folders currently exist. File-specific run instructions will be added as modules are implemented.*

## Course milestones

| Due date (2026) | Deliverable |
| --- | --- |
| October 30 | Project proposal (submitted) |
| **November 20** | **Verification checkpoint:** error-versus-mesh-size or time-step plot with observed order |
| **December 7, 11:59 p.m.** | Final report, slides, and Git repository link |
| December 10 | Presentation and live Q&A defense |

## Reproducibility and attribution

We will document the Python dependencies, input datasets, mesh sizes, run commands, verification tolerances, and the script that produces each figure as the code is developed. Data and software sources will be credited, and the final report will include the course-required disclosure of AI use.

**Course:** NUEN 430, Fall 2026.  
**Note:** This README describes the **revised project plan**, not a claim that the core transport solver has already been implemented or verified.
