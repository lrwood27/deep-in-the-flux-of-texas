"""Starter setup for your absorption-only infinite-medium check.

This is an unfinished template, not a working verification case.
"""

import numpy as np

# Grid: nx cells need nx+1 boundary positions.
nx = 10
ny = 10
hx = 0.1
hy = 0.1

x_edges = np.linspace(0.0, nx * hx, nx + 1)
y_edges = np.linspace(0.0, ny * hy, ny + 1)

# Synthetic homogeneous test data; scattering and fission are zero.
sigma_t = 1.0
Q = 1.0  # External source integrated over all directions.

# TODO 1: Choose directions and quadrature weights (weights should sum to 4*pi).
# TODO 2: Calculate the analytical scalar flux and incoming angular flux.
# TODO 3: Specify the angular source and boundary inflows for your sweep.
# TODO 4: Import and call your completed sweep function.
# TODO 5: Integrate angular flux to scalar flux using the quadrature weights.
# TODO 6: Compare your numerical flux with the analytical answer.
