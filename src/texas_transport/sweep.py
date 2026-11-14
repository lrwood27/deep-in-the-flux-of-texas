"""Starter outline for your Cartesian step-difference sweep."""

import numpy as np


def sweep_cartesian(x_edges, y_edges, directions, sigma_t, angular_source,
                    boundary_inflow):
    """Proposed inputs; feel free to change this interface as you develop it.

    x_edges, y_edges: NumPy arrays of cell boundary positions.
    directions: direction cosines (mu, eta) for each angular direction.
    sigma_t: total macroscopic cross section.
    angular_source: source per unit solid angle.
    boundary_inflow: prescribed incoming angular flux on each domain side.

    Suggested output: psi[j, i, m], where j is y, i is x, and m is direction.
    """
    # TODO 1: Find nx, ny, hx, and hy from the edge arrays.
    nx = len(x_edges) - 1
    ny = len(y_edges) - 1
    
    hx = x_edges[1] - x_edges[0]
    hy = y_edges[1] - y_edges[0]
    # TODO 2: Create an array to hold the angular flux.
    # TODO 3: Loop over the angular directions.
    # TODO 4: Choose the cell traversal order from the signs of mu and eta.
    # TODO 5: For each cell, get its incoming fluxes from neighbors or boundaries.
    # TODO 6: Apply the step-difference cell balance and store the cell flux.
    # TODO 7: Return the angular flux array.
    raise NotImplementedError("Write your transport sweep here.")
