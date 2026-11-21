"""Step-difference transport sweep on a uniform Cartesian grid."""

import numpy as np


def sweep_cartesian(x_edges, y_edges, directions, sigma_t, angular_source,
                    boundary_inflow):
    """Return angular flux on uniformly spaced x and y grids.

    x_edges, y_edges: NumPy arrays of cell boundary positions.
    directions: direction cosines (mu, eta) for each angular direction.
    sigma_t: scalar total macroscopic cross section in inverse grid-length units.
    angular_source: scalar source per unit solid angle.
    boundary_inflow: dictionary of scalar incoming angular fluxes keyed by
        "left", "right", "bottom", or "top". Missing sides have zero inflow.

    Output: psi[j, i, m], where j is y, i is x, and m is direction.
    """
    nx = len(x_edges) - 1
    ny = len(y_edges) - 1
    hx = x_edges[1] - x_edges[0]
    hy = y_edges[1] - y_edges[0]

    psi = np.zeros((ny, nx, len(directions)))
    # Each direction fills its own slice psi[:, :, m].
    for m, (mu, eta) in enumerate(directions):
        # Start at the incoming side so the upstream neighbor is already solved.
        # For negative directions, the stop value -1 is excluded, so we end at 0.
        if mu >= 0:
            x_order = range(nx)
        else:
            x_order = range(nx - 1, -1, -1)
        if eta >= 0:
            y_order = range(ny)
        else:
            y_order = range(ny - 1, -1, -1)

        # Dividing the integrated balance by cell area hx*hy gives these
        # streaming coefficients. They are magnitudes; the signs set the order.
        ax = abs(mu) / hx
        ay = abs(eta) / hy
        denominator = sigma_t + ax + ay
        if denominator <= 0:
            raise ValueError("The cell needs a positive collision or outflow coefficient.")

        for j in y_order:
            for i in x_order:
                # At a boundary, get the prescribed inflow before attempting
                # to index a neighbor. Otherwise, Python's -1 index wraps around.
                if mu > 0:
                    x_in = boundary_inflow.get("left", 0.0) if i == 0 else psi[j, i - 1, m]
                elif mu < 0:
                    x_in = boundary_inflow.get("right", 0.0) if i == nx - 1 else psi[j, i + 1, m]
                else:
                    x_in = 0.0  # This direction carries no flux across x faces.

                if eta > 0:
                    y_in = boundary_inflow.get("bottom", 0.0) if j == 0 else psi[j - 1, i, m]
                elif eta < 0:
                    y_in = boundary_inflow.get("top", 0.0) if j == ny - 1 else psi[j + 1, i, m]
                else:
                    y_in = 0.0  # This direction carries no flux across y faces.

                # Source plus inflow balances collisions plus outflow. Step
                # differencing uses this same cell value on every outgoing face.
                psi[j, i, m] = (angular_source + ax * x_in + ay * y_in) / denominator

    # Every direction is complete before the function returns.
    return psi