"""Run from the repository root: python -m verification.attenuation_refinement."""

import matplotlib.pyplot as plt
import numpy as np

from src.texas_transport.sweep import sweep_cartesian


def main():
    # Send a beam in the +x direction across a unit square. With absorption
    # only, the exact flux is exp(-x) and the outlet at x = 1 has flux exp(-1).
    directions = [(1.0, 0.0)]
    sigma_t = 1.0
    angular_source = 0.0
    boundary_inflow = {"left": 1.0}
    exact_outlet = np.exp(-1.0)

    cell_widths = []
    errors = []
    for n in (10, 20, 40, 80):
        # n cells along each side need n + 1 edge positions.
        edges = np.linspace(0.0, 1.0, n + 1)
        hx = edges[1] - edges[0]  # The cells are squares, so hy = hx.
        psi = sweep_cartesian(
            edges, edges, directions, sigma_t, angular_source, boundary_inflow
        )

        # Take all rows, the last column, and the only direction. Step
        # differencing uses the last cell's flux at its outgoing face.
        outlet_flux = psi[:, -1, 0]
        error = np.max(np.abs(outlet_flux - exact_outlet))
        cell_widths.append(hx)
        errors.append(error)
        print(f"N = {n:2d}, h = {hx:.4f}, outlet error = {error:.6e}")

    # If error is proportional to h**p, the change in log(error) divided
    # by the change in log(h) gives p. Use the two finest grids.
    order = np.log(errors[-2] / errors[-1]) / np.log(cell_widths[-2] / cell_widths[-1])

    plt.figure(figsize=(7, 5))
    plt.loglog(cell_widths, errors, "o-", label="Step differencing")
    plt.xlabel("Cell width h")
    plt.ylabel("Maximum absolute outlet angular-flux error")
    plt.title("Cartesian pure-absorber beam: mesh convergence")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(loc="upper left", title=f"Observed order p = {order:.3f}")
    plt.tight_layout()

    output = "results/cartesian_convergence.png"
    plt.savefig(output, dpi=200)
    plt.close()
    print(f"Observed order p = {order:.3f}. Saved {output}")


if __name__ == "__main__":
    main()
