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
    cell_counts = []
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
        cell_counts.append(n * n)  # Total cells on an n-by-n grid.
        errors.append(error)
        print(f"N = {n:2d}, h = {hx:.4f}, outlet error = {error:.6e}")

    # If error is proportional to h**p, the change in log(error) divided
    # by the change in log(h) gives p. Use the two finest grids.
    order = np.log(errors[-2] / errors[-1]) / np.log(cell_widths[-2] / cell_widths[-1])

    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    ax.loglog(cell_widths, errors, "o-", label="Step differencing")
    ax.set_xlabel("Cell width h")
    ax.set_xticks(cell_widths, labels=[f"{width:g}" for width in cell_widths])
    ax.xaxis.set_ticks([], minor=True)
    ax.set_ylabel("Maximum absolute outlet angular-flux error")
    ax.set_title("Cartesian pure-absorber beam: mesh convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", title=f"Observed order p = {order:.3f}")

    # Place total-cell labels above the corresponding widths. The error
    # curve and the observed order still use cell width h.
    top = ax.secondary_xaxis("top")
    top.set_xticks(cell_widths, labels=[f"{count:,}" for count in cell_counts])
    top.minorticks_off()
    top.set_xlabel("Total cells")

    output = "results/cartesian_convergence.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Observed order p = {order:.3f}. Saved {output}")


if __name__ == "__main__":
    main()
