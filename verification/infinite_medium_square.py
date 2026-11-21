import numpy as np

from src.texas_transport.sweep import sweep_cartesian

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

# Use one (mu, eta) pair in each quadrant. Each pair represents both signs
# of the out-of-plane direction cosine, so the four weights sum to 4*pi.
directions = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]]) / np.sqrt(3)
weights = np.full(4, np.pi)

# In a uniform infinite medium, the source balances absorption everywhere.
# The angular flux is isotropic, so divide scalar flux by 4*pi.
exact_phi = Q / sigma_t
exact_psi = exact_phi / (4 * np.pi)

# Q is integrated over all directions. The sweep needs the source per unit
# solid angle and the exact incoming angular flux on each side of the grid.
angular_source = Q / (4 * np.pi)
boundary_inflow = {
    "left": exact_psi,
    "right": exact_psi,
    "bottom": exact_psi,
    "top": exact_psi,
}

psi = sweep_cartesian(x_edges, y_edges, directions, sigma_t, angular_source, boundary_inflow)

# psi[j, i, m] holds one value per direction m. Summing over the last
# axis gives one scalar flux value phi[j, i] per cell.
phi = np.sum(psi * weights, axis=2)

# Check every cell and direction, allowing for floating-point roundoff.
np.testing.assert_allclose(psi, exact_psi, rtol=1e-12, atol=1e-12)
np.testing.assert_allclose(phi, exact_phi, rtol=1e-12, atol=1e-12)
print(f"PASS: maximum angular flux error = {np.max(np.abs(psi - exact_psi)):.3e}")
print(f"      maximum scalar flux error  = {np.max(np.abs(phi - exact_phi)):.3e}")