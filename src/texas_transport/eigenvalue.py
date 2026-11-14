# eigenvalue.py
# Outline of helpers the transport solver calls to find one global k_eff.
# k_eff is the effective neutron multiplication factor for the whole system.
# The transport solver handles the sweeps and iteration loop
# The eigenvalue calculation excludes external volume sources (and the possible incident beam).
#
# scalar_flux: average scalar neutron flux in each cell.
# nu_sigma_f: average neutrons per fission times the macroscopic fission cross section.
# fission_source: neutron production in each cell, before dividing by k_eff.
# cell_areas: area of each cell; all cell arrays must use the same cell order.
# Void cells remain in the arrays and have zero fission production.
# Old/new inputs refer to consecutive iterations; tolerances are allowed changes.

#################
#### Modules ####
#################

import numpy as np

#################
### Functions ###
#################

def build_cell_id_map(cell_ids):  # Match each Gmsh cell ID to its array position.
    cell_ids = np.asarray(cell_ids)

    ### Checks ###
    if cell_ids.ndim != 1:
        raise ValueError('Dimension of cell_ids array is not valid. Should be a 1D array.')
    
    if not (cell_ids.size >= 1):
        raise ValueError('Erm array should not be empty. Need at least one value.')

    if not np.issubdtype(cell_ids.dtype , np.integer):
        raise ValueError('Cell IDs should be integer values.')


def compute_fission_source(scalar_flux , nu_sigma_f):  # Multiply each cell's flux by its fission production cross section.
    pass


def integrate_fission_source(fission_source , cell_areas):  # Multiply each cell's source by its area, then add the results.
    pass


def update_keff(old_keff , old_total , new_total):  # Multiply the old k by the ratio of new to old total fission production.
    pass


def normalization_factor(fission_source , cell_areas , target_total = 1.0):  # Find one scale factor that makes total production equal target_total.
    pass


def check_eigenvalue_convergence(old_keff , new_keff , old_source , new_source , cell_areas , k_tolerance , source_tolerance): # Check that both k and the normalized fission-source shape have stopped changing.
    pass