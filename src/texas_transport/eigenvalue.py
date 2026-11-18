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

def compute_fission_source(scalar_flux , nu_sigma_f , cell_areas):  # Calculate each cell's fission source and the area-weighted total.
    scalar_flux = np.asarray(scalar_flux , dtype = float)
    nu_sigma_f = np.asarray(nu_sigma_f , dtype = float)
    cell_areas = np.asarray(cell_areas , dtype = float)

    ### Checks ###
    if scalar_flux.ndim != 1 or nu_sigma_f.ndim != 1 or cell_areas.ndim != 1:
        raise ValueError('Inputs should be one dimensional.')

    if scalar_flux.size == 0 or nu_sigma_f.size == 0 or cell_areas.size == 0:
        raise ValueError('Inputs should not be empty.')

    if scalar_flux.shape != nu_sigma_f.shape or scalar_flux.shape != cell_areas.shape:
        raise ValueError('Inputs must have the same array shape.')

    if np.any(scalar_flux < 0):
        raise ValueError('Negative scalar flux is nonphysical.')

    if np.any(nu_sigma_f < 0):
        raise ValueError('Negative fission XS is nonphysical.')

    if np.any(cell_areas <= 0):
        raise ValueError('Negative/zero cell area is nonphysical.')

    # np.isfinite checks each entry; np.all is True only if every entry is finite.
    if not np.all(np.isfinite(scalar_flux)) or not np.all(np.isfinite(nu_sigma_f)) or not np.all(np.isfinite(cell_areas)):
        raise ValueError('Input contains a NaN or Infinite quantity.')

    ### Fission Source Calculation ###
    fission_source = scalar_flux * nu_sigma_f # Production rate per unit volume (this is fission source per cell), before k_eff or angular scaling.

    ### Integration Calculation ###
    fP = fission_source * cell_areas # Area-weighted fission production for each cell.
    total_fiss_source = np.sum(fP)

    return fission_source , total_fiss_source # Return the cell source array and its integrated total.

def update_keff(old_keff , old_total , new_total , target_total = 1.0):  # Return the updated k and the normalization factor.
    ### Checks ###
    for value in (old_keff , old_total , new_total , target_total):
        if np.ndim(value) != 0:
            raise ValueError('Inputs should be single values, not arrays.')

        if np.asarray(value).dtype.kind not in ('i' , 'u' , 'f'):
            raise ValueError('Inputs should be real numeric values.')

        if not np.isfinite(value):
            raise ValueError('Input contains a NaN or Infinite quantity.')

        if value <= 0:
            raise ValueError('Inputs should be greater than zero.')

    ### Calculations ###
    ratio = new_total / old_total
    new_keff = ratio * old_keff

    ### Normalization ###
    normalization_factor = target_total / new_total # Apply this factor to the new flux and fission source in the transport loop.

    return new_keff , normalization_factor

def check_eigenvalue_convergence(old_keff , new_keff , old_source , new_source , cell_areas , k_tolerance , source_tolerance): # Check that both k and the normalized fission-source shape have stopped changing.
    old_source = np.asarray(old_source , dtype = float)
    new_source = np.asarray(new_source , dtype = float)
    cell_areas = np.asarray(cell_areas , dtype = float)

    ### Checks ###
    for value in (old_keff , new_keff , k_tolerance , source_tolerance):
        if np.ndim(value) != 0:
            raise ValueError('k_eff values and tolerances should be simple numeric values.')

        if np.asarray(value).dtype.kind not in ('i' , 'u' , 'f'): # i = regular integer , u = unsigned integer (basically just the magnitude of the number; zero and positive whole numbers.); f = float
            raise ValueError('k_eff values and tolerances should be real numeric values.')

        if not np.isfinite(value):
            raise ValueError('k_eff values and tolerances should be finite.')

    if old_keff <= 0 or new_keff <= 0:
        raise ValueError('Old and new k_eff must be greater than zero.')

    if k_tolerance < 0 or source_tolerance < 0:
        raise ValueError('Tolerances must be zero or greater.')

    for value in (old_source , new_source , cell_areas):
        if value.ndim != 1:
            raise ValueError('Sources and cell areas should be one dimensional arrays.')

        if value.size == 0:
            raise ValueError('Sources and cell areas should not be empty.')

        if value.shape != old_source.shape:
            raise ValueError('Sources and cell areas must have the same array shape.')

        if not np.all(np.isfinite(value)):
            raise ValueError('Sources and cell areas must contain finite numbers.')

    if np.any(old_source < 0) or np.any(new_source < 0):
        raise ValueError('Negative fission source is nonphysical.')

    if np.any(cell_areas <= 0):
        raise ValueError('Negative/zero cell area is nonphysical.')

    # With positive areas, each source needs a positive entry to be normalized.
    if not np.any(old_source > 0) or not np.any(new_source > 0):
        raise ValueError('Old and new sources must each have positive fission production.')

    ### Calculations ###
    old_keff = float(old_keff)
    new_keff = float(new_keff)
    rel_k_change = abs(new_keff - old_keff) / new_keff

    old_source_contribution = old_source * cell_areas # Finds (old) source contribution in each individual cell
    old_source_contribution /= np.sum(old_source_contribution) # Normalizes to one

    new_source_contribution = new_source * cell_areas # Finds (new) source contribution in each individual cell
    new_source_contribution /= np.sum(new_source_contribution) # Normalizes to one

    source_change = np.sum(np.abs(new_source_contribution - old_source_contribution))

    if rel_k_change > k_tolerance or source_change > source_tolerance:
        return False
    else:
        return True
