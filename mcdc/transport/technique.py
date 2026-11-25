import numpy as np
import math

from numba import njit

####

import mcdc.mcdc_get as mcdc_get
import mcdc.mcdc_get.weight_windows as ww_get
import mcdc.numba_types as type_
import mcdc.transport.particle as particle_module
import mcdc.transport.particle_bank as particle_bank_module
import mcdc.transport.physics as physics
import mcdc.transport.rng as rng
import mcdc.transport.util as util

from mcdc.constant import (
    COINCIDENCE_TOLERANCE_TIME,
    EVENT_TIME_CENSUS,
    NEUTRON_REACTION_TOTAL,
    NEUTRON_REACTION_ELASTIC_SCATTERING,
    NEUTRON_REACTION_FISSION,
    PARTICLE_NEUTRON,
    PARTICLE_ELECTRON,
    PARTICLE_PROTON,
    PARTICLE_PHOTON,
    PARTICLE_POSITRON,
)
from mcdc.transport.mesh import get_indices as get_mesh_indices

# ======================================================================================
# Implicit collision
# ======================================================================================


@njit
def implicit_collision(particle_container, flight_data_container, program, data):
    simulation = util.access_simulation(program)
    settings = simulation["settings"]
    particle = particle_container[0]
    particle["alive"] = False

    # Expected neutron production from each reaction
    nuSigmaS = physics.neutron.multigroup.neutron_production_xs(
        NEUTRON_REACTION_ELASTIC_SCATTERING, particle_container, simulation, data
    )
    nuSigmaF = (
        physics.neutron.multigroup.neutron_production_xs(
            NEUTRON_REACTION_FISSION, particle_container, simulation, data
        )
        / simulation["k_eff"]
    )
    production = nuSigmaS + nuSigmaF
    if production == 0.0:
        return

    # Sample the product count and distribute the expected outgoing weight.
    SigmaT = physics.neutron.multigroup.macro_xs(
        NEUTRON_REACTION_TOTAL, particle_container, simulation, data
    )
    config = simulation["technique"]["implicit_collision"]
    fixed_multiplicity = config["relative_density_length"] == 0

    if fixed_multiplicity:
        multiplicity = config["multiplicity"]
        weight_product = particle["w"] * production / (SigmaT * multiplicity)

    else:
        initial = flight_data_container[0]["initial_particle"]
        times = mcdc_get.implicit_collision.time_grid_all(config, data)

        density = mcdc_get.implicit_collision.relative_density_all(config, data)
        density_start = _get_density(initial["t"], times, density)
        density_end = _get_density(particle["t"], times, density)

        target = initial["w"] * (density_end / density_start)

        multiplicity = particle["w"] * (production / SigmaT) / target
        weight_product = target

    N = 1
    if multiplicity != 1.0:
        N = int(math.floor(multiplicity + rng.lcg(particle_container)))

    material = simulation["materials"][particle["material_ID"]]
    mgxs = simulation["neutron_multigroup_data"][material["neutron_multigroup_ID"]]
    particle_container_new = util.local_array(1, type_.particle_data)
    particle_new = particle_container_new[0]

    # Each product samples its own reaction.
    for n in range(N):
        scatter = nuSigmaF == 0.0 or rng.lcg(particle_container) * production < nuSigmaS
        particle_module.copy_as_child(particle_container_new, particle_container)
        particle_new["w"] = weight_product

        prompt = True
        if scatter:
            physics.neutron.multigroup.sample_scattering_neutron(
                particle_container_new, mgxs, simulation, data
            )
        else:
            physics.neutron.multigroup.sample_fission_neutron(
                particle_container_new, particle_container, mgxs, simulation, data
            )
            prompt = particle_new["t"] == particle["t"]

        if not scatter and settings["neutron_eigenvalue_mode"]:
            particle_bank_module.bank_census_particle(particle_container_new, program)
            continue

        # Split or roulette the delayed product to the target weight at emission.
        if not prompt and not fixed_multiplicity:
            t_emit = particle_new["t"]
            if t_emit > settings["time_boundary"] - COINCIDENCE_TOLERANCE_TIME:
                continue
            density_emit = _get_density(t_emit, times, density)
            target_emit = weight_product * (density_emit / density_end)
            survive = split_roulette(particle_container_new, target_emit, program, data)
            if not survive:
                continue

        is_last = n == N - 1
        particle_bank_module.bank_or_revive_particle(
            particle_container, particle_container_new, is_last, program, data
        )


@njit
def _get_density(t, times, density):
    if t <= times[0]:
        return density[0]
    if t >= times[-1]:
        return density[-1]
    i = util.find_bin(t, times)
    return util.linear_interpolation(
        t, times[i], times[i + 1], density[i], density[i + 1]
    )


# ======================================================================================
# Forced decay
# ======================================================================================


@njit
def force_decay(particle_container, decay_rate, time_boundary):
    """Force a finite-rate decay before the time boundary and correct its weight."""
    particle = particle_container[0]
    xi = rng.lcg(particle_container)

    # Uniform delayed emission probability on the remaining physical time domain.
    duration = time_boundary - particle["t"]
    delay = xi * duration
    particle["t"] += delay

    # Weight adjustment
    particle["w"] *= decay_rate * duration * math.exp(-decay_rate * delay)


# ======================================================================================
# Splitting and roulette
# ======================================================================================


@njit
def split_roulette(particle_container, weight_target, program, data):
    """Bank extra copies at the target weight and return whether the particle survives."""
    particle = particle_container[0]

    multiplicity = particle["w"] / weight_target
    particle["w"] = weight_target

    if multiplicity == 1.0:
        return True

    N = int(math.floor(multiplicity + rng.lcg(particle_container)))
    if N == 0:
        return False

    for _ in range(N - 1):
        copy_container = util.local_array(1, type_.particle_data)
        particle_module.copy_as_child(copy_container, particle_container)
        if not particle_bank_module.bank_if_not_active(copy_container, program, data):
            particle_bank_module.bank_active_particle(copy_container, program)
    return True


# ======================================================================================
# Weight Roulette
# ======================================================================================


@njit
def weight_roulette(particle_container, w_threshold, w_target):
    """
    Russian roulette particle if weight is below threshold.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    w_threshold : float
        Lower weight bound triggering roulette.
    w_target : float
        Target weight assigned upon survival.
    """
    particle = particle_container[0]
    if particle["w"] < w_threshold:
        survival_probability = particle["w"] / w_target
        # sample random number to determine survival
        if rng.lcg(particle_container) < survival_probability:
            particle["w"] = w_target
        else:
            particle["alive"] = False


# ======================================================================================
# Global weight Roulette
# ======================================================================================


@njit
def global_weight_roulette(particle_container, simulation):
    """
    Russian roulette particle with the global weight roulette parameters.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    simulation : object
        Simulation state containing global weight roulette parameters.
    """
    w_threshold = simulation["technique"]["global_weight_roulette"]["weight_threshold"]
    w_target = simulation["technique"]["global_weight_roulette"]["weight_target"]
    weight_roulette(particle_container, w_threshold, w_target)


# ======================================================================================
# Weight Windows
# ======================================================================================


@njit
def get_weight_window_object(particle_container, program):
    """
    Get the appropriate weight window object that applies to the provided
    particle.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    program : object
        Program object containing simulation state with weight window objects.

    Returns
    -------
    ww_obj : object
        The weight window object that corresponds to the provided particle.
    """
    simulation = util.access_simulation(program)
    technique = simulation["technique"]
    ptype = particle_container[0]["particle_type"]
    if ptype == PARTICLE_NEUTRON:
        ww_obj = technique["neutron_weight_windows"]
    elif ptype == PARTICLE_ELECTRON or ptype == PARTICLE_POSITRON:
        ww_obj = technique["electron_weight_windows"]
    elif ptype == PARTICLE_PROTON:
        ww_obj = technique["proton_weight_windows"]
    elif ptype == PARTICLE_PHOTON:
        ww_obj = technique["photon_weight_windows"]

    return ww_obj


@njit
def weight_windows(particle_container, program, data):
    """
    Apply weight window splitting and rouletting to a particle.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    program : object
        Program object containing simulation state with weight window and mesh data.
    data : object
        Simulation data for array access.
    """
    [lower, target, upper] = query_weight_window(particle_container, program, data)
    if upper == -1:
        return
    # split
    split_from_weight_window(particle_container, upper, target, lower, program)
    # roulette original particle
    weight_roulette(particle_container, lower, target)


@njit
def query_weight_window(particle_container, program, data):
    """
    Query weight window bounds for the particle.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    program : object
        Program object containing simulation state with weight window and mesh data.
    data : object
        Simulation data for array access.

    Returns
    -------
    lower : float
        Lower weight bound.
    target : float
        Target weight.
    upper : float
        Upper weight bound.
    """
    # grab objects
    ww_obj = get_weight_window_object(particle_container, program)
    indices = get_ww_indices(particle_container, ww_obj, program, data)
    for index in indices:
        if index == -1:
            return (-1.0, -1.0, -1.0)
    lower = ww_get.weights(*indices, 0, ww_obj, data)
    target = ww_get.weights(*indices, 1, ww_obj, data)
    upper = ww_get.weights(*indices, 2, ww_obj, data)
    return lower, target, upper


@njit
def get_ww_indices(particle_container, ww_obj, program, data):
    """
    Get the particle's bin index in each weight-window dimension.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    ww_obj : object
        The weight window object containing index information.
    program : object
        Program object containing simulation state with weight window and mesh data.
    data : object
        Simulation data for array access.

    Returns
    -------
    indices : tuple of int
        Seven bin indices in (time, energy, mu, azimuthal, x, y, z) order.
    """
    simulation = util.access_simulation(program)
    particle = particle_container[0]

    # get time index
    time_bounds = ww_get.time_bounds_all(ww_obj, data)
    it = util.find_bin_with_rules(
        particle["t"], time_bounds, COINCIDENCE_TOLERANCE_TIME, False
    )

    # get energy index
    energy_bounds = ww_get.energy_bounds_all(ww_obj, data)
    energy = particle["E"]
    ie = util.find_bin(energy, energy_bounds)

    # get angular indices
    mu, azimuthal = util.calculate_angles(particle_container, ww_obj["polar_reference"])

    # mu
    mu_bounds = ww_get.mu_bounds_all(ww_obj, data)
    imu = util.find_bin(mu, mu_bounds)

    # azimuthal
    azi_bounds = ww_get.azi_bounds_all(ww_obj, data)
    ia = util.find_bin(azimuthal, azi_bounds)

    # get spatial index
    mesh = simulation["meshes"][ww_obj["mesh_ID"]]
    idx, idy, idz = get_mesh_indices(particle_container, mesh, simulation, data)

    return (it, ie, imu, ia, idx, idy, idz)


@njit
def split_from_weight_window(particle_container, w_upper, w_target, w_lower, program):
    """
    Split a particle if its weight exceeds the threshold.

    Parameters
    ----------
    particle_container : ndarray
        Container holding the particle.
    w_upper : float
        Upper weight bound triggering splitting.
    w_target : float
        Target weight to assign to split particles.
    w_lower : float
        Lower weight bound triggering roulette on residual particle.
    program : object
        Program object containing simulation state with access to active bank.
    """
    particle = particle_container[0]
    weight = particle["w"]
    if weight > w_upper:
        # determine how many to split into
        num_split_to_target = math.floor(weight / w_target)

        # bank target particles
        particle["w"] = w_target
        for _ in range(num_split_to_target - 1):
            container_copy = util.local_array(1, type_.particle)
            particle_module.copy_as_child(container_copy, particle_container)
            if particle["event"] & EVENT_TIME_CENSUS:
                particle_bank_module.bank_census_particle(container_copy, program)
            else:
                particle_bank_module.bank_active_particle(container_copy, program)

        # bank residual particle
        residual_weight = weight - num_split_to_target * w_target
        if residual_weight > 0.0:
            residual_copy = util.local_array(1, type_.particle)
            particle_module.copy_as_child(residual_copy, particle_container)
            residual_copy[0]["w"] = residual_weight
            residual_copy[0]["alive"] = True
            weight_roulette(residual_copy, w_lower, w_target)
            if residual_copy[0]["alive"]:
                if particle["event"] & EVENT_TIME_CENSUS:
                    particle_bank_module.bank_census_particle(residual_copy, program)
                else:
                    particle_bank_module.bank_active_particle(residual_copy, program)


# ======================================================================================
# Population Control
# ======================================================================================


@njit
def population_control(simulation):
    """Uniform Splitting-Roulette technique"""

    bank_census = simulation["bank_census"]
    M = simulation["settings"]["N_particle"]
    bank_source = simulation["bank_source"]

    # Scan the bank
    idx_start, N_local, N = particle_bank_module.bank_scanning(bank_census, simulation)
    idx_end = idx_start + N_local

    # Abort if census bank is empty
    if N == 0:
        return

    # Weight scaling
    ws = float(N) / float(M)

    # Splitting Number
    sn = 1.0 / ws

    P_rec_arr = util.local_array(1, type_.particle_data)
    P_rec = P_rec_arr[0]

    # Perform split-roulette to all particles in local bank
    particle_bank_module.set_bank_size(bank_source, 0)
    for idx in range(N_local):
        # Weight of the surviving particles
        w = bank_census["particle_data"][idx]["w"]
        w_survive = w * ws

        # Determine number of guaranteed splits
        N_split = math.floor(sn)

        # Survive the russian roulette?
        xi = rng.lcg(bank_census["particle_data"][idx : idx + 1])
        if xi < sn - N_split:
            N_split += 1

        # Split the particle
        for i in range(N_split):
            particle_module.copy_as_child(
                P_rec_arr, bank_census["particle_data"][idx : idx + 1]
            )
            # Set weight
            P_rec["w"] = w_survive
            particle_bank_module.bank_source_particle(P_rec_arr, simulation)
