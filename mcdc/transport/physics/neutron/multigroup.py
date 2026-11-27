import math

from numba import njit

####

import mcdc.mcdc_get as mcdc_get
import mcdc.numba_types as type_
import mcdc.transport.geometry.motion as motion
import mcdc.transport.particle as particle_module
import mcdc.transport.particle_bank as particle_bank_module
import mcdc.transport.rng as rng
import mcdc.transport.util as util

from mcdc.constant import (
    PI,
    NEUTRON_MULTIGROUP_ENERGY_MIDPOINT,
    NEUTRON_MULTIGROUP_ENERGY_MIDPOINT_LOG,
    NEUTRON_MULTIGROUP_ENERGY_UNIFORM,
    NEUTRON_MULTIGROUP_ENERGY_UNIFORM_LOG,
    NEUTRON_REACTION_TOTAL,
    NEUTRON_REACTION_CAPTURE,
    NEUTRON_REACTION_ELASTIC_SCATTERING,
    NEUTRON_REACTION_FISSION,
    NEUTRON_REACTION_FISSION_DELAYED,
    NEUTRON_REACTION_FISSION_PROMPT,
)
from mcdc.transport.distribution import sample_isotropic_direction
from mcdc.transport.physics.util import apply_delayed_decay_time, scatter_direction

# ======================================================================================
# Applicability
# ======================================================================================


@njit
def applicable(particle_container, simulation, data):
    particle = particle_container[0]
    material = simulation["materials"][particle["material_ID"]]

    if not material["has_neutron_multigroup"]:
        return False

    if simulation["settings"]["neutron_multigroup"]["hybrid"]:
        mgxs_ID = material["neutron_multigroup_ID"]
        mgxs = simulation["neutron_multigroup_data"][mgxs_ID]

        E = particle["E"]
        E_min = mcdc_get.neutron_multigroup_data.energy_grid(0, mgxs, data)
        E_max = mcdc_get.neutron_multigroup_data.energy_grid_last(mgxs, data)

        if E < E_min or E >= E_max:
            return False

    return True


# ======================================================================================
# Particle attributes
# ======================================================================================


@njit
def particle_speed(particle_container, simulation, data):
    particle = particle_container[0]
    material = simulation["materials"][particle["material_ID"]]

    mgxs_ID = material["neutron_multigroup_ID"]
    mgxs = simulation["neutron_multigroup_data"][mgxs_ID]

    group = _get_energy_group(particle["E"], mgxs, simulation, data)

    return mcdc_get.neutron_multigroup_data.speed(group, mgxs, data)


# ======================================================================================
# Material properties
# ======================================================================================


@njit
def macro_xs(reaction_type, particle_container, simulation, data):
    particle = particle_container[0]
    material = simulation["materials"][particle["material_ID"]]

    mgxs_ID = material["neutron_multigroup_ID"]
    mgxs = simulation["neutron_multigroup_data"][mgxs_ID]

    group = _get_energy_group(particle["E"], mgxs, simulation, data)

    if reaction_type == NEUTRON_REACTION_TOTAL:
        return mcdc_get.neutron_multigroup_data.total(group, mgxs, data)
    elif reaction_type == NEUTRON_REACTION_CAPTURE:
        return mcdc_get.neutron_multigroup_data.capture(group, mgxs, data)
    elif reaction_type == NEUTRON_REACTION_ELASTIC_SCATTERING:
        return mcdc_get.neutron_multigroup_data.scatter(group, mgxs, data)
    elif reaction_type == NEUTRON_REACTION_FISSION:
        return mcdc_get.neutron_multigroup_data.fission(group, mgxs, data)
    return 0.0


@njit
def neutron_production_xs(reaction_type, particle_container, simulation, data):
    particle = particle_container[0]
    material = simulation["materials"][particle["material_ID"]]

    mgxs_ID = material["neutron_multigroup_ID"]
    mgxs = simulation["neutron_multigroup_data"][mgxs_ID]

    group = _get_energy_group(particle["E"], mgxs, simulation, data)

    # Total production
    if reaction_type == NEUTRON_REACTION_TOTAL:
        total = 0.0

        # Scattering production
        nu = mcdc_get.neutron_multigroup_data.nu_s(group, mgxs, data)
        xs = mcdc_get.neutron_multigroup_data.scatter(group, mgxs, data)
        total += nu * xs

        # Fission production
        nu = mcdc_get.neutron_multigroup_data.nu_f(group, mgxs, data)
        xs = mcdc_get.neutron_multigroup_data.fission(group, mgxs, data)
        total += nu * xs
        return total

    # Capture production (none)
    elif reaction_type == NEUTRON_REACTION_CAPTURE:
        return 0.0

    # Scattering production
    elif reaction_type == NEUTRON_REACTION_ELASTIC_SCATTERING:
        nu = mcdc_get.neutron_multigroup_data.nu_s(group, mgxs, data)
        xs = mcdc_get.neutron_multigroup_data.scatter(group, mgxs, data)
        return nu * xs

    # Fission production
    elif reaction_type == NEUTRON_REACTION_FISSION:
        nu = mcdc_get.neutron_multigroup_data.nu_f(group, mgxs, data)
        xs = mcdc_get.neutron_multigroup_data.fission(group, mgxs, data)
        return nu * xs

    # Prompt fission production
    elif reaction_type == NEUTRON_REACTION_FISSION_PROMPT:
        nu = mcdc_get.neutron_multigroup_data.nu_p(group, mgxs, data)
        xs = mcdc_get.neutron_multigroup_data.fission(group, mgxs, data)
        return nu * xs

    # Delayed neutron production
    elif reaction_type == NEUTRON_REACTION_FISSION_DELAYED:
        nu = mcdc_get.neutron_multigroup_data.nu_d_total(group, mgxs, data)
        xs = mcdc_get.neutron_multigroup_data.fission(group, mgxs, data)
        return nu * xs

    # Unsupported default
    return 0.0


# ======================================================================================
# Collision
# ======================================================================================


@njit
def collision(
    particle_container,
    flight_data_container,
    interaction_data_container,
    program,
    data,
):
    simulation = util.access_simulation(program)
    particle = particle_container[0]

    # Get the reaction cross-sections
    SigmaT = macro_xs(NEUTRON_REACTION_TOTAL, particle_container, simulation, data)
    SigmaS = macro_xs(
        NEUTRON_REACTION_ELASTIC_SCATTERING, particle_container, simulation, data
    )
    SigmaF = macro_xs(NEUTRON_REACTION_FISSION, particle_container, simulation, data)

    # Implicit capture
    if simulation["technique"]["implicit_capture"]["active"]:
        SigmaC = macro_xs(
            NEUTRON_REACTION_CAPTURE, particle_container, simulation, data
        )
        particle["w"] *= (SigmaT - SigmaC) / SigmaT
        SigmaT -= SigmaC

    # Sample reaction type and perform the reaction
    xi = rng.lcg(particle_container) * SigmaT
    total = SigmaS
    if total > xi:
        scattering(particle_container, program, data)
    else:
        total += SigmaF
        if total > xi:
            fission(particle_container, program, data)
        else:
            particle["alive"] = False


# ======================================================================================
# Reactions
# ======================================================================================


@njit
def scattering(particle_container, program, data):
    simulation = util.access_simulation(program)
    particle = particle_container[0]
    material = simulation["materials"][particle["material_ID"]]

    # Material attributes
    mgxs_ID = material["neutron_multigroup_ID"]
    mgxs = simulation["neutron_multigroup_data"][mgxs_ID]

    # Kill the current particle
    particle["alive"] = False

    # Adjust production and product weights if weighted emission
    weight_production = 1.0
    weight_product = particle["w"]
    if simulation["technique"]["weighted_emission"]["active"]:
        weight_target = simulation["technique"]["weighted_emission"]["weight_target"]
        weight_production = particle["w"] / weight_target
        weight_product = weight_target

    # Get number of secondaries
    group = _get_energy_group(particle["E"], mgxs, simulation, data)
    nu_s = mcdc_get.neutron_multigroup_data.nu_s(group, mgxs, data)
    N = int(math.floor(weight_production * nu_s + rng.lcg(particle_container)))

    # Set up secondary particle container
    particle_container_new = util.local_array(1, type_.particle_data)
    particle_new = particle_container_new[0]

    # Create the secondaries
    for n in range(N):
        # Set default attributes
        particle_module.copy_as_child(particle_container_new, particle_container)

        # Set weight
        particle_new["w"] = weight_product

        # Set phase-space
        sample_scattering_neutron(particle_container_new, mgxs, simulation, data)

        # Bank or revive
        is_last = n == N - 1
        particle_bank_module.bank_or_revive_particle(
            particle_container, particle_container_new, is_last, program, data
        )


@njit
def fission(particle_container, program, data):
    simulation = util.access_simulation(program)
    particle = particle_container[0]
    material = simulation["materials"][particle["material_ID"]]
    settings = simulation["settings"]

    # Material attributes
    mgxs_ID = material["neutron_multigroup_ID"]
    mgxs = simulation["neutron_multigroup_data"][mgxs_ID]

    # Kill the current particle
    particle["alive"] = False

    # Adjust production and product weights if weighted emission
    weight_production = 1.0
    weight_product = particle["w"]
    if simulation["technique"]["weighted_emission"]["active"]:
        weight_target = simulation["technique"]["weighted_emission"]["weight_target"]
        weight_production = particle["w"] / weight_target
        weight_product = weight_target

    # Fission yields
    group = _get_energy_group(particle["E"], mgxs, simulation, data)
    nu = mcdc_get.neutron_multigroup_data.nu_f(group, mgxs, data)

    # Get number of secondaries
    N = int(
        math.floor(
            weight_production * nu / simulation["k_eff"] + rng.lcg(particle_container)
        )
    )

    # Set up secondary particle container
    particle_container_new = util.local_array(1, type_.particle_data)
    particle_new = particle_container_new[0]

    # Create the secondaries
    for n in range(N):
        # Set default attributes
        particle_module.copy_as_child(particle_container_new, particle_container)

        # Set weight
        particle_new["w"] = weight_product

        # Set phase-space
        sample_fission_neutron(
            particle_container_new, particle_container, mgxs, simulation, data
        )

        # Bank the particle
        if settings["neutron_eigenvalue_mode"]:
            particle_bank_module.bank_census_particle(particle_container_new, program)
        else:
            is_last = n == N - 1
            particle_bank_module.bank_or_revive_particle(
                particle_container, particle_container_new, is_last, program, data
            )


# ======================================================================================
# Reaction products
# ======================================================================================


@njit
def implicit_collision_xs(particle_container, simulation, data):
    """Return total, scattering production, and fission production."""
    SigmaT = macro_xs(NEUTRON_REACTION_TOTAL, particle_container, simulation, data)
    nuSigmaS = neutron_production_xs(
        NEUTRON_REACTION_ELASTIC_SCATTERING, particle_container, simulation, data
    )
    nuSigmaF = (
        neutron_production_xs(
            NEUTRON_REACTION_FISSION, particle_container, simulation, data
        )
        / simulation["k_eff"]
    )
    return SigmaT, nuSigmaS, nuSigmaF


@njit
def sample_implicit_neutron(
    neutron_container, incident_container, threshold, nuSigmaS, simulation, data
):
    """Sample one neutron by production and return whether it came from fission."""
    material = simulation["materials"][incident_container[0]["material_ID"]]
    mgxs = simulation["neutron_multigroup_data"][material["neutron_multigroup_ID"]]
    if threshold < nuSigmaS:
        sample_scattering_neutron(neutron_container, mgxs, simulation, data)
        return False
    sample_fission_neutron(
        neutron_container, incident_container, mgxs, simulation, data
    )
    return True


@njit
def sample_scattering_neutron(particle_container_new, mgxs, simulation, data):
    """Sample one scattered neutron from its incident direction and energy group."""
    particle = particle_container_new[0]
    group = _get_energy_group(particle["E"], mgxs, simulation, data)

    mu0 = 2.0 * rng.lcg(particle_container_new) - 1.0
    azi = 2.0 * PI * rng.lcg(particle_container_new)
    ux, uy, uz = scatter_direction(
        particle["ux"], particle["uy"], particle["uz"], mu0, azi
    )
    particle["ux"] = ux
    particle["uy"] = uy
    particle["uz"] = uz

    spectrum = mcdc_get.neutron_multigroup_data.chi_s_vector(group, mgxs, data)

    # Sample outgoing energy
    xi = rng.lcg(particle_container_new)
    total = 0.0
    for group_out in range(mgxs["G"]):
        total += spectrum[group_out]
        if total > xi:
            break
    particle["E"] = _get_group_energy(
        group_out, particle_container_new, mgxs, simulation, data
    )


@njit
def sample_fission_neutron(
    fission_neutron_container, incident_container, mgxs, simulation, data
):
    """Set the neutron's phase space at emission and return whether it is prompt."""
    fission_neutron = fission_neutron_container[0]
    incident = incident_container[0]

    group = _get_energy_group(incident["E"], mgxs, simulation, data)

    ux, uy, uz = sample_isotropic_direction(fission_neutron_container)
    fission_neutron["ux"] = ux
    fission_neutron["uy"] = uy
    fission_neutron["uz"] = uz

    nu = mcdc_get.neutron_multigroup_data.nu_f(group, mgxs, data)
    total = mcdc_get.neutron_multigroup_data.nu_p(group, mgxs, data)
    xi = rng.lcg(fission_neutron_container) * nu
    prompt = xi < total or mgxs["J"] == 0
    if prompt:
        spectrum = mcdc_get.neutron_multigroup_data.chi_p_vector(group, mgxs, data)
    else:
        nu_d = mcdc_get.neutron_multigroup_data.nu_d_vector(group, mgxs, data)
        for j in range(mgxs["J"]):
            total += nu_d[j]
            if xi < total:
                break
        spectrum = mcdc_get.neutron_multigroup_data.chi_d_vector(j, mgxs, data)
        decay_rate = mcdc_get.neutron_multigroup_data.decay_rate(j, mgxs, data)

    # Sample outgoing energy
    xi = rng.lcg(fission_neutron_container)
    total = 0.0
    for group_out in range(mgxs["G"]):
        total += spectrum[group_out]
        if total > xi:
            break
    fission_neutron["E"] = _get_group_energy(
        group_out, fission_neutron_container, mgxs, simulation, data
    )

    # Sample delayed emission
    if not prompt:
        apply_delayed_decay_time(fission_neutron_container, decay_rate, simulation)

        # Apply produc motion as needed
        motion.move_product(
            fission_neutron_container, incident_container, simulation, data
        )


# ======================================================================================
# Helpers
# ======================================================================================


@njit
def _get_energy_group(E, mgxs, simulation, data):
    if simulation["settings"]["neutron_multigroup"]["hybrid"]:
        offset = mgxs["energy_grid_offset"]
        length = mgxs["energy_grid_length"]
        E_grid = data[offset : offset + length]
        group = util.find_bin_with_rules(E, E_grid, 0.0, False)

    else:
        group = int(E)

    return group


@njit
def _get_group_energy(group, rng_state, mgxs, simulation, data):
    if simulation["settings"]["neutron_multigroup"]["hybrid"]:
        E_low = mcdc_get.neutron_multigroup_data.energy_grid(group, mgxs, data)
        E_high = mcdc_get.neutron_multigroup_data.energy_grid(group + 1, mgxs, data)
        representation = mgxs["energy_representation"]

        if representation == NEUTRON_MULTIGROUP_ENERGY_MIDPOINT:
            energy = 0.5 * (E_low + E_high)
        elif representation == NEUTRON_MULTIGROUP_ENERGY_MIDPOINT_LOG:
            energy = math.sqrt(E_low * E_high)
        elif representation == NEUTRON_MULTIGROUP_ENERGY_UNIFORM:
            xi = rng.lcg(rng_state)
            energy = E_low + xi * (E_high - E_low)
        elif representation == NEUTRON_MULTIGROUP_ENERGY_UNIFORM_LOG:
            xi = rng.lcg(rng_state)
            log_E_low = math.log(E_low)
            energy = math.exp(log_E_low + xi * (math.log(E_high) - log_E_low))
        else:
            # Unreachable
            energy = -1.0

    else:
        energy = float(group)

    return energy
