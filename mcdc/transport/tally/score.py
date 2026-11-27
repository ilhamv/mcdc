from numba import njit

####

import mcdc.mcdc_get as mcdc_get
import mcdc.transport.mesh as mesh_module
import mcdc.transport.physics as physics
import mcdc.transport.util as util

from mcdc.constant import (
    AXIS_T,
    AXIS_X,
    AXIS_Y,
    AXIS_Z,
    COINCIDENCE_TOLERANCE,
    COINCIDENCE_TOLERANCE_TIME,
    INF,
    NEUTRON_REACTION_CAPTURE,
    NEUTRON_REACTION_FISSION,
    NEUTRON_REACTION_TOTAL,
    SCORE_FLUX,
    SCORE_DENSITY,
    SCORE_COLLISION,
    SCORE_CAPTURE,
    SCORE_FISSION,
    SCORE_CURRENT_NET,
    SCORE_CURRENT_IN,
    SCORE_CURRENT_OUT,
    SCORE_ENERGY_DEPOSITION,
)
from mcdc.transport.geometry.interface import check_cell
from mcdc.transport.geometry.surface import get_normal_component
from mcdc.transport.tally.filter import get_filter_indices

# ======================================================================================
# Surface-crossing
# ======================================================================================


@njit
def surface_crossing(
    particle_container,
    surface,
    tally,
    simulation,
    data,
):
    particle = particle_container[0]
    sub_ID = tally["sub_ID"]
    surface_crossing_tally = simulation["surface_crossing_tallies"][sub_ID]

    # Get filter indices
    i_mu, i_azi, i_energy, i_time = get_filter_indices(
        particle_container,
        tally,
        data,
    )

    # No score if outside non-changing phase-space bins
    if i_mu == -1 or i_azi == -1 or i_energy == -1 or i_time == -1:
        return

    # Tally index
    idx_base = (
        tally["bin_offset"]
        + i_mu * tally["stride_mu"]
        + i_azi * tally["stride_azi"]
        + i_energy * tally["stride_energy"]
        + i_time * tally["stride_time"]
    )

    # Flux
    speed = physics.particle_speed(particle_container, simulation, data)
    mu = get_normal_component(particle_container, speed, surface, data)
    flux = particle["w"] / abs(mu)

    # Non-cell-filtered score
    if not surface_crossing_tally["cell_filtered"]:
        for i_score in range(tally["scores_length"]):
            score_type = mcdc_get.tally.scores(i_score, tally, data)

            score = 0.0
            if score_type == SCORE_CURRENT_NET:
                score = flux * mu
            elif score_type == SCORE_CURRENT_IN:
                if mu < 0.0:
                    score = particle["w"]
            elif score_type == SCORE_CURRENT_OUT:
                if mu > 0.0:
                    score = particle["w"]

            util.atomic_add(data, idx_base + i_score, score)
        return

    # Cell-filtered score
    previous_cell_ID = particle["cell_ID"]
    filter_cell_ID = surface_crossing_tally["cell_filter_ID"]
    filter_cell = simulation["cells"][filter_cell_ID]
    was_in_filter_cell = previous_cell_ID == filter_cell_ID
    now_in_filter_cell = check_cell(particle_container, filter_cell, simulation, data)
    entered_filter_cell = not was_in_filter_cell and now_in_filter_cell
    exited_filter_cell = was_in_filter_cell and not now_in_filter_cell
    for i_score in range(tally["scores_length"]):
        score_type = mcdc_get.tally.scores(i_score, tally, data)

        score = 0.0
        if score_type == SCORE_CURRENT_NET:
            if exited_filter_cell:
                score = particle["w"]
            elif entered_filter_cell:
                score = -particle["w"]
        elif score_type == SCORE_CURRENT_IN:
            if entered_filter_cell:
                score = particle["w"]
        elif score_type == SCORE_CURRENT_OUT:
            if exited_filter_cell:
                score = particle["w"]

        util.atomic_add(data, idx_base + i_score, score)


# ======================================================================================
# Interaction
# ======================================================================================


@njit
def interaction(interaction_data_container, tally, simulation, data):
    """Score a discrete or condensed contribution using its saved incoming state."""
    interaction_data = interaction_data_container[0]
    sub_ID = tally["sub_ID"]
    interaction_tally = simulation["interaction_tallies"][sub_ID]

    # Get filter indices
    i_mu, i_azi, i_energy, i_time = get_filter_indices(
        interaction_data_container["incident_particle"], tally, data
    )

    # No score if outside non-changing phase-space bins
    if i_mu == -1 or i_azi == -1 or i_energy == -1 or i_time == -1:
        return

    # Mesh tally indices if needed
    i_x, i_y, i_z = 0, 0, 0
    if interaction_tally["mesh_filtered"]:
        mesh = simulation["meshes"][interaction_tally["mesh_filter_ID"]]
        i_x, i_y, i_z = mesh_module.get_indices(
            interaction_data_container["incident_particle"], mesh, simulation, data
        )

        # No score outside mesh bins
        if i_x == -1 or i_y == -1 or i_z == -1:
            return

    # Tally index
    idx_base = (
        tally["bin_offset"]
        + i_mu * tally["stride_mu"]
        + i_azi * tally["stride_azi"]
        + i_energy * tally["stride_energy"]
        + i_time * tally["stride_time"]
    )
    if interaction_tally["mesh_filtered"]:
        idx_base += (
            +i_x * interaction_tally["mesh_stride_x"]
            + i_y * interaction_tally["mesh_stride_y"]
            + i_z * interaction_tally["mesh_stride_z"]
        )

    # Score
    for i_score in range(tally["scores_length"]):
        score_type = mcdc_get.tally.scores(i_score, tally, data)
        score = 0.0
        if score_type == SCORE_ENERGY_DEPOSITION:
            score = interaction_data["energy_deposition"]
        util.atomic_add(data, idx_base + i_score, score)


# ======================================================================================
# Tracklength tally
# ======================================================================================


@njit
def tracklength(particle_container, distance, tally, simulation, data):
    particle = particle_container[0]
    sub_ID = tally["sub_ID"]
    tracklength_tally = simulation["tracklength_tallies"][sub_ID]

    # Get filter indices
    i_mu, i_azi, i_energy, i_time = get_filter_indices(
        particle_container,
        tally,
        data,
    )

    # No score if outside non-changing phase-space bins
    if i_mu == -1 or i_azi == -1 or i_energy == -1:
        return

    # Particle/track properties
    x = particle["x"]
    y = particle["y"]
    z = particle["z"]
    t = particle["t"]
    ux = particle["ux"]
    uy = particle["uy"]
    uz = particle["uz"]
    ut = 1.0 / physics.particle_speed(particle_container, simulation, data)
    x_final = x + ux * distance
    y_final = y + uy * distance
    z_final = z + uz * distance
    t_final = t + ut * distance

    # No score if particle does not cross the time bins. Tracks shorter than
    # the coincidence tolerance still score: electron lifetimes are below it.
    t_min = mcdc_get.tally.time(0, tally, data)
    t_max = mcdc_get.tally.time_last(tally, data)
    if t_final <= t_min or t > t_max - COINCIDENCE_TOLERANCE_TIME:
        return

    # Get the appropriate time index if the filter starts in the future
    if t < t_min + COINCIDENCE_TOLERANCE_TIME:
        i_time = 0

    # ==============================
    # Mesh tally preparation [START]
    #   - Get mesh bin indices
    #   - Return if it's outside mesh grid

    # Mesh axis indices
    i_x, i_y, i_z = 0, 0, 0
    if tracklength_tally["mesh_filtered"]:
        mesh = simulation["meshes"][tracklength_tally["mesh_filter_ID"]]

        # Mesh axis indices
        i_x, i_y, i_z = mesh_module.get_indices(
            particle_container, mesh, simulation, data
        )

        # No score if particle does not cross the mesh bins
        # Also get the appropriate index if needed
        x_min = mesh_module.get_x(0, mesh, simulation, data)
        x_max = mesh_module.get_x(mesh["Nx"], mesh, simulation, data)
        if ux > 0.0:
            if (
                x_final < x_min + COINCIDENCE_TOLERANCE
                or x > x_max - COINCIDENCE_TOLERANCE
            ):
                return
            if x < x_min + COINCIDENCE_TOLERANCE:
                i_x = 0
        else:
            if (
                x < x_min + COINCIDENCE_TOLERANCE
                or x_final > x_max - COINCIDENCE_TOLERANCE
            ):
                return
            if x > x_max - COINCIDENCE_TOLERANCE:
                i_x = mesh["Nx"]
        #
        y_min = mesh_module.get_y(0, mesh, simulation, data)
        y_max = mesh_module.get_y(mesh["Ny"], mesh, simulation, data)
        if uy > 0.0:
            if (
                y_final < y_min + COINCIDENCE_TOLERANCE
                or y > y_max - COINCIDENCE_TOLERANCE
            ):
                return
            if y < y_min + COINCIDENCE_TOLERANCE:
                i_y = 0
        else:
            if (
                y < y_min + COINCIDENCE_TOLERANCE
                or y_final > y_max - COINCIDENCE_TOLERANCE
            ):
                return
            if y > y_max - COINCIDENCE_TOLERANCE:
                i_y = mesh["Ny"]
        #
        z_min = mesh_module.get_z(0, mesh, simulation, data)
        z_max = mesh_module.get_z(mesh["Nz"], mesh, simulation, data)
        if uz > 0.0:
            if (
                z_final < z_min + COINCIDENCE_TOLERANCE
                or z > z_max - COINCIDENCE_TOLERANCE
            ):
                return
            if z < z_min + COINCIDENCE_TOLERANCE:
                i_z = 0
        else:
            if (
                z < z_min + COINCIDENCE_TOLERANCE
                or z_final > z_max - COINCIDENCE_TOLERANCE
            ):
                return
            if z > z_max - COINCIDENCE_TOLERANCE:
                i_z = mesh["Nz"]

    # Mesh tally preparation [END]
    # ============================

    # Tally base index
    idx_base = (
        tally["bin_offset"]
        + i_mu * tally["stride_mu"]
        + i_azi * tally["stride_azi"]
        + i_energy * tally["stride_energy"]
        + i_time * tally["stride_time"]
    )
    if tracklength_tally["mesh_filtered"]:
        idx_base += (
            i_x * tracklength_tally["mesh_stride_x"]
            + i_y * tracklength_tally["mesh_stride_y"]
            + i_z * tracklength_tally["mesh_stride_z"]
        )

    # Sweep through the distance
    distance_swept = 0.0
    while distance_swept < distance - COINCIDENCE_TOLERANCE:
        # The next time grid
        t_next = mcdc_get.tally.time(i_time + 1, tally, data)

        # Get the distance to score in this segment
        if t_final < t_next - COINCIDENCE_TOLERANCE_TIME:
            distance_scored = distance - distance_swept
        else:
            distance_scored = (t_next - t) / ut

        # ===========================================
        # Mesh tally grid crossing evaluation [START]
        #   - Determine smaller distance and which
        #     axis is crossed

        axis_crossed = AXIS_T
        if tracklength_tally["mesh_filtered"]:
            mesh = simulation["meshes"][tracklength_tally["mesh_filter_ID"]]

            # x-direction
            if ux == 0.0:
                dx = INF
            else:
                if ux > 0.0:
                    x_next = mesh_module.get_x(i_x + 1, mesh, simulation, data)
                    x_next = min(x_next, x_final)
                else:
                    x_next = mesh_module.get_x(i_x, mesh, simulation, data)
                    x_next = max(x_next, x_final)
                dx = (x_next - x) / ux
            if dx <= distance_scored:
                axis_crossed = AXIS_X
                distance_scored = dx

            # y-direction
            if uy == 0.0:
                dy = INF
            else:
                if uy > 0.0:
                    y_next = mesh_module.get_y(i_y + 1, mesh, simulation, data)
                    y_next = min(y_next, y_final)
                else:
                    y_next = mesh_module.get_y(i_y, mesh, simulation, data)
                    y_next = max(y_next, y_final)
                dy = (y_next - y) / uy
            if dy <= distance_scored:
                axis_crossed = AXIS_Y
                distance_scored = dy

            # z-direction
            if uz == 0.0:
                dz = INF
            else:
                if uz > 0.0:
                    z_next = mesh_module.get_z(i_z + 1, mesh, simulation, data)
                    z_next = min(z_next, z_final)
                else:
                    z_next = mesh_module.get_z(i_z, mesh, simulation, data)
                    z_next = max(z_next, z_final)
                dz = (z_next - z) / uz
            if dz <= distance_scored:
                axis_crossed = AXIS_Z
                distance_scored = dz

        # Mesh tally grid crossing evaluation [END]
        # =========================================

        # Score
        flux = distance_scored * particle["w"]
        for i_score in range(tally["scores_length"]):
            score_type = mcdc_get.tally.scores(i_score, tally, data)
            score = 0.0
            if score_type == SCORE_FLUX:
                score = flux
            elif score_type == SCORE_DENSITY:
                speed = physics.particle_speed(particle_container, simulation, data)
                score = flux / speed
            elif score_type == SCORE_COLLISION:
                score = flux * physics.macro_xs(
                    NEUTRON_REACTION_TOTAL, particle_container, simulation, data
                )
            elif score_type == SCORE_CAPTURE:
                score = flux * physics.macro_xs(
                    NEUTRON_REACTION_CAPTURE, particle_container, simulation, data
                )
            elif score_type == SCORE_FISSION:
                score = flux * physics.macro_xs(
                    NEUTRON_REACTION_FISSION, particle_container, simulation, data
                )
            util.atomic_add(data, idx_base + i_score, score)

        # Accumulate distance swept
        distance_swept += distance_scored

        # Move the 4D position
        if tracklength_tally["mesh_filtered"]:
            x += distance_scored * ux
            y += distance_scored * uy
            z += distance_scored * uz
        t += distance_scored * ut

        # Increment index and check if out of bounds
        if axis_crossed == AXIS_T:
            i_time += 1
            idx_base += tally["stride_time"]
            if i_time == tally["time_length"] - 1:
                return
        elif tracklength_tally["mesh_filtered"]:
            mesh = simulation["meshes"][tracklength_tally["mesh_filter_ID"]]
            if axis_crossed == AXIS_X:
                if ux > 0.0:
                    i_x += 1
                    if i_x == mesh["Nx"]:
                        return
                    idx_base += tracklength_tally["mesh_stride_x"]
                else:
                    i_x -= 1
                    if i_x == -1:
                        return
                    idx_base -= tracklength_tally["mesh_stride_x"]
            elif axis_crossed == AXIS_Y:
                if uy > 0.0:
                    i_y += 1
                    if i_y == mesh["Ny"]:
                        return
                    idx_base += tracklength_tally["mesh_stride_y"]
                else:
                    i_y -= 1
                    if i_y == -1:
                        return
                    idx_base -= tracklength_tally["mesh_stride_y"]
            elif axis_crossed == AXIS_Z:
                if uz > 0.0:
                    i_z += 1
                    if i_z == mesh["Nz"]:
                        return
                    idx_base += tracklength_tally["mesh_stride_z"]
                else:
                    i_z -= 1
                    if i_z == -1:
                        return
                    idx_base -= tracklength_tally["mesh_stride_z"]


# =============================================================================
# Eigenvalue tally
# =============================================================================


@njit
def eigenvalue_tally(particle_container, distance, simulation, data):
    # Weighted track length
    particle = particle_container[0]
    flux = distance * particle["w"]

    # Fission production for the cycle eigenvalue
    nuSigmaF = physics.neutron_production_xs(
        NEUTRON_REACTION_FISSION, particle_container, simulation, data
    )

    util.atomic_add(simulation["eigenvalue_tally_nuSigmaF"], 0, flux * nuSigmaF)

    # Score densities during active cycles
    if not simulation["cycle_active"]:
        return

    # ==================================================================================
    # Neutron density
    # ==================================================================================

    v = physics.particle_speed(particle_container, simulation, data)
    n_density = flux / v
    util.atomic_add(simulation["eigenvalue_tally_n"], 0, n_density)

    # Maximum neutron density
    util.atomic_max(simulation["eigenvalue_tally_n_max"], 0, n_density)

    # ==================================================================================
    # Delayed neutron production and precursor density
    # ==================================================================================

    if nuSigmaF == 0.0:
        return

    # Calculate the precursor density for this flight
    material = simulation["materials"][particle["material_ID"]]
    E = particle["E"]
    if physics.neutron.multigroup.applicable(particle_container, simulation, data):
        C_density = _calculate_precursor_density_multigroup(
            flux, E, material, simulation, data
        )
    else:
        C_density = _calculate_precursor_density(flux, E, material, simulation, data)

    # Total precursor density for this flight
    util.atomic_add(simulation["eigenvalue_tally_C"], 0, C_density)
    util.atomic_max(simulation["eigenvalue_tally_C_max"], 0, C_density)


@njit
def _calculate_precursor_density(flux, E, material, simulation, data):
    """Score native precursor groups and return their total density contribution."""
    C_density = 0.0
    for i in range(material["N_nuclide"]):
        nuclide_ID = mcdc_get.material.nuclide_IDs(i, material, data)
        nuclide = simulation["nuclides"][nuclide_ID]
        density = mcdc_get.material.nuclide_densities(i, material, data)

        # Select fission reactions that produce delayed neutrons
        for r in range(nuclide["N_neutron_fission_reaction"]):
            reaction_ID = mcdc_get.nuclide.neutron_fission_reaction_IDs(
                r, nuclide, data
            )
            reaction = simulation["neutron_reactions"][reaction_ID]
            fission = simulation["neutron_fission_reactions"][reaction["sub_ID"]]
            if fission["N_delayed_precursor"] == 0:
                continue

            # Delayed production from this reaction
            xs = physics.neutron.native.reaction_micro_xs(E, reaction, nuclide, data)
            nu_d = physics.neutron.native.neutron_fission_delayed_multiplicity(
                E, fission, simulation, data
            )
            production = flux * density * xs * nu_d

            # Partition production using the reaction’s group fractions
            for j in range(fission["N_delayed_precursor"]):
                fraction = mcdc_get.neutron_fission_reaction.delayed_fractions(
                    j, fission, data
                )
                decay = mcdc_get.neutron_fission_reaction.delayed_decay_rates(
                    j, fission, data
                )
                C_density += _score_precursor_group(
                    j, production * fraction, decay, simulation, data
                )
    return C_density


@njit
def _calculate_precursor_density_multigroup(flux, E, material, simulation, data):
    """Score multigroup precursor groups and return their total density contribution."""
    # Fission contribution from the incident energy group
    mgxs = simulation["neutron_multigroup_data"][material["neutron_multigroup_ID"]]
    group = physics.neutron.multigroup._get_energy_group(E, mgxs, simulation, data)
    fissions = flux * mcdc_get.neutron_multigroup_data.fission(group, mgxs, data)

    # Score production and density by precursor group
    C_density = 0.0
    for j in range(mgxs["J"]):
        production = fissions * mcdc_get.neutron_multigroup_data.nu_d(
            group, j, mgxs, data
        )
        decay = mcdc_get.neutron_multigroup_data.decay_rate(j, mgxs, data)
        C_density += _score_precursor_group(j, production, decay, simulation, data)
    return C_density


@njit
def _score_precursor_group(j, production, decay, simulation, data):
    C = production / decay

    # Score the group at its index in the global tallies
    util.atomic_add(
        data, simulation["eigenvalue_tally_nuSigmaF_delayed_offset"] + j, production
    )
    util.atomic_add(data, simulation["eigenvalue_tally_C_group_offset"] + j, C)
    return C
