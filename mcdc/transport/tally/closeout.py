import math
import numpy as np

from numba import literal_unroll, njit, objmode
from mpi4py import MPI

####

import mcdc.config as config
import mcdc.mcdc_get as mcdc_get
import mcdc.mcdc_set as mcdc_set
import mcdc.output as output_module
import mcdc.transport.particle_bank as particle_bank_module

from mcdc.constant import (
    GYRATION_RADIUS_ALL,
    GYRATION_RADIUS_INFINITE_X,
    GYRATION_RADIUS_INFINITE_Y,
    GYRATION_RADIUS_INFINITE_Z,
    GYRATION_RADIUS_ONLY_X,
    GYRATION_RADIUS_ONLY_Y,
    GYRATION_RADIUS_ONLY_Z,
)

# ======================================================================================
# Tally score preparation
# ======================================================================================


@njit
def reduce(simulation, data):
    """Normalize and reduce all tally scores to the master rank."""
    for tally in simulation["tallies"]:
        _reduce(tally, simulation, data)


@njit
def _reduce(tally, simulation, data):
    """Normalize one tally's scores and sum them across MPI ranks."""
    N = tally["bin_length"]
    start = tally["bin_offset"]
    if N == 0:
        return
    chunks = min(N, simulation["settings"]["tally_reduction_chunks"])

    # Normalize
    N_particle = simulation["settings"]["N_particle"]
    for i in range(N):
        data[start + i] /= N_particle

    # MPI reduce
    master = simulation["mpi_master"]
    with objmode():
        chunk_size, remainder = divmod(N, chunks)
        gpu_backed = config.target == "gpu" and config.gpu_state_storage != "separate"
        if gpu_backed:
            # Reuse one chunk-sized host buffer for GPU-backed data.
            buff = np.zeros(chunk_size + (remainder > 0))

        offset = start
        for chunk in range(chunks):
            length = chunk_size + (chunk < remainder)
            end = offset + length
            scores = data[offset:end]
            if gpu_backed:
                received = buff[:length]
                MPI.COMM_WORLD.Reduce(scores, received, MPI.SUM, 0)
                if master:
                    data[offset:end] = received
            elif master:
                MPI.COMM_WORLD.Reduce(MPI.IN_PLACE, scores, MPI.SUM, 0)
            else:
                MPI.COMM_WORLD.Reduce(scores, None, MPI.SUM, 0)
            offset = end


@njit
def reset_scores(simulation, data):
    """Reset scores for all tallies."""
    for tally in simulation["tallies"]:
        _reset_scores(tally, data)


@njit
def _reset_scores(tally, data):
    """Reset scores of a tally."""
    start = tally["bin_offset"]
    for i in range(tally["bin_length"]):
        data[start + i] = 0.0


# ======================================================================================
# Tally statistics lifecycle
# ======================================================================================


@njit
def accumulate_statistics_and_reset_scores(simulation, data):
    """Accumulate sample statistics and reset scores on a rank owning moments."""
    # Batch/cycle samples: only the master updates moments after score reduction.
    # History samples: each rank updates its own moments after each source history.
    # Update sample counter
    simulation["N_tally_sample"] += 1
    N_sample = simulation["N_tally_sample"]

    for tally in simulation["tallies"]:
        # Update the running mean and sum of squared deviations.
        N_bin = tally["bin_length"]
        offset_bin = tally["bin_offset"]
        offset_mean = tally["bin_mean_offset"]
        offset_sum_squared_deviations = tally["bin_sum_squared_deviations_offset"]
        values = data[offset_bin : offset_bin + N_bin]
        mean = data[offset_mean : offset_mean + N_bin]
        sum_squared_deviations = data[
            offset_sum_squared_deviations : offset_sum_squared_deviations + N_bin
        ]
        _update_moments(values, mean, sum_squared_deviations, N_sample)

        # Reset score bins
        _reset_scores(tally, data)


@njit
def finalize(simulation, data):
    """Finalize all tally statistics and store their combined effective variance."""
    relative_variance = 0.0
    nonzero_bins = 0

    for tally in simulation["tallies"]:
        subtotal, count = _finalize(tally, simulation, data)
        relative_variance += subtotal
        nonzero_bins += count

    if nonzero_bins == 0:
        simulation["effective_variance"] = np.nan
    else:
        simulation["effective_variance"] = relative_variance / nonzero_bins


@njit
def _finalize(tally, simulation, data):
    """Finalize one tally's statistics and return its relative-variance sum and count."""
    N_bin = tally["bin_length"]

    # Reuse the buffers for the running mean and sum of squared deviations.
    offset_mean = tally["bin_mean_offset"]
    offset_sum_squared_deviations = tally["bin_sum_squared_deviations_offset"]
    mean = data[offset_mean : offset_mean + N_bin]
    sum_squared_deviations = data[
        offset_sum_squared_deviations : offset_sum_squared_deviations + N_bin
    ]

    N_sample = simulation["N_tally_sample"]

    # Only history-based statistics need rank-local moments merged at finalization.
    # Batch/cycle moments are already accumulated on the master rank.
    if simulation["history_based_statistics"]:
        N_sample = _reduce_moments(mean, sum_squared_deviations, N_sample, simulation)

    # All ranks must finish any reductions before workers can return.
    if not simulation["mpi_master"]:
        return 0.0, 0

    return _finalize_statistics(mean, sum_squared_deviations, N_sample)


@njit
def finalize_census(simulation, data):
    """Finalize all census tally statistics and store their combined effective variance."""
    relative_variance = 0.0
    nonzero_bins = 0

    for tally in simulation["tallies"]:
        subtotal, count = _finalize_census(tally, simulation, data)
        relative_variance += subtotal
        nonzero_bins += count

    if nonzero_bins == 0:
        simulation["effective_variance"] = np.nan
    else:
        simulation["effective_variance"] = relative_variance / nonzero_bins


@njit
def _finalize_census(tally, simulation, data):
    """Finalize one census tally's statistics and return its relative-variance sum and count."""
    # Census tallies use batches as independent samples.
    N_sample = simulation["settings"]["N_batch"]
    if not simulation["mpi_master"] or N_sample < 2:
        return 0.0, 0

    relative_variance = 0.0
    nonzero_bins = 0

    for score in range(tally["scores_length"]):
        for census in range(simulation["settings"]["N_census"] - 1):
            # The running mean and sum of squared deviations
            N_bin = tally["bin_length"] // tally["scores_length"]
            mean = np.zeros(N_bin)
            sum_squared_deviations = np.zeros(N_bin)
            for batch in range(N_sample):
                with objmode(values="float64[:]"):
                    values = output_module.read_census_score(
                        simulation, data, tally, score, batch, census
                    )
                _update_moments(values, mean, sum_squared_deviations, batch + 1)

            subtotal, count = _finalize_statistics(
                mean, sum_squared_deviations, N_sample
            )
            relative_variance += subtotal
            nonzero_bins += count

    return relative_variance, nonzero_bins


@njit
def reset_statistics(simulation, data):
    """Reset the sample count and accumulated statistics for all tallies."""
    simulation["N_tally_sample"] = 0
    for tally in simulation["tallies"]:
        _reset_statistics(tally, data)


@njit
def _reset_statistics(tally, data):
    """Reset one tally's mean and squared-deviation buffers."""
    N_bin = tally["bin_length"]
    offset_mean = tally["bin_mean_offset"]
    offset_sum_squared_deviations = tally["bin_sum_squared_deviations_offset"]

    for i in range(N_bin):
        data[offset_mean + i] = 0.0
        data[offset_sum_squared_deviations + i] = 0.0


# ======================================================================================
# Shared statistical calculations (history, batch, and cycle samples)
# ======================================================================================


@njit
def _update_moments(values, mean, sum_squared_deviations, N_sample):
    """Update the running mean and squared deviations using Welford's algorithm."""
    for i in range(len(values)):
        delta = values[i] - mean[i]
        mean[i] += delta / N_sample
        sum_squared_deviations[i] += delta * (values[i] - mean[i])


@njit
def _finalize_statistics(mean, sum_squared_deviations, N_sample):
    """Store standard errors and return the relative-variance sum and bin count."""
    relative_variance = 0.0
    nonzero_bins = 0

    for i in range(len(mean)):
        if N_sample > 1:
            # Convert squared deviations into standard error.
            radicand = sum_squared_deviations[i] / (N_sample - 1)
            radicand = radicand / N_sample
        else:
            radicand = 0.0

        # Clamp negative variance caused by round-off error.
        radicand = max(radicand, 0.0)
        sum_squared_deviations[i] = math.sqrt(radicand)

        # Accumulate squared relative errors for nonzero means.
        if mean[i] != 0.0:
            relative_error = sum_squared_deviations[i] / mean[i]
            relative_variance += relative_error * relative_error
            nonzero_bins += 1

    return relative_variance, nonzero_bins


# ======================================================================================
# History-based statistics only: parallel moment merging
# ======================================================================================
# Single-batch fixed-source CPU runs accumulate moments independently on each rank.
# These helpers merge those moments at finalization; batch/cycle statistics do not
# use this path because their scores are reduced before the master updates moments.
# Time-census and GPU fixed-source runs require batch-based statistics.


@njit
def _reduce_moments(mean, sum_squared_deviations, N_sample, simulation):
    """Merge rank-local history statistics to rank zero through a binary tree."""
    rank = simulation["mpi_rank"]
    size = simulation["mpi_size"]
    if size == 1:
        return N_sample

    count = np.zeros(1, dtype=np.int64)
    # Reuse the receive buffers at each tree level; never gather all rank data.
    N_receive = len(mean) if rank % 2 == 0 and rank + 1 < size else 0
    other_mean = np.empty(N_receive)
    other_sum_squared_deviations = np.empty(N_receive)
    stride = 1
    while stride < size:
        if rank % (2 * stride) == 0:
            source = rank + stride
            if source < size:
                # Isolate object-mode MPI calls from the tree's branches to avoid
                # Numba lowering errors.
                _receive_moments(
                    count, other_mean, other_sum_squared_deviations, source
                )
                N_sample = _merge_moments(
                    mean,
                    sum_squared_deviations,
                    N_sample,
                    other_mean,
                    other_sum_squared_deviations,
                    count[0],
                )
        else:
            destination = rank - stride
            count[0] = N_sample

            # Isolate object-mode MPI calls from the tree's branches to avoid
            # Numba lowering errors.
            _send_moments(count, mean, sum_squared_deviations, destination)
            break
        stride *= 2
    return N_sample


@njit
def _merge_moments(
    mean,
    sum_squared_deviations,
    N_sample,
    other_mean,
    other_sum_squared_deviations,
    N_other,
):
    """Merge two history-sample groups using their counts and centered moments."""
    if N_other == 0:
        return N_sample
    if N_sample == 0:
        mean[:] = other_mean
        sum_squared_deviations[:] = other_sum_squared_deviations
        return N_other

    N_total = N_sample + N_other
    fraction = N_other / N_total
    for i in range(len(mean)):
        delta = other_mean[i] - mean[i]
        mean[i] += delta * fraction
        sum_squared_deviations[i] += other_sum_squared_deviations[i] + delta * delta * (
            N_sample * fraction
        )
    return N_total


@njit
def _receive_moments(count, mean, sum_squared_deviations, source):
    """Receive a rank's history count, mean, and sum of squared deviations."""
    with objmode():
        MPI.COMM_WORLD.Recv(count, source=source, tag=0)
        MPI.COMM_WORLD.Recv(mean, source=source, tag=1)
        MPI.COMM_WORLD.Recv(sum_squared_deviations, source=source, tag=2)


@njit
def _send_moments(count, mean, sum_squared_deviations, destination):
    """Send a rank's history count, mean, and sum of squared deviations."""
    with objmode():
        MPI.COMM_WORLD.Send(count, dest=destination, tag=0)
        MPI.COMM_WORLD.Send(mean, dest=destination, tag=1)
        MPI.COMM_WORLD.Send(sum_squared_deviations, dest=destination, tag=2)


# ======================================================================================
# Eigenvalue statistics and diagnostics
# ======================================================================================


@njit
def eigenvalue_cycle(simulation, data):
    """Close out one eigenvalue cycle and update global statistics and diagnostics."""
    idx_cycle = simulation["idx_cycle"]
    N_particle = simulation["settings"]["N_particle"]

    # Local precursor group scores
    tally_delayed = mcdc_get.simulation.eigenvalue_tally_nuSigmaF_delayed_all(
        simulation, data
    )
    tally_C_group = mcdc_get.simulation.eigenvalue_tally_C_group_all(simulation, data)

    # MPI reduction buffers
    buff_nuSigmaF = np.zeros(1, np.float64)
    buff_n = np.zeros(1, np.float64)
    buff_C = np.zeros(1, np.float64)
    buff_delayed = np.zeros_like(tally_delayed)
    buff_C_group = np.zeros_like(tally_C_group)

    # Reduce total fission production every cycle
    with objmode():
        MPI.COMM_WORLD.Allreduce(
            np.array(simulation["eigenvalue_tally_nuSigmaF"]), buff_nuSigmaF, MPI.SUM
        )

        # Reduce density and precursor group tallies during active cycles
        if simulation["cycle_active"]:
            MPI.COMM_WORLD.Reduce(
                np.array(simulation["eigenvalue_tally_n"]), buff_n, MPI.SUM, 0
            )
            MPI.COMM_WORLD.Reduce(
                np.array(simulation["eigenvalue_tally_C"]), buff_C, MPI.SUM, 0
            )
            MPI.COMM_WORLD.Reduce(tally_delayed, buff_delayed, MPI.SUM, 0)
            MPI.COMM_WORLD.Reduce(tally_C_group, buff_C_group, MPI.SUM, 0)

    # Update k_eff on every rank for the next cycle
    simulation["k_eff"] = buff_nuSigmaF[0] / N_particle
    if simulation["mpi_master"]:
        mcdc_set.simulation.k_cycle(
            idx_cycle, simulation, data, value=simulation["k_eff"]
        )

    # =====================================================================
    # Active-cycle statistics
    # =====================================================================

    if simulation["cycle_active"] and simulation["mpi_master"]:
        # Normalize densities per source particle
        tally_n = buff_n[0] / N_particle
        tally_C = buff_C[0] / N_particle

        # Accumulate eigenvalue and density moments
        simulation["k_avg"] += simulation["k_eff"]
        simulation["k_sdv"] += simulation["k_eff"] * simulation["k_eff"]
        simulation["n_avg"] += tally_n
        simulation["n_sdv"] += tally_n * tally_n
        simulation["C_avg"] += tally_C
        simulation["C_sdv"] += tally_C * tally_C

        # Mean generation time from neutron population and fission production
        N = 1 + simulation["idx_cycle"] - simulation["settings"]["N_inactive"]
        production = buff_nuSigmaF[0]
        generation_time = buff_n[0] / production
        _update_ratio_moments(
            generation_time,
            N,
            0,
            simulation["generation_time_avg"],
            simulation["generation_time_sdv"],
        )

        # Sample each precursor group ratio once per active cycle
        beta_avg = mcdc_get.simulation.beta_avg_all(simulation, data)
        beta_sdv = mcdc_get.simulation.beta_sdv_all(simulation, data)
        lambda_avg = mcdc_get.simulation.lambda_avg_all(simulation, data)
        lambda_sdv = mcdc_get.simulation.lambda_sdv_all(simulation, data)
        for j in range(len(buff_delayed)):
            delayed = buff_delayed[j]
            beta = delayed / production
            decay = delayed / buff_C_group[j]
            _update_ratio_moments(beta, N, j, beta_avg, beta_sdv)
            _update_ratio_moments(decay, N, j, lambda_avg, lambda_sdv)

        # Running eigenvalue mean and standard error
        simulation["k_avg_running"] = simulation["k_avg"] / N
        if N == 1:
            simulation["k_sdv_running"] = 0.0
        else:
            simulation["k_sdv_running"] = math.sqrt(
                (simulation["k_sdv"] / N - simulation["k_avg_running"] ** 2) / (N - 1)
            )

    # Reset accumulators
    simulation["eigenvalue_tally_nuSigmaF"][0] = 0.0
    simulation["eigenvalue_tally_n"][0] = 0.0
    simulation["eigenvalue_tally_C"][0] = 0.0
    tally_C_group[:] = 0.0
    tally_delayed[:] = 0.0

    # =====================================================================
    # Gyration radius
    # =====================================================================

    if simulation["settings"]["use_gyration_radius"]:
        # Center of mass
        N_local = particle_bank_module.get_bank_size(simulation["bank_census"])
        total_local = np.zeros(4, np.float64)  # [x,y,z,W]
        total = np.zeros(4, np.float64)
        for i in range(N_local):
            P = simulation["bank_census"]["particle_data"][i]
            total_local[0] += P["x"] * P["w"]
            total_local[1] += P["y"] * P["w"]
            total_local[2] += P["z"] * P["w"]
            total_local[3] += P["w"]
        # MPI Allreduce
        with objmode():
            MPI.COMM_WORLD.Allreduce(total_local, total, MPI.SUM)
        # COM
        W = total[3]
        com_x = total[0] / W
        com_y = total[1] / W
        com_z = total[2] / W

        # Distance RMS
        rms_local = np.zeros(1, np.float64)
        rms = np.zeros(1, np.float64)
        gr_type = simulation["settings"]["gyration_radius_type"]
        if gr_type == GYRATION_RADIUS_ALL:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += (
                    (P["x"] - com_x) ** 2
                    + (P["y"] - com_y) ** 2
                    + (P["z"] - com_z) ** 2
                ) * P["w"]
        elif gr_type == GYRATION_RADIUS_INFINITE_X:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += ((P["y"] - com_y) ** 2 + (P["z"] - com_z) ** 2) * P["w"]
        elif gr_type == GYRATION_RADIUS_INFINITE_Y:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += ((P["x"] - com_x) ** 2 + (P["z"] - com_z) ** 2) * P["w"]
        elif gr_type == GYRATION_RADIUS_INFINITE_Z:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += ((P["x"] - com_x) ** 2 + (P["y"] - com_y) ** 2) * P["w"]
        elif gr_type == GYRATION_RADIUS_ONLY_X:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += ((P["x"] - com_x) ** 2) * P["w"]
        elif gr_type == GYRATION_RADIUS_ONLY_Y:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += ((P["y"] - com_y) ** 2) * P["w"]
        elif gr_type == GYRATION_RADIUS_ONLY_Z:
            for i in range(N_local):
                P = simulation["bank_census"]["particle_data"][i]
                rms_local[0] += ((P["z"] - com_z) ** 2) * P["w"]

        # Reduce squared distances to the rank reporting the gyration radius
        with objmode():
            MPI.COMM_WORLD.Reduce(rms_local, rms, MPI.SUM, 0)
        if simulation["mpi_master"]:
            radius = math.sqrt(rms[0] / W)
            mcdc_set.simulation.gyration_radius(
                idx_cycle, simulation, data, value=radius
            )


@njit
def eigenvalue_simulation(simulation, data):
    """Finalize global tally means and standard errors over active cycles."""
    # Reduce the largest single-flight density scores over all active cycles
    local_max = np.array(
        [
            simulation["eigenvalue_tally_n_max"][0],
            simulation["eigenvalue_tally_C_max"][0],
        ]
    )
    density_max = np.zeros(2, np.float64)
    with objmode():
        MPI.COMM_WORLD.Reduce(local_max, density_max, MPI.MAX, 0)
    if not simulation["mpi_master"]:
        return
    simulation["n_max"] = density_max[0]
    simulation["C_max"] = density_max[1]

    # Neutron and precursor density statistics
    N = simulation["settings"]["N_active"]
    simulation["n_avg"] /= N
    simulation["C_avg"] /= N
    if N > 1:
        simulation["n_sdv"] = math.sqrt(
            (simulation["n_sdv"] / N - simulation["n_avg"] ** 2) / (N - 1)
        )
        simulation["C_sdv"] = math.sqrt(
            max(simulation["C_sdv"] / N - simulation["C_avg"] ** 2, 0.0) / (N - 1)
        )
    else:
        simulation["n_sdv"] = 0.0
        simulation["C_sdv"] = 0.0

    # Standard errors of the active-cycle ratios
    if N > 1:
        simulation["generation_time_sdv"][0] = math.sqrt(
            simulation["generation_time_sdv"][0] / (N * (N - 1))
        )
        for sdv in (
            mcdc_get.simulation.beta_sdv_all(simulation, data),
            mcdc_get.simulation.lambda_sdv_all(simulation, data),
        ):
            sdv[:] = np.sqrt(sdv / (N * (N - 1)))


@njit
def _update_ratio_moments(value, N, j, mean, M2):
    """Update the running mean and sum of squared deviations for a cycle ratio."""
    delta = value - mean[j]
    mean[j] += delta / N
    M2[j] += delta * (value - mean[j])
