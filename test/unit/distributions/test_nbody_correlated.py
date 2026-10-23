import math

import h5py
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

import mcdc.transport.distribution as dist
from mcdc.object_.distribution import DistributionNBody
from mcdc.object_.hdf5_reader import read_distribution


def test_nbody_sample_correlated(mock_rng_sequence, prepare_simulation):
    # MCNP Theory & User Manual §2.4.3.5.4.13 (Law 66: N-body Phase Space Distribution)
    distribution = DistributionNBody(
        values=[2.0, 4.0, 6.0],
        probabilities=[1.0, 1.0, 1.0],
    )
    simulation_container, data = prepare_simulation(objects=[distribution])
    simulation = simulation_container[0]
    distribution_base = simulation["distributions"][distribution.ID]

    # First value samples energy, second value samples isotropic cosine.
    xi1, xi2 = 0.2, 0.75
    mock_rng = mock_rng_sequence(xi1, xi2)

    sampled_E, sampled_mu = dist.sample_correlated_distribution(
        2.0,
        distribution_base,
        mock_rng,
        simulation,
        data,
    )

    # The current implementation samples energy from the tabulated distribution and
    # samples the cosine isotropically. This test is therefore checking the current
    # reduced implementation, not reconstructing the full Law 66 rejection sampler
    # from Eq. (2.103) through Eq. (2.106).
    # The constant PDF is normalized to 0.25 over [2, 6], so inverse-CDF sampling
    # in the first bin gives E_out = 2 + xi_1 / 0.25.
    # For the angular part, MCNP Eq. (2.107) gives mu = 2 * xi_10 - 1 for isotropic
    # center-of-mass sampling.
    expected_E = 2.0 + xi1 / 0.25
    expected_mu = 2.0 * xi2 - 1.0

    np.testing.assert_allclose(sampled_E, expected_E, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(sampled_mu, expected_mu, rtol=0.0, atol=1e-12)


@pytest.mark.parametrize(
    "law, name",
    [(1, "histogram"), (2, "linear"), (3, "semilog-x"), (4, "semilog-y"), (5, "log")],
)
@pytest.mark.parametrize("from_hdf5", [False, True])
def test_nbody_interpolation_quantiles(
    law, name, from_hdf5, tmp_path, mock_rng_sequence, prepare_simulation
):
    """Check every density law against independently integrated quantiles."""
    if from_hdf5:
        with h5py.File(tmp_path / "nbody.h5", "w") as group:
            group.attrs["type"] = "n-body"
            group["value"] = [1.0, 4.0]
            group["pdf"] = [0.2, 1.1]
            group.create_group("interpolation")
            group["interpolation/laws"] = np.asarray([name], dtype=h5py.string_dtype())
            group["interpolation/boundaries"] = [2]
            distribution = read_distribution(group)
    else:
        distribution = DistributionNBody([1.0, 4.0], [0.2, 1.1], interpolations=law)
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    base = simulation["distributions"][distribution.ID]

    def density(x):
        fraction = math.log(x) / math.log(4.0) if law in (3, 5) else (x - 1.0) / 3.0
        if law == 1:
            return 0.2
        if law in (4, 5):
            return 0.2 * (1.1 / 0.2) ** fraction
        return 0.2 + fraction * 0.9

    integral = quad(density, 1.0, 4.0)[0]
    for probability in (0.01, 0.25, 0.9):
        expected = brentq(
            lambda x: quad(density, 1.0, x)[0] / integral - probability, 1.0, 4.0
        )
        energy, cosine = dist.sample_correlated_distribution(
            10.0, base, mock_rng_sequence(probability, 0.75), simulation, data
        )
        assert energy == pytest.approx(expected, abs=1e-10)
        assert cosine == pytest.approx(0.5)


def test_nbody_multiple_interpolation_regions(mock_rng_sequence, prepare_simulation):
    """A histogram interval and a logarithmic region share one normalized CDF."""
    distribution = DistributionNBody(
        [1.0, 2.0, 2.0, 4.0],
        [1.0, 1.0, 2.0, 4.0],
        interpolations=[1, 5],
        interpolation_boundaries=[2, 4],
    )
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    base = simulation["distributions"][distribution.ID]
    # The first region has mass 1; integrating p(x) = x from 2 to 4 adds 6.
    for probability, expected in [(0.1, 1.7), (0.5, 3.0)]:
        energy, cosine = dist.sample_correlated_distribution(
            10.0, base, mock_rng_sequence(probability, 0.75), simulation, data
        )
        assert energy == pytest.approx(expected, abs=1e-10)
        assert cosine == pytest.approx(0.5)
