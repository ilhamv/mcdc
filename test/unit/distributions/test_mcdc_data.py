"""Interpolation and conditional-sampling checks for mcdc-data representations."""

import math

import h5py
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from mcdc.object_.hdf5_reader import read_distribution
from mcdc.object_.distribution import DistributionTabulated, DistributionMultiTable
from mcdc.transport import distribution as sampling


@pytest.mark.parametrize("law", [1, 2, 3, 4, 5])
def test_density_interpolation_quantiles(law, mock_rng_sequence, prepare_simulation):
    distribution = DistributionTabulated([1.0, 4.0], pdf=[0.2, 1.1], interpolations=law)
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    table = simulation["tabulated_distributions"][distribution.sub_ID]

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
        sampled = sampling.sample_tabulated(
            table, mock_rng_sequence(probability), simulation, data
        )
        assert sampled == pytest.approx(expected, abs=1.0e-10)


def test_density_regions_and_repeated_coordinates(
    mock_rng_sequence, prepare_simulation
):
    distribution = DistributionTabulated(
        [1.0, 2.0, 2.0, 4.0],
        pdf=[1.0, 1.0, 2.0, 4.0],
        interpolations=[1, 5],
        interpolation_boundaries=[2, 4],
    )
    # The duplicate value is a density jump: its interval carries zero mass.
    np.testing.assert_allclose(
        distribution.pdf.aux[0], [0.0, 1.0 / 7.0, 1.0 / 7.0, 1.0]
    )
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    table = simulation["tabulated_distributions"][distribution.sub_ID]
    # Integral from 2 to x is (x**2 - 4)/2; total density integral is 7.
    sampled = sampling.sample_tabulated(table, mock_rng_sequence(0.5), simulation, data)
    assert sampled == pytest.approx(3.0)


def test_histogram_incident_energy_and_endpoints(mock_rng_sequence, prepare_simulation):
    distribution = DistributionMultiTable(
        [1.0, 3.0],
        tables=[
            DistributionTabulated([10.0, 20.0], pdf=[1.0, 1.0]),
            DistributionTabulated([100.0, 200.0], pdf=[1.0, 1.0]),
        ],
        interpolations=1,
    )
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    table = simulation["multi_table_distributions"][distribution.sub_ID]
    for energy, expected in [
        (0.0, 12.5),
        (1.0, 12.5),
        (2.0, 12.5),
        (3.0, 125.0),
        (4.0, 125.0),
    ]:
        # Interior histogram selection consumes the selection variate first.
        values = (0.9, 0.25) if 1.0 < energy < 3.0 else (0.25,)
        result = sampling._sample_multi_table(
            energy, mock_rng_sequence(*values), table, simulation, data, True
        )
        assert result == pytest.approx(expected)


@pytest.mark.parametrize("kind", ["evaporation", "maxwellian"])
@pytest.mark.parametrize("temperature_kind", ["constant", "polynomial"])
def test_analytic_temperature_data(
    tmp_path, kind, temperature_kind, mock_rng_sequence, prepare_simulation
):
    with h5py.File(tmp_path / "spectrum.h5", "w") as group:
        group.attrs["type"] = kind
        group["restriction-energy"] = 0.0
        temperature = group.create_group("temperature")
        temperature.attrs["type"] = temperature_kind
        temperature["value" if temperature_kind == "constant" else "coefficient"] = (
            1.0 if temperature_kind == "constant" else [1.0]
        )
        distribution = read_distribution(group)
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    sample = sampling.sample_distribution(
        10.0,
        simulation["distributions"][distribution.ID],
        mock_rng_sequence(0.5, 0.5, 0.5),
        simulation,
        data,
    )
    assert 0.0 < sample < 10.0


@pytest.mark.parametrize("kind", ["kalbach-mann", "energy-angle-tabulated"])
def test_correlated_histogram_tables(
    tmp_path, kind, mock_rng_sequence, prepare_simulation
):
    """Preserve histogram densities, conditional angles, and parameter Data."""

    def rule(group, size):
        group["interpolation/laws"] = np.array(["histogram"], dtype=h5py.string_dtype())
        group["interpolation/boundaries"] = [size - 1]

    with h5py.File(tmp_path / "spectrum.h5", "w") as group:
        group.attrs["type"] = kind
        group["energy"] = [10.0, 20.0]
        for index, factor in enumerate([1.0, 2.0]):
            child = group.create_group(f"tables/{index}")
            outgoing = child.create_group("outgoing-energy")
            outgoing.attrs["type"] = "tabulated"
            outgoing["value"] = np.array([1.0, 2.0, 4.0]) * factor
            outgoing["pdf"] = np.array([0.2, 0.4, 0.0]) / factor
            rule(outgoing, 3)
            if kind == "kalbach-mann":
                child.create_group("precompound-factor").attrs["type"] = "constant"
                child["precompound-factor/value"] = 1.0
                child.create_group("angular-slope").attrs["type"] = "polynomial"
                child["angular-slope/coefficient"] = [0.0]
            else:
                for j in range(3):
                    angular = child.create_group(f"cosine-distributions/{j}")
                    angular.attrs["type"] = "tabulated"
                    angular["value"] = [-1.0, 0.0] if j < 2 else [0.0, 1.0]
                    angular["pdf"] = [1.0, 1.0]
                    rule(angular, 2)
        distribution = read_distribution(group)
    np.testing.assert_allclose(distribution.cdf, [0.0, 0.2, 1.0] * 2)
    np.testing.assert_allclose(distribution.pdf, [0.2, 0.4, 0.0, 0.1, 0.2, 0.0])
    container, data = prepare_simulation(objects=[distribution])
    simulation = container[0]
    base = simulation["distributions"][distribution.ID]
    for energy, expected in [(5.0, 3.5), (10.0, 3.5), (20.0, 7.0), (25.0, 7.0)]:
        E, mu = sampling.sample_correlated_distribution_with_scale(
            energy, base, mock_rng_sequence(0.1, 0.8, 0.75, 0.75), simulation, data
        )
        assert E == pytest.approx(expected)
        assert mu == pytest.approx(0.5 if kind == "kalbach-mann" else -0.25)
