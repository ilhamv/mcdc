import numpy as np
import pytest

import mcdc
import mcdc.numba_types as type_
from mcdc.constant import PARTICLE_NEUTRON
from mcdc.transport.technique import implicit_collision


def test_implicit_collision_configuration():
    technique = mcdc.Simulation().technique
    assert not technique.implicit_collision.active

    technique.implicit_collision(multiplicity=0.5)
    assert technique.implicit_collision.active
    assert technique.implicit_collision.multiplicity == 0.5

    technique.implicit_collision(relative_density=[3.0, 1.0], time=[0.5, 2.0])
    np.testing.assert_array_equal(technique.implicit_collision.time_grid, [0.5, 2.0])

    technique.branchless_collision()
    assert technique.implicit_collision.relative_density.size == 0
    assert technique.implicit_collision.active
    assert technique.implicit_collision.multiplicity == 1.0

    technique.branchless_collision(active=False)
    assert not technique.implicit_collision.active


def test_time_dependent_implicit_collision(prepare_simulation):
    material = mcdc.Material.multigroup(scatter=[[1.0]], capture=[1.0])

    def configure(simulation):
        simulation.technique.implicit_collision([3.0, 1.0], time=[0.5, 2.0])

    container, data = prepare_simulation(
        cells=[mcdc.Cell(fill=material)], configure=configure
    )
    simulation = container[0]
    particles = np.zeros(1, dtype=type_.particle)
    particle = particles[0]
    particle["material_ID"] = 0
    particle["rng_seed"] = 1
    particle["uz"] = 1.0
    particle["w"] = 8.0  # Mean neutron yield of 1/2 gives outgoing weight 4.
    particle["t"] = 1.25  # Interpolated density = 2.
    flight = np.zeros(1, dtype=type_.flight_data)
    initial = flight[0]["initial_particle"]
    initial["t"] = 0.25  # Clamped density = 3.
    initial["w"] = 3.0  # Target child weight = 3 * (2 / 3) = 2.
    implicit_collision(particles, flight, simulation, data)

    assert particle["alive"]
    assert simulation["bank_active"]["size"][0] == 1
    for child in (particle, simulation["bank_active"]["particle_data"][0]):
        assert child["w"] == 2.0


def test_relative_density_flight_starts_after_surface_crossing(prepare_simulation):
    from mcdc.transport.simulation import step_particle

    surface = mcdc.Surface.PlaneX(x=1.5)
    void = mcdc.Material.multigroup(capture=[0.0], speed=[1.0])
    scatterer = mcdc.Material.multigroup(scatter=[[1.0]], speed=[1.0])
    cells = [
        mcdc.Cell(region=-surface, fill=void),
        mcdc.Cell(region=+surface, fill=scatterer),
    ]
    time_grid = [0.0, 1.5, 2.5]
    density = [8.0, 2.0, 1.0]

    def configure(simulation):
        simulation.technique.implicit_collision(
            relative_density=density, time=time_grid
        )

    container, data = prepare_simulation(cells=cells, configure=configure)
    particles = np.zeros(1, dtype=type_.particle)
    particle = particles[0]
    particle["cell_ID"] = particle["material_ID"] = particle["surface_ID"] = -1
    particle["particle_type"] = PARTICLE_NEUTRON
    particle["rng_seed"] = 1
    particle["alive"] = True
    particle["ux"] = 1.0
    particle["w"] = 4.0

    step_particle(particles, container[0], data)
    assert particle["t"] == pytest.approx(1.5)
    step_particle(particles, container[0], data)
    assert particle["alive"]
    assert particle["t"] > 1.5
    assert particle["w"] == pytest.approx(
        4.0 * np.interp(particle["t"], time_grid, density) / 2.0
    )
