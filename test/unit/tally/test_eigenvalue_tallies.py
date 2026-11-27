import numpy as np
import pytest

import mcdc
import mcdc.mcdc_get as mcdc_get
from mcdc.transport.tally import closeout, score


def test_eigenvalue_statistics_use_active_cycle_ratios(
    prepare_simulation, crossing_particle
):
    def material(nu_p, nu_d, decay):
        return mcdc.Material.multigroup(
            fission=[1.0], nu_p=[nu_p], nu_d=[[nu_d]], decay_rate=[decay], speed=[2.0]
        )

    materials = [material(3.0, 1.0, 2.0), material(1.0, 1.0, 4.0)]
    cells = [mcdc.Cell(fill=m) for m in materials]

    def configure(sim):
        sim.settings.N_particle = 1
        sim.settings.set_eigenmode(N_inactive=1, N_active=2)

    container, data = prepare_simulation(cells=cells, configure=configure)
    sim = container[0]
    particle = crossing_particle(-1, 0.0, 1.0, w=1.0)
    for cycle, (m, distance) in enumerate(
        [(materials[0], 100.0), (materials[0], 1.0), (materials[1], 4.0)]
    ):
        sim["idx_cycle"] = cycle
        sim["cycle_active"] = cycle > 0
        particle[0]["material_ID"] = m.ID
        score.eigenvalue_tally(particle, distance, sim, data)
        closeout.eigenvalue_cycle(sim, data)
    closeout.eigenvalue_simulation(sim, data)

    # Active-cycle generation times are (1/2)/4 and (4/2)/8 seconds.
    np.testing.assert_allclose(
        [sim["generation_time_avg"][0], sim["generation_time_sdv"][0]], [0.1875, 0.0625]
    )

    # Active cycles give beta = (1/4, 1/2), lambda = (2, 4), C = (1/2, 1).
    np.testing.assert_allclose(
        [
            mcdc_get.simulation.beta_avg(0, sim, data),
            mcdc_get.simulation.beta_sdv(0, sim, data),
        ],
        [0.375, 0.125],
    )
    np.testing.assert_allclose(
        [
            mcdc_get.simulation.lambda_avg(0, sim, data),
            mcdc_get.simulation.lambda_sdv(0, sim, data),
        ],
        [3.0, 1.0],
    )
    np.testing.assert_allclose(
        [sim["C_avg"], sim["C_sdv"], sim["C_max"]], [0.75, 0.25, 1.0]
    )
    assert sim["beta_avg_length"] == sim["lambda_avg_length"] == 1


@pytest.mark.parametrize("group_counts", [(0, 0), (0, 8), (8, 8), (6, 8)])
def test_precursor_group_counts(group_counts, capsys):
    cells = []
    for count in group_counts:
        yields = np.ones((count, 1)) if count else None
        material = mcdc.Material.multigroup(fission=[1.0], nu_p=[1.0], nu_d=yields)
        cells.append(mcdc.Cell(fill=material))
    sim = mcdc.Simulation()
    sim.set_model(cells)
    sim.settings.set_eigenmode(N_inactive=1, N_active=2)

    if group_counts == (6, 8):
        with pytest.raises(SystemExit):
            sim.compile()
        assert "Inconsistent delayed precursor group counts" in capsys.readouterr().out
    else:
        sim.compile()
        assert len(sim.beta_avg) == len(sim.lambda_avg) == max(group_counts)
