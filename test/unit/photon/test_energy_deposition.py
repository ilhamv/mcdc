"""Photon energy deposition, branch by branch.

Electrons are not transported, so every charged product deposits at the
collision site. The deposition each branch writes into
``interaction_data["energy_deposition"]`` is therefore fixed by what the branch
does *not* carry away:

===================  ==================================================
Coherent             0 -- elastic, nothing transferred
Incoherent           ``E - E_scattered``, the recoil electron
Photoelectric        ``E - E_fluorescence``
Pair production      ``E - 2 m_e c^2``, the pair's kinetic energy
===================  ==================================================

Ported from the pre-refactor ``test_energy_deposition.py``.
"""

import pytest

import mcdc.transport.physics.photon.native as native

from mcdc.constant import ELECTRON_MASS, PHOTON_CUTOFF_ENERGY

from conftest import make_interaction_data, make_photon

PAIR_THRESHOLD = 2.0 * ELECTRON_MASS


def _element(simulation):
    return simulation["elements"][0]


# ======================================================================================
# The cutoff
# ======================================================================================


def test_below_cutoff_the_photon_deposits_everything(photon_model):
    simulation, data = photon_model()
    E = PHOTON_CUTOFF_ENERGY * 0.5
    container = make_photon(E, seed=3)
    interaction = make_interaction_data()

    native.collision(container, interaction, simulation, data)

    assert interaction[0]["energy_deposition"] == pytest.approx(E)
    assert not container[0]["alive"]
    assert container[0]["E"] == 0.0


# ======================================================================================
# Conservation across a whole collision
# ======================================================================================


@pytest.mark.parametrize("E", [1.0e3, 5.0e5, 1.0e6, 5.0e6, 1.0e7])
def test_collision_never_deposits_more_than_it_was_given(photon_model, E):
    """Deposition is bounded by the incident energy, at every energy and seed.

    This is the invariant that catches a sign error or a unit error in any
    branch: whatever the sampled channel, the site cannot absorb more than
    arrived, and it cannot absorb a negative amount.
    """
    simulation, data = photon_model()
    for seed in range(1, 25):
        container = make_photon(E, seed=seed)
        interaction = make_interaction_data()

        native.collision(container, interaction, simulation, data)

        deposited = interaction[0]["energy_deposition"]
        assert deposited >= 0.0
        assert deposited <= E * (1.0 + 1e-12)
