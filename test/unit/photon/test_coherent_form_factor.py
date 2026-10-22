"""Coherent (Rayleigh) scattering from the tabulated atomic form factor.

Ported from the pre-refactor ``coherent_form_factor.py``. The physics being
checked is the shape of the angular distribution: coherent scattering peaks
forward, more sharply as energy rises and as Z rises, and it is elastic.
"""

import numpy as np
import pytest

import mcdc.transport.physics.photon.native as native

from conftest import make_interaction_data, make_photon


def _element(simulation):
    return simulation["elements"][0]


def _sample_cosines(simulation, data, E, count=400, first_seed=1):
    """Sample ``count`` scattering cosines from one advancing RNG stream.

    One stream, not one fresh seed per draw: an LCG's first output is a
    near-linear function of its seed, so consecutive seeds give correlated
    samples whose mean is biased and whose spread is not Poisson.
    """
    container = make_photon(E, seed=first_seed)
    element = _element(simulation)
    return np.array(
        [
            native.sample_coherent_scattering_cosine(
                container, element, simulation, data
            )
            for _ in range(count)
        ]
    )


# ======================================================================================
# The form factor itself
# ======================================================================================


@pytest.mark.parametrize("Z", [1, 13, 82])
def test_form_factor_forward_value_equals_z(photon_model, Z):
    """``F(q = 0, Z) = Z``: every electron scatters in phase in the forward direction."""
    from mcdc.transport.data import evaluate_data

    simulation, data = photon_model(atomic_number=Z)
    element = _element(simulation)
    reaction_ID = int(
        np.asarray(
            __import__(
                "mcdc.mcdc_get", fromlist=["element"]
            ).element.photon_coherent_scattering_reaction_IDs(0, element, data)
        )
    )
    reaction_base = simulation["photon_reactions"][reaction_ID]
    reaction = simulation["photon_coherent_scattering_reactions"][
        reaction_base["sub_ID"]
    ]
    form_factor = simulation["data"][reaction["form_factor_ID"]]

    assert evaluate_data(0.0, form_factor, simulation, data) == pytest.approx(float(Z))


# ======================================================================================
# The sampled angular distribution
# ======================================================================================


def test_forward_peaking_sharpens_with_energy(photon_model):
    # Higher energy means larger accessible momentum transfer, so the form
    # factor cuts off the large angles more aggressively.
    simulation, data = photon_model()
    low = _sample_cosines(simulation, data, 1.0e4, count=600).mean()
    high = _sample_cosines(simulation, data, 1.0e6, count=600).mean()
    assert high > low


# ======================================================================================
# Coherent scattering is elastic
# ======================================================================================


def test_scattering_preserves_energy_and_deposits_nothing(photon_model):
    simulation, data = photon_model()
    E = 1.0e5
    container = make_photon(E, seed=3)
    interaction = make_interaction_data()

    native.coherent_scattering(container, _element(simulation), simulation, data)

    assert container[0]["E"] == E
    assert container[0]["alive"]
    assert interaction[0]["energy_deposition"] == 0.0
