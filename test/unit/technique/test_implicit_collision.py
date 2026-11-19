import mcdc


def test_implicit_collision_configuration():
    technique = mcdc.Simulation().technique
    assert not technique.implicit_collision.active

    technique.implicit_collision(multiplicity=0.5)
    assert technique.implicit_collision.active
    assert technique.implicit_collision.multiplicity == 0.5

    technique.branchless_collision()
    assert technique.implicit_collision.active
    assert technique.implicit_collision.multiplicity == 1.0

    technique.branchless_collision(active=False)
    assert not technique.implicit_collision.active
