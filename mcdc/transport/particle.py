from numba import njit

####

import mcdc.transport.physics as physics
import mcdc.transport.rng as rng


@njit
def move(particle_container, distance, simulation, data):
    """Move the particle along its direction and advance time by its flight duration."""
    particle = particle_container[0]
    ut = 1.0 / physics.particle_speed(particle_container, simulation, data)

    particle["x"] += particle["ux"] * distance
    particle["y"] += particle["uy"] * distance
    particle["z"] += particle["uz"] * distance
    particle["t"] += ut * distance


@njit
def reset_IDs(particle_container):
    """Clear the particle's cell, material, and surface IDs."""
    particle = particle_container[0]
    particle["cell_ID"] = -1
    particle["material_ID"] = -1
    particle["surface_ID"] = -1


@njit
def copy(target_particle_container, source_particle_container):
    """Copy all ParticleData fields, including the RNG seed, from source to target.

    WARNING: This leaves source and target with identical RNG seeds. When creating
    a child of the source, use ``copy_as_child`` to split the child's seed and
    advance the parent's RNG stream.
    """
    target_particle = target_particle_container[0]
    source_particle = source_particle_container[0]

    target_particle["x"] = source_particle["x"]
    target_particle["y"] = source_particle["y"]
    target_particle["z"] = source_particle["z"]
    target_particle["t"] = source_particle["t"]
    target_particle["ux"] = source_particle["ux"]
    target_particle["uy"] = source_particle["uy"]
    target_particle["uz"] = source_particle["uz"]
    target_particle["E"] = source_particle["E"]
    target_particle["w"] = source_particle["w"]
    target_particle["particle_type"] = source_particle["particle_type"]
    target_particle["rng_seed"] = source_particle["rng_seed"]


@njit
def copy_as_child(child_particle_container, parent_particle_container):
    """Copy the parent into the child, split the child's RNG seed, and advance the parent's RNG."""
    parent_particle = parent_particle_container[0]
    child_particle = child_particle_container[0]

    copy(child_particle_container, parent_particle_container)

    # Set child RNG seed based of the parent
    parent_seed = parent_particle["rng_seed"]
    child_particle["rng_seed"] = rng.split_seed(parent_seed, rng.SEED_SPLIT_PARTICLE)

    # Evolve parent seed
    rng.lcg(parent_particle_container)


@njit
def revive_from_secondary(target_particle_container, source_particle_container):
    """Resume the target particle with the secondary's state and the target's RNG stream."""
    target = target_particle_container[0]
    seed = target["rng_seed"]
    copy(target_particle_container, source_particle_container)
    target["rng_seed"] = seed
    target["alive"] = True
