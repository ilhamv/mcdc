import mcdc
import numpy as np
import os

os.environ["MCDC_LIB"] = "../mcdc-regression_test_data/"

# Create MC/DC simulation
simulation = mcdc.Simulation("Pincell")

# Material
# Minimal fuel/moderator composition keeps continuous-energy data preparation short.
fuel = mcdc.Material(
    number_density={
        "U235": 0.0001654509603995036,
    }
)
moderator = mcdc.Material(
    number_density={
        "H1": 0.05129627050184732,
    }
)

# Geometry
cylinder = mcdc.Surface.CylinderZ(radius=0.45720)
pitch = 1.25984
x0 = mcdc.Surface.PlaneX(x=-pitch / 2, boundary_condition="reflective")
x1 = mcdc.Surface.PlaneX(x=pitch / 2, boundary_condition="reflective")
y0 = mcdc.Surface.PlaneY(y=-pitch / 2, boundary_condition="reflective")
y1 = mcdc.Surface.PlaneY(y=pitch / 2, boundary_condition="reflective")
#
fuel_cell = mcdc.Cell(-cylinder, fill=fuel)
moderator_cell = mcdc.Cell(+x0 & -x1 & +y0 & -y1 & +cylinder, fill=moderator)
simulation.set_model([fuel_cell, moderator_cell])

# Source
source = mcdc.Source(position=[0.0, 0.0, 0.0], isotropic=True, time=0.0, energy=14.1e6)
simulation.set_sources([source])

# Setting
simulation.settings.N_particle = 2
simulation.settings.N_batch = 2
simulation.settings.time_boundary = 1e-8
simulation.settings.active_bank_buffer = 1000

# Tally
t_grid = np.insert(np.logspace(-9, -4, 200), 0, 0.0)
e_min, e_max = 1e-5, 20.0e6
groups = 500
energies = np.logspace(np.log10(e_min), np.log10(e_max), groups + 1)

tally = mcdc.Tally(scores=["flux"], time=t_grid, energy=energies)
simulation.set_tallies([tally])

simulation.run()
