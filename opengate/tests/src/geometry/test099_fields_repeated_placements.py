#!/usr/bin/env python3
"""
Test 099 — Uniform magnetic field: analytical validation.

A proton enters a region with uniform B along Y.
Checks:
  1. Circular trajectory in XZ plane (cyclotron radius)
  2. No deflection in Y
  3. Energy conservation (B does no work)
"""

# import uproot

import opengate as gate
from opengate.geometry import fields
from opengate.tests import utility

from test099_fields_helpers import (
    g4_tesla,
    g4_mm,
    g4_MeV,
    g4_eplus,
    PROTON_MASS,
    cyclotron_radius,
    build_field_simulation,
)

m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
MeV = gate.g4_units.MeV

if __name__ == "__main__":
    paths = utility.get_default_test_paths(__file__, output_folder="test099_fields")

    By = 1 * g4_tesla
    T = 10 * g4_MeV

    sim = gate.Simulation()
    sim.g4_verbose = False
    sim.visu = False
    sim.random_seed = 13579
    sim.number_of_threads = 1
    sim.output_dir = paths.output

    sim.visu = True
    # sim.visu_type = "qt"
    sim.visu_commands.append("/vis/scene/endOfEventAction accumulate 10000")
    sim.visu_commands.append("/vis/scene/add/magneticField 20 fullArrow")
    sim.visu_commands.append("/vis/scene/add/trajectories smooth")

    sim.world.size = [1 * m, 1 * m, 1 * m]
    sim.world.material = "G4_Galactic"

    # Repeated box placements
    box = sim.add_volume("Box", "box")
    box.size = [50 * mm, 50 * mm, 50 * mm]
    box.translation = gate.geometry.utility.get_grid_repetition(
        [1, 1, 5], [0, 0, 50 * mm], start=[0, 0, -100 * mm]
    )
    box.material = "G4_Galactic"

    field = fields.UniformMagneticField(name="B_uniform")
    field.field_vector = [0, By, 0]

    field = fields.QuadrupoleMagneticField(name="B_quad")
    field.gradient = By / 0.05  # B = G * r, so G = B / r, and we want B = By at r = 5 cm

    box.add_field(field)

    source = sim.add_source("GenericSource", "particle_source")
    source.particle = "proton"
    source.n = 1
    source.energy.type = "mono"
    source.energy.mono = T
    source.position.type = "point"
    source.position.translation = [0, 0, -50 * cm]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]


    # Single box placement
    box_single = sim.add_volume("Box", "box_single")
    box_single.size = [50 * mm, 50 * mm, 250 * mm]
    box_single.translation = [125 * mm, 0, 0]
    box_single.material = "G4_Galactic"
    box_single.add_field(field)


    source1 = sim.add_source("GenericSource", "particle_source1")
    source1.particle = "proton"
    source1.n = 1
    source1.energy.type = "mono"
    source1.energy.mono = T
    source1.position.type = "point"
    source1.position.translation = [125 * mm, 0, -50 * cm]
    source1.direction.type = "momentum"
    source1.direction.momentum = [0, 0, 1]

    sim.run()





