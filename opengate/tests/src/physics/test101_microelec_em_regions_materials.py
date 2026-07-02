#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# MicroElec has no valid surface/cross-section treatment for a material that
# lacks MicroElec data unless that material is vacuum-like (its boundary work
# function can then be safely forced to 0). A MicroElec region surrounded by a
# dense, unsupported material (e.g. G4_WATER) therefore cannot be initialized
# and must be rejected with a clear error rather than silently assumed to be
# vacuum. This test verifies that rejection.

import opengate as gate
from opengate.tests import utility

TARGET_NAME = "target"


def build_sim(world_material):
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    nm = gate.g4_units.nm
    MeV = gate.g4_units.MeV

    sim.volume_manager.add_material_database("opengate/contrib/GateMaterials.db")

    sim.world.size = [10 * cm, 10 * cm, 10 * cm]
    sim.world.material = world_material

    target = sim.add_volume("Box", TARGET_NAME)
    target.size = [5 * cm, 5 * cm, 1 * mm]
    target.material = "G4_Al"

    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    sim.physics_manager.set_microelec_em_physics(TARGET_NAME)

    region = sim.physics_manager.regions[f"{TARGET_NAME}_region"]
    region.production_cuts.electron = 1e3 * nm
    region.production_cuts.gamma = 1e3 * nm
    region.production_cuts.positron = 1e3 * nm
    region.production_cuts.proton = 1e3 * nm

    source = sim.add_source("GenericSource", "source")
    source.particle = "proton"
    source.energy.mono = 70 * MeV
    source.position.type = "point"
    source.position.translation = [0, 0, -4 * cm]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.n = 1

    sim.add_actor("SimulationStatisticsActor", "stats")
    return sim


if __name__ == "__main__":
    # A MicroElec region surrounded by dense, unsupported G4_WATER must be
    # rejected during initialization.
    sim = build_sim("G4_WATER")

    is_ok = False
    try:
        sim.run(start_new_process=False)
        print(
            "  [FAIL] expected a fatal for a MicroElec region in dense G4_WATER, "
            "but the simulation initialized successfully"
        )
    except Exception as e:
        message = str(e)
        expected = "not vacuum-like" in message and "G4_WATER" in message
        status = "PASS" if expected else "FAIL"
        print(f"  [{status}] dense unsupported surrounder G4_WATER is rejected")
        print(f"    error: {message}")
        is_ok = expected

    utility.test_ok(is_ok)
