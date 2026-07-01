#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

TARGET_NAME = "target"


def microelec_hook(simulation_engine):
    eV = g4.G4UnitDefinition.GetValueOf("eV")
    results = {}
    for (
        volume_name,
        volume,
    ) in simulation_engine.simulation.volume_manager.volumes.items():
        results[volume_name] = g4.check_em_model_in_volume(
            volume.g4_logical_volume,
            "e-",
            "e-_G4MicroElecElastic",
            1000.0 * eV,
        )
    simulation_engine.user_hook_log.append(results)


def check_microelec_models(sim):
    model_checks = sim.user_hook_log[0]
    target_model = model_checks.get(TARGET_NAME)
    world_model = model_checks.get("world")

    checks = [
        (
            f"MicroElec elastic model active in {TARGET_NAME}",
            target_model,
            "contains MicroElec",
        ),
        (
            "World volume uses DummyModel (not MicroElec)",
            world_model,
            "not MicroElec",
        ),
    ]

    is_ok = True
    for label, actual, expected in checks:
        if expected == "contains MicroElec":
            passed = actual is not None and "MicroElec" in str(actual)
        elif expected == "not MicroElec":
            passed = actual is None or "MicroElec" not in str(actual)
        else:
            passed = actual == expected
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}")
        print(f"    actual  : {actual}")
        print(f"    expected: {expected}")
        is_ok = is_ok and passed

    return is_ok


VISU = False

if __name__ == "__main__":
    sim = gate.Simulation()
    sim.g4_verbose = True
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    if VISU:
        sim.visu = True
        sim.visu_type = "qt"
        sim.visu_commands.append("/vis/scene/endOfEventAction accumulate")
        sim.visu_commands.append("/vis/scene/add/trajectories smooth")
        sim.visu_commands.append("/vis/scene/add/magneticField 15 fullArrow")

    m = gate.g4_units.m
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    nm = gate.g4_units.nm
    MeV = gate.g4_units.MeV

    sim.world.size = [10 * cm, 10 * cm, 10 * cm]
    sim.world.material = "Vacuum"

    sim.volume_manager.add_material_database("opengate/contrib/GateMaterials.db")

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
    sim.user_hook_after_run = microelec_hook

    sim.run(start_new_process=False)

    is_ok = check_microelec_models(sim)
    utility.test_ok(is_ok)
