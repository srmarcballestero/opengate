#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

TARGET_AUTO_NAME = "target_auto"
TARGET_EXPLICIT_NAME = "target_explicit"
EXPLICIT_REGION_NAME = "microelec_region_explicit"


def check_microelec_regions(sim):
    hook_output = sim.user_hook_log[0]
    microelec_regions = hook_output["microelec_regions"]
    volume_regions = hook_output["volume_regions"]
    model_checks = sim.user_hook_log[1]

    target_auto_region = f"{TARGET_AUTO_NAME}_region"

    print("Checking configured MicroElec regions:")
    print(microelec_regions)
    print("Checking volume-to-region association:")
    print(volume_regions)

    checks = [
        (
            f"MicroElec region registered for {target_auto_region}",
            target_auto_region in microelec_regions,
            True,
        ),
        (
            f"MicroElec region registered for {EXPLICIT_REGION_NAME}",
            EXPLICIT_REGION_NAME in microelec_regions,
            True,
        ),
        (
            f"Region attached to {TARGET_AUTO_NAME}",
            volume_regions.get(TARGET_AUTO_NAME),
            target_auto_region,
        ),
        (
            f"Region attached to {TARGET_EXPLICIT_NAME}",
            volume_regions.get(TARGET_EXPLICIT_NAME),
            EXPLICIT_REGION_NAME,
        ),
        (
            "World region is not a MicroElec region",
            "DefaultRegionForTheWorld" in microelec_regions,
            False,
        ),
        (
            f"MicroElec elastic model active in {TARGET_AUTO_NAME}",
            model_checks.get(TARGET_AUTO_NAME),
            "contains MicroElec",
        ),
        (
            f"World volume uses standard model (not MicroElec)",
            model_checks.get("world"),
            "not MicroElec",
        ),
    ]

    is_ok = True
    print("Detailed MicroElec region checks:")
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

    print("Diagnostic EM-model checks per volume:")
    for volume_name, model_name in model_checks.items():
        print(f"  {volume_name}: {model_name}")

    return is_ok


def microelec_hook(simulation_engine):
    em = simulation_engine.physics_engine.g4_em_parameters
    microelec_regions = set(str(r) for r in em.RegionsMicroElec())

    volume_regions = {}
    for volume_name, volume in simulation_engine.simulation.volume_manager.volumes.items():
        region_name = None
        if volume.g4_region is not None:
            region_name = volume.g4_region.GetName()
        volume_regions[volume_name] = region_name

    simulation_engine.user_hook_log.append(
        {
            "microelec_regions": microelec_regions,
            "volume_regions": volume_regions,
        }
    )

    # Check which EM model is active for e- elastic in each volume
    eV = g4.G4UnitDefinition.GetValueOf("eV")
    model_checks = {}
    for volume_name, volume in simulation_engine.simulation.volume_manager.volumes.items():
        model_checks[volume_name] = g4.check_em_model_in_volume(
            volume.g4_logical_volume,
            "e-",
            "e-_G4MicroElecElastic",
            1000.0 * eV,
        )
    simulation_engine.user_hook_log.append(model_checks)


if __name__ == "__main__":
    sim = gate.Simulation()

    sim.g4_verbose = True
    sim.visu = True
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    m = gate.g4_units.m
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    MeV = gate.g4_units.MeV

    sim.world.size = [2 * m, 2 * m, 2 * m]
    sim.world.material = "G4_AIR"

    # Auto-region target (Si)
    target_auto = sim.add_volume("Box", TARGET_AUTO_NAME)
    target_auto.size = [2 * cm, 2 * cm, 2 * cm]
    target_auto.translation = [-4 * cm, 0, 0]
    target_auto.material = "G4_Si"

    # Explicit-region target (Si)
    target_explicit = sim.add_volume("Box", TARGET_EXPLICIT_NAME)
    target_explicit.size = [2 * cm, 2 * cm, 2 * cm]
    target_explicit.translation = [4 * cm, 0, 0]
    target_explicit.material = "G4_Si"

    # Standard EM base list — MicroElec overrides its models in the target regions
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"

    # Auto-region via volume name convenience method
    sim.physics_manager.set_microelec_em_physics(TARGET_AUTO_NAME, "MicroElec")

    # Explicit region
    explicit_region = sim.physics_manager.add_region(EXPLICIT_REGION_NAME)
    explicit_region.associate_volume(target_explicit)
    sim.physics_manager.set_microelec_em_physics_in_region(EXPLICIT_REGION_NAME, "MicroElec")

    source = sim.add_source("GenericSource", "source")
    source.particle = "e-"
    source.energy.mono = 1 * MeV
    source.position.type = "point"
    source.position.translation = [0, 0, 0]
    source.direction.type = "iso"
    source.n = 1

    sim.add_actor("SimulationStatisticsActor", "stats")

    sim.user_hook_after_init = microelec_hook

    sim.run(start_new_process=False)

    is_ok = check_microelec_regions(sim)
    utility.test_ok(is_ok)
