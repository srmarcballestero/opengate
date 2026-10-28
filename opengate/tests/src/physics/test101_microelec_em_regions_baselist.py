#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Standard models selected in a MicroElec region above the thresholds, for the
# opt3 and opt4 base lists:
#   - proton hIoni: Bragg below 2 MeV, BetheBloch above
#   - e- eIoni:     opt4: Penelope below 100 keV, MollerBhabha above
#                   opt3: MollerBhabha
# (e- msc is not checked: check_em_model_in_volume does not handle
# G4VMultipleScattering.)

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

TARGET_NAME = "target"

# Below the 2 MeV Bragg/BetheBloch crossover: Bragg is used in [P_MAX, 2 MeV].
P_MAX_MeV = 0.5

_MeV_in_eV = 1e6

# (particle, process, energy_eV, expected model-name substring for opt4 / opt3)
PROBE_POINTS = [
    ("proton", "hIoni", 1.0 * _MeV_in_eV, "Bragg", "Bragg"),
    ("proton", "hIoni", 5.0 * _MeV_in_eV, "BetheBloch", "BetheBloch"),
    ("e-", "eIoni", 0.05 * _MeV_in_eV, "PenIoni", "MollerBhabha"),
    ("e-", "eIoni", 0.5 * _MeV_in_eV, "MollerBhabha", "MollerBhabha"),
]


def microelec_hook(simulation_engine):
    eV_unit = g4.G4UnitDefinition.GetValueOf("eV")
    target = simulation_engine.simulation.volume_manager.volumes[TARGET_NAME]
    simulation_engine.user_hook_log["models"] = [
        g4.check_em_model_in_volume(
            target.g4_logical_volume, particle, process, energy_eV * eV_unit
        )
        for particle, process, energy_eV, *_ in PROBE_POINTS
    ]


def check_baselist(sim, expected_idx, base_list):
    models = sim.user_hook_log["models"]
    is_ok = True
    for (particle, process, energy_eV, *expected), model in zip(PROBE_POINTS, models):
        want = expected[expected_idx]
        passed = model is not None and want in str(model)
        label = f"[{base_list}] {process} @ {energy_eV / _MeV_in_eV:g} MeV"
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}")
        print(f"    model   : {model}")
        print(f"    expected: contains '{want}'")
        is_ok = is_ok and passed
    return is_ok


def run(base_list, physics_list, expected_idx):
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

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

    sim.physics_manager.physics_list_name = physics_list
    sim.physics_manager.set_microelec_em_physics(
        TARGET_NAME,
        proton_threshold=P_MAX_MeV * MeV,
        base_list=base_list,
    )

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
    source.number_of_primaries = 1

    sim.add_actor("SimulationStatisticsActor", "stats")
    sim.user_hook_after_run = microelec_hook
    sim.run(start_new_process=True)

    return check_baselist(sim, expected_idx, base_list)


if __name__ == "__main__":
    ok_opt4 = run("opt4", "G4EmStandardPhysics_option4", expected_idx=0)
    ok_opt3 = run("opt3", "G4EmStandardPhysics_option3", expected_idx=1)
    utility.test_ok(ok_opt4 and ok_opt3)
