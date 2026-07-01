#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Verify that, above the MicroElec thresholds, a MicroElec region reproduces the
# selected standard EM base list (opt3 or opt4), mirroring TOPAS
# TsEmMicroElecPhysics. We probe which *standard* model each process selects in
# the target region:
#   - proton hIoni: Bragg below 2 MeV, BetheBloch above (both base lists)
#   - e- eIoni:     opt4 -> Penelope below 100 keV then MollerBhabha;
#                   opt3 -> MollerBhabha everywhere
# (The e- msc Urban/GoudsmitSaunderson difference is not probed here: msc is a
# G4VMultipleScattering, which the check_em_model_in_volume helper does not
# resolve. The eIoni Penelope-vs-MollerBhabha split already distinguishes the
# two base lists.)

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

TARGET_NAME = "target"

# Lower the proton threshold below the 2 MeV Bragg/BetheBloch crossover so a
# Bragg window [P_MAX, 2 MeV] opens inside the region.
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
    results = {}
    for particle, process, energy_eV, *_ in PROBE_POINTS:
        # These standard models are active above threshold, so the plain
        # (name-selected) model query is what we want here.
        results[(process, energy_eV)] = g4.check_em_model_in_volume(
            target.g4_logical_volume, particle, process, energy_eV * eV_unit
        )
    simulation_engine.user_hook_log.append(results)


def check_baselist(sim, expected_idx, base_list):
    results = sim.user_hook_log[0]
    is_ok = True
    for particle, process, energy_eV, *expected in PROBE_POINTS:
        want = expected[expected_idx]
        model = results.get((process, energy_eV))
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
    source.n = 1

    sim.add_actor("SimulationStatisticsActor", "stats")
    sim.user_hook_after_run = microelec_hook
    # A fresh process per base list so both configurations can run in one script.
    sim.run(start_new_process=True)

    return check_baselist(sim, expected_idx, base_list)


if __name__ == "__main__":
    ok_opt4 = run("opt4", "G4EmStandardPhysics_option4", expected_idx=0)
    ok_opt3 = run("opt3", "G4EmStandardPhysics_option3", expected_idx=1)
    utility.test_ok(ok_opt4 and ok_opt3)
