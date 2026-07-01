#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Verify the two MicroElec thresholds (electron / proton) and the automatic
# capping of each MicroElec model at its data-validity ceiling:
#   e- elastic   valid to 500 keV
#   e- inelastic valid to 10 keV   <- caps below the electron threshold here
#   proton inelastic valid to 10 MeV
# With electron_threshold = 50 keV the elastic model reaches 50 keV, but the
# inelastic model is still capped at its 10 keV ceiling; with
# proton_threshold = 5 MeV the proton inelastic model reaches 5 MeV.

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

TARGET_NAME = "target"

ELECTRON_THRESHOLD_keV = 50.0  # elastic reaches this; inelastic capped at 10 keV
PROTON_THRESHOLD_MeV = 5.0

_MeV_in_eV = 1e6

# (particle, process, energy_eV, expect MicroElec active in target)
PROBE_POINTS = [
    # elastic: active up to the 50 keV threshold (< 500 keV validity)
    ("e-", "e-_G4MicroElecElastic", 0.030 * _MeV_in_eV, True),
    ("e-", "e-_G4MicroElecElastic", 0.100 * _MeV_in_eV, False),
    # inelastic: capped at its 10 keV ceiling, below the 50 keV threshold
    ("e-", "e-_G4MicroElecInelastic", 0.005 * _MeV_in_eV, True),
    ("e-", "e-_G4MicroElecInelastic", 0.030 * _MeV_in_eV, False),
    # proton inelastic: active up to the 5 MeV threshold
    ("proton", "p_G4MicroElecInelastic", 1.0 * _MeV_in_eV, True),
    ("proton", "p_G4MicroElecInelastic", 8.0 * _MeV_in_eV, False),
]


def microelec_hook(simulation_engine):
    eV_unit = g4.G4UnitDefinition.GetValueOf("eV")
    target = simulation_engine.simulation.volume_manager.volumes[TARGET_NAME]
    results = {}
    for particle, process, energy_eV, _ in PROBE_POINTS:
        # An energy limit deactivates the MicroElec model without changing which
        # model is name-selected, so query the *active* model.
        results[(process, energy_eV)] = g4.check_active_em_model_in_volume(
            target.g4_logical_volume, particle, process, energy_eV * eV_unit
        )
    simulation_engine.user_hook_log.append(results)


def check_microelec_limits(sim):
    results = sim.user_hook_log[0]
    is_ok = True
    for particle, process, energy_eV, expect_microelec in PROBE_POINTS:
        model = results.get((process, energy_eV))
        is_microelec = model is not None and "MicroElec" in str(model)
        passed = is_microelec == expect_microelec
        side = "active" if expect_microelec else "inactive"
        label = f"{process} @ {energy_eV / _MeV_in_eV:g} MeV ({side})"
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}")
        print(f"    model   : {model}")
        print(f"    expected: {'MicroElec' if expect_microelec else 'not MicroElec'}")
        is_ok = is_ok and passed
    return is_ok


VISU = False

if __name__ == "__main__":
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    nm = gate.g4_units.nm
    keV = gate.g4_units.keV
    MeV = gate.g4_units.MeV

    sim.world.size = [10 * cm, 10 * cm, 10 * cm]
    sim.world.material = "Vacuum"

    sim.volume_manager.add_material_database("opengate/contrib/GateMaterials.db")

    target = sim.add_volume("Box", TARGET_NAME)
    target.size = [5 * cm, 5 * cm, 1 * mm]
    target.material = "G4_Al"

    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    # Activate MicroElec in the target with custom handoff thresholds.
    sim.physics_manager.set_microelec_em_physics(
        TARGET_NAME,
        electron_threshold=ELECTRON_THRESHOLD_keV * keV,
        proton_threshold=PROTON_THRESHOLD_MeV * MeV,
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
    # Material-cuts couples are only assigned during the run, so query the models
    # after the run (not after init).
    sim.user_hook_after_run = microelec_hook

    sim.run(start_new_process=False)

    is_ok = check_microelec_limits(sim)
    utility.test_ok(is_ok)
