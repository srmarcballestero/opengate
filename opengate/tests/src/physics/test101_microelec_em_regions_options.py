#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# MicroElec options low_energy_ionisation and activate_{proton,alpha,ion}_models.
# Region "on" keeps the defaults; region "off" deactivates every heavy species
# and the low-energy ionisation model. The test checks the model each process
# selects in each volume.

import functools

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

_keV = 1e3  # in eV
_MeV = 1e6  # in eV


def add_target(sim, name, material, z_mm):
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    nm = gate.g4_units.nm
    vol = sim.add_volume("Box", name)
    vol.size = [2 * cm, 2 * cm, 1 * mm]
    vol.material = material
    vol.translation = [0, 0, z_mm * mm]
    return vol


def set_cuts(sim, volume_name):
    nm = gate.g4_units.nm
    region = sim.physics_manager.regions[f"{volume_name}_region"]
    region.production_cuts.electron = 1e3 * nm
    region.production_cuts.gamma = 1e3 * nm
    region.production_cuts.positron = 1e3 * nm
    region.production_cuts.proton = 1e3 * nm


def build_sim(physics_list):
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    cm = gate.g4_units.cm
    MeV = gate.g4_units.MeV

    sim.world.size = [10 * cm, 10 * cm, 10 * cm]
    sim.world.material = "Vacuum"
    sim.volume_manager.add_material_database("opengate/contrib/GateMaterials.db")
    sim.physics_manager.physics_list_name = physics_list

    source = sim.add_source("GenericSource", "source")
    source.particle = "proton"
    source.energy.mono = 70 * MeV
    source.position.type = "point"
    source.position.translation = [0, 0, -4 * cm]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.number_of_primaries = 1
    sim.add_actor("SimulationStatisticsActor", "stats")
    return sim


def probe_hook(probes, simulation_engine):
    eV = g4.G4UnitDefinition.GetValueOf("eV")
    volumes = simulation_engine.simulation.volume_manager.volumes
    simulation_engine.user_hook_log["models"] = [
        g4.check_em_model_in_volume(
            volumes[volume].g4_logical_volume, particle, process, energy_eV * eV
        )
        for volume, particle, process, energy_eV, _ in probes
    ]


def make_hook(probes):
    # picklable, for start_new_process=True
    return functools.partial(probe_hook, probes)


def check(sim, probes, label):
    models = sim.user_hook_log["models"]
    is_ok = True
    for (volume, particle, process, energy_eV, expected), model in zip(probes, models):
        if expected is None:
            passed = model is None
        elif expected.startswith("!"):
            passed = model is not None and expected[1:] not in str(model)
        else:
            passed = model is not None and expected in str(model)
        status = "PASS" if passed else "FAIL"
        print(
            f"  [{status}] [{label}] {volume}: {particle} {process} @ "
            f"{energy_eV / _keV:g} keV -> {model} (expected: {expected})"
        )
        is_ok = is_ok and passed
    return is_ok


def run_two_regions():
    """opt4 base list, one region with everything on, one with all off."""
    sim = build_sim("G4EmStandardPhysics_option4")
    add_target(sim, "on", "G4_Al", -10)
    add_target(sim, "off", "G4_Si", 10)
    sim.physics_manager.set_microelec_em_physics("on", base_list="opt4")
    sim.physics_manager.set_microelec_em_physics(
        "off",
        base_list="opt4",
        low_energy_ionisation="none",
        activate_proton_models=False,
        activate_alpha_models=False,
        activate_ion_models=False,
    )
    set_cuts(sim, "on")
    set_cuts(sim, "off")

    probes = [
        # opt4 base list outside the regions.
        ("world", "e-", "eIoni", 50 * _keV, "PenIoni"),
        ("on", "e-", "eIoni", 50 * _keV, "PenIoni"),
        ("on", "e-", "eIoni", 500 * _keV, "MollerBhabha"),
        ("off", "e-", "eIoni", 50 * _keV, "MollerBhabha"),
        # Electrons are always MicroElec.
        ("on", "e-", "e-_G4MicroElecInelastic", 500, "MicroElec"),
        ("off", "e-", "e-_G4MicroElecInelastic", 500, "MicroElec"),
        # Heavy species: MicroElec in "on", dummy model in "off".
        ("on", "proton", "p_G4MicroElecInelastic", 1 * _MeV, "MicroElec"),
        ("off", "proton", "p_G4MicroElecInelastic", 1 * _MeV, "!MicroElec"),
        ("on", "alpha", "alpha_G4MicroElecInelastic", 1 * _MeV, "MicroElec"),
        ("off", "alpha", "alpha_G4MicroElecInelastic", 1 * _MeV, "!MicroElec"),
        ("on", "GenericIon", "ion_G4MicroElecInelastic", 1 * _MeV, "MicroElec"),
        ("off", "GenericIon", "ion_G4MicroElecInelastic", 1 * _MeV, "!MicroElec"),
        ("off", "proton", "hIoni", 1 * _MeV, "Bragg"),
    ]
    sim.user_hook_after_run = make_hook(probes)
    sim.run(start_new_process=True)
    return check(sim, probes, "opt4 two regions")


def run_reference_list_all_off():
    """*_EMZ reference list, no heavy species activated: their MicroElec
    processes are not registered."""
    sim = build_sim("QGSP_BIC_EMZ")
    add_target(sim, "target", "G4_Al", 0)
    sim.physics_manager.set_microelec_em_physics(
        "target",
        base_list="opt4",
        activate_proton_models=False,
        activate_alpha_models=False,
        activate_ion_models=False,
    )
    set_cuts(sim, "target")

    probes = [
        ("world", "e-", "eIoni", 50 * _keV, "PenIoni"),
        ("target", "e-", "eIoni", 50 * _keV, "PenIoni"),
        ("target", "proton", "p_G4MicroElecInelastic", 1 * _MeV, None),
        ("target", "alpha", "alpha_G4MicroElecInelastic", 1 * _MeV, None),
        ("target", "GenericIon", "ion_G4MicroElecInelastic", 1 * _MeV, None),
    ]
    sim.user_hook_after_run = make_hook(probes)
    sim.run(start_new_process=True)
    return check(sim, probes, "QGSP_BIC_EMZ all off")


def run_opt3_inconsistent():
    """opt3 with low_energy_ionisation='penelope' is fatal."""
    sim = build_sim("G4EmStandardPhysics_option3")
    add_target(sim, "target", "G4_Al", 0)
    sim.physics_manager.set_microelec_em_physics(
        "target", base_list="opt3", low_energy_ionisation="penelope"
    )
    try:
        sim.run(start_new_process=False)
    except Exception as e:  # GATE's fatal() raises a plain Exception
        passed = "not consistent with base_list='opt3'" in str(e)
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] [opt3 + penelope] rejected: {str(e).splitlines()[0]}")
        return passed
    print("  [FAIL] [opt3 + penelope] was not rejected")
    return False


if __name__ == "__main__":
    ok = run_two_regions()
    ok = run_reference_list_all_off() and ok
    ok = run_opt3_inconsistent() and ok
    utility.test_ok(ok)
