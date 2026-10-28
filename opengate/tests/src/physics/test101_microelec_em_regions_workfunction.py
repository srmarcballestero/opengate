#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Surface work-function override: a custom "Brass" world (no MicroElec data,
# Cu 0.7 / Zn 0.3, 8.55 g/cm3, I = 324.4 eV) around a G4_Cu MicroElec region,
# with the Brass work function set to Cu's. The simulation initializes and
# Data_Brass.dat carries the supplied work function.

from pathlib import Path

import opengate as gate
import opengate_core as g4
from opengate.tests import utility

TARGET_NAME = "target"
BULK_MATERIAL = "Brass"  # custom-named, dense (8.55 g/cm3), no MicroElec data
SKIN_MATERIAL = "G4_Cu"  # database material, in the MicroElec region
WORK_FUNCTION_EV = 4.2  # Cu work function (eV): makes the Cu/Brass step zero

# Data files are written for both name forms: first three characters erased
# ("ss") and "G4_" stripped ("Brass"). G4MicroElecSurface reads the latter.
BULK_TOKEN_ERASE3 = BULK_MATERIAL[3:]
BULK_TOKEN_SURFACE = (
    BULK_MATERIAL[3:] if BULK_MATERIAL.startswith("G4_") else BULK_MATERIAL
)


def build_sim(with_override):
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    nm = gate.g4_units.nm
    eV = gate.g4_units.eV
    gcm3 = gate.g4_units.g_cm3
    MeV = gate.g4_units.MeV

    sim.volume_manager.material_database.add_material_weights(
        BULK_MATERIAL, ["Cu", "Zn"], [0.7, 0.3], 8.550 * gcm3
    )
    sim.physics_manager.material_ionisation_potential[BULK_MATERIAL] = 324.4 * eV

    sim.world.size = [10 * cm, 10 * cm, 10 * cm]
    sim.world.material = BULK_MATERIAL

    target = sim.add_volume("Box", TARGET_NAME)
    target.size = [5 * cm, 5 * cm, 1 * mm]
    target.material = SKIN_MATERIAL

    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    sim.physics_manager.set_microelec_em_physics(TARGET_NAME)

    if with_override:
        sim.physics_manager.set_microelec_surface_work_function(
            BULK_MATERIAL, WORK_FUNCTION_EV * eV
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
    return sim


def read_work_function_ev(structure_file):
    """Return the WorkFunction value (in eV) stored in a MicroElec structure file."""
    for line in structure_file.read_text().splitlines():
        tokens = line.split()
        if len(tokens) >= 4 and tokens[1] == "WorkFunction" and tokens[2] == "eV":
            return float(tokens[3])
    return None


def cleanup_generated_files(microelec_dir):
    for token in {BULK_TOKEN_ERASE3, BULK_TOKEN_SURFACE}:
        for f in microelec_dir.glob(f"**/*_{token}.dat"):
            f.unlink()


if __name__ == "__main__":
    microelec_dir = Path(g4.get_g4_data_paths()["G4LEDATA"]) / "microelec"
    cleanup_generated_files(microelec_dir)  # start from a clean state

    init_ok = False
    wf_ok = False
    sim = build_sim(with_override=True)
    try:
        sim.run(start_new_process=False)
        init_ok = True
        print(f"  [PASS] {BULK_MATERIAL} accepted with an explicit work function")
    except Exception as e:
        print(f"  [FAIL] override run did not initialize: {e}")

    if init_ok:
        structure_file = microelec_dir / "Structure" / f"Data_{BULK_TOKEN_SURFACE}.dat"
        if not structure_file.exists():
            print(f"  [FAIL] expected {structure_file} to be written")
        else:
            value = read_work_function_ev(structure_file)
            wf_ok = value is not None and abs(value - WORK_FUNCTION_EV) < 1e-6
            status = "PASS" if wf_ok else "FAIL"
            print(
                f"  [{status}] Data_{BULK_TOKEN_SURFACE}.dat WorkFunction = {value} "
                f"eV (expected {WORK_FUNCTION_EV} eV)"
            )

    cleanup_generated_files(microelec_dir)  # do not pollute the shared G4 data dir

    utility.test_ok(init_ok and wf_ok)
