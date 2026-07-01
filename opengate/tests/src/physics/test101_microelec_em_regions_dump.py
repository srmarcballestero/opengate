#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Reproduce some microelec configurations and dump the geant4 output to files,
# to be reviewed manually.

import os
import sys
import contextlib

import opengate as gate
from opengate.tests import utility

TARGET_NAME = "target"
CONFIGS = [
    # name         physics_list              base   e-thr  penelope
    ("opt3", "G4EmStandardPhysics_option3", "opt3", None, None),
    ("opt3_1keV", "G4EmStandardPhysics_option3", "opt3", 1.0, None),
    ("opt3_10keV", "G4EmStandardPhysics_option3", "opt3", 10.0, None),
    ("opt4", "G4EmStandardPhysics_option4", "opt4", None, None),
    ("opt4_1keV", "G4EmStandardPhysics_option4", "opt4", 1.0, None),
    ("opt4_10keV", "G4EmStandardPhysics_option4", "opt4", 10.0, None),
    ("opt4_nopen_1keV", "G4EmStandardPhysics_option4", "opt4", 1.0, False),
]


@contextlib.contextmanager
def redirect_all_output_to(path):
    """Redirect OS-level stdout/stderr (fd 1 & 2) to `path`.

    Redirecting at the file-descriptor level captures the C++ Geant4 output as
    well as Python prints, and is inherited by the child process spawned by
    sim.run(start_new_process=True).
    """
    sys.stdout.flush()
    sys.stderr.flush()
    saved_out = os.dup(1)
    saved_err = os.dup(2)
    with open(path, "w") as f:
        try:
            os.dup2(f.fileno(), 1)
            os.dup2(f.fileno(), 2)
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_out, 1)
            os.dup2(saved_err, 2)
            os.close(saved_out)
            os.close(saved_err)


def build_sim(physics_list, base_list, thr_keV, use_penelope):
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456

    # Dump the Geant4 EM output.
    sim.g4_verbose = True
    sim.g4_verbose_level_tracking = -1  # keep the per-step tracking dump off
    sim.g4_commands_before_init.append("/process/em/verbose 1")
    sim.g4_commands_after_init.append("/process/em/printParameters")

    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    nm = gate.g4_units.nm
    eV = gate.g4_units.eV
    keV = gate.g4_units.keV
    MeV = gate.g4_units.MeV
    GeV = gate.g4_units.GeV

    sim.world.size = [10 * cm, 10 * cm, 10 * cm]
    sim.world.material = "Vacuum"

    sim.volume_manager.add_material_database("opengate/contrib/GateMaterials.db")

    # G4_Al cylinder, RMax = 1 mm, half-length = 1 mm (TOPAS TsCylinder).
    target = sim.add_volume("Tubs", TARGET_NAME)
    target.rmin = 0
    target.rmax = 1 * mm
    target.dz = 1 * mm
    target.material = "G4_Al"

    sim.physics_manager.physics_list_name = physics_list

    sim.physics_manager.energy_range_min = 0.1 * eV
    sim.physics_manager.energy_range_max = 10 * GeV

    if thr_keV is not None:
        kwargs = {
            "electron_threshold": thr_keV * keV,
            "base_list": base_list,
        }
        if use_penelope is not None:
            kwargs["use_penelope"] = use_penelope
        sim.physics_manager.set_microelec_em_physics(TARGET_NAME, **kwargs)

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
    paths = utility.get_default_test_paths(__file__)

    for name, physics_list, base_list, thr_keV, use_penelope in CONFIGS:
        dump_file = paths.output / f"microelec_dump_{name}.log"
        print(f"[{name}] -> {dump_file}")
        sim = build_sim(physics_list, base_list, thr_keV, use_penelope)
        with redirect_all_output_to(dump_file):
            # write a header into the dump file itself
            print("=" * 78)
            print(f"CONFIG: {name}")
            print(f"  physics_list        = {physics_list}")
            print(f"  microelec base_list = {base_list}")
            print(
                f"  electron_threshold  = {thr_keV} keV"
                + ("  (baseline, no MicroElec)" if thr_keV is None else "")
            )
            print(f"  use_penelope        = {use_penelope}")
            print("=" * 78)
            sys.stdout.flush()
            # fresh process per config: physics list / EM parameters are global.
            sim.run(start_new_process=True)

    print(f"\nWrote {len(CONFIGS)} dump files to {paths.output}")
