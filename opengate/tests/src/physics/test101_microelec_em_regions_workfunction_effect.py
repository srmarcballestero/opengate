#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Verify that Geant4 actually CONSUMES the MicroElec surface work-function
# override -- not merely that GATE writes the data file (that is checked by
# test101_microelec_em_regions_workfunction.py).
#
# G4MicroElecSurface applies the work-function difference across an interface as
# a barrier: an electron crossing Cu -> Brass has WF_Cu - WF_Brass removed from
# its kinetic energy, and that energy is NOT deposited anywhere. So the total
# energy deposited in the geometry decreases as the barrier grows. Cu's tabulated
# work function is 4.2 eV, so:
#   - Brass override = 4.2 eV  -> barrier 0     -> ~all primary energy deposited
#   - Brass override = 0.02 eV -> barrier ~4.2  -> less energy deposited
#
# We fire 50 eV electrons (well above the barrier, so they cross) from a thin Cu
# slab (MicroElec region) into a surrounding Brass block, with IDENTICAL primary
# histories (same seed) for the two runs, so the ONLY difference is the Brass
# work function. A monotonic, significant change in the deposited energy proves
# the surface process reads our per-material value; if it ignored the override,
# the two totals would be identical.

from pathlib import Path

import itk
import opengate as gate
from opengate.tests import utility

paths = utility.get_default_test_paths(
    __file__, output_folder="test101_microelec_wf_effect"
)

E0_eV = 50.0  # primary electron energy (above the ~4.2 eV barrier)
N_PRIMARIES = 20000
SEED = 20260707


def build_sim(brass_wf_eV, tag):
    sim = gate.Simulation()
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = SEED  # same histories for both runs
    sim.progress_bar = False
    sim.output_dir = paths.output

    nm = gate.g4_units.nm
    um = gate.g4_units.um
    eV = gate.g4_units.eV
    gcm3 = gate.g4_units.g_cm3

    sim.volume_manager.add_material_database("opengate/contrib/GateMaterials.db")
    sim.volume_manager.material_database.add_material_weights(
        "Brass", ["Cu", "Zn"], [0.7, 0.3], 8.550 * gcm3
    )
    sim.physics_manager.material_ionisation_potential["Brass"] = 324.4 * eV

    sim.world.size = [1 * um, 1 * um, 1 * um]
    sim.world.material = "Vacuum"

    # A thin Cu slab (the MicroElec region) embedded in a Brass block (the bulk).
    brass = sim.add_volume("Box", "brass")
    brass.size = [400 * nm, 400 * nm, 400 * nm]
    brass.material = "Brass"

    cu = sim.add_volume("Box", "cu")
    cu.mother = "brass"
    cu.size = [200 * nm, 200 * nm, 10 * nm]
    cu.material = "G4_Cu"

    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    sim.physics_manager.set_microelec_em_physics("cu")
    sim.physics_manager.set_microelec_surface_work_function("Brass", brass_wf_eV * eV)

    region = sim.physics_manager.regions["cu_region"]
    for c in ("electron", "gamma", "positron", "proton"):
        setattr(region.production_cuts, c, 0.1 * nm)

    # Electrons born just inside the +z face of Cu, aimed at the Cu/Brass boundary.
    src = sim.add_source("GenericSource", "src")
    src.particle = "e-"
    src.energy.mono = E0_eV * eV
    src.position.type = "box"
    src.position.size = [50 * nm, 50 * nm, 0.5 * nm]
    src.position.translation = [0, 0, 4.5 * nm]  # 0.5 nm from the +z Cu face (z=5nm)
    src.direction.type = "momentum"
    src.direction.momentum = [0, 0, 1]
    src.n = N_PRIMARIES

    # Total energy deposited: a single voxel spanning the whole world. The barrier
    # energy removed at the interface never appears here, so this shrinks with the
    # barrier.
    dose = sim.add_actor("DoseActor", "dose_total")
    dose.attached_to = "world"
    dose.size = [1, 1, 1]
    dose.spacing = [1 * um, 1 * um, 1 * um]
    dose.output_filename = f"edep_{tag}.mhd"

    sim.add_actor("SimulationStatisticsActor", "stats")
    return sim, dose


def run_mean_edep_eV(brass_wf_eV, tag):
    sim, dose = build_sim(brass_wf_eV, tag)
    sim.run(start_new_process=True)
    total_MeV = float(
        itk.array_from_image(itk.imread(str(dose.edep.get_output_path()))).sum()
    )
    return total_MeV / N_PRIMARIES * 1e6  # mean eV deposited per primary


if __name__ == "__main__":
    # Barrier ~0: Brass work function set equal to Cu's (4.2 eV).
    mean_no_barrier = run_mean_edep_eV(4.2, "matched")
    # Barrier ~4.2 eV: Brass work function driven to nearly zero.
    mean_barrier = run_mean_edep_eV(0.02, "lowered")

    drop = mean_no_barrier - mean_barrier
    print(f"  E0 = {E0_eV} eV, {N_PRIMARIES} primaries, identical histories")
    print(f"  mean edep, Brass WF = 4.2 eV (barrier 0)   : {mean_no_barrier:.3f} eV")
    print(f"  mean edep, Brass WF = 0.02 eV (barrier ~4.2): {mean_barrier:.3f} eV")
    print(f"  drop attributable to the work-function barrier: {drop:.3f} eV")

    # With the barrier removed, almost all the 50 eV is deposited.
    retained_ok = mean_no_barrier > 49.0
    # Raising the barrier measurably removes energy at the interface.
    barrier_ok = mean_barrier < 48.0
    # The effect is significant and in the expected direction: only the Brass work
    # function changed between the two runs, so any change proves Geant4 uses it.
    effect_ok = drop > 1.0

    for label, ok in [
        ("barrier-0 run retains ~all energy", retained_ok),
        ("raising the barrier lowers deposited energy", barrier_ok),
        ("effect is significant (Geant4 consumes the override)", effect_ok),
    ]:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")

    utility.test_ok(retained_ok and barrier_ok and effect_ok)
