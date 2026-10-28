/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateG4MicroElecPhysics_h
#define GateG4MicroElecPhysics_h

#include "CLHEP/Units/SystemOfUnits.h"
#include "G4String.hh"
#include "G4VPhysicsConstructor.hh"
#include <map>
#include <vector>

// Standard EM list whose models a MicroElec region uses above the thresholds.
enum class GateMicroElecBaseList { opt3, opt4 };

// e- ionisation model between the MicroElec electron threshold and 100 keV.
// "automatic": Penelope for opt4, none (MollerBhabha only) for opt3.
enum class GateMicroElecLowEnergyIonisation { automatic, penelope, none };

// MicroElec track-structure physics in selected regions, on top of a standard
// EM physics list (opt3 or opt4).
class GateG4MicroElecPhysics : public G4VPhysicsConstructor {
public:
  explicit GateG4MicroElecPhysics(G4int ver = 1);
  ~GateG4MicroElecPhysics() override = default;

  void AddRegion(const G4String &regionName);

  // Energies below which the MicroElec models are used, capped at the MicroElec
  // data limits.
  void SetRegionElectronThreshold(const G4String &regionName, G4double value);
  void SetRegionProtonThreshold(const G4String &regionName, G4double value);

  // "opt3" or "opt4".
  void SetRegionBaseList(const G4String &regionName, const G4String &baseList);

  // "penelope", "none" or "automatic". opt3 accepts only "none" and
  // "automatic".
  void SetRegionLowEnergyIonisation(const G4String &regionName,
                                    const G4String &model);

  // Per-species MicroElec switches (electrons always use MicroElec). A
  // deactivated species keeps the base-list models inside the region.
  void SetRegionActivateProtonModels(const G4String &regionName, G4bool value);
  void SetRegionActivateAlphaModels(const G4String &regionName, G4bool value);
  void SetRegionActivateIonModels(const G4String &regionName, G4bool value);

  void ConstructParticle() override;
  void ConstructProcess() override;

private:
  struct RegionConfig {
    G4double electronThreshold = 1.0 * CLHEP::keV;
    G4double protonThreshold = 2.0 * CLHEP::MeV;
    GateMicroElecBaseList baseList = GateMicroElecBaseList::opt4;
    GateMicroElecLowEnergyIonisation lowEnergyIonisation =
        GateMicroElecLowEnergyIonisation::automatic;
    G4bool activateProtonModels = true;
    G4bool activateAlphaModels = true;
    G4bool activateIonModels = true;
  };

  RegionConfig &configFor(const G4String &regionName);

  // Resolves "automatic"; G4Exception if inconsistent with the base list.
  GateMicroElecLowEnergyIonisation
  ResolveLowEnergyIonisation(const G4String &regionName,
                             const RegionConfig &cfg) const;

  std::vector<G4String> fRegions;
  std::map<G4String, RegionConfig> fConfigs;
};

#endif
