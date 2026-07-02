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

// Standard EM base list whose models the region reproduces above the MicroElec
// thresholds. Mirrors TOPAS Ph/MicroElec/BaseList.
enum class GateMicroElecBaseList { opt3, opt4 };

// MicroElec track-structure physics overlaid on a standard EM base list
// (em_opt3 or em_opt4).
class GateG4MicroElecPhysics : public G4VPhysicsConstructor {
public:
  explicit GateG4MicroElecPhysics(G4int ver = 1);
  ~GateG4MicroElecPhysics() override = default;

  void AddRegion(const G4String &regionName);

  // MicroElec -> base-list handoff energies (Geant4 internal units). Each
  // MicroElec model is additionally capped at its data-validity ceiling.
  void SetRegionElectronThreshold(const G4String &regionName, G4double value);
  void SetRegionProtonThreshold(const G4String &regionName, G4double value);

  // Standard EM base list ("opt3" or "opt4") reproduced above the thresholds.
  void SetRegionBaseList(const G4String &regionName, const G4String &baseList);

  // opt4 only: whether e- ionisation uses G4PenelopeIonisationModel below
  // 100 keV (as in em_opt4). Ignored for opt3.
  void SetRegionUsePenelope(const G4String &regionName, G4bool usePenelope);

  void ConstructParticle() override;
  void ConstructProcess() override;

private:
  // Per-region configuration bundle (defaults mirror TOPAS).
  struct RegionConfig {
    G4double electronThreshold = 1.0 * CLHEP::keV;
    G4double protonThreshold = 2.0 * CLHEP::MeV;
    GateMicroElecBaseList baseList = GateMicroElecBaseList::opt4;
    G4bool usePenelope = true;
  };

  RegionConfig &configFor(const G4String &regionName);

  std::vector<G4String> fRegions;
  std::map<G4String, RegionConfig> fConfigs;
};

#endif
