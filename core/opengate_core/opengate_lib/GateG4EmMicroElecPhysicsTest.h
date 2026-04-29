/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateG4EmMicroElecPhysicsTest_h
#define GateG4EmMicroElecPhysicsTest_h

#include "G4String.hh"
#include "G4VPhysicsConstructor.hh"
#include <vector>

class GateG4EmMicroElecPhysicsTest : public G4VPhysicsConstructor {
public:
  explicit GateG4EmMicroElecPhysicsTest(G4int ver = 1);
  ~GateG4EmMicroElecPhysicsTest() override = default;

  void AddRegion(const G4String& regionName);

  void ConstructParticle() override;
  void ConstructProcess() override;

private:
  std::vector<G4String> fRegions;
};

#endif