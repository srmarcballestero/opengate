/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateMagneticField_h
#define GateMagneticField_h

#include "GateField.h"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <vector>

// Local-coordinate wrapper for any G4MagneticField.
//
// Geant4 built-in fields (G4QuadrupoleMagField, …) evaluate GetFieldValue in
// world coordinates — their "centre" is the world origin.  This wrapper converts
// the incoming world-space query point to the local frame of the physical
// volume(s) the field is attached to, delegates to the inner field, then rotates
// the result back to world coordinates.
//
// For repeated placements the correct copy is identified by finding the one
// whose local origin is closest to the query point (see GateField for the
// geometric justification).
//
// Inheriting from G4MagneticField (not G4Field) lets the object be passed to
// G4Mag_UsualEqRhs.  GateField (no G4 base) provides the transform logic;
// there is no diamond-inheritance issue.
class GateMagneticField : public G4MagneticField, protected GateField {
public:
  // inner – wrapped built-in field (not owned; caller must keep it alive)
  GateMagneticField(
    G4MagneticField*              inner,
    std::vector<G4ThreeVector>    translations,
    std::vector<G4RotationMatrix> rotations
  );

  void GetFieldValue(const G4double Point[4], G4double* Bfield) const override;

private:
  G4MagneticField* m_inner;  // not owned
};

#endif // GateMagneticField_h
