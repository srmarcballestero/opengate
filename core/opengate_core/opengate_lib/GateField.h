/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateField_h
#define GateField_h

#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <vector>

// Non-G4 helper base: stores local-to-world transforms for one or more physical
// placements of a volume and provides the coordinate conversions shared by
// GateMagneticField (and future GateEMField).
//
// Deliberately does NOT inherit from G4Field so that concrete subclasses can
// inherit from the specific G4 type they need (G4MagneticField,
// G4ElectroMagneticField, …) without triggering diamond inheritance.
//
// Placement identification strategy: pick the transform whose local origin is
// closest to the query point in world space.  Geant4's field interface does not
// expose copy-number or touchable information, so navigator lookup is the only
// alternative — but it is prohibitively expensive inside GetFieldValue (called
// hundreds of times per step by the chord finder).  The closest-centre approach
// is exact for non-overlapping placements: if a particle is inside copy i it is
// necessarily closer to copy i's centre than to any other copy's centre, which
// is a geometric consequence of the non-overlap requirement already enforced by
// Geant4's geometry checker.
class GateField {
public:
  GateField(
    std::vector<G4ThreeVector>    translations,
    std::vector<G4RotationMatrix> rotations
  );

protected:
  // Return the local coordinates of worldPoint in the closest placement,
  // and set outTransform to point to that placement's transform.
  // Always succeeds (falls back to the first transform for a single placement).
  G4ThreeVector closestPlacement(
    const G4ThreeVector&      worldPoint,
    const G4AffineTransform*& outTransform
  ) const;

  // Rotate a field vector from local to world using the given transform.
  static G4ThreeVector rotateToWorld(
    const G4ThreeVector&     localField,
    const G4AffineTransform& transform
  );

  std::vector<G4AffineTransform> m_transforms;
};

#endif // GateField_h
