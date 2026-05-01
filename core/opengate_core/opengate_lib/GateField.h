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

class G4VSolid;

// Non-G4 helper base: stores local-to-world transforms for one or more physical
// placements of a logical volume and provides the coordinate conversions shared
// by GateMagneticField (and future GateEMField).
//
// Deliberately does NOT inherit from G4Field so that concrete subclasses can
// inherit from the specific G4 type they need (G4MagneticField,
// G4ElectroMagneticField, …) without triggering diamond inheritance.
//
// Placement identification: for each candidate placement, transform the world
// query point into the placement's local frame and ask the shared G4VSolid
// whether that local point lies inside.  This is exact for any solid
// (convex or not, axis-aligned or rotated) — the only general test that
// works once Geant4's chord finder strips touchable information from the
// query (GetFieldValue receives only (x,y,z,t)).
//
// Why not closest-centre? It is wrong for rotated, elongated, or non-convex
// placements (e.g. a long ellipsoid rotated 90° next to its sibling: a point
// inside one copy can be closer to the other copy's centre).  Inside() is
// the only criterion that always agrees with Geant4's own definition of
// "which volume contains this point".
//
// Robust fallback: at chord-finder sub-step positions slightly outside the
// volume (a normal occurrence near boundaries) no copy reports kInside.  In
// that case we pick the placement whose local origin is closest to the query
// point — this keeps the field defined and continuous across boundaries
// without returning a wildly wrong frame.
class GateField {
public:
  GateField(const G4VSolid *solid, std::vector<G4ThreeVector> translations,
            std::vector<G4RotationMatrix> rotations);

protected:
  // Return the local coordinates of worldPoint in the placement that contains
  // it, and set outTransform to point to that placement's transform.  Falls
  // back to the closest-centre placement if no copy reports kInside.
  G4ThreeVector
  findContainingPlacement(const G4ThreeVector &worldPoint,
                          const G4AffineTransform *&outTransform) const;

  // Rotate a field vector from local to world using the given transform.
  static G4ThreeVector rotateToWorld(const G4ThreeVector &localField,
                                     const G4AffineTransform &transform);

  const G4VSolid *m_solid; // not owned (Geant4 owns solids)
  std::vector<G4AffineTransform> m_transforms;
};

#endif // GateField_h
