/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateMagneticField.h"

GateMagneticField::GateMagneticField(
    G4MagneticField*              inner,
    std::vector<G4ThreeVector>    translations,
    std::vector<G4RotationMatrix> rotations
) : G4MagneticField(),
    GateField(std::move(translations), std::move(rotations)),
    m_inner(inner)
{}

void GateMagneticField::GetFieldValue(
    const G4double Point[4],
    G4double*      Bfield
) const {
  const G4ThreeVector worldPoint(Point[0], Point[1], Point[2]);
  const G4AffineTransform* transform = nullptr;
  const G4ThreeVector localPoint = closestPlacement(worldPoint, transform);

  const G4double localPos[4] = {
    localPoint.x(), localPoint.y(), localPoint.z(), Point[3]
  };
  G4double localB[3] = {0.0, 0.0, 0.0};
  m_inner->GetFieldValue(localPos, localB);

  const G4ThreeVector worldB =
    rotateToWorld({localB[0], localB[1], localB[2]}, *transform);
  Bfield[0] = worldB.x();
  Bfield[1] = worldB.y();
  Bfield[2] = worldB.z();
}