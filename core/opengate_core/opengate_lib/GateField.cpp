/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateField.h"

#include <cmath>
#include <stdexcept>

GateField::GateField(
    std::vector<G4ThreeVector>    translations,
    std::vector<G4RotationMatrix> rotations
)
{
  if (translations.size() != rotations.size() || translations.empty()) {
    throw std::invalid_argument(
      "GateField: translations and rotations must be non-empty and have the same size"
    );
  }
  m_transforms.reserve(translations.size());
  for (std::size_t i = 0; i < translations.size(); ++i) {
    m_transforms.emplace_back(rotations[i], translations[i]);
  }
}

G4ThreeVector GateField::closestPlacement(
    const G4ThreeVector&      worldPoint,
    const G4AffineTransform*& outTransform
) const {
  // For a single placement this is trivially correct.
  // For repeated placements, the closest centre is the right copy for any
  // non-overlapping geometry (a Geant4 validity requirement).
  G4ThreeVector bestLocal = m_transforms[0].InverseTransformPoint(worldPoint);
  outTransform = &m_transforms[0];
  double minDist2 = bestLocal.mag2();

  for (std::size_t i = 1; i < m_transforms.size(); ++i) {
    G4ThreeVector lp = m_transforms[i].InverseTransformPoint(worldPoint);
    double d2 = lp.mag2();
    if (d2 < minDist2) {
      minDist2 = d2;
      bestLocal = lp;
      outTransform = &m_transforms[i];
    }
  }
  return bestLocal;
}

G4ThreeVector GateField::rotateToWorld(
    const G4ThreeVector&     localField,
    const G4AffineTransform& transform
) {
  return transform.TransformAxis(localField);
}
