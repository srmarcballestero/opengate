/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateField.h"

#include "G4VSolid.hh"

#include <stdexcept>

GateField::GateField(const G4VSolid *solid,
                     std::vector<G4ThreeVector> translations,
                     std::vector<G4RotationMatrix> rotations)
    : m_solid(solid) {
  if (solid == nullptr) {
    throw std::invalid_argument("GateField: solid must not be null");
  }
  if (translations.size() != rotations.size() || translations.empty()) {
    throw std::invalid_argument("GateField: translations and rotations must be "
                                "non-empty and have the same size");
  }
  m_transforms.reserve(translations.size());
  for (std::size_t i = 0; i < translations.size(); ++i) {
    m_transforms.emplace_back(rotations[i], translations[i]);
  }
}

G4ThreeVector GateField::findContainingPlacement(
    const G4ThreeVector &worldPoint,
    const G4AffineTransform *&outTransform) const {
  // First pass: pick the placement whose solid actually contains worldPoint.
  // Track the closest centre as a fallback for queries that fall slightly
  // outside every copy (chord-finder sub-steps near boundaries).
  std::size_t fallbackIdx = 0;
  G4ThreeVector fallbackLocal =
      m_transforms[0].InverseTransformPoint(worldPoint);
  double minDist2 = fallbackLocal.mag2();

  for (std::size_t i = 0; i < m_transforms.size(); ++i) {
    const G4ThreeVector lp =
        (i == 0) ? fallbackLocal
                 : m_transforms[i].InverseTransformPoint(worldPoint);

    if (m_solid->Inside(lp) != kOutside) {
      outTransform = &m_transforms[i];
      return lp;
    }

    if (i > 0) {
      const double d2 = lp.mag2();
      if (d2 < minDist2) {
        minDist2 = d2;
        fallbackIdx = i;
        fallbackLocal = lp;
      }
    }
  }

  outTransform = &m_transforms[fallbackIdx];
  return fallbackLocal;
}

G4ThreeVector GateField::rotateToWorld(const G4ThreeVector &localField,
                                       const G4AffineTransform &transform) {
  return transform.TransformAxis(localField);
}
