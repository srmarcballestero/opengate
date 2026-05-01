/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "G4MagneticField.hh"
#include "G4VSolid.hh"
#include "GateMagneticField.h"

void init_GateMagneticField(py::module &m) {

  py::class_<GateMagneticField, G4MagneticField,
             std::unique_ptr<GateMagneticField, py::nodelete>>(
      m, "GateMagneticField")

      .def(py::init([](G4MagneticField *inner, const G4VSolid *solid,
                       std::vector<G4ThreeVector> translations,
                       std::vector<G4RotationMatrix> rotations) {
             return new GateMagneticField(inner, solid, translations,
                                          rotations);
           }),
           py::arg("inner_field"), py::arg("solid"), py::arg("translations"),
           py::arg("rotations"));
}
