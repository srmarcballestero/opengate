/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GateMagneticField.h"
#include "G4MagneticField.hh"

void init_GateMagneticField(py::module &m) {

  py::class_<GateMagneticField, G4MagneticField,
             std::unique_ptr<GateMagneticField, py::nodelete>>(
      m, "GateMagneticField")

      .def(
        py::init(
          [](
            G4MagneticField*              inner,
            std::vector<G4ThreeVector>    translations,
            std::vector<G4RotationMatrix> rotations
          ) {
            return new GateMagneticField(inner, translations, rotations);
          }),
          py::arg("inner_field"),
          py::arg("translations"),
          py::arg("rotations")
      );
}
