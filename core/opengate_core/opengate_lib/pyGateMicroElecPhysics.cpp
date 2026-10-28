/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateMicroElecPhysics.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_GateMicroElecPhysics(py::module &m) {
  // Owned and deleted by the G4VModularPhysicsList it is registered in.
  py::class_<GateMicroElecPhysics, G4VPhysicsConstructor,
             std::unique_ptr<GateMicroElecPhysics, py::nodelete>>(
      m, "GateMicroElecPhysics")
      .def(py::init<G4int>())
      .def("AddRegion", &GateMicroElecPhysics::AddRegion)
      .def("SetRegionElectronThreshold",
           &GateMicroElecPhysics::SetRegionElectronThreshold)
      .def("SetRegionProtonThreshold",
           &GateMicroElecPhysics::SetRegionProtonThreshold)
      .def("SetRegionBaseList", &GateMicroElecPhysics::SetRegionBaseList)
      .def("SetRegionLowEnergyIonisation",
           &GateMicroElecPhysics::SetRegionLowEnergyIonisation)
      .def("SetRegionActivateProtonModels",
           &GateMicroElecPhysics::SetRegionActivateProtonModels)
      .def("SetRegionActivateAlphaModels",
           &GateMicroElecPhysics::SetRegionActivateAlphaModels)
      .def("SetRegionActivateIonModels",
           &GateMicroElecPhysics::SetRegionActivateIonModels);
}
