/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

// MicroElec track-structure physics in selected regions. Below the per-species
// thresholds the MicroElec models are used; above them, the region uses the
// models of the base list (em_opt3 or em_opt4).
//
// Only the MicroElec processes and the region model overrides are added; the
// simulation physics list must provide the matching base list.
//
// Each MicroElec model ends at min(threshold, data limit) and the base-list
// model above it starts at that same energy. MicroElec data limits:
//   e- elastic           : 500 keV (FIXME: electron threshold capped at 10 keV)
//   e- inelastic         : 10 keV
//   p/alpha/ion inelastic: 100 eV - 10 MeV
//
// Base-list model crossovers (Geant4 11.4.0):
//   opt4:
//     electron msc:          GoudsmitSaunderson <= 100 MeV, WentzelVI > 100 MeV
//     electron ionisation:   Penelope <= 100 keV, MollerBhabha > 100 keV
//   opt3:
//     electron msc:          Urban
//     electron ionisation:   MollerBhabha
//   all:
//     proton ionisation:     Bragg <= 2 MeV, BetheBloch > 2 MeV
//     alpha ionisation:      Bragg <= 7.9452 MeV, BetheBloch > 7.9452 MeV
//     ion ionisation:        LindhardSorensen

#include "GateG4MicroElecPhysics.h"

#include "G4EmConfigurator.hh"
#include "G4EmParameters.hh"
#include "G4EmStandUtil.hh"
#include "G4LossTableManager.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>

// particles
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

// standard models
#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4DummyModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4LindhardSorensenIonModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"

// MicroElec processes and models
#include "G4MicroElecElastic.hh"
#include "G4MicroElecElasticModel_new.hh"
#include "G4MicroElecInelastic.hh"
#include "G4MicroElecInelasticModel_new.hh"
#include "G4MicroElecLOPhononModel.hh"
#include "G4MicroElecLOPhononScattering.hh"
#include "G4MicroElecSurface.hh"

namespace {
// MicroElec model limits.
constexpr G4double kElectronMin = 0.1 * CLHEP::eV;
constexpr G4double kElectronElasticMax = 500.0 * CLHEP::keV;
constexpr G4double kElectronInelasticMax = 10.0 * CLHEP::keV;
constexpr G4double kElectronPhononMax = 10.0 * CLHEP::MeV;
constexpr G4double kHadronMin = 100.0 * CLHEP::eV;
constexpr G4double kHadronInelasticMax = 10.0 * CLHEP::MeV;
// Base-list model crossovers.
constexpr G4double kMscWentzelCrossover =
    100.0 * CLHEP::MeV;                                     // opt4 GS/WentzelVI
constexpr G4double kLowEIoniHighLimit = 100.0 * CLHEP::keV; // opt4 low-E/MB
constexpr G4double kBraggHighLimit = 2.0 * CLHEP::MeV;      // proton Bragg/BB
constexpr G4double kAlphaCrossover = 7.9452 * CLHEP::MeV;   // alpha ionIoni
constexpr G4double kIonCrossover = 10.0 * CLHEP::MeV;       // ion ionIoni
constexpr G4double kTopEnergy = 10.0 * CLHEP::TeV;

const char *ToString(GateMicroElecBaseList bl) {
  switch (bl) {
  case GateMicroElecBaseList::opt3:
    return "opt3";
  case GateMicroElecBaseList::opt4:
    return "opt4";
  }
  return "?";
}

const char *ToString(GateMicroElecLowEnergyIonisation m) {
  switch (m) {
  case GateMicroElecLowEnergyIonisation::automatic:
    return "automatic";
  case GateMicroElecLowEnergyIonisation::penelope:
    return "penelope";
  case GateMicroElecLowEnergyIonisation::none:
    return "none";
  }
  return "?";
}
} // namespace

// G4MicroElecLOPhononScattering uses subtype fLowEnergyElastic (51), which is
// already taken by G4MicroElecElastic.
class GateMicroElecLOPhononProcess : public G4MicroElecLOPhononScattering {
public:
  explicit GateMicroElecLOPhononProcess(
      const G4String &name = "e-_G4MicroElecLOPhonon")
      : G4MicroElecLOPhononScattering(name) {
    SetProcessSubType(fLowEnergyVibrationalExcitation);
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GateG4MicroElecPhysics::GateG4MicroElecPhysics(G4int ver)
    : G4VPhysicsConstructor("GateG4MicroElecPhysics", ver) {
  // Set at construction so that G4EmParameters set later (e.g. by the GATE
  // physics manager) take precedence.
  G4EmParameters *param = G4EmParameters::Instance();
  param->SetBuildCSDARange(true);
  param->SetMinEnergy(0.1 * CLHEP::eV);
  param->SetMaxEnergy(kTopEnergy);
  param->SetLowestElectronEnergy(0.0);
  param->SetNumberOfBinsPerDecade(20);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetFluo(true);
  param->SetAuger(true);
  param->SetPixe(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GateG4MicroElecPhysics::RegionConfig &
GateG4MicroElecPhysics::configFor(const G4String &regionName) {
  return fConfigs[regionName]; // default-constructs on first access
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::AddRegion(const G4String &regionName) {
  if (std::find(fRegions.begin(), fRegions.end(), regionName) ==
      fRegions.end()) {
    fRegions.push_back(regionName);
  }
  configFor(regionName); // ensure a config exists
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::SetRegionElectronThreshold(
    const G4String &regionName, G4double value) {
  configFor(regionName).electronThreshold = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::SetRegionProtonThreshold(
    const G4String &regionName, G4double value) {
  configFor(regionName).protonThreshold = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::SetRegionBaseList(const G4String &regionName,
                                               const G4String &baseList) {
  G4String bl = baseList;
  G4StrUtil::to_lower(bl);
  if (bl == "opt3") {
    configFor(regionName).baseList = GateMicroElecBaseList::opt3;
  } else if (bl == "opt4") {
    configFor(regionName).baseList = GateMicroElecBaseList::opt4;
  } else {
    G4String msg = "Unknown MicroElec base list '" + baseList +
                   "'. Allowed values are 'opt3' and 'opt4'.";
    G4Exception("GateG4MicroElecPhysics::SetRegionBaseList", "gate0002",
                FatalException, msg.c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::SetRegionLowEnergyIonisation(
    const G4String &regionName, const G4String &model) {
  G4String m = model;
  G4StrUtil::to_lower(m);
  auto &cfg = configFor(regionName);
  if (m == "automatic") {
    cfg.lowEnergyIonisation = GateMicroElecLowEnergyIonisation::automatic;
  } else if (m == "penelope") {
    cfg.lowEnergyIonisation = GateMicroElecLowEnergyIonisation::penelope;
  } else if (m == "none") {
    cfg.lowEnergyIonisation = GateMicroElecLowEnergyIonisation::none;
  } else {
    G4String msg = "Unknown MicroElec low-energy ionisation model '" + model +
                   "'. Allowed values are 'penelope', 'none' and 'automatic'.";
    G4Exception("GateG4MicroElecPhysics::SetRegionLowEnergyIonisation",
                "gate0003", FatalException, msg.c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::SetRegionActivateProtonModels(
    const G4String &regionName, G4bool value) {
  configFor(regionName).activateProtonModels = value;
}

void GateG4MicroElecPhysics::SetRegionActivateAlphaModels(
    const G4String &regionName, G4bool value) {
  configFor(regionName).activateAlphaModels = value;
}

void GateG4MicroElecPhysics::SetRegionActivateIonModels(
    const G4String &regionName, G4bool value) {
  configFor(regionName).activateIonModels = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GateMicroElecLowEnergyIonisation
GateG4MicroElecPhysics::ResolveLowEnergyIonisation(
    const G4String &regionName, const RegionConfig &cfg) const {
  using LEI = GateMicroElecLowEnergyIonisation;
  if (cfg.lowEnergyIonisation == LEI::automatic) {
    switch (cfg.baseList) {
    case GateMicroElecBaseList::opt4:
      return LEI::penelope;
    case GateMicroElecBaseList::opt3:
      return LEI::none;
    }
  }
  // em_opt3 uses MollerBhabha over the whole range.
  if (cfg.baseList == GateMicroElecBaseList::opt3 &&
      cfg.lowEnergyIonisation != LEI::none) {
    G4String msg = "MicroElec region '" + regionName +
                   "': low-energy ionisation '" +
                   ToString(cfg.lowEnergyIonisation) +
                   "' is not consistent with base list 'opt3', which uses "
                   "MollerBhabha over the whole range.";
    G4Exception("GateG4MicroElecPhysics::ConstructProcess", "gate0004",
                FatalException, msg.c_str());
  }
  return cfg.lowEnergyIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::ConstructParticle() {
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Proton::ProtonDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4MicroElecPhysics::ConstructProcess() {
  if (fRegions.empty())
    return;

  // Processes are global, models per region: a heavy-particle MicroElec
  // process is registered only if at least one region activates it.
  G4bool anyProton = false;
  G4bool anyAlpha = false;
  G4bool anyIon = false;
  for (const auto &regionName : fRegions) {
    const RegionConfig &cfg = fConfigs.at(regionName);
    anyProton = anyProton || cfg.activateProtonModels;
    anyAlpha = anyAlpha || cfg.activateAlphaModels;
    anyIon = anyIon || cfg.activateIonModels;
  }

  auto *eProcMgr = G4Electron::Electron()->GetProcessManager();

  // MicroElec processes carry a dummy model outside the regions.
  auto *eElastic = new G4MicroElecElastic("e-_G4MicroElecElastic");
  eElastic->SetEmModel(new G4DummyModel());
  eProcMgr->AddDiscreteProcess(eElastic);

  auto *eInelastic = new G4MicroElecInelastic("e-_G4MicroElecInelastic");
  eInelastic->SetEmModel(new G4DummyModel());
  eProcMgr->AddDiscreteProcess(eInelastic);

  auto *ePhonon = new GateMicroElecLOPhononProcess();
  ePhonon->SetEmModel(new G4DummyModel());
  eProcMgr->AddDiscreteProcess(ePhonon);

  auto *eSurface = new G4MicroElecSurface("e-_G4MicroElecSurface");
  eSurface->SetProcessManager(eProcMgr);
  eProcMgr->AddDiscreteProcess(eSurface);

  if (anyProton) {
    auto *pInelastic = new G4MicroElecInelastic("p_G4MicroElecInelastic");
    pInelastic->SetEmModel(new G4DummyModel());
    G4Proton::Proton()->GetProcessManager()->AddDiscreteProcess(pInelastic);
  }

  if (anyAlpha) {
    auto *alphaInelastic =
        new G4MicroElecInelastic("alpha_G4MicroElecInelastic");
    alphaInelastic->SetEmModel(new G4DummyModel());
    G4Alpha::Alpha()->GetProcessManager()->AddDiscreteProcess(alphaInelastic);
  }

  if (anyIon) {
    auto *iInelastic = new G4MicroElecInelastic("ion_G4MicroElecInelastic");
    iInelastic->SetEmModel(new G4DummyModel());
    G4GenericIon::GenericIon()->GetProcessManager()->AddDiscreteProcess(
        iInelastic);
  }

  // Per-region model injection via G4EmConfigurator.
  G4EmConfigurator *em_config =
      G4LossTableManager::Instance()->EmConfigurator();

  for (const auto &regionName : fRegions) {
    const RegionConfig &cfg = fConfigs.at(regionName);
    const GateMicroElecLowEnergyIonisation lowEIoni =
        ResolveLowEnergyIonisation(regionName, cfg);

    const G4double eElasticMax =
        std::min(cfg.electronThreshold, kElectronElasticMax);
    const G4double eInelasticMax =
        std::min(cfg.electronThreshold, kElectronInelasticMax);
    const G4double pMax = std::min(cfg.protonThreshold, kHadronInelasticMax);

    if (verboseLevel > 0) {
      G4cout << "GateG4MicroElecPhysics: region '" << regionName
             << "': base list = " << ToString(cfg.baseList)
             << ", electron threshold = " << cfg.electronThreshold / CLHEP::keV
             << " keV, proton threshold = " << cfg.protonThreshold / CLHEP::MeV
             << " MeV, low-energy ionisation = " << ToString(lowEIoni)
             << ", protons = " << (cfg.activateProtonModels ? "on" : "off")
             << ", alphas = " << (cfg.activateAlphaModels ? "on" : "off")
             << ", ions = " << (cfg.activateIonModels ? "on" : "off") << G4endl;
    }

    // ================================================================
    // ELECTRONS (always active)
    // ================================================================

    // ---- multiple scattering above eElasticMax ----
    if (cfg.baseList == GateMicroElecBaseList::opt3) {
      auto *mscUrban = new G4UrbanMscModel();
      mscUrban->SetActivationLowEnergyLimit(eElasticMax);
      em_config->SetExtraEmModel("e-", "msc", mscUrban, regionName, eElasticMax,
                                 kTopEnergy);
    } else { // opt4
      auto *mscGS = new G4GoudsmitSaundersonMscModel();
      mscGS->SetActivationLowEnergyLimit(eElasticMax);
      em_config->SetExtraEmModel("e-", "msc", mscGS, regionName, eElasticMax,
                                 kMscWentzelCrossover);

      auto *mscWentzel = new G4WentzelVIModel();
      mscWentzel->SetActivationLowEnergyLimit(kMscWentzelCrossover);
      em_config->SetExtraEmModel("e-", "msc", mscWentzel, regionName,
                                 kMscWentzelCrossover, kTopEnergy);
    }

    // ---- ionisation above eInelasticMax ----
    if (lowEIoni != GateMicroElecLowEnergyIonisation::none &&
        eInelasticMax < kLowEIoniHighLimit) {
      auto *pen = new G4PenelopeIonisationModel();
      pen->SetActivationLowEnergyLimit(eInelasticMax);
      em_config->SetExtraEmModel("e-", "eIoni", pen, regionName, eInelasticMax,
                                 kLowEIoniHighLimit,
                                 G4EmStandUtil::ModelOfFluctuations());

      auto *mb = new G4MollerBhabhaModel();
      mb->SetActivationLowEnergyLimit(kLowEIoniHighLimit);
      em_config->SetExtraEmModel("e-", "eIoni", mb, regionName,
                                 kLowEIoniHighLimit, kTopEnergy,
                                 G4EmStandUtil::ModelOfFluctuations());
    } else {
      auto *mb = new G4MollerBhabhaModel();
      mb->SetActivationLowEnergyLimit(eInelasticMax);
      em_config->SetExtraEmModel("e-", "eIoni", mb, regionName, eInelasticMax,
                                 kTopEnergy,
                                 G4EmStandUtil::ModelOfFluctuations());
    }

    // ---- MicroElec electron models (below threshold) ----
    em_config->SetExtraEmModel("e-", "e-_G4MicroElecElastic",
                               new G4MicroElecElasticModel_new(), regionName,
                               kElectronMin, eElasticMax);

    em_config->SetExtraEmModel("e-", "e-_G4MicroElecInelastic",
                               new G4MicroElecInelasticModel_new(), regionName,
                               kElectronMin, eInelasticMax);

    em_config->SetExtraEmModel("e-", "e-_G4MicroElecLOPhonon",
                               new G4MicroElecLOPhononModel(), regionName,
                               kElectronMin, kElectronPhononMax);

    // ================================================================
    // PROTONS: Bragg (<= 2 MeV) + BetheBloch above; MicroElec inelastic below.
    // ================================================================
    if (cfg.activateProtonModels) {
      if (pMax < kBraggHighLimit) {
        auto *bragg = new G4BraggModel();
        bragg->SetActivationLowEnergyLimit(pMax);
        em_config->SetExtraEmModel("proton", "hIoni", bragg, regionName, pMax,
                                   kBraggHighLimit,
                                   G4EmStandUtil::ModelOfFluctuations());

        auto *bb = new G4BetheBlochModel();
        bb->SetActivationLowEnergyLimit(kBraggHighLimit);
        em_config->SetExtraEmModel("proton", "hIoni", bb, regionName,
                                   kBraggHighLimit, kTopEnergy,
                                   G4EmStandUtil::ModelOfFluctuations());
      } else {
        auto *bb = new G4BetheBlochModel();
        bb->SetActivationLowEnergyLimit(pMax);
        em_config->SetExtraEmModel("proton", "hIoni", bb, regionName, pMax,
                                   kTopEnergy,
                                   G4EmStandUtil::ModelOfFluctuations());
      }

      em_config->SetExtraEmModel("proton", "p_G4MicroElecInelastic",
                                 new G4MicroElecInelasticModel_new(),
                                 regionName, kHadronMin, pMax);
    }

    // ================================================================
    // ALPHA: BetheBloch above 7.9452 MeV; MicroElec inelastic below.
    // ================================================================
    if (cfg.activateAlphaModels) {
      auto *alphaBB = new G4BetheBlochModel();
      alphaBB->SetActivationLowEnergyLimit(kAlphaCrossover);
      em_config->SetExtraEmModel("alpha", "ionIoni", alphaBB, regionName,
                                 kAlphaCrossover, kTopEnergy,
                                 G4EmStandUtil::ModelOfFluctuations(true));

      em_config->SetExtraEmModel("alpha", "alpha_G4MicroElecInelastic",
                                 new G4MicroElecInelasticModel_new(),
                                 regionName, kHadronMin, kAlphaCrossover);
    }

    // ================================================================
    // GENERIC IONS: LindhardSorensen above 10 MeV; MicroElec inelastic below.
    // ================================================================
    if (cfg.activateIonModels) {
      auto *ionLS = new G4LindhardSorensenIonModel();
      ionLS->SetActivationLowEnergyLimit(kIonCrossover);
      em_config->SetExtraEmModel("GenericIon", "ionIoni", ionLS, regionName,
                                 kIonCrossover, kTopEnergy,
                                 G4EmStandUtil::ModelOfFluctuations(true));

      em_config->SetExtraEmModel("GenericIon", "ion_G4MicroElecInelastic",
                                 new G4MicroElecInelasticModel_new(),
                                 regionName, kHadronMin, kIonCrossover);
    }
  }
}
