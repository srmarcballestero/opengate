/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

// MicroElec track-structure physics for user-selected regions, designed to
// overlay a standard EM base list (em_opt3 or em_opt4). Above the per-species
// thresholds the region reproduces the base list's models; below the thresholds
// the MicroElec models take over.
//
// The overlaps with the base list are managed automatically. Each MicroElec
// model is capped at its data-validity ceiling, and the base-list model just
// above it is activated at exactly that energy, so there is neither a gap nor a
// double-count. MicroElec data-validity ceilings:
//   e- elastic   : 500 keV // FIXME: this is not yet taken into account, 10 keV
//   is used globally. e- inelastic : 10 keV p/alpha/ion inelastic : 100 eV - 10
//   MeV
//
// Base-list model crossovers (Geant4 11.4.0):
//   opt4:
//     electron msc:          GoudsmitSaunderson <= 100 MeV, WentzelVI > 100 MeV
//     electron ionisation:   Penelope <= 100 keV, MollerBhabha > 100 keV
//     proton ionisation:     Bragg <= 2 MeV, BetheBloch > 2 MeV
//     alpha ionisation:      Bragg <= 7.9452 MeV, BetheBloch > 7.9452 MeV
//     ion ionisation:        LindhardSorensen
//
//   opt3:
//     electron msc:          Urban
//     electron ionisation:   MollerBhabha
//     proton ionisation:     Bragg <= 2 MeV, BetheBloch > 2 MeV
//     alpha ionisation:      Bragg <= 7.9452 MeV, BetheBloch > 7.9452 MeV
//     ion ionisation:        LindhardSorensen

#include "GateG4EmMicroElecPhysicsTest.h"

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
// MicroElec data-validity ceilings.
constexpr G4double kElectronMin = 0.1 * CLHEP::eV;
constexpr G4double kElectronElasticMax = 500.0 * CLHEP::keV;
constexpr G4double kElectronInelasticMax = 10.0 * CLHEP::keV;
constexpr G4double kElectronPhononMax = 10.0 * CLHEP::MeV;
constexpr G4double kHadronMin = 100.0 * CLHEP::eV;
constexpr G4double kHadronInelasticMax = 10.0 * CLHEP::MeV;
// Base-list model crossovers.
constexpr G4double kMscWentzelCrossover =
    100.0 * CLHEP::MeV;                                     // opt4 GS/WentzelVI
constexpr G4double kPenelopeHighLimit = 100.0 * CLHEP::keV; // opt4 Penelope/MB
constexpr G4double kBraggHighLimit = 2.0 * CLHEP::MeV;      // proton Bragg/BB
constexpr G4double kAlphaCrossover = 7.9452 * CLHEP::MeV;   // alpha ionIoni
constexpr G4double kIonCrossover = 10.0 * CLHEP::MeV;       // ion ionIoni
constexpr G4double kTopEnergy = 10.0 * CLHEP::TeV;
} // namespace

// G4MicroElecLOPhononScattering sets subtype fLowEnergyElastic (51), which
// collides with G4MicroElecElastic.  Reassign to
// fLowEnergyVibrationalExcitation (54).
class GateMicroElecLOPhononProcess : public G4MicroElecLOPhononScattering {
public:
  explicit GateMicroElecLOPhononProcess(
      const G4String &name = "e-_G4MicroElecLOPhonon")
      : G4MicroElecLOPhononScattering(name) {
    SetProcessSubType(fLowEnergyVibrationalExcitation);
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GateG4EmMicroElecPhysicsTest::GateG4EmMicroElecPhysicsTest(G4int ver)
    : G4VPhysicsConstructor("GateG4EmMicroElecPhysicsTest", ver) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GateG4EmMicroElecPhysicsTest::RegionConfig &
GateG4EmMicroElecPhysicsTest::configFor(const G4String &regionName) {
  return fConfigs[regionName]; // default-constructs on first access
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::AddRegion(const G4String &regionName) {
  fRegions.push_back(regionName);
  configFor(regionName); // ensure a config exists
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::SetRegionElectronThreshold(
    const G4String &regionName, G4double value) {
  configFor(regionName).electronThreshold = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::SetRegionProtonThreshold(
    const G4String &regionName, G4double value) {
  configFor(regionName).protonThreshold = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::SetRegionBaseList(const G4String &regionName,
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
    G4Exception("GateG4EmMicroElecPhysicsTest::SetRegionBaseList", "gate0002",
                FatalException, msg.c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::SetRegionUsePenelope(
    const G4String &regionName, G4bool usePenelope) {
  configFor(regionName).usePenelope = usePenelope;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::ConstructParticle() {
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Proton::ProtonDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::ConstructProcess() {
  if (fRegions.empty())
    return;

  // ---- EM parameters: allow MicroElec to track down to sub-keV energies ----
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

  auto *elec = G4Electron::Electron();
  auto *prot = G4Proton::Proton();
  auto *alph = G4Alpha::Alpha();
  auto *gion = G4GenericIon::GenericIon();

  auto *eProcMgr = elec->GetProcessManager();
  auto *pProcMgr = prot->GetProcessManager();
  auto *alphaProcMgr = alph->GetProcessManager();
  auto *ionProcMgr = gion->GetProcessManager();

  // Register MicroElec processes with dummy models to ensure they are in the
  // process manager. The actual models are injected per-region below.
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

  auto *pInelastic = new G4MicroElecInelastic("p_G4MicroElecInelastic");
  pInelastic->SetEmModel(new G4DummyModel());
  pProcMgr->AddDiscreteProcess(pInelastic);

  auto *alphaInelastic = new G4MicroElecInelastic("alpha_G4MicroElecInelastic");
  alphaInelastic->SetEmModel(new G4DummyModel());
  alphaProcMgr->AddDiscreteProcess(alphaInelastic);

  auto *iInelastic = new G4MicroElecInelastic("ion_G4MicroElecInelastic");
  iInelastic->SetEmModel(new G4DummyModel());
  ionProcMgr->AddDiscreteProcess(iInelastic);

  // Per-region model injection via G4EmConfigurator.
  G4EmConfigurator *em_config =
      G4LossTableManager::Instance()->EmConfigurator();

  for (const auto &regionName : fRegions) {
    const RegionConfig &cfg = fConfigs.at(regionName);

    const G4double eElasticMax =
        std::min(cfg.electronThreshold, kElectronElasticMax);
    const G4double eInelasticMax =
        std::min(cfg.electronThreshold, kElectronInelasticMax);
    const G4double pMax = std::min(cfg.protonThreshold, kHadronInelasticMax);

    // ================================================================
    // ELECTRONS
    // ================================================================

    // ---- multiple scattering (reproduces the base list above eElasticMax)
    if (cfg.baseList == GateMicroElecBaseList::opt3) {
      auto *mscUrban = new G4UrbanMscModel();
      mscUrban->SetActivationLowEnergyLimit(eElasticMax);
      em_config->SetExtraEmModel("e-", "msc", mscUrban, regionName, eElasticMax,
                                 kTopEnergy);
    } else {
      auto *mscGS = new G4GoudsmitSaundersonMscModel();
      mscGS->SetActivationLowEnergyLimit(eElasticMax);
      em_config->SetExtraEmModel("e-", "msc", mscGS, regionName, eElasticMax,
                                 kMscWentzelCrossover);

      auto *mscWentzel = new G4WentzelVIModel();
      mscWentzel->SetActivationLowEnergyLimit(kMscWentzelCrossover);
      em_config->SetExtraEmModel("e-", "msc", mscWentzel, regionName,
                                 kMscWentzelCrossover, kTopEnergy);
    }

    // ---- ionisation (reproduces the base list above eInelasticMax) ----
    if (cfg.baseList == GateMicroElecBaseList::opt4 && cfg.usePenelope &&
        eInelasticMax < kPenelopeHighLimit) {
      auto *pen = new G4PenelopeIonisationModel();
      pen->SetActivationLowEnergyLimit(eInelasticMax);
      em_config->SetExtraEmModel("e-", "eIoni", pen, regionName, eInelasticMax,
                                 kPenelopeHighLimit,
                                 G4EmStandUtil::ModelOfFluctuations());

      auto *mb = new G4MollerBhabhaModel();
      mb->SetActivationLowEnergyLimit(kPenelopeHighLimit);
      em_config->SetExtraEmModel("e-", "eIoni", mb, regionName,
                                 kPenelopeHighLimit, kTopEnergy,
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
                               new G4MicroElecInelasticModel_new(), regionName,
                               kHadronMin, pMax);

    // ================================================================
    // ALPHA: BetheBloch above 7.9452 MeV; MicroElec inelastic below.
    // ================================================================
    auto *alphaBB = new G4BetheBlochModel();
    alphaBB->SetActivationLowEnergyLimit(kAlphaCrossover);
    em_config->SetExtraEmModel("alpha", "ionIoni", alphaBB, regionName,
                               kAlphaCrossover, kTopEnergy,
                               G4EmStandUtil::ModelOfFluctuations(true));

    em_config->SetExtraEmModel("alpha", "alpha_G4MicroElecInelastic",
                               new G4MicroElecInelasticModel_new(), regionName,
                               kHadronMin, kAlphaCrossover);

    // ================================================================
    // GENERIC IONS: LindhardSorensen above 10 MeV; MicroElec inelastic below.
    // ================================================================
    auto *ionLS = new G4LindhardSorensenIonModel();
    ionLS->SetActivationLowEnergyLimit(kIonCrossover);
    em_config->SetExtraEmModel("GenericIon", "ionIoni", ionLS, regionName,
                               kIonCrossover, kTopEnergy,
                               G4EmStandUtil::ModelOfFluctuations(true));

    em_config->SetExtraEmModel("GenericIon", "ion_G4MicroElecInelastic",
                               new G4MicroElecInelasticModel_new(), regionName,
                               kHadronMin, kIonCrossover);
  }
}
