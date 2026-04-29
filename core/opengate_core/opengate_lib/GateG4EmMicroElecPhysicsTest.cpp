/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateG4EmMicroElecPhysicsTest.h"

#include "G4EmParameters.hh"
#include "G4EmUtility.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PhysListUtil.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"

// particles
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

// standard models used to deactivate global physics in MicroElec regions
#include "G4BetheBlochModel.hh"
#include "G4DummyModel.hh"
#include "G4EmProcessSubType.hh"
#include "G4IonFluctuations.hh"
#include "G4LowECapture.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UrbanMscModel.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VMultipleScattering.hh"

// MicroElec processes and models
#include "G4MicroElecElastic.hh"
#include "G4MicroElecElasticModel_new.hh"
#include "G4MicroElecInelastic.hh"
#include "G4MicroElecInelasticModel_new.hh"
#include "G4MicroElecLOPhononModel.hh"
#include "G4MicroElecLOPhononScattering.hh"
#include "G4MicroElecSurface.hh"

#include <set>

// Energy limits matching the validity range of the _new cross-section data
static const G4double eMicroElecMin   =   0.1 * CLHEP::eV;
static const G4double eMicroElecElMax = 100.0 * CLHEP::MeV;  // elastic
static const G4double eMicroElecInMax =  10.0 * CLHEP::MeV;  // inelastic / phonon
static const G4double pMicroElecMin   = 100.0 * CLHEP::eV;
static const G4double pMicroElecMax   =  10.0 * CLHEP::MeV;
static const G4double ionMicroElecMax =  10.0 * CLHEP::MeV;
static const G4double eCaptureThresh  =  16.7 * CLHEP::eV;

// ---------------------------------------------------------------------------
// Anonymous-namespace helpers
// ---------------------------------------------------------------------------
namespace {

// ---- Material filter -------------------------------------------------
//
// G4MicroElecElasticModel_new::Initialise() and G4MicroElecInelasticModel_new::
// Initialise() iterate the GLOBAL G4ProductionCutsTable and call
// G4MicroElecMaterialStructure for every material in the geometry.  That
// constructor throws FatalException for materials without data files (e.g.
// G4_AIR).  Both loops already contain:
//   if (material->GetName() == "Vacuum") continue;
//
// Fix: temporarily replace unsupported material-cuts couples with a "Vacuum"
// sentinel before calling parent::Initialise(), then restore.
// G4MaterialCutsCouple::SetMaterial() is a documented public API used by
// G4RunManager during geometry updates, so const_cast is intentional and safe
// during single-threaded physics initialization.

// Material-name stems (G4_ prefix stripped) that have MicroElec data files.
// Source: G4EMLOW microelec/Structure/Data_*.dat
static const std::set<G4String> kSupportedStems = {
    "Ag", "Al", "ALUMINUM_OXIDE", "Au", "Be", "BORON_NITRIDE",
    "C",  "Cu", "Fe",  "Ge",  "KAPTON", "Ni", "Si",
    "SILICON_DIOXIDE", "Ti", "TITANIUM_NITRIDE", "W",
};

static G4Material* VacuumSentinel() {
    G4Material* vac = G4Material::GetMaterial("Vacuum", /*warn=*/false);
    if (!vac) {
        G4Element* H = G4NistManager::Instance()->FindOrBuildElement("H");
        vac = new G4Material("Vacuum", 1e-25 * CLHEP::g / CLHEP::cm3, 1, kStateGas);
        vac->AddElement(H, 1);
    }
    return vac;
}

static bool IsSupportedMaterial(const G4Material* mat) {
    const G4String& name = mat->GetName();
    if (name.size() <= 3 || name[0] != 'G' || name[1] != '4' || name[2] != '_')
        return false;
    return kSupportedStems.count(name.substr(3)) > 0;
}

// RAII guard: swaps unsupported couples to the Vacuum sentinel for the
// duration of the _new model Initialise() call, then restores them.
struct MaterialFilterGuard {
    std::vector<std::pair<G4MaterialCutsCouple*, const G4Material*>> saved;

    MaterialFilterGuard() {
        G4Material* sentinel = VacuumSentinel();
        auto* table = G4ProductionCutsTable::GetProductionCutsTable();
        G4int n = (G4int)table->GetTableSize();
        for (G4int i = 0; i < n; ++i) {
            auto* couple = const_cast<G4MaterialCutsCouple*>(
                table->GetMaterialCutsCouple(i));
            const G4Material* mat = couple->GetMaterial();
            if (mat == sentinel) continue;
            if (!IsSupportedMaterial(mat)) {
                saved.emplace_back(couple, mat);
                couple->SetMaterial(sentinel);
            }
        }
    }

    ~MaterialFilterGuard() {
        for (auto& [couple, mat] : saved)
            couple->SetMaterial(mat);
    }
};

// Filtered wrappers: override Initialise() to hide unsupported material-cuts
// couples, delegate to parent (which skips "Vacuum"), then restore.
class FilteredElasticModel_new : public G4MicroElecElasticModel_new {
public:
    using G4MicroElecElasticModel_new::G4MicroElecElasticModel_new;
    void Initialise(const G4ParticleDefinition* p, const G4DataVector& v) override {
        MaterialFilterGuard guard;
        G4MicroElecElasticModel_new::Initialise(p, v);
    }
};

class FilteredInelasticModel_new : public G4MicroElecInelasticModel_new {
public:
    using G4MicroElecInelasticModel_new::G4MicroElecInelasticModel_new;
    void Initialise(const G4ParticleDefinition* p, const G4DataVector& v) override {
        MaterialFilterGuard guard;
        G4MicroElecInelasticModel_new::Initialise(p, v);
    }
};

// ---- Phonon process subtype fix --------------------------------------
//
// G4MicroElecLOPhononScattering sets subtype fLowEnergyElastic (51), which
// collides with G4MicroElecElastic.  G4PhysicsListHelper::RegisterProcess()
// rejects the second registration with the same subtype.
//
// Fix: subclass and reassign to fLowEnergyVibrationalExcitation (54), which
// is the correct semantic category for LO-phonon scattering.
class GateMicroElecLOPhononProcess : public G4MicroElecLOPhononScattering {
public:
    explicit GateMicroElecLOPhononProcess(const G4String& name = "e-_G4MicroElecLOPhonon")
        : G4MicroElecLOPhononScattering(name)
    {
        SetProcessSubType(fLowEnergyVibrationalExcitation);
    }
};

} // anonymous namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GateG4EmMicroElecPhysicsTest::GateG4EmMicroElecPhysicsTest(G4int ver)
    : G4VPhysicsConstructor("GateG4EmMicroElecPhysicsTest", ver) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::AddRegion(const G4String& regionName)
{
  fRegions.push_back(regionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Proton::ProtonDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GateG4EmMicroElecPhysicsTest::ConstructProcess()
{
  if (fRegions.empty()) return;

  auto* elec = G4Electron::Electron();
  auto* prot = G4Proton::Proton();
  auto* gion = G4GenericIon::GenericIon();

  auto* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // ------------------------------------------------------------------
  // Build processes with DummyModel as global default (inert outside
  // MicroElec regions); per-region real models are injected below.
  // G4PhysicsListHelper handles ordering and deduplication by subtype.
  // G4MicroElecSurface and G4LowECapture are non-EM; register directly.
  // ------------------------------------------------------------------

  auto* eElastic = new G4MicroElecElastic("e-_G4MicroElecElastic");
  eElastic->SetEmModel(new G4DummyModel());
  ph->RegisterProcess(eElastic, elec);

  auto* eInelastic = new G4MicroElecInelastic("e-_G4MicroElecInelastic");
  eInelastic->SetEmModel(new G4DummyModel());
  ph->RegisterProcess(eInelastic, elec);

  // GateMicroElecLOPhononProcess uses fLowEnergyVibrationalExcitation so the
  // helper does not reject it as a duplicate of the elastic process.
  auto* ePhonon = new GateMicroElecLOPhononProcess();
  ePhonon->SetEmModel(new G4DummyModel());
  ph->RegisterProcess(ePhonon, elec);

  // G4MicroElecSurface: plain G4VDiscreteProcess — register directly.
  auto* eProcMgr = elec->GetProcessManager();
  auto* eSurface = new G4MicroElecSurface("e-_G4MicroElecSurface");
  eSurface->SetProcessManager(eProcMgr);
  eProcMgr->AddDiscreteProcess(eSurface);

  // G4LowECapture: plain G4VDiscreteProcess — register directly.
  eProcMgr->AddDiscreteProcess(new G4LowECapture(eCaptureThresh));

  auto* pInelastic = new G4MicroElecInelastic("p_G4MicroElecInelastic");
  pInelastic->SetEmModel(new G4DummyModel());
  ph->RegisterProcess(pInelastic, prot);

  auto* iInelastic = new G4MicroElecInelastic("ion_G4MicroElecInelastic");
  iInelastic->SetEmModel(new G4DummyModel());
  ph->RegisterProcess(iInelastic, gion);

  // ------------------------------------------------------------------
  // Per-region configuration:
  //   1) Push standard models above the MicroElec range via
  //      SetActivationLowEnergyLimit + AddEmModel(-2, ..., region)
  //   2) Inject _new MicroElec models at index -1 (highest priority)
  //      using FilteredElasticModel_new / FilteredInelasticModel_new,
  //      which hide unsupported geometry materials from Initialise()
  //      so no FatalException is thrown for e.g. G4_AIR
  // ------------------------------------------------------------------

  const G4double emax = G4EmParameters::Instance()->MaxKinEnergy();

  for (const auto& regionName : fRegions) {
    const G4Region* reg = G4EmUtility::FindRegion(regionName);
    if (nullptr == reg) continue;

    // ---- electron ------------------------------------------------

    {
      auto* p = G4PhysListUtil::FindProcess(elec, fMultipleScattering);
      auto* msc = dynamic_cast<G4VMultipleScattering*>(p);
      if (msc) {
        auto* mod = new G4UrbanMscModel();
        mod->SetActivationLowEnergyLimit(eMicroElecElMax);
        mod->SetHighEnergyLimit(emax);
        msc->AddEmModel(-2, mod, reg);
      }
    }

    {
      auto* p = G4PhysListUtil::FindProcess(elec, fIonisation);
      auto* ioni = dynamic_cast<G4VEnergyLossProcess*>(p);
      if (ioni) {
        auto* mod = new G4MollerBhabhaModel();
        mod->SetActivationLowEnergyLimit(eMicroElecInMax);
        mod->SetHighEnergyLimit(emax);
        ioni->AddEmModel(-2, mod, new G4UniversalFluctuation(), reg);
      }
    }

    {
      auto* mod = new FilteredElasticModel_new();
      mod->SetLowEnergyLimit(eMicroElecMin);
      mod->SetHighEnergyLimit(eMicroElecElMax);
      eElastic->AddEmModel(-1, mod, reg);
    }

    {
      auto* mod = new FilteredInelasticModel_new();
      mod->SetLowEnergyLimit(eMicroElecMin);
      mod->SetHighEnergyLimit(eMicroElecInMax);
      eInelastic->AddEmModel(-1, mod, reg);
    }

    {
      auto* mod = new G4MicroElecLOPhononModel();
      mod->SetLowEnergyLimit(eMicroElecMin);
      mod->SetHighEnergyLimit(eMicroElecInMax);
      ePhonon->AddEmModel(-1, mod, reg);
    }

    // ---- proton --------------------------------------------------

    {
      auto* p = G4PhysListUtil::FindProcess(prot, fIonisation);
      auto* ioni = dynamic_cast<G4VEnergyLossProcess*>(p);
      if (ioni) {
        auto* mod = new G4BetheBlochModel();
        mod->SetActivationLowEnergyLimit(pMicroElecMax);
        mod->SetHighEnergyLimit(emax);
        ioni->AddEmModel(-2, mod, new G4IonFluctuations(), reg);
      }
    }

    {
      auto* mod = new FilteredInelasticModel_new();
      mod->SetLowEnergyLimit(pMicroElecMin);
      mod->SetHighEnergyLimit(pMicroElecMax);
      pInelastic->AddEmModel(-1, mod, reg);
    }

    // ---- GenericIon ----------------------------------------------

    {
      auto* p = G4PhysListUtil::FindProcess(gion, fIonisation);
      auto* ioni = dynamic_cast<G4VEnergyLossProcess*>(p);
      if (ioni) {
        auto* mod = new G4BetheBlochModel();
        mod->SetActivationLowEnergyLimit(ionMicroElecMax);
        mod->SetHighEnergyLimit(emax);
        ioni->AddEmModel(-2, mod, new G4IonFluctuations(), reg);
      }
    }

    {
      auto* mod = new FilteredInelasticModel_new();
      mod->SetLowEnergyLimit(0.0);
      mod->SetHighEnergyLimit(ionMicroElecMax);
      iInelastic->AddEmModel(-1, mod, reg);
    }
  }
}
