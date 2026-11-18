/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// lQPD_Nb -- G4CMP-TFR (BogoliubovQP physics) version
//
// Port of G4CMP/examples/RISQTutorial/src/RISQTutorialDetectorConstruction_lQPD_Nb.cc
//
// This file REPLACES QuasiparticleDetectorConstruction.cc when the example is
// configured with   cmake -DQUASIPARTICLE_GEOMETRY=lQPD_Nb ...   (see
// CMakeLists.txt).  Both files implement the same class, so only one of them
// is ever compiled.
//
// ============================== WHAT CHANGED ================================
// The geometry (every box/tube, its dimensions and its position on the chip)
// is the same as in the G4CMP version.  What changed is the PHYSICS MODEL of
// the metal films:
//
//  G4CMP / Kaplan version                   TFR version (this file)
//  ---------------------------------------  ---------------------------------
//  Metal volumes are only absorbing         Metal volumes are real superconducting
//  "electrodes": a phonon that reaches a    thin films. Each one gets its own
//  film surface is absorbed by              G4LatticePhysical with gap, Teff, Dn,
//  G4CMPPhononElectrode (KaplanQP) with a   QP trapping time, polycrystalline
//  film thickness/absorption/lifetime       scattering length ...
//  parameter table (AttachPhononSensor_Al/  Phonons are TRANSMITTED into the film
//  _Nb, G4CMPPhononElectrode, filmThickness by the TFR boundary process (or
//  in nm only a model parameter).           absorbed/reflected according to the
//                                           surface pAbsProb/pReflProb), break
//                                           Cooper pairs (scPairBreaking) and the
//                                           resulting BogoliubovQPs diffuse in XY.
//  One surface per metal (AlSurfProp,       One G4CMPSurfaceProperty per interface
//  NbSurfProp) applied Si->metal only.      (Si/Al, Si/Nb, Al/Vac, Nb/Vac, Al/Al,
//                                           Nb/Nb, Al/Nb) applied in BOTH
//                                           directions, as required for phonon
//                                           and QP transmission.
//  filmThicknessAl/Nb = Kaplan parameter.   filmThicknessAl/Nb = REAL thickness of
//                                           the films (see "Thickness" below).
//
// Rules from the quasiparticle README that this file follows:
//  * QPs only diffuse in the world XY plane and keep their z, so every
//    top-side thin film has the SAME thickness and the SAME mid-plane z.
//  * Every superconducting volume has its own G4LatticePhysical.
//  * Every pair of touching volumes needs a G4CMPLogicalBorderSurface in both
//    directions (otherwise QPs/phonons cannot be handled at that boundary).
//
// Macro variables ported from RISQTutorial (read from
// QuasiparticleConfigManager, set with /g4cmp/... before /run/initialize):
//   numsensors              chip half-width = 0.5 cm / sqrt(numsensors)
//   Al                      true: patterned film is Al; false: Nb (the
//                           junction/absorber "sensor" is Al in both cases)
//   filmThicknessAl/Nb      film thickness in nm (<=0: default 200 nm)
//   pAbsProbSideWallSi      phonon absorption prob. on the Si side walls
//   pAbsProbPolishedWallSi  phonon absorption prob. on the polished Si faces
//
// Thickness: all top-side films share one thickness, equal to
// filmThicknessAl if Al=true, otherwise filmThicknessNb (the Al sensor takes
// the same thickness so that QPs can cross from Nb into Al at the same z).
// The Nb backplane is on the other face of the chip, is its own film plane,
// and always uses filmThicknessNb.

#include "QuasiparticleDetectorConstruction.hh"
#include "QuasiparticleDetectorParameters.hh"
#include "QuasiparticleConfigManager.hh"
#include "QuasiparticleSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace QuasiparticleDetectorParameters;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Same member list/ordering as the resonator geometry (QuasiparticleDetector
// Construction.cc) because both files implement the same class.
QuasiparticleDetectorConstruction::QuasiparticleDetectorConstruction()
  : fLiquidHelium(0), fSilicon(0), fAluminum(0), fTungsten(0), fCopper(0),
    fNiobium(0), fWorldPhys(0), topSurfProp(0), topSurfProp2(0), vacSurfProp(0),
    wallSurfProp(0), copTopSurfProp(0), copTopSurfProp2(0), alNbSurfProp(0),
    fConstructed(false), fVacVacInterface(0), fAlVacInterface(0), fAlAlInterface(0),
    fSiAlInterface(0), fSiCuInterface(0), fSiVacInterface(0), fCuVacInterface(0),
    fAir(0), fVacuum(0), fSiNbInterface(0), fNbVacInterface(0), fNbNbInterface(0),
    fAlNbInterface(0), fSiSideWallInterface(0), fSiPolishedWallInterface(0),
    fSuperconductorSensitivity(0) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Only the surface properties that SetupGeometry() creates are deleted here
// (the sensitive detector is owned by the G4SDManager)
QuasiparticleDetectorConstruction::~QuasiparticleDetectorConstruction() {
  delete fSiAlInterface;
  delete fSiNbInterface;
  delete fAlVacInterface;
  delete fNbVacInterface;
  delete fAlAlInterface;
  delete fNbNbInterface;
  delete fAlNbInterface;
  delete fSiSideWallInterface;
  delete fSiPolishedWallInterface;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* QuasiparticleDetectorConstruction::Construct()
{
  if (fConstructed) {
    if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
      // Run manager hasn't cleaned volume stores. This code shouldn't execute
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    // Have to completely remove all lattices to avoid warning on reconstruction
    G4LatticeManager::GetLatticeManager()->Reset();
    // Clear all LogicalSurfaces
    // NOTE: No need to redefine the G4CMPSurfaceProperties
    G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void QuasiparticleDetectorConstruction::DefineMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();

  fAir = nistManager->FindOrBuildMaterial("G4_AIR");

  // TFR: look the vacuum up first so that re-constructing the geometry does
  // not try to create a second material called "VACUUM"
  fVacuum = G4Material::GetMaterial("VACUUM", false);
  if (!fVacuum) {
    fVacuum = new G4Material("VACUUM",
          1.,
		1.008*CLHEP::g/CLHEP::mole,
		1.0e-25*CLHEP::g/CLHEP::cm3,
		kStateGas,
		0.01*CLHEP::kelvin,
	    3.0e-18*pascal);
  }
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void QuasiparticleDetectorConstruction::SetupGeometry()
{
  //
  // TFR: read the macro variables (QuasiparticleConfigManager) ----------------
  //
  const G4bool useAl = QuasiparticleConfigManager::GetAl();

  // numsensors: the G4CMP version divided by sqrt(0) when the macro did not
  // set it (default -1 was mapped to 0). Fall back to 1 sensor instead.
  G4double numsensors = QuasiparticleConfigManager::Getnumsensors();
  if (numsensors <= 0.) {
    G4cout << "QuasiparticleDetectorConstruction (lQPD_Nb): /g4cmp/numsensors "
           << "not set (or <= 0); using numsensors = 1" << G4endl;
    numsensors = 1.;
  }

  // Film thicknesses (macro values are in nm, <= 0 means "not set").
  // Top-side films (pattern + Al sensor): one common thickness.
  // Backplane: always the Nb value.
  G4double filmThicknessNb_nm = QuasiparticleConfigManager::GetfilmThicknessNb();
  G4double filmThicknessTop_nm = useAl ? QuasiparticleConfigManager::GetfilmThicknessAl()
                                       : filmThicknessNb_nm;
  const G4double filmThicknessTop = (filmThicknessTop_nm > 0.) ? filmThicknessTop_nm*nm
                                                               : dp_lQPDDefaultFilmThickness;
  const G4double filmThicknessBack = (filmThicknessNb_nm > 0.) ? filmThicknessNb_nm*nm
                                                               : dp_lQPDDefaultFilmThickness;

  // The G4CMP version used the constant half-thickness  0.1*um  for all metal
  // volumes ("thickness" below is a HALF length, as G4Box/G4Tubs expect)
  const G4double thickness = 0.5*filmThicknessTop;
  const G4double nbBackplaneHalfZ = 0.5*filmThicknessBack;

  // Wall losses (the G4CMP version took -1 = unset as 0)
  G4double pAbsProbSideWallSi = 0;//.01
  if (QuasiparticleConfigManager::GetpAbsProbSideWallSi() != -1.0) pAbsProbSideWallSi = QuasiparticleConfigManager::GetpAbsProbSideWallSi();
  G4double pAbsProbPolishedWallSi = 0; //.0025
  if (QuasiparticleConfigManager::GetpAbsProbPolishedWallSi() != -1.0) pAbsProbPolishedWallSi = QuasiparticleConfigManager::GetpAbsProbPolishedWallSi();

  G4cout << "SR--lQPD_Nb (TFR): patterned film = " << (useAl ? "Al" : "Nb")
         << ", top film thickness = " << filmThicknessTop/nm << " nm"
         << ", Nb backplane thickness = " << filmThicknessBack/nm << " nm"
         << ", numsensors = " << numsensors
         << ", pAbsProbSideWallSi = " << pAbsProbSideWallSi
         << ", pAbsProbPolishedWallSi = " << pAbsProbPolishedWallSi << G4endl;

  // Set to true to have Geant4 check every placement for overlaps. The
  // junction pieces junct1/junct2/junct3 intentionally overlap slightly (the
  // original G4CMP geometry), so this prints warnings; it is off by default
  // because it also slows down start-up (~160 volumes).
  const G4bool checkOverlaps = false;

  //
  // TFR: surface properties for every interface ------------------------------
  // (created once; NOTE: like in the G4CMP version the wall-loss values are
  // therefore latched the first time the geometry is built)
  //
  // Each G4CMPSurfaceProperty is (name, qAbs, qRefl, eMinK, hMinK,
  //   pAbs, pRefl, pSpec, pMinK, qpAbs, qpRefl)
  // A phonon that is neither absorbed (pAbs) nor reflected (pRefl) at the
  // boundary is TRANSMITTED; likewise for BogoliubovQPs with qpAbs/qpRefl
  // (a QP below the gap of the next superconductor is always reflected).
  //
  if (!fConstructed) {
    const G4double GHz = 1e9 * hertz;

    // Reflection coefficients used by the TFR resonator example for all of its
    // interfaces ...
    const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
    const std::vector<G4double> diffCoeffs = {};
    const std::vector<G4double> specCoeffs = {};
    // ... and the (all-diffuse) coefficients that the G4CMP lQPD_Nb version used
    // for the bare silicon walls
    const std::vector<G4double> wallAnhCoeffs = {0, 0, 0, 0, 0, 0};
    const std::vector<G4double> wallDiffCoeffs = {1, 0, 0, 0, 0, 0};
    const std::vector<G4double> wallSpecCoeffs = {0, 0, 0, 0, 0, 0};

    const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external

    // Helper so that each interface below is a single, readable line
    auto makeInterface = [&](const G4String& name, G4double pAbs, G4double pRefl,
                             G4double qpAbs, G4double qpRefl, G4bool siWall) {
      auto* surf = new G4CMPSurfaceProperty(name, 0.0, 1.0, 0.0, 0.0,
                                            pAbs, pRefl, 0.0, 0.0,
                                            qpAbs, qpRefl);
      surf->AddScatteringProperties(anhCutoff, reflCutoff,
                                    siWall ? wallAnhCoeffs : anhCoeffs,
                                    siWall ? wallDiffCoeffs : diffCoeffs,
                                    siWall ? wallSpecCoeffs : specCoeffs,
                                    GHz, GHz, GHz);
      return surf;
    };

    // Silicon <-> film: replaces AlSurfProp/NbSurfProp (+ Kaplan electrode).
    // Phonon transmission is 1 - pAbs - pRefl (see QuasiparticleDetectorParameters.hh).
    // QPs live only in the film, so they always reflect off the silicon.
    fSiAlInterface = makeInterface("SiAlSurf", dp_lQPD_pAbsProbSiFilm,
                                   dp_lQPD_pReflProbSiFilm, 0.0, 1.0, false);
    fSiNbInterface = makeInterface("SiNbSurf", dp_lQPD_pAbsProbSiFilm,
                                   dp_lQPD_pReflProbSiFilm, 0.0, 1.0, false);

    // Film <-> vacuum (top surface and side walls of the film): all reflection
    fAlVacInterface = makeInterface("AlVacSurf", 0.0, 1.0, 0.0, 1.0, false);
    fNbVacInterface = makeInterface("NbVacSurf", 0.0, 1.0, 0.0, 1.0, false);

    // Film <-> film (touching volumes): everything transmits (a QP is only
    // blocked from Al -> Nb below the Nb gap, which the QP process handles)
    fAlAlInterface = makeInterface("AlAlSurf", 0.0, 0.0, 0.0, 0.0, false);
    fNbNbInterface = makeInterface("NbNbSurf", 0.0, 0.0, 0.0, 0.0, false);
    fAlNbInterface = makeInterface("AlNbSurf", 0.0, 0.0, 0.0, 0.0, false);

    // Bare silicon: replaces sidewallSurfProp/polishedwallSurfProp.
    // The absorption probabilities come from the macro; everything that is
    // not absorbed is reflected (diffusely, like in the G4CMP version).
    fSiSideWallInterface = makeInterface("SideWallSurf", pAbsProbSideWallSi,
                                         1.0, 0.0, 1.0, true);
    fSiPolishedWallInterface = makeInterface("polishedWallSurf",
                                             pAbsProbPolishedWallSi, 1.0,
                                             0.0, 1.0, true);
  }

  //
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",8.*cm,8.*cm,8.*cm); // half (-16,16)
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fVacuum,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0); // physical placement

  //
  // Silicon box - this is the volume in which we will propagate phonons
  //
  const G4double siHalfX = 0.5*cm / std::sqrt(numsensors);
  const G4double siHalfY = 0.5*cm / std::sqrt(numsensors);
  const G4double siHalfZ = 0.05*mm;
  G4VSolid* fSiliconSolid = new G4Box("fSiliconSolid", siHalfX, siHalfY, siHalfZ);
  G4LogicalVolume* fSiliconLogical = new G4LogicalVolume(fSiliconSolid,fSilicon,"fSiliconLogical");
  G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),fSiliconLogical,"fSiliconPhysical", worldLogical,false,0);
  // placing physical volume at center of world logical

  //
  //Silicon lattice information
  //

  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");

  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
  SiPhysical->SetMillerOrientation(1,0,0); // how crystal is oriented, also 4 coord orient. (online calculator)
  LM->RegisterLattice(SiPhys, SiPhysical); // connects physical lattice to volume

  // NOTE:  Above registration can also be done in single step:
  // G4LatticlePhysical* SiPhysical = LM->LoadLattice(SiPhys, "Si");

  //
  // TFR: logical lattices of the two superconductors. The "Al" and "Nb"
  // directories of CrystalMaps hold the sc_tau0_qp / sc_tau0_ph values (the
  // numbers that were the Kaplan "phononLifetime" in the G4CMP version).
  //
  G4LatticeLogical* AlLogical = LM->LoadLattice(fAluminum, "Al");
  G4LatticeLogical* NbLogical = LM->LoadLattice(fNiobium, "Nb");

  //
  // Air boxes touching the 4 side faces of the Si crystal (+X, -X, +Y, -Y)
  //
  const G4double airSideThickness = 0.5*cm; // Distance the air extends away from the chip

  // X boxes cover the +X and -X faces of the Si crystal
  G4VSolid* airSideSolidX = new G4Box("airSideSolidX", airSideThickness/2, siHalfY, siHalfZ);
  G4LogicalVolume* airSideLogicalX = new G4LogicalVolume(airSideSolidX, fAir, "airSideLogicalX");
  G4VPhysicalVolume* airSideRightXPhys = new G4PVPlacement(0, G4ThreeVector(siHalfX + airSideThickness/2, 0, 0), airSideLogicalX, "airSideRightX", worldLogical, false, 0);
  G4VPhysicalVolume* airSideLeftXPhys = new G4PVPlacement(0, G4ThreeVector(-siHalfX - airSideThickness/2, 0, 0), airSideLogicalX, "airSideLeftX", worldLogical, false, 0);

  // Y boxes cover the +Y and -Y faces of the Si crystal
  G4VSolid* airSideSolidY = new G4Box("airSideSolidY", siHalfX, airSideThickness/2, siHalfZ);
  G4LogicalVolume* airSideLogicalY = new G4LogicalVolume(airSideSolidY, fAir, "airSideLogicalY");
  G4VPhysicalVolume* airSideRightYPhys = new G4PVPlacement(0, G4ThreeVector(0, siHalfY + airSideThickness/2, 0), airSideLogicalY, "airSideRightY", worldLogical, false, 0);
  G4VPhysicalVolume* airSideLeftYPhys = new G4PVPlacement(0, G4ThreeVector(0, -siHalfY - airSideThickness/2, 0), airSideLogicalY, "airSideLeftY", worldLogical, false, 0);

  //
  // TFR: bookkeeping for the superconducting films ----------------------------
  //
  // Every film volume is registered with RegisterFilm() right after it is
  // placed. RegisterFilm()
  //   1. gives the volume its own G4LatticePhysical with the superconducting
  //      parameters of its metal,
  //   2. creates the two-way Si <-> film and film <-> world border surfaces,
  //   3. remembers the XY bounding box, which is used after all volumes have
  //      been placed to create the film <-> film border surfaces between
  //      volumes that touch.
  //
  struct FilmVolume {
    G4VPhysicalVolume* phys;
    G4bool isAl;			// true: aluminum, false: niobium
    G4double xmin, xmax, ymin, ymax;	// Bounding box in the world XY plane
  };
  std::vector<FilmVolume> films;

  auto RegisterFilm = [&](G4VPhysicalVolume* phys, G4bool isAl,
                          G4double xmin, G4double xmax,
                          G4double ymin, G4double ymax) {
    // 1. Superconducting lattice (one physical lattice per volume).
    //    Arguments: logical lattice, polycrystalline scattering length, Delta0,
    //    Teff, normal-state diffusion constant, QP local trapping time
    G4LatticePhysical* filmLattice = isAl ?
      new G4LatticePhysical(AlLogical, dp_polycryElScatMFP_Al, dp_scDelta0_Al,
                            dp_scTeff_Al, dp_scDn_Al, dp_scTauQPTrap_Al) :
      new G4LatticePhysical(NbLogical, dp_polycryElScatMFP_Nb, dp_scDelta0_Nb,
                            dp_scTeff_Nb, dp_scDn_Nb, dp_scTauQPTrap_Nb);
    filmLattice->SetMillerOrientation(1,0,0);
    LM->RegisterLattice(phys, filmLattice);

    // 2. Border surfaces with the silicon below and the world above/around
    G4CMPSurfaceProperty* siSurf = isAl ? fSiAlInterface : fSiNbInterface;
    G4CMPSurfaceProperty* vacSurf = isAl ? fAlVacInterface : fNbVacInterface;
    const G4String& n = phys->GetName();
    new G4CMPLogicalBorderSurface("border_Si_" + n, SiPhys, phys, siSurf);
    new G4CMPLogicalBorderSurface("border_" + n + "_Si", phys, SiPhys, siSurf);
    new G4CMPLogicalBorderSurface("border_world_" + n, fWorldPhys, phys, vacSurf);
    new G4CMPLogicalBorderSurface("border_" + n + "_world", phys, fWorldPhys, vacSurf);

    // 3. Remember the bounding box for the film <-> film surfaces
    films.push_back({phys, isAl, xmin, xmax, ymin, ymax});
  };

  // Convenience for the (many) G4Box films: centre and half lengths in XY
  auto RegisterBoxFilm = [&](G4VPhysicalVolume* phys, G4bool isAl,
                             G4double cx, G4double cy, G4double hx, G4double hy) {
    RegisterFilm(phys, isAl, cx - hx, cx + hx, cy - hy, cy + hy);
  };

  // Material of the patterned film (feedline, ground planes, capacitor,
  // inductor ...) is selected by /g4cmp/Al. The "sensor" volumes (junction +
  // absorber, names ending in "shunt") are always aluminum.
  G4Material* patternMaterial = useAl ? fAluminum : fNiobium;

  //
  // Metal films on the +Z face of the chip. This is where phonon hits are
  // registered in the G4CMP version. (Every film below is placed on
  // z = siHalfZ + thickness, i.e. it sits directly on the silicon.)
  //

  // Aluminum (or Nb, see above) feedline
  const G4double alFeedlineHalfX = siHalfX;
  const G4double alFeedlineHalfY = 1.5*um;
  // "thickness" (half film thickness) is defined above from the macro values

  // Outer feedline dimensions (100um width)
  const G4double alUPGhalfy = 50.0*um;
  const G4double alLPGhalfy = 2.0*um;
  const G4double feedlineGap = 2.0*um;

  G4Box* feedline = new G4Box("feedlineCenter", alFeedlineHalfX, alFeedlineHalfY, thickness);
  G4Box* uppergroundplane  = new G4Box("uppergroundplane",  alFeedlineHalfX, alUPGhalfy, thickness);
  G4Box* lowergroundplane  = new G4Box("lowergroundplane",  alFeedlineHalfX, alLPGhalfy, thickness);

  G4LogicalVolume* alFLlogical = new G4LogicalVolume(feedline,patternMaterial,"alFLlogical"); // logical feedline
  G4LogicalVolume* alUGPlogical = new G4LogicalVolume(uppergroundplane,patternMaterial,"alUGPlogical"); // logical feedline
  G4LogicalVolume* alLGPlogical = new G4LogicalVolume(lowergroundplane,patternMaterial,"alLGPlogical"); // logical feedline

  G4VPhysicalVolume* alFLphysical = new G4PVPlacement(
    0, G4ThreeVector(0.,0., siHalfZ + thickness), alFLlogical, "alFLphysical",
    worldLogical, false, 0, checkOverlaps); // physical feedline
  RegisterBoxFilm(alFLphysical, useAl, 0., 0., alFeedlineHalfX, alFeedlineHalfY);

  G4double yUGPoffset = alFeedlineHalfY + feedlineGap + alUPGhalfy;
  G4double yLGPoffset = alFeedlineHalfY + feedlineGap + alLPGhalfy;

  G4VPhysicalVolume* alUGPphysical = new G4PVPlacement(
    0, G4ThreeVector(0., yUGPoffset, siHalfZ + thickness), alUGPlogical, "alUGPphysical",
    worldLogical, false, 0, checkOverlaps); // physical UPG
  RegisterBoxFilm(alUGPphysical, useAl, 0., yUGPoffset, alFeedlineHalfX, alUPGhalfy);
  G4VPhysicalVolume* alLGPphysical = new G4PVPlacement(
    0, G4ThreeVector(0., -yLGPoffset, siHalfZ + thickness), alLGPlogical, "alLGPphysical",
    worldLogical, false, 0, checkOverlaps); // physical LGP
  RegisterBoxFilm(alLGPphysical, useAl, 0., -yLGPoffset, alFeedlineHalfX, alLPGhalfy);


  // QPD

  // coupling capacitor
  const G4double couplCapy = 3.0*um;
  const G4double couplCapx = 300*um;
  const G4double couplCapgap = 2.0*um;

  G4Box* couplCap = new G4Box("couplCap", couplCapx, couplCapy, thickness);
  G4LogicalVolume* couplCaplogical = new G4LogicalVolume(couplCap,patternMaterial,"couplCaplogical");

  G4double ycouplcapoffset = yLGPoffset + alLPGhalfy + couplCapgap+ couplCapy;

  G4VPhysicalVolume* couplCapphysical = new G4PVPlacement(
    0, G4ThreeVector(0., -ycouplcapoffset, siHalfZ + thickness), couplCaplogical, "couplCapphysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(couplCapphysical, useAl, 0., -ycouplcapoffset, couplCapx, couplCapy);

  // Placeholder removed for Side Ground Planes

  // connector
  const G4double connectory = 5.0*um;
  const G4double connectorx = 3.0*um;

  G4Box* connector = new G4Box("connector", connectorx, connectory, thickness);
  G4LogicalVolume* connectorlogical = new G4LogicalVolume(connector,patternMaterial,"connectorlogical");

  G4double yconnectoroffset = ycouplcapoffset + couplCapy + connectory;

  G4VPhysicalVolume* connectorphysical = new G4PVPlacement(
    0, G4ThreeVector(0., -yconnectoroffset, siHalfZ + thickness), connectorlogical, "connectorphysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(connectorphysical, useAl, 0., -yconnectoroffset, connectorx, connectory);

  // cap ind connector
  const G4double capindconnectory = 3.0*um;
  const G4double capindconnectorx = 200.5*um;

  G4Box* capindconnector = new G4Box("capindconnector", capindconnectorx, capindconnectory, thickness);
  G4LogicalVolume* capindconnectorlogical = new G4LogicalVolume(capindconnector,patternMaterial,"capindconnectorlogical");

  G4double ycapindconnectoroffset = yconnectoroffset + connectory + capindconnectory;
  G4double xcapindconnectoroffset = -32.5*um;

  G4VPhysicalVolume* capindconnectorphysical = new G4PVPlacement(
    0, G4ThreeVector(xcapindconnectoroffset, -ycapindconnectoroffset, siHalfZ + thickness), capindconnectorlogical, "capindconnectorphysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(capindconnectorphysical, useAl, xcapindconnectoroffset, -ycapindconnectoroffset, capindconnectorx, capindconnectory);

  // indvert1
  const G4double indcenter = 0 - 208*um;

  const G4double indhorzx = 23*um;
  const G4double indhorzy = 1.0*um;

  const G4double indvertx = 1.0*um;
  const G4double indverty = 8.0*um;

  G4Box* indvert1 = new G4Box("indvert1", indvertx, indverty, thickness);
  G4LogicalVolume* indvert1logical = new G4LogicalVolume(indvert1, patternMaterial, "indvert1logical");

  G4double xindvert1offset = indcenter - indhorzx - indvertx;
  G4double yindvert1offset = ycapindconnectoroffset + capindconnectory + indverty;

  G4VPhysicalVolume* indvert1physical = new G4PVPlacement(
    0, G4ThreeVector(xindvert1offset, -yindvert1offset, siHalfZ + thickness), indvert1logical, "indvert1physical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(indvert1physical, useAl, xindvert1offset, -yindvert1offset, indvertx, indverty);

  // Generate meandering inductor pattern from indhorz1 to indhorz40
  // (TFR: the volumes are registered with RegisterBoxFilm inside the loop, so
  //  inductorPhys is no longer needed to attach border surfaces afterwards)
  std::vector<G4LogicalVolume*> inductorLog;

  G4double currentYVertOffset = yindvert1offset;
  G4double yIndHorzLastPos = 0.0;

  // Total Inductor Height (Top of indhorz1 to Bottom of indhorz40):
  // Pitch = 2 * (indverty - indhorzy) = 2 * (8.0 - 1.0) = 14.0 um
  // Distance = 39 * Pitch + 2 * indhorzy
  //          = 39 * 14.0 + 2 * 1.0 = 548.0 um

  for (int i = 1; i <= 40; ++i) {
    // indhorz_i
    G4String hName = "indhorz" + std::to_string(i);
    G4Box* solidH = new G4Box(hName, indhorzx, indhorzy, thickness);
    G4LogicalVolume* logicH = new G4LogicalVolume(solidH, patternMaterial, hName + "logical");
    inductorLog.push_back(logicH);

    G4double yHorzOffset = currentYVertOffset + indverty - indhorzy;
    if (i == 40) yIndHorzLastPos = yHorzOffset;

    G4VPhysicalVolume* physH = new G4PVPlacement(
      0, G4ThreeVector(indcenter, -yHorzOffset, siHalfZ + thickness),
      logicH, hName + "physical", worldLogical, false, 0, checkOverlaps);
    RegisterBoxFilm(physH, useAl, indcenter, -yHorzOffset, indhorzx, indhorzy);

    // indvert_{i+1}
    if (i < 40) {
      int nextVertIdx = i + 1;
      G4String vName = "indvert" + std::to_string(nextVertIdx);
      G4Box* solidV = new G4Box(vName, indvertx, indverty, thickness);
      G4LogicalVolume* logicV = new G4LogicalVolume(solidV, patternMaterial, vName + "logical");
      inductorLog.push_back(logicV);

      G4double xVertOffset = (nextVertIdx % 2 == 0) ?
                             (indcenter + indhorzx + indvertx) : // Even -> Right
                             (indcenter - indhorzx - indvertx);  // Odd -> Left

      G4double yVertOffset = yHorzOffset - indhorzy + indverty;
      currentYVertOffset = yVertOffset;

      G4VPhysicalVolume* physV = new G4PVPlacement(
        0, G4ThreeVector(xVertOffset, -yVertOffset, siHalfZ + thickness),
        logicV, vName + "physical", worldLogical, false, 0, checkOverlaps);
      RegisterBoxFilm(physV, useAl, xVertOffset, -yVertOffset, indvertx, indverty);
    }
  }

  // leftcapwall
  const G4double leftcapwally = 284*um;
  const G4double leftcapwallx = 1*um;

  G4Box* leftcapwall = new G4Box("leftcapwall", leftcapwallx, leftcapwally, thickness);
  G4LogicalVolume* leftcapwalllogical = new G4LogicalVolume(leftcapwall, patternMaterial, "leftcapwalllogical");

  G4double xleftcapwalloffset = 0 - 161*um;
  G4double yleftcapwalloffset = ycapindconnectoroffset + 305*um;

  G4VPhysicalVolume* leftcapwallphysical = new G4PVPlacement(
    0, G4ThreeVector(xleftcapwalloffset, -yleftcapwalloffset, siHalfZ + thickness), leftcapwalllogical, "leftcapwallphysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(leftcapwallphysical, useAl, xleftcapwalloffset, -yleftcapwalloffset, leftcapwallx, leftcapwally);

  // rightcapwall
  const G4double rightcapwally = 291*um;
  const G4double rightcapwallx = 1*um;

  G4Box* rightcapwall = new G4Box("rightcapwall", rightcapwallx, rightcapwally, thickness);
  G4LogicalVolume* rightcapwalllogical = new G4LogicalVolume(rightcapwall, patternMaterial, "rightcapwalllogical");

  G4double xrightcapwalloffset = 0 + 167*um;
  G4double yrightcapwalloffset = ycapindconnectoroffset + 294*um;

  G4VPhysicalVolume* rightcapwallphysical = new G4PVPlacement(
    0, G4ThreeVector(xrightcapwalloffset, -yrightcapwalloffset, siHalfZ + thickness), rightcapwalllogical, "rightcapwallphysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(rightcapwallphysical, useAl, xrightcapwalloffset, -yrightcapwalloffset, rightcapwallx, rightcapwally);

  // Interdigitated Capacitor (IDC) Fingers: caphorz1 to caphorz24
  const G4double xcapcenter = 0 + 0;
  const G4double xcapoffset = 6*um;
  const G4double caphorzy = 1*um;
  const G4double caphorzx = 160*um;
  const G4double captocapgap = 8*um;
  const G4double ycapstart = ycapindconnectoroffset + 22*um;
  // Total IDC height (top of 1 to bottom of 24) = 2*caphorzy + 23*captocapgap = 186.0 um

  std::vector<G4LogicalVolume*> idcLog;

  for (int i = 1; i <= 57; ++i) {
    G4String name = "caphorz" + std::to_string(i);
    G4Box* solid = new G4Box(name, caphorzx, caphorzy, thickness);
    G4LogicalVolume* logic = new G4LogicalVolume(solid, patternMaterial, name + "logical");
    idcLog.push_back(logic);

    // Odd indices shift left, Even indices shift right
    G4double xpos = (i % 2 != 0) ? (xcapcenter) : (xcapcenter + xcapoffset);
    G4double ypos = ycapstart + (i - 1) * captocapgap;

    G4VPhysicalVolume* phys = new G4PVPlacement(
      0, G4ThreeVector(xpos, -ypos, siHalfZ + thickness),
      logic, name + "physical", worldLogical, false, 0, checkOverlaps);
    RegisterBoxFilm(phys, useAl, xpos, -ypos, caphorzx, caphorzy);
  }

  // botcapindconnector1
  const G4double botcapindconnector1x = 1*um;
  const G4double botcapindconnector1y = 8*um;

  G4Box* botcapindconnector1 = new G4Box("botcapindconnector1", botcapindconnector1x, botcapindconnector1y, thickness);
  G4LogicalVolume* botcapindconnector1logical = new G4LogicalVolume(botcapindconnector1, patternMaterial, "botcapindconnector1logical");

  G4double xbotcapindconnector1offset = xindvert1offset;
  G4double ybotcapindconnector1offset = yIndHorzLastPos + 7*um;

  G4VPhysicalVolume* botcapindconnector1physical = new G4PVPlacement(
    0, G4ThreeVector(xbotcapindconnector1offset, -ybotcapindconnector1offset, siHalfZ + thickness),
    botcapindconnector1logical, "botcapindconnector1physical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(botcapindconnector1physical, useAl, xbotcapindconnector1offset, -ybotcapindconnector1offset, botcapindconnector1x, botcapindconnector1y);

  // botcapindconnector2
  // Connects botcapindconnector1 and leftcapwall
  const G4double botcapindconnector2x = 34.5*um; // 69um / 2
  const G4double botcapindconnector2y = 1.0*um;  // 2um / 2

  G4Box* botcapindconnector2 = new G4Box("botcapindconnector2", botcapindconnector2x, botcapindconnector2y, thickness);
  G4LogicalVolume* botcapindconnector2logical = new G4LogicalVolume(botcapindconnector2, patternMaterial, "botcapindconnector2logical");

  G4double xbotcapindconnector2offset = (xleftcapwalloffset + xbotcapindconnector1offset) / 2.0;
  // Align bottom edge with botcapindconnector1 bottom edge to ensure connection
  G4double ybotcapindconnector2offset = (ybotcapindconnector1offset + botcapindconnector1y) - botcapindconnector2y;

  G4VPhysicalVolume* botcapindconnector2physical = new G4PVPlacement(
    0, G4ThreeVector(xbotcapindconnector2offset, -ybotcapindconnector2offset, siHalfZ + thickness),
    botcapindconnector2logical, "botcapindconnector2physical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(botcapindconnector2physical, useAl, xbotcapindconnector2offset, -ybotcapindconnector2offset, botcapindconnector2x, botcapindconnector2y);

  // capjunctconnect
  const G4double capjunctconnectx = 150*um;
  const G4double capjunctconnecty = 1.0*um;

  G4Box* capjunctconnect = new G4Box("capjunctconnect", capjunctconnectx, capjunctconnecty, thickness);
  G4LogicalVolume* capjunctconnectlogical = new G4LogicalVolume(capjunctconnect, patternMaterial, "capjunctconnectlogical");

  // Right align with rightcapwall: (xrightcapwalloffset + rightcapwallx) - capjunctconnectx
  // Top align with rightcapwall bottom: (yrightcapwalloffset + rightcapwally) + capjunctconnecty
  G4double xcapjunctconnectoffset = (xrightcapwalloffset + rightcapwallx) - capjunctconnectx;
  G4double ycapjunctconnectoffset = (yrightcapwalloffset + rightcapwally) + capjunctconnecty;

  G4VPhysicalVolume* capjunctconnectphysical = new G4PVPlacement(
    0, G4ThreeVector(xcapjunctconnectoffset, -ycapjunctconnectoffset, siHalfZ + thickness),
    capjunctconnectlogical, "capjunctconnectphysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(capjunctconnectphysical, useAl, xcapjunctconnectoffset, -ycapjunctconnectoffset, capjunctconnectx, capjunctconnecty);

  // ---- Sensor: junction + absorber ----
  // TFR: these volumes (names ending in "shunt") are ALWAYS aluminum, even
  // when /g4cmp/Al false makes the rest of the pattern Nb. This matches the
  // G4CMP version, where only the pattern switched to NbSurfProp.
  // They are the volumes QuasiparticleSensitivity looks for.

  // junct1
  const G4double junct1x = 5.0*um; // 10um / 2
  const G4double junct1y = 0.05*um; // 0.1um / 2

  const G4double capjunctconnectjunct1gap = 0.25*um;

  G4Box* junct1 = new G4Box("junct1", junct1x, junct1y, thickness);
  G4LogicalVolume* junct1logical = new G4LogicalVolume(junct1, fAluminum, "junct1logical");

  G4double xjunct1offset = 2.0*um;
  // 0.25um below capjunctconnect
  G4double yjunct1offset = ycapjunctconnectoffset + capjunctconnecty + capjunctconnectjunct1gap + junct1y;

  G4VPhysicalVolume* junct1physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct1offset, -yjunct1offset, siHalfZ + thickness),
    junct1logical, "junct1physicalshunt",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(junct1physical, true, xjunct1offset, -yjunct1offset, junct1x, junct1y);

  // junct2
  const G4double junct2x = 0.05*um; // 0.1um / 2
  const G4double junct2y = 0.5*um;  // 1.0um / 2

  G4Box* junct2 = new G4Box("junct2", junct2x, junct2y, thickness);
  G4LogicalVolume* junct2logical = new G4LogicalVolume(junct2, fAluminum, "junct2logical");

  G4double xjunct2offset = xjunct1offset;
  // Connects to center of junct1 and extends downward
  G4double yjunct2offset = yjunct1offset + junct2y;

  G4VPhysicalVolume* junct2physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct2offset, -yjunct2offset, siHalfZ + thickness),
    junct2logical, "junct2physicalshunt",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(junct2physical, true, xjunct2offset, -yjunct2offset, junct2x, junct2y);

  // junct3
  const G4double junct3x = 1.0*um; // 2.0um / 2
  const G4double junct3y = 0.05*um; // 0.1um / 2

  G4Box* junct3 = new G4Box("junct3", junct3x, junct3y, thickness);
  G4LogicalVolume* junct3logical = new G4LogicalVolume(junct3, fAluminum, "junct3logical");

  // Left side at x=0 -> Center at +1.0um
  G4double xjunct3offset = 1.0*um;
  // Intersect junct2 0.2um above bottom of junct2
  // Bottom of junct2 is at -(yjunct2offset + junct2y)
  // Target Y is -(yjunct2offset + junct2y - 0.2*um)
  G4double yjunct3offset = yjunct2offset + junct2y - 0.2*um;

  G4VPhysicalVolume* junct3physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct3offset, -yjunct3offset, siHalfZ + thickness),
    junct3logical, "junct3physicalshunt",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(junct3physical, true, xjunct3offset, -yjunct3offset, junct3x, junct3y);

  // junct4
  const G4double junct4x = 0.05*um; // 0.1um / 2
  const G4double junct4y = 1.0*um;  // 2.0um / 2

  G4Box* junct4 = new G4Box("junct4", junct4x, junct4y, thickness);
  G4LogicalVolume* junct4logical = new G4LogicalVolume(junct4, fAluminum, "junct4logical");

  // Left align with junct3: (xjunct3offset - junct3x) + junct4x
  G4double xjunct4offset = (xjunct3offset - junct3x) + junct4x;
  // Top connects to bottom of junct3: (yjunct3offset + junct3y) + junct4y
  G4double yjunct4offset = yjunct3offset + junct3y + junct4y;

  G4VPhysicalVolume* junct4physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct4offset, -yjunct4offset, siHalfZ + thickness),
    junct4logical, "junct4physicalshunt",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(junct4physical, true, xjunct4offset, -yjunct4offset, junct4x, junct4y);

  // junct5
  const G4double junct5x = 0.5*um; // 1.0um / 2
  const G4double junct5y = 1.0*um; // 2.0um / 2

  G4Box* junct5 = new G4Box("junct5", junct5x, junct5y, thickness);
  G4LogicalVolume* junct5logical = new G4LogicalVolume(junct5, fAluminum, "junct5logical");

  // Centered with junct4
  G4double xjunct5offset = xjunct4offset;
  // Top connects to bottom of junct4
  G4double yjunct5offset = yjunct4offset + junct4y + junct5y;

  G4VPhysicalVolume* junct5physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct5offset, -yjunct5offset, siHalfZ + thickness),
    junct5logical, "junct5physicalshunt",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(junct5physical, true, xjunct5offset, -yjunct5offset, junct5x, junct5y);

  // absorber
  const G4double absorberRadius = 350.0*um;
  // Semicircle extending downwards (180 to 360 degrees)
  G4Tubs* absorber = new G4Tubs("absorber", 0.0, absorberRadius, thickness, 180.0*deg, 180.0*deg);
  G4LogicalVolume* absorberlogical = new G4LogicalVolume(absorber, fAluminum, "absorberlogical");

  G4double xabsorberoffset = xjunct5offset;
  // Flat face connects to bottom of junct5
  G4double yabsorberoffset = yjunct5offset + junct5y;

  G4VPhysicalVolume* absorberphysical = new G4PVPlacement(
    0, G4ThreeVector(xabsorberoffset, -yabsorberoffset, siHalfZ + thickness),
    absorberlogical, "absorberphysicalshunt",
    worldLogical, false, 0, checkOverlaps);
  // Half disk (phi 180..360 deg) occupies local y <= 0, i.e. world y below -yabsorberoffset
  RegisterFilm(absorberphysical, true,
               xabsorberoffset - absorberRadius, xabsorberoffset + absorberRadius,
               -yabsorberoffset - absorberRadius, -yabsorberoffset);

  // Side Ground Planes
  const G4double sideGroundx = 10.0*um; // 20um width / 2

  // They extend from the bottom of the lowergroundplane to the top flat edge of the bridge
  G4double yLGPBottom = yLGPoffset + alLPGhalfy;
  G4double ySideGroundTop = yLGPBottom;
  // The bridge top flat edge is at yabsorberoffset (same y as the absorber flat edge)
  G4double ySideGroundBottom = yabsorberoffset;

  G4double sideGroundy = (ySideGroundBottom - ySideGroundTop) / 2.0;
  // Note: Y coordinates are negative in placement
  G4double ySideGroundCenter = ySideGroundTop + sideGroundy;

  G4Box* sideGround = new G4Box("sideGround", sideGroundx, sideGroundy, thickness);
  G4LogicalVolume* sideGroundlogical = new G4LogicalVolume(sideGround, patternMaterial, "sideGroundlogical");

  // Centers of the side ground planes
  // They wrap around the absorber, so their inner edge aligns with the absorber outer edge.
  // Left: center = xabsorberoffset - absorberRadius - sideGroundx
  // Right: center = xabsorberoffset + absorberRadius + sideGroundx
  G4double xLeftGround = xabsorberoffset - absorberRadius - sideGroundx;
  G4double xRightGround = xabsorberoffset + absorberRadius + sideGroundx;

  G4VPhysicalVolume* leftSideGroundPhysical = new G4PVPlacement(
    0, G4ThreeVector(xLeftGround, -ySideGroundCenter, siHalfZ + thickness),
    sideGroundlogical, "leftSideGroundPhysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(leftSideGroundPhysical, useAl, xLeftGround, -ySideGroundCenter, sideGroundx, sideGroundy);

  G4VPhysicalVolume* rightSideGroundPhysical = new G4PVPlacement(
    0, G4ThreeVector(xRightGround, -ySideGroundCenter, siHalfZ + thickness),
    sideGroundlogical, "rightSideGroundPhysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterBoxFilm(rightSideGroundPhysical, useAl, xRightGround, -ySideGroundCenter, sideGroundx, sideGroundy);

  // Bottom bridge wrapping around the absorber
  const G4double bridgeInnerRadius = absorberRadius;
  const G4double bridgeOuterRadius = absorberRadius + 2.0 * sideGroundx; // 20.0um wide

  G4Tubs* sideGroundBridge = new G4Tubs("sideGroundBridge",
                                        bridgeInnerRadius,
                                        bridgeOuterRadius,
                                        thickness,
                                        180.0*deg, 180.0*deg);
  G4LogicalVolume* sideGroundBridgeLogical =
    new G4LogicalVolume(sideGroundBridge, patternMaterial, "sideGroundBridgeLogical");

  G4VPhysicalVolume* sideGroundBridgePhysical = new G4PVPlacement(
    0, G4ThreeVector(xabsorberoffset, -yabsorberoffset, siHalfZ + thickness),
    sideGroundBridgeLogical, "sideGroundBridgePhysical",
    worldLogical, false, 0, checkOverlaps);
  RegisterFilm(sideGroundBridgePhysical, useAl,
               xabsorberoffset - bridgeOuterRadius, xabsorberoffset + bridgeOuterRadius,
               -yabsorberoffset - bridgeOuterRadius, -yabsorberoffset);

  //
  // Niobium backplane: covers the entire backside (-Z) face of the chip,
  // opposite the patterned metal on the front (+Z) face. Scales with siHalfX/siHalfY.
  // TFR: thickness comes from /g4cmp/filmThicknessNb (nbBackplaneHalfZ above;
  // default 200 nm as before). It is its own film plane, so it is not part of
  // the film <-> film pairing below (it never touches the top-side films); it
  // gets a Nb lattice and Si/Nb and Nb/world border surfaces.
  //
  G4Box* nbBackplaneSolid = new G4Box("nbBackplaneSolid", siHalfX, siHalfY, nbBackplaneHalfZ);
  G4LogicalVolume* nbBackplaneLogical = new G4LogicalVolume(nbBackplaneSolid, fNiobium, "nbBackplaneLogical");

  G4VPhysicalVolume* nbBackplanePhysical = new G4PVPlacement(
    0, G4ThreeVector(0., 0., -(siHalfZ + nbBackplaneHalfZ)), nbBackplaneLogical, "nbBackplanePhysical",
    worldLogical, false, 0, checkOverlaps);

  G4LatticePhysical* nbBackplaneLattice =
    new G4LatticePhysical(NbLogical, dp_polycryElScatMFP_Nb, dp_scDelta0_Nb,
                          dp_scTeff_Nb, dp_scDn_Nb, dp_scTauQPTrap_Nb);
  nbBackplaneLattice->SetMillerOrientation(1,0,0);
  LM->RegisterLattice(nbBackplanePhysical, nbBackplaneLattice);

  // ---- End of Geometry definition ----

  //
  // detector -- Note : "sensitive detector" is attached to Silicon crystal
  // want a phonon sensitive detector, attached to Si crystal
  // TFR: QuasiparticleSensitivity writes the primary + hit files and defines
  // a hit as a phonon entering (or being absorbed on) a "shunt" volume.
  // The detector is created once (it owns the output files) and re-attached
  // to the new silicon volume if the geometry is rebuilt.
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity) {
    fSuperconductorSensitivity = new QuasiparticleSensitivity("PhononElectrode");
    SDman->AddNewDetector(fSuperconductorSensitivity);
  }
  fSiliconLogical->SetSensitiveDetector(fSuperconductorSensitivity);

  //
  // TFR: film <-> film border surfaces --------------------------------------
  //
  // Phonons and QPs can cross from one film volume into the next one, but only
  // if a border surface exists in BOTH directions. The G4CMP version had no
  // such surfaces (the films were independent absorbers). Rather than listing
  // every touching pair by hand (the inductor and the interdigitated capacitor
  // alone are ~100 volumes), create a pair of surfaces for every two films
  // whose XY bounding boxes touch or overlap (within 1 nm). Bounding boxes of
  // the half-disk absorber/bridge are larger than the shapes, so a few pairs
  // that never really touch also get a surface; that is harmless because a
  // border surface is only used when a track actually crosses between the
  // two volumes.
  //
  const G4double touchTolerance = 1.*nm;
  size_t nFilmPairs = 0;
  for (size_t i = 0; i < films.size(); ++i) {
    for (size_t j = i+1; j < films.size(); ++j) {
      const FilmVolume& a = films[i];
      const FilmVolume& b = films[j];
      const G4bool touching =
        a.xmin <= b.xmax + touchTolerance && b.xmin <= a.xmax + touchTolerance &&
        a.ymin <= b.ymax + touchTolerance && b.ymin <= a.ymax + touchTolerance;
      if (!touching) continue;

      G4CMPSurfaceProperty* filmFilm =
        (a.isAl && b.isAl) ? fAlAlInterface :
        (!a.isAl && !b.isAl) ? fNbNbInterface : fAlNbInterface;

      const G4String& na = a.phys->GetName();
      const G4String& nb = b.phys->GetName();
      new G4CMPLogicalBorderSurface("border_" + na + "_" + nb, a.phys, b.phys, filmFilm);
      new G4CMPLogicalBorderSurface("border_" + nb + "_" + na, b.phys, a.phys, filmFilm);
      ++nFilmPairs;
    }
  }
  G4cout << "SR--lQPD_Nb (TFR): " << films.size() << " top-side film volumes, "
         << nFilmPairs << " touching film/film pairs" << G4endl;

  // walls
  // TFR: same volume pairs as the G4CMP version, but with the new Si wall
  // surfaces (absorption from the macro, everything else reflected)
  new G4CMPLogicalBorderSurface("SiToWorld", SiPhys, fWorldPhys, fSiPolishedWallInterface);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideRightXPhys, fSiSideWallInterface);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideLeftXPhys, fSiSideWallInterface);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideRightYPhys, fSiSideWallInterface);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideLeftYPhys, fSiSideWallInterface);

  // backplane (two-way Si <-> Nb, and Nb <-> world, as for the top films)
  new G4CMPLogicalBorderSurface("border_Si_nbBackplane", SiPhys, nbBackplanePhysical, fSiNbInterface);
  new G4CMPLogicalBorderSurface("border_nbBackplane_Si", nbBackplanePhysical, SiPhys, fSiNbInterface);
  new G4CMPLogicalBorderSurface("border_world_nbBackplane", fWorldPhys, nbBackplanePhysical, fNbVacInterface);
  new G4CMPLogicalBorderSurface("border_nbBackplane_world", nbBackplanePhysical, fWorldPhys, fNbVacInterface);

  //
  // Visualization attributes
  //
  // World remains invisible
  G4VisAttributes* wrldVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.1, 0.1));
  wrldVis->SetVisibility(false);
  worldLogical->SetVisAttributes(wrldVis);

  // Silicon crystal: light gray, solid
  G4VisAttributes* siVis = new G4VisAttributes(G4Colour(0.85, 0.85, 0.85, 0.4));
  siVis->SetVisibility(true);
  fSiliconLogical->SetVisAttributes(siVis);

  // Air boxes: light blue, semi-transparent
  G4VisAttributes* airVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.1));
  airVis->SetVisibility(true);
  airSideLogicalX->SetVisAttributes(airVis);
  airSideLogicalY->SetVisAttributes(airVis);

  // Patterned film: red if aluminum (as in the G4CMP version), cyan if niobium
  // (the color the resonator example uses for its niobium ground plane)
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  alVis->SetVisibility(true);
  G4VisAttributes* patternVis = useAl ? alVis
                                      : new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.5));
  patternVis->SetVisibility(true);
  alUGPlogical->SetVisAttributes(patternVis);
  alLGPlogical->SetVisAttributes(patternVis);
  alFLlogical->SetVisAttributes(patternVis);
  connectorlogical->SetVisAttributes(patternVis);
  couplCaplogical->SetVisAttributes(patternVis);
  sideGroundlogical->SetVisAttributes(patternVis);
  sideGroundBridgeLogical->SetVisAttributes(patternVis);
  capindconnectorlogical->SetVisAttributes(patternVis);
  indvert1logical->SetVisAttributes(patternVis);
  for (auto log : inductorLog) {
    log->SetVisAttributes(patternVis);
  }
  leftcapwalllogical->SetVisAttributes(patternVis);
  rightcapwalllogical->SetVisAttributes(patternVis);
  for (auto log : idcLog) {
    log->SetVisAttributes(patternVis);
  }
  botcapindconnector1logical->SetVisAttributes(patternVis);
  botcapindconnector2logical->SetVisAttributes(patternVis);
  capjunctconnectlogical->SetVisAttributes(patternVis);

  // Sensor (junction + absorber): always aluminum, so always red
  junct1logical->SetVisAttributes(alVis);
  junct2logical->SetVisAttributes(alVis);
  junct3logical->SetVisAttributes(alVis);
  junct4logical->SetVisAttributes(alVis);
  junct5logical->SetVisAttributes(alVis);
  absorberlogical->SetVisAttributes(alVis);

  // Niobium backplane: distinct color from the front-side aluminum/niobium pattern
  G4VisAttributes* nbVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5));
  nbVis->SetVisibility(true);
  nbBackplaneLogical->SetVisAttributes(nbVis);

}
