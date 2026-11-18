/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "QuasiparticleSensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPUtils.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "QuasiparticleConfigManager.hh"
#include <fstream>
#include <iostream>


// Both output files are opened here (at detector construction, i.e. at
// /run/initialize), so /g4cmp/primFileName and /g4cmp/hitsFileName must be
// issued BEFORE /run/initialize.
QuasiparticleSensitivity::QuasiparticleSensitivity(G4String name) :
  G4CMPElectrodeSensitivity(name), primaryFileName(""), hitFileName("") {
  SetHitOutputFile(QuasiparticleConfigManager::GetHitOutput());
  SetPrimaryOutputFile(QuasiparticleConfigManager::GetPrimaryOutput());
}

QuasiparticleSensitivity::~QuasiparticleSensitivity() {
  //Close file and check: primaries
  if (primaryOutput.is_open()) primaryOutput.close();
  if (!primaryOutput.good()) {
    G4cerr << "Error closing primary output file, " << primaryFileName << ".\n"
           << "Expect bad things like loss of data.";
  }

  //Close file and check: hits
  if (hitOutput.is_open()) hitOutput.close();
  if (!hitOutput.good()) {
    G4cerr << "Error closing hit output file, " << hitFileName << ".\n"
           << "Expect bad things like loss of data.";
  }
}

// New event: let the base class create the hit collection, then forget which
// tracks have already entered the sensor in the previous event
void QuasiparticleSensitivity::Initialize(G4HCofThisEvent* HCE) {
  G4CMPElectrodeSensitivity::Initialize(HCE);
  fEnteredShuntTracks.clear();
}

void QuasiparticleSensitivity::EndOfEvent(G4HCofThisEvent* HCE) {
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  auto* hitCol = static_cast<G4CMPElectrodeHitsCollection*>(HCE->GetHC(HCID));
  std::vector<G4CMPElectrodeHit*>* hitVec = hitCol->GetVector();

  G4RunManager* runMan = G4RunManager::GetRunManager();

  // Primary output: one line per event (whitespace separated, mm and ns), in
  // the order  runID eventID particleName energy[eV] X Y Z T
  if (primaryOutput.good()) {
    primaryOutput << runMan->GetCurrentRun()->GetRunID() << " "
		  << runMan->GetCurrentEvent()->GetEventID() << " "
		  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName() << " "
		  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy()/eV << " "
		  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetX0()/mm << " "
		  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetY0()/mm << " "
		  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetZ0()/mm << " "
		  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetT0()/ns << "\n";
  }

  // Hit output: one line per hit (whitespace separated, mm and ns)
  if (hitOutput.good()) {
    for (G4CMPElectrodeHit* hit : *hitVec) {
      hitOutput << runMan->GetCurrentRun()->GetRunID() << ' '
                << runMan->GetCurrentEvent()->GetEventID() << ' '
                << hit->GetTrackID() << ' '
                << hit->GetParticleName() << ' '
                << hit->GetStartEnergy()/eV << ' '
                << hit->GetStartPosition().getX()/mm << ' '
                << hit->GetStartPosition().getY()/mm << ' '
                << hit->GetStartPosition().getZ()/mm << ' '
                << hit->GetStartTime()/ns << ' '
                << hit->GetEnergyDeposit()/eV << ' '
                << hit->GetWeight() << ' '
                << hit->GetFinalPosition().getX()/mm << ' '
                << hit->GetFinalPosition().getY()/mm << ' '
                << hit->GetFinalPosition().getZ()/mm << ' '
                << hit->GetFinalTime()/ns << '\n';
    }
  }
}

void QuasiparticleSensitivity::SetHitOutputFile(const G4String &fn) {
  if (hitFileName != fn) {
    if (hitOutput.is_open()) hitOutput.close();
    hitFileName = fn;
    hitOutput.open(hitFileName, std::ios_base::trunc);
    if (!hitOutput.good()) {
      G4ExceptionDescription msg;
      msg << "Error opening hit output file " << hitFileName;
      G4Exception("QuasiparticleSensitivity::SetHitOutputFile", "PhonSense003",
                  FatalException, msg);
      hitOutput.close();
    } else {
      // 15 columns, matching what is written in EndOfEvent()
      hitOutput << "Run ID,Event ID,Track ID,Particle Name,Start Energy [eV],"
		<< "Start X [mm],Start Y [mm],Start Z [mm],Start Time [ns],"
		<< "Energy Deposited [eV],Track Weight,End X [mm],End Y [mm],End Z [mm],"
		<< "Final Time [ns]\n";
    }
  }
}

void QuasiparticleSensitivity::SetPrimaryOutputFile(const G4String &fn) {
  if (primaryFileName != fn) {
    if (primaryOutput.is_open()) primaryOutput.close();
    primaryFileName = fn;
    primaryOutput.open(primaryFileName, std::ios_base::trunc);
    if (!primaryOutput.good()) {
      G4ExceptionDescription msg;
      msg << "Error opening output file " << primaryFileName;
      G4Exception("QuasiparticleSensitivity::SetPrimaryOutputFile", "PhonSense003",
                  FatalException, msg);
      primaryOutput.close();
    } else {
      primaryOutput << "Run ID,Event ID,Particle Name,Start Energy [eV],"
		    << "Start X [mm],Start Y [mm],Start Z [mm],Start Time [ns]\n";
    }
  }
}


// Hit definition for the TFR lQPD_Nb device.
//
// In the G4CMP (Kaplan) version, a hit was a phonon KILLED at the boundary of
// a "shunt" volume with a non-zero non-ionizing energy deposit, because the
// Kaplan electrode absorbed the phonon at the surface. In the TFR physics the
// Si <-> superconductor boundary process decides (using the surface's pAbs,
// pRefl and transmission probabilities) whether the phonon ENTERS the film;
// a transmitted phonon is not killed and goes on to break Cooper pairs. So a
// hit is now either
//   (1) the old criterion (phonon killed on a boundary of a shunt volume with
//       energy deposited -- still possible if the Si/film pAbsProb is > 0), or
//   (2) a phonon that crosses a boundary from the silicon (the volume this
//       detector is attached to) INTO a "shunt" volume, i.e. the TFR
//       boundary process let it enter the superconducting sensor.
// Case (2) is counted only the first time a given track enters the sensor in
// an event, and its "Energy Deposited" is the phonon energy at entry.
G4bool QuasiparticleSensitivity::IsHit(const G4Step* step,
                                const G4TouchableHistory*) const {
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4ParticleDefinition* particle = track->GetDefinition();

  G4bool correctParticle = particle == G4PhononLong::Definition() ||
                           particle == G4PhononTransFast::Definition() ||
                           particle == G4PhononTransSlow::Definition();
  if (!correctParticle) return false;

  // The phonon has to be at a volume boundary and the volume it lands in
  // (the post-step volume) must be one of the sensor ("shunt") volumes.
  // post-step volume is null if the phonon leaves the world.
  if (postStepPoint->GetStepStatus() != fGeomBoundary) return false;
  G4VPhysicalVolume* postPV = postStepPoint->GetPhysicalVolume();
  G4bool landedOnTargetSurface =
    postPV && (postPV->GetName().find("shunt") != std::string::npos);
  if (!landedOnTargetSurface) return false;

  // Option (1): old criterion, phonon absorbed (killed with deposit) here
  G4bool absorbed = track->GetTrackStatus() == fStopAndKill &&
                    step->GetNonIonizingEnergyDeposit() > 0.;
  if (absorbed) return true;

  // Option (2): phonon entered the sensor (once per track and event).
  // Geant4 reports the NEXT volume as the post-step volume even when the
  // boundary process reflects the phonon back into the silicon, so the
  // post-step volume alone cannot tell "transmitted" from "reflected".
  // The boundary process has already updated the post-step momentum direction,
  // so look 1 nm ahead along it: if we are in the shunt volume the phonon was
  // transmitted; if we are back in the silicon it was reflected.
  if (track->GetTrackStatus() == fStopAndKill) return false;
  if (fEnteredShuntTracks.find(track->GetTrackID()) !=
      fEnteredShuntTracks.end()) return false;

  G4ThreeVector probe = postStepPoint->GetPosition() +
                        1.*nm * postStepPoint->GetMomentumDirection();
  return G4CMP::GetVolumeAtPoint(probe) == postPV;
}


// Same as G4CMPElectrodeSensitivity::ProcessHits(), plus bookkeeping for
// option (2) above.
G4bool QuasiparticleSensitivity::ProcessHits(G4Step* aStep,
                                             G4TouchableHistory* ROhist) {
  if (IsHit(aStep, ROhist)) {
    auto hit = new G4CMPElectrodeHit;
    G4CMP::FillHit(aStep, hit); // Mutates hit

    // A phonon that was transmitted into the film deposits nothing at the
    // boundary: record the energy it carried into the sensor instead
    if (hit->GetEnergyDeposit() <= 0.) {
      hit->SetEnergyDeposit(aStep->GetPostStepPoint()->GetKineticEnergy());
    }

    fEnteredShuntTracks.insert(aStep->GetTrack()->GetTrackID());
    hitsCollection->insert(hit);
  }

  return true;
}
