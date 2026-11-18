/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  QuasiparticleConfigManager.cc
//
// Description:	Singleton container class for user configuration of G4CMP
//		quasiparticle example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see QuasiparticleConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "QuasiparticleConfigManager.hh"
#include "QuasiparticleConfigMessenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>


// Constructor and Singleton Initializer

QuasiparticleConfigManager* QuasiparticleConfigManager::theInstance = 0;

QuasiparticleConfigManager* QuasiparticleConfigManager::Instance() {
  if (!theInstance) theInstance = new QuasiparticleConfigManager;
  return theInstance;
}

QuasiparticleConfigManager::QuasiparticleConfigManager()
  : Hit_file(getenv("G4CMP_HIT_FILE")?getenv("G4CMP_HIT_FILE"):"phonon_hits.txt"),
    Primary_file(getenv("G4CMP_PRIMARY_FILE")?getenv("G4CMP_PRIMARY_FILE"):"phonon_primary.txt"),
    filmThicknessAl(-1.0), filmThicknessNb(-1.0), Al(true),
    pAbsProbSideWallSi(-1.0), pAbsProbPolishedWallSi(-1.0), numsensors(-1.0),
    StepOutput(true),	// Default keeps the original TFR behavior (step file written)
    messenger(new QuasiparticleConfigMessenger(this)) {;}

QuasiparticleConfigManager::~QuasiparticleConfigManager() {
  delete messenger; messenger=0;
}


// Trigger rebuild of geometry if parameters change

void QuasiparticleConfigManager::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

// Setters used by the Messenger.  Unlike SetHitOutput()/SetPrimaryOutput()
// (inline in the header), these do NOT call UpdateGeometry(): the macros set
// them before /run/initialize, so the geometry is built once with the final
// values.

void QuasiparticleConfigManager::SetfilmThicknessAl(G4double filmThicknessAl) {
  Instance()->filmThicknessAl = filmThicknessAl;
}

void QuasiparticleConfigManager::SetfilmThicknessNb(G4double filmThicknessNb) {
  Instance()->filmThicknessNb = filmThicknessNb;
}

void QuasiparticleConfigManager::SetAl(bool Al) {
  Instance()->Al = Al;
}

void QuasiparticleConfigManager::SetpAbsProbSideWallSi(G4double pAbs) {
  Instance()->pAbsProbSideWallSi = pAbs;
}

void QuasiparticleConfigManager::SetpAbsProbPolishedWallSi(G4double pAbs) {
  Instance()->pAbsProbPolishedWallSi = pAbs;
}

void QuasiparticleConfigManager::Setnumsensors(G4double num) {
  Instance()->numsensors = num;
}

void QuasiparticleConfigManager::SetStepOutput(bool enable) {
  Instance()->StepOutput = enable;
}

double QuasiparticleConfigManager::GetfilmThicknessAl() { return Instance()->filmThicknessAl; }
double QuasiparticleConfigManager::GetfilmThicknessNb() { return Instance()->filmThicknessNb; }
bool QuasiparticleConfigManager::GetAl() { return Instance()->Al; }
double QuasiparticleConfigManager::GetpAbsProbSideWallSi() { return Instance()->pAbsProbSideWallSi; }
double QuasiparticleConfigManager::GetpAbsProbPolishedWallSi() { return Instance()->pAbsProbPolishedWallSi; }
double QuasiparticleConfigManager::Getnumsensors() { return Instance()->numsensors; }
bool QuasiparticleConfigManager::GetStepOutput() { return Instance()->StepOutput; }
