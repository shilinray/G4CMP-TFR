/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  QuasiparticleConfigMessenger.cc
//
// Description:	Macro command defitions to set user configuration in
//		QuasiparticleConfigManager.
//
// 20170816  Michael Kelsey

#include "QuasiparticleConfigMessenger.hh"
#include "QuasiparticleConfigManager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"


// Constructor and destructor

QuasiparticleConfigMessenger::
QuasiparticleConfigMessenger(QuasiparticleConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), hitsCmd(0) {
  hitsCmd =
    CreateCommand<G4UIcmdWithAString>("HitsFile",
                                      "Set filename for output of phonon hit locations");

  // ---- Commands ported from RISQTutorialConfigMessenger (lQPD_Nb macros) ----
  // File names.  "hitsFileName" does the same thing as the original TFR
  // "HitsFile" command; it is kept so that the RISQ macros/scripts (which
  // use sed to fill in "primFileName <file>" and "hitsFileName <file>") work
  // unchanged.
  primFileNameCmd = CreateCommand<G4UIcmdWithAString>("primFileName",
			      "Set text file name for primary particle information");
  hitsFileNameCmd = CreateCommand<G4UIcmdWithAString>("hitsFileName",
			      "Set text file name for phonon hit information");

  // Film thicknesses [nm] (geometric thickness in the TFR version)
  filmThicknessAlCmd = CreateCommand<G4UIcmdWithADouble>("filmThicknessAl",
            "Set the Al film thickness in nm (<=0 uses the default)");
  filmThicknessNbCmd = CreateCommand<G4UIcmdWithADouble>("filmThicknessNb",
            "Set the Nb film thickness in nm (<=0 uses the default)");

  // true: patterned film is Al; false: patterned film is Nb
  AlCmd = CreateCommand<G4UIcmdWithABool>("Al",
            "Set the boolean parameter for Al (true: Al film, false: Nb film)");

  // Phonon absorption probabilities on the bare Si surfaces
  pAbsProbSideWallSiCmd = CreateCommand<G4UIcmdWithADouble>("pAbsProbSideWallSi",
            "Set the phonon absorption probability on the Si side walls");
  pAbsProbPolishedWallSiCmd = CreateCommand<G4UIcmdWithADouble>("pAbsProbPolishedWallSi",
            "Set the phonon absorption probability on the polished Si faces");

  // Number of sensors tiled over a 1 cm x 1 cm area (scales chip size)
  numsensorsCmd = CreateCommand<G4UIcmdWithADouble>("numsensors",
            "Set the number of sensors (chip half-width = 0.5 cm / sqrt(numsensors))");

  // New in the TFR port: switch off the very large QP step output file
  stepOutputCmd = CreateCommand<G4UIcmdWithABool>("stepOutput",
            "Enable/disable QuasiparticleStepInformationFile.txt (default true)");
}


QuasiparticleConfigMessenger::~QuasiparticleConfigMessenger() {
  delete hitsCmd; hitsCmd=0;
  delete primFileNameCmd; primFileNameCmd=0;
  delete hitsFileNameCmd; hitsFileNameCmd=0;
  delete filmThicknessAlCmd; filmThicknessAlCmd=0;
  delete filmThicknessNbCmd; filmThicknessNbCmd=0;
  delete AlCmd; AlCmd=0;
  delete pAbsProbSideWallSiCmd; pAbsProbSideWallSiCmd=0;
  delete pAbsProbPolishedWallSiCmd; pAbsProbPolishedWallSiCmd=0;
  delete numsensorsCmd; numsensorsCmd=0;
  delete stepOutputCmd; stepOutputCmd=0;
}


// Parse user input and add to configuration

void QuasiparticleConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == hitsCmd) theManager->SetHitOutput(value);

  // Ported RISQTutorial commands
  if (cmd == primFileNameCmd) theManager->SetPrimaryOutput(value);
  if (cmd == hitsFileNameCmd) theManager->SetHitOutput(value);
  if (cmd == filmThicknessAlCmd) theManager->SetfilmThicknessAl(filmThicknessAlCmd->GetNewDoubleValue(value));
  if (cmd == filmThicknessNbCmd) theManager->SetfilmThicknessNb(filmThicknessNbCmd->GetNewDoubleValue(value));
  if (cmd == AlCmd) theManager->SetAl(AlCmd->GetNewBoolValue(value));
  if (cmd == pAbsProbSideWallSiCmd) theManager->SetpAbsProbSideWallSi(pAbsProbSideWallSiCmd->GetNewDoubleValue(value));
  if (cmd == pAbsProbPolishedWallSiCmd) theManager->SetpAbsProbPolishedWallSi(pAbsProbPolishedWallSiCmd->GetNewDoubleValue(value));
  if (cmd == numsensorsCmd) theManager->Setnumsensors(numsensorsCmd->GetNewDoubleValue(value));
  if (cmd == stepOutputCmd) theManager->SetStepOutput(stepOutputCmd->GetNewBoolValue(value));
}
