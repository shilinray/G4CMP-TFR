/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef QuasiparticleConfigMessenger_hh
#define QuasiparticleConfigMessenger_hh 1

// $Id$
// File:  QuasiparticleConfigMessenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		QuasiparticleConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class QuasiparticleConfigManager;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcommand;


class QuasiparticleConfigMessenger : public G4UImessenger {
public:
  QuasiparticleConfigMessenger(QuasiparticleConfigManager* theData);
  virtual ~QuasiparticleConfigMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  QuasiparticleConfigManager* theManager;
  G4UIcmdWithAString* hitsCmd;		// /g4cmp/HitsFile (original TFR command)

  // Commands ported from RISQTutorialConfigMessenger (lQPD_Nb macros).
  // All live in the /g4cmp/ directory, e.g. "/g4cmp/numsensors 4".
  G4UIcmdWithAString* primFileNameCmd;	// /g4cmp/primFileName
  G4UIcmdWithAString* hitsFileNameCmd;	// /g4cmp/hitsFileName
  G4UIcmdWithADouble* filmThicknessAlCmd;	// /g4cmp/filmThicknessAl [nm]
  G4UIcmdWithADouble* filmThicknessNbCmd;	// /g4cmp/filmThicknessNb [nm]
  G4UIcmdWithABool* AlCmd;			// /g4cmp/Al true|false
  G4UIcmdWithADouble* pAbsProbSideWallSiCmd;	// /g4cmp/pAbsProbSideWallSi
  G4UIcmdWithADouble* pAbsProbPolishedWallSiCmd;// /g4cmp/pAbsProbPolishedWallSi
  G4UIcmdWithADouble* numsensorsCmd;		// /g4cmp/numsensors
  G4UIcmdWithABool* stepOutputCmd;		// /g4cmp/stepOutput true|false (new in TFR)

private:
  QuasiparticleConfigMessenger(const QuasiparticleConfigMessenger&);	// Copying is forbidden
  QuasiparticleConfigMessenger& operator=(const QuasiparticleConfigMessenger&);
};

#endif /* QuasiparticleConfigMessenger_hh */
