/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef QuasiparticleConfigManager_hh
#define QuasiparticleConfigManager_hh 1

// $Id$
// File:  QuasiparticleConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		quasiparticle example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see QuasiparticleConfigMessenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.
//
// TFR port of the RISQTutorial lQPD_Nb configuration variables.  The values
// below are set from macros (see QuasiparticleConfigMessenger) BEFORE
// /run/initialize, and are read when the detector is constructed
// (QuasiparticleDetectorConstruction_lQPD_Nb.cc) and when the sensitive
// detector opens its output files (QuasiparticleSensitivity.cc).
//   - Primary_file             : second output file (one line per primary)
//   - filmThicknessAl/Nb       : film thicknesses in nm (see NOTE below)
//   - Al                       : true -> Al patterned film, false -> Nb
//   - pAbsProbSideWallSi       : phonon absorption prob. on the Si side walls
//   - pAbsProbPolishedWallSi   : phonon absorption prob. on polished Si faces
//   - numsensors               : number of QPD sensors tiled on a 1 cm x 1 cm
//                                area (scales the chip size)
//   - StepOutput               : enable/disable the (very large) QP step file
// A value of -1 means "not set by the user" for the numeric variables, so the
// detector construction can fall back to its default.
//
// NOTE: in the G4CMP (Kaplan) version filmThickness* was only a Kaplan model
// parameter.  In the TFR version it sets the real geometric thickness of the
// thin films, because BogoliubovQPs live inside real superconducting volumes.

#include "globals.hh"

class QuasiparticleConfigMessenger;


class QuasiparticleConfigManager {
public:
  ~QuasiparticleConfigManager();	// Must be public for end-of-job cleanup
  static QuasiparticleConfigManager* Instance();   // Only needed by static accessors

  // Access current values
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }
  static const G4String& GetPrimaryOutput()  { return Instance()->Primary_file; }
  static double GetfilmThicknessAl();
  static double GetfilmThicknessNb();
  static bool GetAl();
  static double GetpAbsProbSideWallSi();
  static double GetpAbsProbPolishedWallSi();
  static double Getnumsensors();
  static bool GetStepOutput();

  // Change values (e.g., via Messenger)
  static void SetHitOutput(const G4String& name)
  { Instance()->Hit_file=name; UpdateGeometry(); }
  static void SetPrimaryOutput(const G4String& name)
  { Instance()->Primary_file=name; UpdateGeometry(); }
  static void SetfilmThicknessAl(G4double filmThicknessAl);
  static void SetfilmThicknessNb(G4double filmThicknessNb);
  static void SetAl(bool Al);
  static void SetpAbsProbSideWallSi(G4double pAbs);
  static void SetpAbsProbPolishedWallSi(G4double pAbs);
  static void Setnumsensors(G4double num);
  static void SetStepOutput(bool enable);

  static void UpdateGeometry();

private:
  QuasiparticleConfigManager();		// Singleton: only constructed on request
  QuasiparticleConfigManager(const QuasiparticleConfigManager&) = delete;
  QuasiparticleConfigManager(QuasiparticleConfigManager&&) = delete;
  QuasiparticleConfigManager& operator=(const QuasiparticleConfigManager&) = delete;
  QuasiparticleConfigManager& operator=(QuasiparticleConfigManager&&) = delete;

  static QuasiparticleConfigManager* theInstance;

private:
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)
  G4String Primary_file;	// Output file of primaries ($G4CMP_PRIMARY_FILE)
  G4double filmThicknessAl;	// [nm]; -1 = unset
  G4double filmThicknessNb;	// [nm]; -1 = unset
  bool Al;			// true: Al patterned film; false: Nb patterned film
  G4double pAbsProbSideWallSi;		// -1 = unset
  G4double pAbsProbPolishedWallSi;	// -1 = unset
  G4double numsensors;		// -1 = unset
  bool StepOutput;		// write QuasiparticleStepInformationFile.txt?

  QuasiparticleConfigMessenger* messenger;
};

#endif	/* QuasiparticleConfigManager_hh */
