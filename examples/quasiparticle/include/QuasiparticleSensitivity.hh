/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef QuasiparticleSensitivity_h
#define QuasiparticleSensitivity_h 1

// TFR port of RISQTutorialSensitivity (see QuasiparticleSensitivity.cc for
// the definition of a "hit" and for the output file formats).
//   - Writes TWO text files (as in RISQTutorial): one with the primary
//     particle(s) of each event and one with the hits.
//   - File names come from QuasiparticleConfigManager (/g4cmp/primFileName and
//     /g4cmp/hitsFileName), positions are written in mm so that the migrated
//     AnalysisTools/RISQTutorialAnalysis_lQPD_Nb.cc can read them unchanged.

#include "G4CMPElectrodeSensitivity.hh"
#include <fstream>
#include <set>

class QuasiparticleSensitivity final : public G4CMPElectrodeSensitivity {
public:
  QuasiparticleSensitivity(G4String name);
  virtual ~QuasiparticleSensitivity();
  // No copies
  QuasiparticleSensitivity(const QuasiparticleSensitivity&) = delete;
  QuasiparticleSensitivity& operator=(const QuasiparticleSensitivity&) = delete;
  /* Move is disabled for now because old versions of GCC can't move ofstream
  // Move OK
  QuasiparticleSensitivity(QuasiparticleSensitivity&&);
  QuasiparticleSensitivity& operator=(QuasiparticleSensitivity&&);
  */
  QuasiparticleSensitivity(QuasiparticleSensitivity&&) = delete;
  QuasiparticleSensitivity& operator=(QuasiparticleSensitivity&&) = delete;

  virtual void Initialize(G4HCofThisEvent*) override;
  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetHitOutputFile(const G4String& fn);
  void SetPrimaryOutputFile(const G4String& fn);

protected:
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream primaryOutput;
  std::ofstream hitOutput;
  G4String primaryFileName;
  G4String hitFileName;

  // Tracks that have already been recorded entering a "shunt" volume in the
  // current event. Cleared in Initialize(); used so that a phonon that enters,
  // leaves and re-enters the sensor is counted (and its energy summed) once.
  std::set<G4int> fEnteredShuntTracks;
};

#endif
