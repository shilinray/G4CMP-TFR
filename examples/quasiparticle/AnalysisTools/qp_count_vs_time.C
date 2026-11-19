//---------------------------------------------------------
//
// qp_count_vs_time.C
//
// ROOT macro that counts the BogoliubovQPs present vs time in the volumes
// whose name contains a given string (default "shunt": junct1..junct5 and the
// absorber of the lQPD_Nb geometry), from the QP step file
// QuasiparticleStepInformationFile.txt (written by QuasiparticleSteppingAction
// when /g4cmp/stepOutput is true, which is the default).
//
// Every line of the step file is one QP step:
//   run event track particle  preX preY preZ preT preE preKE
//                             postX postY postZ postT postE postKE
//                             nReflections process preVolume postVolume
// (mm, ns, eV). A QP is counted as present in its pre-step volume for the
// whole time span [preT, postT] of the step. The number of QPs at a time t is
// therefore the number of steps with preT <= t < postT in that volume.
//
// Energy collection: a QP counts as collected the first time one of its steps
// in a volume containing the search string ends with its kinetic energy at the
// gap: postKE <= (1 + gapTolerance) * gapEnergy. Its energy at that step (about
// the gap) is added to the collected energy; every QP track is counted once
// (a QP that is already at the gap when it enters such a volume counts there).
// The energy collected per event is the sum over the QPs, averaged over the
// events in the file, and the collection efficiency is that divided by the
// initial energy (the photon energy). The energy collection is not limited to
// the time window of the plot.
//
// Usage:
//   root -l -b -q 'qp_count_vs_time.C("QuasiparticleStepInformationFile.txt")'
//   root -l 'qp_count_vs_time.C("QuasiparticleStepInformationFile.txt", "shunt", 0, 10000, 1000)'
//   root -l 'qp_count_at_time("QuasiparticleStepInformationFile.txt", 500.)'
//
// Notes:
//  * Counts are averaged over the events found in the step file (events in
//    which no QP was created are not in the file and are not counted).
//  * The step file has no track weight column, so every QP counts as 1.
//  * To get a density, divide by the volume (the file has no volumes), or pass
//    the total volume of the selected volumes in um^3 as the totalVolume_um3
//    argument of qp_count_vs_time().
//  * The defaults for the gap (176 ueV, Al), the tolerance (10%, about
//    k_B Teff / gap for Teff = 0.2 K) and the initial energy (5.78993 meV, a
//    1.4 THz photon) match the lQPD_Nb photon macro; change them for other runs.
//
//---------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"

// One QP step, the columns of the step file that this macro needs
struct QPCVStep {
  int run, event, track;
  double preT, postT;	// ns
  double postKE;	// eV, QP energy at the end of the step
  std::string preVolume;
};

// Read the next QP step; returns false at the end of the file
bool QPCVReadStep(std::ifstream& in, QPCVStep& s) {
  std::string particle, process, postVolume;
  double x, y, z, e, ke, nRefl;
  if (!(in >> s.run >> s.event >> s.track >> particle)) return false;
  in >> x >> y >> z >> s.preT >> e >> ke;
  in >> x >> y >> z >> s.postT >> e >> s.postKE;
  in >> nRefl >> process >> s.preVolume >> postVolume;
  return bool(in) && particle.find("BogoliubovQP") != std::string::npos;
}


// Plot the number of QPs vs time in the volumes whose name contains
// volumeContains: one histogram per volume plus the total. The bin content is
// the time-averaged number of QPs present during that bin, per event.
// Also prints the energy collected by QPs that reach the gap in those volumes
// (see the top of the file) and the collection efficiency, and returns them
// through energyCollected_eV and efficiency if those pointers are given.
// Returns the histogram of the total (nullptr if the file cannot be read).
TH1D* qp_count_vs_time(const char* stepFile = "QuasiparticleStepInformationFile.txt",
                       const char* volumeContains = "shunt",
                       double tMin_ns = 0., double tMax_ns = 10000.,
                       int nBins = 1000, bool perVolume = true,
                       bool logY = false,
                       const char* outFile = "qp_count_vs_time.png",
                       double totalVolume_um3 = 0.,
                       double gapEnergy_eV = 176e-6,
                       double gapTolerance = 0.1,
                       double initialEnergy_eV = 5.78993e-3,
                       double* energyCollected_eV = nullptr,
                       double* efficiency = nullptr) {
  std::ifstream in(stepFile);
  if (!in.is_open()) {
    std::cerr << "qp_count_vs_time: cannot open " << stepFile << std::endl;
    return nullptr;
  }
  if (nBins < 1 || tMax_ns <= tMin_ns) {
    std::cerr << "qp_count_vs_time: bad time range or number of bins" << std::endl;
    return nullptr;
  }

  const double binWidth = (tMax_ns - tMin_ns) / nBins;
  const double gapThreshold = (1. + gapTolerance) * gapEnergy_eV;
  std::map<std::string, std::vector<double> > sums;	// volume -> bin sums
  std::map<std::string, double> collectedByVolume;	// eV, summed over events
  std::set<std::tuple<int,int,int> > collectedTracks;
  double collectedTotal = 0.;				// eV, summed over events
  std::set<std::pair<int,int> > events;
  QPCVStep s;
  while (in.good()) {
    if (!QPCVReadStep(in, s)) {
      if (in.eof()) break;
      in.clear();		// skip a line that is not a QP step
      in.ignore(100000, '\n');
      continue;
    }
    events.insert(std::make_pair(s.run, s.event));
    if (s.preVolume.find(volumeContains) == std::string::npos) continue;

    // First time this QP is at the gap in a selected volume: collect it
    if (s.postKE <= gapThreshold &&
        collectedTracks.insert(std::make_tuple(s.run, s.event, s.track)).second) {
      collectedTotal += s.postKE;
      collectedByVolume[s.preVolume] += s.postKE;
    }

    double t0 = std::max(s.preT, tMin_ns);
    double t1 = std::min(s.postT, tMax_ns);
    if (t1 <= t0) continue;

    std::vector<double>& v = sums[s.preVolume];
    if (v.empty()) v.assign(nBins, 0.);
    int first = std::max(0, (int)std::floor((t0 - tMin_ns) / binWidth));
    int last = std::min(nBins - 1, (int)std::floor((t1 - tMin_ns) / binWidth));
    for (int b = first; b <= last; ++b) {
      double lo = std::max(t0, tMin_ns + b*binWidth);
      double hi = std::min(t1, tMin_ns + (b+1)*binWidth);
      if (hi > lo) v[b] += (hi - lo) / binWidth;
    }
  }

  const double norm = events.empty() ? 1. : 1. / events.size();
  const double volumeScale = (totalVolume_um3 > 0.) ? 1. / totalVolume_um3 : 1.;
  std::cout << "qp_count_vs_time: " << events.size() << " events, "
            << sums.size() << " volumes containing \"" << volumeContains
            << "\"" << std::endl;

  // Energy collected per event by the QPs that reached the gap
  const double collected_eV = collectedTotal * norm;
  const double eff = (initialEnergy_eV > 0.) ? collected_eV / initialEnergy_eV : 0.;
  std::cout << "QPs that reached the gap (E <= " << gapThreshold * 1e6
            << " ueV) in volumes containing \"" << volumeContains << "\": "
            << collectedTracks.size() << " in " << events.size() << " events"
            << std::endl;
  for (std::map<std::string, double>::const_iterator it = collectedByVolume.begin();
       it != collectedByVolume.end(); ++it)
    std::cout << "  " << it->first << "  " << it->second * norm * 1e6
              << " ueV per event" << std::endl;
  std::cout << "Energy collected: " << collected_eV * 1e6 << " ueV per event"
            << std::endl;
  std::cout << "Initial energy:   " << initialEnergy_eV * 1e6 << " ueV"
            << std::endl;
  std::cout << "Collection efficiency: " << eff * 100. << " %" << std::endl;
  if (energyCollected_eV) *energyCollected_eV = collected_eV;
  if (efficiency) *efficiency = eff;

  TH1D* total = new TH1D("qpCountTotal",
                         Form("QPs in volumes containing \"%s\";time [ns];%s",
                              volumeContains,
                              totalVolume_um3 > 0. ? "QP density [1/#mum^{3}]"
                                                   : "number of QPs"),
                         nBins, tMin_ns, tMax_ns);
  total->SetLineColor(kBlack);
  total->SetLineWidth(3);
  total->SetStats(0);

  std::vector<TH1D*> perVol;
  const int palette[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2,
                         kViolet, kPink+7};
  const int nPalette = sizeof(palette) / sizeof(palette[0]);
  int iColor = 0;
  for (std::map<std::string, std::vector<double> >::const_iterator it = sums.begin();
       it != sums.end(); ++it) {
    TH1D* h = new TH1D(Form("qpCount_%s", it->first.c_str()), it->first.c_str(),
                       nBins, tMin_ns, tMax_ns);
    h->SetStats(0);
    h->SetLineColor(palette[iColor++ % nPalette]);
    h->SetLineWidth(2);
    for (int b = 0; b < nBins; ++b) {
      double n = it->second[b] * norm;
      h->SetBinContent(b + 1, n);
      total->AddBinContent(b + 1, n);
    }
    perVol.push_back(h);
  }
  total->Scale(volumeScale);

  TCanvas* c = new TCanvas("cQPCount", "QP count vs time", 900, 600);
  c->SetLogy(logY);
  double ymax = total->GetMaximum();
  total->SetMaximum(ymax > 0. ? 1.3 * ymax : 1.);
  if (logY) {
    double ymin = 0.;
    for (int b = 1; b <= nBins; ++b) {
      double v = total->GetBinContent(b);
      if (v > 0. && (ymin == 0. || v < ymin)) ymin = v;
    }
    total->SetMinimum(ymin > 0. ? 0.3 * ymin : 1e-3);
  } else {
    total->SetMinimum(0.);
  }
  total->Draw("HIST");
  TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89);
  leg->AddEntry(total, "total", "l");
  if (perVolume) {
    for (size_t i = 0; i < perVol.size(); ++i) {
      perVol[i]->Scale(volumeScale);
      perVol[i]->Draw("HIST SAME");
      leg->AddEntry(perVol[i], perVol[i]->GetTitle(), "l");
    }
  }
  leg->Draw();
  TLatex* text = new TLatex(0.15, 0.85,
    Form("E_{collected} = %.4g #mueV, efficiency = %.3g%%", collected_eV * 1e6,
         eff * 100.));
  text->SetNDC();
  text->SetTextSize(0.035);
  text->Draw();
  if (outFile && outFile[0]) c->SaveAs(outFile);
  return total;
}


// Number of QPs present at the time t_ns in each volume whose name contains
// volumeContains (steps with preT <= t < postT), averaged over the events in
// the file. Prints the counts and returns the total.
double qp_count_at_time(const char* stepFile, double t_ns,
                        const char* volumeContains = "shunt") {
  std::ifstream in(stepFile);
  if (!in.is_open()) {
    std::cerr << "qp_count_at_time: cannot open " << stepFile << std::endl;
    return -1.;
  }

  std::map<std::string, double> counts;
  std::set<std::pair<int,int> > events;
  QPCVStep s;
  while (in.good()) {
    if (!QPCVReadStep(in, s)) {
      if (in.eof()) break;
      in.clear();
      in.ignore(100000, '\n');
      continue;
    }
    events.insert(std::make_pair(s.run, s.event));
    if (s.preVolume.find(volumeContains) == std::string::npos) continue;
    if (s.preT <= t_ns && t_ns < s.postT) counts[s.preVolume] += 1.;
  }

  const double norm = events.empty() ? 1. : 1. / events.size();
  double total = 0.;
  std::cout << "QPs at t = " << t_ns << " ns (" << events.size()
            << " events):" << std::endl;
  for (std::map<std::string, double>::const_iterator it = counts.begin();
       it != counts.end(); ++it) {
    std::cout << "  " << it->first << "  " << it->second * norm << std::endl;
    total += it->second * norm;
  }
  std::cout << "  total  " << total << std::endl;
  return total;
}
