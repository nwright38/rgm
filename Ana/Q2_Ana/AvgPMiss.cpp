// AvgPMiss.cpp
//
// Standalone cross-check: computes the per-event average pMiss within each
// of Main_Figs_Binned.cpp's pMiss bins (bE_pmiss), for e'p and e'pp events
// separately, directly from raw per-event kinematics. This sidesteps the
// Python side's approach in plot_Emiss_pmiss.py (which estimates the same
// average from the pMiss_ep_note/pMiss_epp_note 50-bin histograms) so the
// two can be compared against each other when the E_miss threshold line
// position looks off.
//
// Same physics inputs, nominal cuts, and MC reweighting as
// Main_Figs_Binned.cpp's nominal selection -- no toy systematics here,
// since this is a numeric cross-check, not a plot.
//
// Usage: ./AvgPMiss isMC A outputfile.root inputfiles.hipo

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "reweighter.h"
#include "Corrections.h"

using namespace std;
using namespace clas12;

// ---- Same physics constants as Main_Figs_Binned.cpp ----
auto db = TDatabasePDG::Instance();
const double beam_E = 5.98636;
const double mass_p = db->GetParticle(2212)->Mass();
const double mD = 1.8756;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2 * me;
const double mN = 0.938272;  // matches SRC_Cuts.cpp's mN

// ---- Same pMiss bin edges as Main_Figs_Binned.cpp's bE_pmiss ----
vector<double> bE_pmiss = {0.4, 0.55, 0.7, 0.85, 1.0};

// ---- Same nominal cut thresholds as Main_Figs_Binned.cpp's
// CutVariation::Nominal() ----
const double xB_lower = 1.2;
const double Q2_lower = 1.5;
const double mMiss_lower = 0.65;
const double mMiss_upper = 1.10;
const double kMiss_lower = 0.3;
const double pLead_lower = 1.0;
const double tLead_upper = 37;
const double pRecoil_lower = 0.3;

int findBin(double x, const vector<double>& edges) {
  for (size_t i = 0; i + 1 < edges.size(); i++) {
    if (x >= edges[i] && x < edges[i + 1]) return (int)i;
  }
  return -1;
}

void Usage() {
  cerr << "Usage: ./AvgPMiss isMC A outputfile.root inputfiles.hipo \n";
}

int main(int argc, char** argv) {
  if (argc < 4) {
    Usage();
    return -1;
  }

  int isMC = atoi(argv[1]);
  char* uType = "AV18";
  if (isMC > 0) {
    if (isMC == 1) uType = "AV18";
    if (isMC == 2) uType = "AV4";
    if (isMC == 3) uType = "N2LO10";
    if (isMC == 4) uType = "N2LO12";
    if (isMC == 5) uType = "NV";
    isMC = 1;
  }

  int nucleus_A = atoi(argv[2]);
  TString outFile = argv[3];
  cout << "Output file " << outFile << endl;

  clas12ana clasAna;
  clasAna.printParams();

  clas12root::HipoChain chain;
  for (int k = 4; k < argc; k++) {
    cout << "Input file " << argv[k] << endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12 = chain.GetC12Reader();
  auto& c12 = chain.C12ref();

  int Z = 2, N = 2;
  if (isMC) {
    Z = nucleus_A / 2;
    N = nucleus_A / 2;
  }
  reweighter newWeight(beam_E, Z, N, kelly, uType, .15);

  TLorentzVector beam(0, 0, beam_E, beam_E);
  TLorentzVector nucleus_ptr(0, 0, 0, m_4He);
  TLorentzVector deut_ptr(0, 0, 0, mD);
  TLorentzVector el(0, 0, 0, db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0, 0, 0, db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0, 0, 0, db->GetParticle(2212)->Mass());

  int nBins = (int)bE_pmiss.size() - 1;
  vector<double> sumW_ep(nBins, 0.0), sumWPmiss_ep(nBins, 0.0);
  vector<double> sumW_epp(nBins, 0.0), sumWPmiss_epp(nBins, 0.0);
  vector<long> nEvents_ep(nBins, 0), nEvents_epp(nBins, 0);

  int ctr = 0;
  int max_ev = 10000000;
  while (chain.Next() && ctr < max_ev) {
    ctr++;
    if (ctr % 100000 == 0) cout << "Event " << ctr << " of " << max_ev << "\r" << flush;

    double wep = 1, wepp = 1, original_weight = 1;
    if (isMC) {
      original_weight = c12->mcevent()->getWeight();
      wep = original_weight * newWeight.get_weight_ep(c12->mcparts());
      wepp = original_weight * newWeight.get_weight_epp(c12->mcparts());
    }

    clasAna.Run(c12);
    auto electrons = clasAna.getByPid(11);
    if (electrons.size() != 1) continue;

    GetLorentzVector_Corrected(el, electrons[0], isMC);
    TLorentzVector q = beam - el;
    double Q2 = -q.M2();
    double xB = Q2 / (2 * mass_p * (beam.E() - el.E()));

    clasAna.getLeadRecoilSRC(beam, deut_ptr, el);
    auto lead = clasAna.getLeadSRC();
    auto recoil = clasAna.getRecoilSRC();
    if (lead.size() != 1) continue;
    GetLorentzVector_Corrected(lead_ptr, lead[0], isMC);
    if (lead[0]->getRegion() != FD) continue;

    TLorentzVector miss = q + deut_ptr - lead_ptr;
    double mmiss = sqrt(miss.M2());
    double pmiss = miss.P();
    double mom_lead = lead_ptr.P();
    double theta_lead = lead_ptr.Theta() * 180 / M_PI;

    TLorentzVector miss_LC = lead_ptr - q;
    TVector3 u = q.Vect().Unit();
    double pmm = miss_LC.E() - miss_LC.Vect().Dot(u);
    double pmp = miss_LC.Vect().Perp(u);
    double kmiss = sqrt(mN * mN * ((pmp * pmp + mN * mN) / (pmm * (2 * mN - pmm))) - mN * mN);

    bool passep = true;
    if (xB < xB_lower) passep = false;
    if (Q2 < Q2_lower) passep = false;
    if (Q2 > 5.0) passep = false;
    if (mmiss < mMiss_lower) passep = false;
    if (mmiss > mMiss_upper) passep = false;
    if (kmiss < kMiss_lower) passep = false;
    if (mom_lead < pLead_lower) passep = false;
    if (theta_lead > tLead_upper) passep = false;
    if (!passep) continue;

    int b = findBin(pmiss, bE_pmiss);
    if (b >= 0) {
      sumW_ep[b] += wep;
      sumWPmiss_ep[b] += wep * pmiss;
      nEvents_ep[b]++;
    }

    bool rec = (recoil.size() == 1);
    if (!rec) continue;
    GetLorentzVector_Corrected(recoil_ptr, recoil[0], isMC);
    double pRecoil = recoil_ptr.P();
    if (pRecoil < pRecoil_lower) continue;

    if (b >= 0) {
      sumW_epp[b] += wepp;
      sumWPmiss_epp[b] += wepp * pmiss;
      nEvents_epp[b]++;
    }
  }
  cout << endl;

  cout << "\npMiss bin averages (e,e'p):\n";
  for (int i = 0; i < nBins; i++) {
    double avg = sumW_ep[i] > 0 ? sumWPmiss_ep[i] / sumW_ep[i] : 0.0;
    cout << "  [" << bE_pmiss[i] << ", " << bE_pmiss[i + 1] << "): "
         << "avg pMiss = " << avg << " GeV, N = " << nEvents_ep[i] << "\n";
  }

  cout << "\npMiss bin averages (e,e'pp):\n";
  for (int i = 0; i < nBins; i++) {
    double avg = sumW_epp[i] > 0 ? sumWPmiss_epp[i] / sumW_epp[i] : 0.0;
    cout << "  [" << bE_pmiss[i] << ", " << bE_pmiss[i + 1] << "): "
         << "avg pMiss = " << avg << " GeV, N = " << nEvents_epp[i] << "\n";
  }

  TFile* f = new TFile(outFile, "RECREATE");
  f->cd();
  TTree* t = new TTree("avgPMiss", "Average pMiss per bin");
  string sel;
  int binIdx;
  double binLo, binHi, avgPmiss;
  long nEv;
  t->Branch("selection", &sel);
  t->Branch("bin", &binIdx, "bin/I");
  t->Branch("lo", &binLo, "lo/D");
  t->Branch("hi", &binHi, "hi/D");
  t->Branch("avg_pmiss", &avgPmiss, "avg_pmiss/D");
  t->Branch("n_events", &nEv, "n_events/L");

  for (int i = 0; i < nBins; i++) {
    sel = "ep";
    binIdx = i;
    binLo = bE_pmiss[i];
    binHi = bE_pmiss[i + 1];
    avgPmiss = sumW_ep[i] > 0 ? sumWPmiss_ep[i] / sumW_ep[i] : 0.0;
    nEv = nEvents_ep[i];
    t->Fill();
  }
  for (int i = 0; i < nBins; i++) {
    sel = "epp";
    binIdx = i;
    binLo = bE_pmiss[i];
    binHi = bE_pmiss[i + 1];
    avgPmiss = sumW_epp[i] > 0 ? sumWPmiss_epp[i] / sumW_epp[i] : 0.0;
    nEv = nEvents_epp[i];
    t->Fill();
  }
  t->Write();
  f->Close();

  return 0;
}
