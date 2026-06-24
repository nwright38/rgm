// Generic Q2/pMiss/kMiss binned-plot framework. Same physics inputs as
// Main_Figs_Sys_Err.cpp (which is left untouched), but the histogram
// bookkeeping is data-driven (see buildFillTasks()) instead of one hand
// written struct field + loop block per plot, and the output is a flat
// ROOT TTree (diffTable / integratedTable) instead of one-off TGraphErrors.
//
// Two known issues found while porting the original file's logic were fixed
// here by construction rather than silently inherited -- see the comments
// next to CutVariation::Randomized() and buildFillTasks()'s "_ep_" tasks
// below for what changed and why.

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "many_plots.h"
#include "reweighter.h"
#include "TGraphErrors.h"
#include "Corrections.h"
#include "TRandom3.h"

#include "BinnedHistStore.h"

using namespace std;
using namespace clas12;

// ---- Same physics constants as Main_Figs_Sys_Err.cpp ----
const double c = 29.9792458;
auto db = TDatabasePDG::Instance();
const double beam_E = 5.98636;
const double mass_n = db->GetParticle(2112)->Mass();
const double mass_p = db->GetParticle(2212)->Mass();
const double mass_pi = db->GetParticle(-211)->Mass();
const double mD = 1.8756;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2 * me;

// ---- Bin edges: single source of truth, used to build Axis objects below ----
vector<double> bE_Q2 = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
vector<double> bE_pmiss = {0.4, 0.55, 0.7, 0.85, 1.0};
vector<double> bE_kmiss = {0.3, 0.45, 0.6, 0.75, 0.9};
vector<double> bE_pmiss_long = {0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.85, 0.95, 1.05, 1.2};
vector<double> bE_kmiss_long = {0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};

Axis Q2Axis{"Q2", bE_Q2};
Axis pMissAxis{"pMiss", bE_pmiss};
Axis kMissAxis{"kMiss", bE_kmiss};

// ---- Per-event physics quantities, computed once per event ----
struct EventKinematics {
  double qSq = 0, x = 0, pL = 0, tL = 0, pR = 0;
  double px = 0, py = 0, pz = 0;
  double mM = 0, pM = 0, kM = 0;
  double E0 = 0, E1 = 0, E2 = 0;
  bool passep = false, passepp = false;
};

// Same kinematic calculation as runEvent() in Main_Figs_Sys_Err.cpp, just
// returning one struct instead of writing through 14 reference parameters
// (the original 19-parameter-by-reference signature is the easiest place to
// accidentally pass an argument in the wrong slot).
EventKinematics computeEventKinematics(const std::unique_ptr<clas12::clas12reader>& c12,
                                        clas12ana& clasAna, bool isMC, int& ctr) {
  EventKinematics ek;

  TLorentzVector beam(0, 0, beam_E, beam_E);
  TLorentzVector nucleus_ptr(0, 0, 0, m_4He);
  TLorentzVector deut_ptr(0, 0, 0, mD);
  TLorentzVector el(0, 0, 0, db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0, 0, 0, db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0, 0, 0, db->GetParticle(2212)->Mass());

  ctr++;
  if ((ctr % 100000) == 0) cerr << "\n" << ctr / 100000 << " hundred thousand completed";
  if ((ctr % 10000) == 0) cerr << ".";

  clasAna.Run(c12);
  auto electrons = clasAna.getByPid(11);
  if (electrons.size() != 1) return ek;

  GetLorentzVector_Corrected(el, electrons[0], isMC);
  TLorentzVector q = beam - el;
  double Q2 = -q.M2();
  double omega = q.E();
  double xB = Q2 / (2 * mass_p * (beam.E() - el.E()));

  clasAna.getLeadRecoilSRC(beam, deut_ptr, el);
  auto lead = clasAna.getLeadSRC();
  auto recoil = clasAna.getRecoilSRC();

  if (lead.size() != 1) return ek;
  GetLorentzVector_Corrected(lead_ptr, lead[0], isMC);
  if (lead[0]->getRegion() != FD) return ek;

  TLorentzVector miss = q + deut_ptr - lead_ptr;
  double mmiss2 = miss.M2();
  double mmiss = sqrt(mmiss2);
  TVector3 miss_neg = -miss.Vect();
  double mom_miss = miss.P();

  double E0miss = sqrt(mom_miss * mom_miss + mN * mN) - mN;
  double mom_lead = lead_ptr.P();
  double theta_lead = lead_ptr.Theta() * 180 / M_PI;
  double phi_lead = lead_ptr.Phi() * 180 / M_PI;
  double EP = lead_ptr.E();
  double EB = omega + nucleus_ptr.M() - EP;
  double TB = EB - sqrt(EB * EB - mom_miss * mom_miss);
  double TP = EP - sqrt(EP * EP - mom_lead * mom_lead);
  double E1miss = omega - TP - TB;

  TLorentzVector miss_LC = lead_ptr - q;
  TVector3 u = q.Vect().Unit();
  double pmm = miss_LC.E() - miss_LC.Vect().Dot(u);
  double pmp = miss_LC.Vect().Perp(u);
  double pmiss = miss.P();
  double kmiss = sqrt(mN * mN * ((pmp * pmp + mN * mN) / (pmm * (2 * mN - pmm))) - mN * mN);

  bool rec = (recoil.size() == 1);

  ek.passep = true;
  ek.qSq = Q2;
  ek.x = xB;
  ek.pL = mom_lead;
  ek.tL = theta_lead;
  ek.mM = mmiss;
  ek.pM = pmiss;
  ek.kM = kmiss;
  ek.E0 = E0miss;
  ek.E1 = E1miss;

  if (!rec) return ek;
  GetLorentzVector_Corrected(recoil_ptr, recoil[0], isMC);

  double TP2 = recoil_ptr.E() - recoil_ptr.M();
  TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr;
  double TB2 = miss_Am2.E() - miss_Am2.M();
  double E2miss = q.E() - TP - TP2 - TB2;

  TVector3 v_rec = recoil_ptr.Vect();
  TVector3 v_cm = miss_neg + v_rec;
  TVector3 vz = miss_neg.Unit();
  TVector3 vy = miss_neg.Cross(q.Vect()).Unit();
  TVector3 vx = vz.Cross(vy).Unit();

  ek.passepp = true;
  ek.pR = v_rec.Mag();
  ek.px = v_cm.Dot(vx);
  ek.py = v_cm.Dot(vy);
  ek.pz = v_cm.Dot(vz);
  ek.E2 = E2miss;

  return ek;
}

// Same cut thresholds as CutRandom() in the original, but restructured so a
// toy's randomized thresholds are generated ONCE (before the event loop) and
// then applied identically to every event for that toy.
//
// Original bug this fixes: the old CutRandom(..., randomize=true) called
// `new TRandom3(0)` *inside the per-event, per-toy loop* (100x per event).
// TRandom3(0) seeds from the clock, so that line (a) leaked 100 TRandom3
// objects per event, and (b) re-rolled the toy's cut thresholds on every
// single event instead of holding them fixed across the dataset -- so the
// resulting "systematic error" was measuring per-event threshold noise, not
// the effect of a single coherent alternate cut choice (which is what a
// cut-variation systematic is supposed to estimate).
struct CutVariation {
  double xB_lower = 1.2;
  double Q2_lower = 1.5;
  double mMiss_lower = 0.65;
  double mMiss_upper = 1.10;
  double kMiss_lower = 0.3;
  double pLead_lower = 1.0;
  double tLead_upper = 37;
  double pRecoil_lower = 0.3;

  static CutVariation Nominal() { return CutVariation{}; }

  static CutVariation Randomized(TRandom3& rng) {
    CutVariation v;
    v.xB_lower += rng.Gaus(0.0, 0.01);
    v.Q2_lower += rng.Gaus(0.0, 0.01);
    v.mMiss_lower += rng.Gaus(0.0, 0.03);
    v.mMiss_upper += rng.Gaus(0.0, 0.03);
    v.kMiss_lower += rng.Gaus(0.0, 0.03);
    v.pLead_lower += rng.Gaus(0.0, 0.03);
    v.tLead_upper += rng.Gaus(0.0, 1.0);
    v.pRecoil_lower += rng.Gaus(0.0, 0.045);
    return v;
  }

  void apply(const EventKinematics& ek, bool& passep, bool& passepp) const {
    passep = ek.passep;
    passepp = ek.passepp;
    if (!passep) return;

    if (ek.x < xB_lower) passep = false;
    if (ek.qSq < Q2_lower) passep = false;
    if (ek.qSq > 5.0) passep = false;
    if (ek.mM < mMiss_lower) passep = false;
    if (ek.mM > mMiss_upper) passep = false;
    if (ek.kM < kMiss_lower) passep = false;
    if (ek.pL < pLead_lower) passep = false;
    if (ek.tL > tLead_upper) passep = false;
    if (!passep) return;

    if (ek.pR < pRecoil_lower) passepp = false;
  }
};

// One FillTask per *category* of plot (e.g. "all E1miss-in-pMiss-bins,
// split by Q2"), not one per individual histogram -- the cartesian product
// over each task's selector axes is what used to be a hand-written
// HistGroup struct field plus a [4] or [4][7] loop block.
vector<FillTask<EventKinematics>> buildFillTasks() {
  vector<FillTask<EventKinematics>> tasks;

  auto passEP = [](const EventKinematics& ek) { return ek.passep; };
  auto passEPP = [](const EventKinematics& ek) { return ek.passepp; };

  // Plain Q2/pMiss/kMiss yields, no selector binning.
  tasks.push_back({"Q2_ep", Selection::EP, passEP, {}, {},
                    [](const EventKinematics& ek) { return ek.qSq; }, 50, 1.5, 5.0, {}});
  tasks.push_back({"Q2_epp", Selection::EPP, passEPP, {}, {},
                    [](const EventKinematics& ek) { return ek.qSq; }, 50, 1.5, 5.0, {}});
  tasks.push_back({"pMiss_ep", Selection::EP, passEP, {}, {},
                    [](const EventKinematics& ek) { return ek.pM; }, 0, 0, 0, bE_pmiss_long});
  tasks.push_back({"pMiss_epp", Selection::EPP, passEPP, {}, {},
                    [](const EventKinematics& ek) { return ek.pM; }, 0, 0, 0, bE_pmiss_long});
  tasks.push_back({"kMiss_ep", Selection::EP, passEP, {}, {},
                    [](const EventKinematics& ek) { return ek.kM; }, 0, 0, 0, bE_kmiss_long});
  tasks.push_back({"kMiss_epp", Selection::EPP, passEPP, {}, {},
                    [](const EventKinematics& ek) { return ek.kM; }, 0, 0, 0, bE_kmiss_long});

  // Q2 yield, selected (rebinned) by its own Q2 bin -- matches the original
  // h_Q2_ep_SRC_Q2 / h_Q2_epp_SRC_Q2 (used to get each Q2 bin's mean Q2).
  tasks.push_back({"Q2_ep_SRC_Q2", Selection::EP, passEP, {&Q2Axis},
                    {[](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.qSq; }, 100, 1.5, 5.0, {}});
  tasks.push_back({"Q2_epp_SRC_Q2", Selection::EPP, passEPP, {&Q2Axis},
                    {[](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.qSq; }, 100, 1.5, 5.0, {}});

  // Q2 yield, selected by pMiss / kMiss bin.
  tasks.push_back({"Q2_ep_SRC_pmiss", Selection::EP, passEP, {&pMissAxis},
                    {[](const EventKinematics& ek) { return ek.pM; }},
                    [](const EventKinematics& ek) { return ek.qSq; }, 0, 0, 0, bE_Q2});
  tasks.push_back({"Q2_epp_SRC_pmiss", Selection::EPP, passEPP, {&pMissAxis},
                    {[](const EventKinematics& ek) { return ek.pM; }},
                    [](const EventKinematics& ek) { return ek.qSq; }, 0, 0, 0, bE_Q2});
  tasks.push_back({"Q2_ep_SRC_kmiss", Selection::EP, passEP, {&kMissAxis},
                    {[](const EventKinematics& ek) { return ek.kM; }},
                    [](const EventKinematics& ek) { return ek.qSq; }, 0, 0, 0, bE_Q2});
  tasks.push_back({"Q2_epp_SRC_kmiss", Selection::EPP, passEPP, {&kMissAxis},
                    {[](const EventKinematics& ek) { return ek.kM; }},
                    [](const EventKinematics& ek) { return ek.qSq; }, 0, 0, 0, bE_Q2});

  // E_miss shapes, selected by pMiss bin.
  tasks.push_back({"E0miss_ep_SRC_pmiss", Selection::EP, passEP, {&pMissAxis},
                    {[](const EventKinematics& ek) { return ek.pM; }},
                    [](const EventKinematics& ek) { return ek.E0; }, 30, -0.1, 0.85, {}});
  tasks.push_back({"E1miss_ep_SRC_pmiss", Selection::EP, passEP, {&pMissAxis},
                    {[](const EventKinematics& ek) { return ek.pM; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.25, 0.55, {}});
  tasks.push_back({"E1miss_epp_SRC_pmiss", Selection::EPP, passEPP, {&pMissAxis},
                    {[](const EventKinematics& ek) { return ek.pM; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.2, 0.5, {}});
  tasks.push_back({"E2miss_epp_SRC_pmiss", Selection::EPP, passEPP, {&pMissAxis},
                    {[](const EventKinematics& ek) { return ek.pM; }},
                    [](const EventKinematics& ek) { return ek.E2; }, 15, -0.25, 0.5, {}});

  // E_miss shapes, selected by kMiss bin.
  tasks.push_back({"E0miss_ep_SRC_kmiss", Selection::EP, passEP, {&kMissAxis},
                    {[](const EventKinematics& ek) { return ek.kM; }},
                    [](const EventKinematics& ek) { return ek.E0; }, 30, -0.1, 0.85, {}});
  tasks.push_back({"E1miss_ep_SRC_kmiss", Selection::EP, passEP, {&kMissAxis},
                    {[](const EventKinematics& ek) { return ek.kM; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.25, 0.55, {}});
  tasks.push_back({"E1miss_epp_SRC_kmiss", Selection::EPP, passEPP, {&kMissAxis},
                    {[](const EventKinematics& ek) { return ek.kM; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.2, 0.5, {}});
  tasks.push_back({"E2miss_epp_SRC_kmiss", Selection::EPP, passEPP, {&kMissAxis},
                    {[](const EventKinematics& ek) { return ek.kM; }},
                    [](const EventKinematics& ek) { return ek.E2; }, 15, -0.25, 0.5, {}});

  // E_miss shapes, selected by pMiss bin AND Q2 bin together.
  //
  // NOTE on a discrepancy found vs. the original file: in
  // Main_Figs_Sys_Err.cpp's fillUpHistGroup(), the "_ep_"-named histograms
  // h_Q2_ep_SRC_Q2, h_E1miss_ep_SRC_pmiss_Q2 and h_E1miss_ep_SRC_kmiss_Q2 are
  // filled with the `wepp` weight, not `wep` -- every sibling "_ep_"
  // histogram elsewhere in that function uses `wep`. That looks like an
  // unintentional copy/paste slip (it's not flagged or commented as
  // intentional anywhere) rather than a deliberate choice, so this new file
  // uses `wep` (passed automatically here via Selection::EP) for these
  // tasks instead of reproducing it. Flagging this explicitly: if there was
  // a reason for the original wepp weighting, tell me and I'll switch these
  // back to match.
  tasks.push_back({"E1miss_ep_SRC_pmiss_Q2", Selection::EP, passEP, {&pMissAxis, &Q2Axis},
                    {[](const EventKinematics& ek) { return ek.pM; },
                     [](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.25, 0.55, {}});
  tasks.push_back({"E1miss_epp_SRC_pmiss_Q2", Selection::EPP, passEPP, {&pMissAxis, &Q2Axis},
                    {[](const EventKinematics& ek) { return ek.pM; },
                     [](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.2, 0.5, {}});
  tasks.push_back({"E2miss_epp_SRC_pmiss_Q2", Selection::EPP, passEPP, {&pMissAxis, &Q2Axis},
                    {[](const EventKinematics& ek) { return ek.pM; },
                     [](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.E2; }, 15, -0.25, 0.5, {}});

  tasks.push_back({"E1miss_ep_SRC_kmiss_Q2", Selection::EP, passEP, {&kMissAxis, &Q2Axis},
                    {[](const EventKinematics& ek) { return ek.kM; },
                     [](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.25, 0.55, {}});
  tasks.push_back({"E1miss_epp_SRC_kmiss_Q2", Selection::EPP, passEPP, {&kMissAxis, &Q2Axis},
                    {[](const EventKinematics& ek) { return ek.kM; },
                     [](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.E1; }, 15, -0.2, 0.5, {}});
  tasks.push_back({"E2miss_epp_SRC_kmiss_Q2", Selection::EPP, passEPP, {&kMissAxis, &Q2Axis},
                    {[](const EventKinematics& ek) { return ek.kM; },
                     [](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.E2; }, 15, -0.25, 0.5, {}});

  // Center-of-mass momentum components, no selector and selected by Q2 bin.
  tasks.push_back({"pcmx_epp", Selection::EPP, passEPP, {}, {},
                    [](const EventKinematics& ek) { return ek.px; }, 15, -0.75, 0.75, {}});
  tasks.push_back({"pcmy_epp", Selection::EPP, passEPP, {}, {},
                    [](const EventKinematics& ek) { return ek.py; }, 15, -0.75, 0.75, {}});
  tasks.push_back({"pcmz_epp", Selection::EPP, passEPP, {}, {},
                    [](const EventKinematics& ek) { return ek.pz; }, 15, -0.75, 0.75, {}});

  tasks.push_back({"pcmx_epp_SRC_Q2", Selection::EPP, passEPP, {&Q2Axis},
                    {[](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.px; }, 15, -0.75, 0.75, {}});
  tasks.push_back({"pcmy_epp_SRC_Q2", Selection::EPP, passEPP, {&Q2Axis},
                    {[](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.py; }, 15, -0.75, 0.75, {}});
  // Original uses a different (asymmetric) range here than plain pcmz_epp --
  // kept faithfully, since p_cm_z is along the boosted q-direction and an
  // asymmetric range there looks like a deliberate physics choice.
  tasks.push_back({"pcmz_epp_SRC_Q2", Selection::EPP, passEPP, {&Q2Axis},
                    {[](const EventKinematics& ek) { return ek.qSq; }},
                    [](const EventKinematics& ek) { return ek.pz; }, 15, -0.5, 1.0, {}});

  return tasks;
}

void Usage() {
  std::cerr << "Usage: ./Main_Figs_Binned isMC A outputfile.root inputfiles.hipo \n\n\n";
}

int main(int argc, char** argv) {
  if (argc < 4) {
    Usage();
    return -1;
  }

  // Detach every histogram we create from ROOT's directory bookkeeping, so
  // booking 101 stores (nominal + 100 toys) with identically-named
  // histograms in memory never collides with ROOT's "Replacing existing
  // TH1D" duplicate-name behavior. We place histograms into TDirectories
  // explicitly (HistStore::writeAll) only at the very end, when writing out.
  TH1::AddDirectory(kFALSE);

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

  const int nToys = 100;
  TRandom3 toyRng(1234);  // fixed seed: deterministic, reproducible across runs
  CutVariation nominalCut = CutVariation::Nominal();
  vector<CutVariation> toyCuts(nToys);
  vector<reweighter> toyWeighters;
  toyWeighters.reserve(nToys);
  for (int i = 0; i < nToys; i++) {
    toyCuts[i] = CutVariation::Randomized(toyRng);
    reweighter w(beam_E, Z, N, kelly, uType, .15);
    w.randomize_Config();
    toyWeighters.push_back(w);
  }

  vector<FillTask<EventKinematics>> tasks = buildFillTasks();

  HistStore<EventKinematics> nominalStore("nominal");
  nominalStore.book(tasks);
  vector<HistStore<EventKinematics>> toyStores;
  toyStores.reserve(nToys);
  for (int i = 0; i < nToys; i++) {
    char label[32];
    snprintf(label, sizeof(label), "toy_%03d", i);
    toyStores.emplace_back(label);
    toyStores.back().book(tasks);
  }

  int ctr = 0;
  int max_ev = 10000000;
  while (chain.Next() && ctr < max_ev) {
    if (ctr % 1000 == 0) cout << "Event " << ctr << " of " << max_ev << ". " << '\t' << (double)ctr / max_ev * 100 << "%\r" << flush;

    double wep = 1, wepp = 1, original_weight = 1;
    if (isMC) {
      original_weight = c12->mcevent()->getWeight();
      wep = original_weight * newWeight.get_weight_ep(c12->mcparts());
      wepp = original_weight * newWeight.get_weight_epp(c12->mcparts());
    }

    EventKinematics ek = computeEventKinematics(c12, clasAna, false, ctr);

    bool passep_nom, passepp_nom;
    nominalCut.apply(ek, passep_nom, passepp_nom);
    EventKinematics ekNom = ek;
    ekNom.passep = passep_nom;
    ekNom.passepp = passepp_nom;
    nominalStore.fill(tasks, ekNom, wep, wepp);

    for (int i = 0; i < nToys; i++) {
      double wep_sys = wep, wepp_sys = wepp;
      if (isMC) {
        wep_sys = original_weight * toyWeighters[i].get_weight_ep(c12->mcparts());
        wepp_sys = original_weight * toyWeighters[i].get_weight_epp(c12->mcparts());
      }
      bool passep_i, passepp_i;
      toyCuts[i].apply(ek, passep_i, passepp_i);
      EventKinematics ekToy = ek;
      ekToy.passep = passep_i;
      ekToy.passepp = passepp_i;
      toyStores[i].fill(tasks, ekToy, wep_sys, wepp_sys);
    }
  }

  // ---- Output ----
  TFile* f = new TFile(outFile, "RECREATE");
  f->cd();

  TTree* diffTree = new TTree("diffTable", "Differential bin table");
  string d_task, d_sel;
  vector<int> d_axisBin;
  vector<double> d_axisCenter;
  int d_valueBin;
  double d_valueCenter, d_count, d_statErr, d_sysErr;
  diffTree->Branch("task_name", &d_task);
  diffTree->Branch("selection", &d_sel);
  diffTree->Branch("axis_bin", &d_axisBin);
  diffTree->Branch("axis_center", &d_axisCenter);
  diffTree->Branch("value_bin", &d_valueBin, "value_bin/I");
  diffTree->Branch("value_center", &d_valueCenter, "value_center/D");
  diffTree->Branch("count", &d_count, "count/D");
  diffTree->Branch("stat_error", &d_statErr, "stat_error/D");
  diffTree->Branch("sys_error", &d_sysErr, "sys_error/D");

  TTree* intTree = new TTree("integratedTable", "Integrated (collapsed-axis) table");
  string i_task, i_sel, i_pattern;
  vector<int> i_axisBin;
  vector<double> i_axisCenter;
  double i_count, i_statErr, i_sysErr;
  intTree->Branch("task_name", &i_task);
  intTree->Branch("selection", &i_sel);
  intTree->Branch("pattern", &i_pattern);
  intTree->Branch("axis_bin", &i_axisBin);
  intTree->Branch("axis_center", &i_axisCenter);
  intTree->Branch("count", &i_count, "count/D");
  intTree->Branch("stat_error", &i_statErr, "stat_error/D");
  intTree->Branch("sys_error", &i_sysErr, "sys_error/D");

  for (size_t t = 0; t < tasks.size(); t++) {
    auto diffRows = buildDiffRows(tasks[t], (int)t, nominalStore, toyStores);
    for (auto& r : diffRows) {
      d_task = r.task_name;
      d_sel = r.selection;
      d_axisBin = r.axis_bin;
      d_axisCenter = r.axis_center;
      d_valueBin = r.value_bin;
      d_valueCenter = r.value_center;
      d_count = r.count;
      d_statErr = r.stat_error;
      d_sysErr = r.sys_error;
      diffTree->Fill();
    }

    // Integrated quantities: collapse each selector axis individually (e.g.
    // "integrate out pMiss within each Q2 bin") plus collapse all of them
    // (grand total). Add more entries to collapsePatterns to support more
    // integration patterns without restructuring anything.
    // "collapse:none" always runs first: it sums each selector bin's
    // histogram over its own full range with no axes collapsed (e.g. total
    // yield per pMiss bin, before integrating across pMiss bins too) --
    // useful on its own, and free since buildIntegratedRows always computes
    // a full-range integral per histogram regardless of pattern.
    vector<pair<vector<int>, string>> collapsePatterns;
    collapsePatterns.push_back({{}, "collapse:none"});
    for (size_t a = 0; a < tasks[t].binAxes.size(); a++) {
      collapsePatterns.push_back({{(int)a}, "collapse:" + tasks[t].binAxes[a]->name});
    }
    if (tasks[t].binAxes.size() > 1) {
      vector<int> allAxes;
      for (size_t a = 0; a < tasks[t].binAxes.size(); a++) allAxes.push_back((int)a);
      collapsePatterns.push_back({allAxes, "collapse:all"});
    }

    for (auto& pat : collapsePatterns) {
      auto intRows = buildIntegratedRows(tasks[t], (int)t, nominalStore, toyStores,
                                          pat.first, pat.second);
      for (auto& r : intRows) {
        i_task = r.task_name;
        i_sel = r.selection;
        i_pattern = r.pattern;
        i_axisBin = r.axis_bin;
        i_axisCenter = r.axis_center;
        i_count = r.count;
        i_statErr = r.stat_error;
        i_sysErr = r.sys_error;
        intTree->Fill();
      }
    }
  }

  intTree->Write();

  // Ratio tasks (epp/ep), same toy-spread-of-the-ratio approach as the
  // original getGraphWithError() but generalized and non-destructive.
  auto findTaskIdx = [&](const string& name) -> int {
    for (size_t t = 0; t < tasks.size(); t++) if (tasks[t].name == name) return (int)t;
    return -1;
  };
  vector<pair<string, string>> ratioSpecs = {
      {"pMiss_epp", "pMiss_ep"},
      {"kMiss_epp", "kMiss_ep"},
      {"Q2_epp_SRC_pmiss", "Q2_ep_SRC_pmiss"},
      {"Q2_epp_SRC_kmiss", "Q2_ep_SRC_kmiss"},
  };
  for (auto& spec : ratioSpecs) {
    int numIdx = findTaskIdx(spec.first);
    int denIdx = findTaskIdx(spec.second);
    if (numIdx < 0 || denIdx < 0) continue;
    string ratioName = spec.first + "_over_" + spec.second;
    auto rows = buildRatioDiffRows(tasks[numIdx], numIdx, denIdx, ratioName, nominalStore, toyStores);
    for (auto& r : rows) {
      d_task = r.task_name;
      d_sel = r.selection;
      d_axisBin = r.axis_bin;
      d_axisCenter = r.axis_center;
      d_valueBin = r.value_bin;
      d_valueCenter = r.value_center;
      d_count = r.count;
      d_statErr = r.stat_error;
      d_sysErr = r.sys_error;
      diffTree->Fill();
    }
  }
  diffTree->Write();

  // Underlying histograms (nominal + every toy), so non-summable quantities
  // (e.g. a Gaussian fit to an E_miss shape) can still be recomputed later
  // at any merged binning by hadd-ing the relevant histograms together.
  TDirectory* histDir = f->mkdir("hists");
  histDir->cd();
  nominalStore.writeAll();
  for (auto& toy : toyStores) toy.writeAll();
  f->cd();

  f->Close();
  return 0;
}
