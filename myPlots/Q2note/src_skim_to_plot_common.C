#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TRandom3.h>
#include <TTree.h>

namespace src_skim_plot {

enum SelectionType { SEL_EP, SEL_EPP };
enum LeadMode { LEAD_FD, LEAD_CD, LEAD_BOTH };

struct AnalysisOptions {
  LeadMode leadMode = LEAD_FD;
};

struct CutConfig {
  bool useXBMin = true;
  double xB_min = 1.2;
  double xB_min_sigma = 0.01;

  bool useQ2Min = true;
  double Q2_min = 1.5;
  double Q2_min_sigma = 0.01;

  bool useQ2Max = true;
  double Q2_max = 5.0;
  double Q2_max_sigma = 0.0;

  bool useMMissMin = true;
  double mMiss_min = 0.65;
  double mMiss_min_sigma = 0.03;
  bool useMMissMax = true;
  double mMiss_max = 1.10;
  double mMiss_max_sigma = 0.03;

  bool useKMissMin = true;
  double kMiss_min = 0.3;
  double kMiss_min_sigma = 0.03;

  bool useKMissMax = true;
  double kMiss_max = 1.0;
  double kMiss_max_sigma = 0.03;

  bool useLeadPMin = true;
  double pLead_min = 1.0;
  double pLead_min_sigma = 0.03;

  bool usePoverQMax = true;
  double pOverQ_max = 0.95;
  double pOverQ_max_sigma = 0.03;

  bool useLeadThetaByRegion = true;
  double tLead_FD_max = 37.0;
  double tLead_FD_max_sigma = 1.0;
  double tLead_CD_min = 45.0;
  double tLead_CD_min_sigma = 1.0;

  bool useRecoilPMin = true;
  double pRecoil_min = 0.3;
  double pRecoil_min_sigma = 0.045;

};

CutConfig DefaultCutConfig() {
  return CutConfig();
}

struct EventRow {
  double qSq = 0.0;
  double x = 0.0;
  double pL = 0.0;
  double tL = 0.0;
  double phiL = 0.0;
  double pR = 0.0;
  double tR = 0.0;
  double phiR = 0.0;
  double pE = 0.0;
  double tE = 0.0;
  double phiE = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double mM = 0.0;
  double pM = 0.0;
  double kM = 0.0;
  double E0 = 0.0;
  double E1 = 0.0;
  double E1adj = 0.0;
  double E2 = 0.0;
  double qMag = 0.0;
  double thPMissQ = 0.0;
  double thPLeadQ = 0.0;
  double thPRecPMiss = 0.0;
  double thPRecQ = 0.0;
  double thPLeadPrec = 0.0;
  double pRel = 0.0;
  double pCM = 0.0;
  int leadRegion = 0;
  bool passep = false;
  bool passepp = false;
};

struct CutVariation {
  // Realized cut values for one variation (nominal or one toy).
  // Edit cut definitions and toy spreads in CutConfig only.
  bool useXBMin = true;
  double xB_lower = 1.2;

  bool useQ2Min = true;
  double Q2_lower = 1.5;

  bool useQ2Max = true;
  double Q2_upper = 5.0;

  bool useMMissMin = true;
  double mMiss_lower = 0.65;
  bool useMMissMax = true;
  double mMiss_upper = 1.10;

  bool useKMissMin = true;
  double kMiss_lower = 0.3;
  bool useKMissMax = true;
  double kMiss_upper = 1.0;

  bool useLeadPMin = true;
  double pLead_lower = 1.0;

  bool usePoverQMax = true;
  double pOverQ_upper = 0.95;

  bool useLeadThetaByRegion = true;
  double tLead_upper = 37.0;
  double tLead_lower = 45.0;

  bool useRecoilPMin = true;
  double pRecoil_lower = 0.3;

  static CutVariation Nominal(const CutConfig& cfg) {
    CutVariation v;
    v.useXBMin = cfg.useXBMin;
    v.xB_lower = cfg.xB_min;
    v.useQ2Min = cfg.useQ2Min;
    v.Q2_lower = cfg.Q2_min;
    v.useQ2Max = cfg.useQ2Max;
    v.Q2_upper = cfg.Q2_max;
    v.useMMissMin = cfg.useMMissMin;
    v.mMiss_lower = cfg.mMiss_min;
    v.useMMissMax = cfg.useMMissMax;
    v.mMiss_upper = cfg.mMiss_max;
    v.useKMissMin = cfg.useKMissMin;
    v.kMiss_lower = cfg.kMiss_min;
    v.useKMissMax = cfg.useKMissMax;
    v.kMiss_upper = cfg.kMiss_max;
    v.useLeadPMin = cfg.useLeadPMin;
    v.pLead_lower = cfg.pLead_min;
    v.usePoverQMax = cfg.usePoverQMax;
    v.pOverQ_upper = cfg.pOverQ_max;
    v.useLeadThetaByRegion = cfg.useLeadThetaByRegion;
    v.tLead_upper = cfg.tLead_FD_max;
    v.tLead_lower = cfg.tLead_CD_min;
    v.useRecoilPMin = cfg.useRecoilPMin;
    v.pRecoil_lower = cfg.pRecoil_min;
    return v;
  }

  static CutVariation Randomized(const CutConfig& cfg,
                                 unsigned int seed,
                                 int toyIndex,
                                 bool varyAngularCut) {
    // Stable per-toy randomization with deterministic seed.
    TRandom3 rng(seed + 7919U * (unsigned int)(toyIndex + 1));
    CutVariation v = Nominal(cfg);

    auto smear = [&](double sigma) {
      if (sigma <= 0.0) return 0.0;
      return rng.Gaus(0.0, sigma);
    };

    v.xB_lower += smear(cfg.xB_min_sigma);
    v.Q2_lower += smear(cfg.Q2_min_sigma);
    v.Q2_upper += smear(cfg.Q2_max_sigma);
    v.mMiss_lower += smear(cfg.mMiss_min_sigma);
    v.mMiss_upper += smear(cfg.mMiss_max_sigma);
    v.kMiss_lower += smear(cfg.kMiss_min_sigma);
    v.kMiss_upper += smear(cfg.kMiss_max_sigma);
    v.pLead_lower += smear(cfg.pLead_min_sigma);
    v.pOverQ_upper += smear(cfg.pOverQ_max_sigma);
    if (varyAngularCut) {
      v.tLead_upper += smear(cfg.tLead_FD_max_sigma);
      v.tLead_lower += smear(cfg.tLead_CD_min_sigma);
    }
    v.pRecoil_lower += smear(cfg.pRecoil_min_sigma);
    return v;
  }

  void apply(const EventRow& ek, bool& passep, bool& passepp,
             const AnalysisOptions& opts) const {
    passep = ek.passep;
    passepp = ek.passepp;
    if (!passep) return;

    bool passLeadRegionAndTheta = false;
    if (opts.leadMode == LEAD_FD) {
      passLeadRegionAndTheta = (ek.leadRegion == 1 && ek.tL < tLead_upper);
    } else if (opts.leadMode == LEAD_CD) {
      passLeadRegionAndTheta = (ek.leadRegion == 2 && ek.tL > tLead_lower);
    } else {
      passLeadRegionAndTheta =
          (ek.leadRegion == 1 && ek.tL < tLead_upper) ||
          (ek.leadRegion == 2 && ek.tL > tLead_lower);
    }

    if (useXBMin && ek.x < xB_lower) passep = false;
    if (useQ2Min && ek.qSq < Q2_lower) passep = false;
    if (useQ2Max && ek.qSq > Q2_upper) passep = false;
    if (useMMissMin && ek.mM < mMiss_lower) passep = false;
    if (useMMissMax && ek.mM > mMiss_upper) passep = false;
    if (useKMissMin && ek.kM < kMiss_lower) passep = false;
    if (useKMissMax && ek.kM > kMiss_upper) passep = false;
    if (useLeadPMin && ek.pL < pLead_lower) passep = false;
    if (usePoverQMax) {
      double pOverQ = (ek.qMag > 0.0) ? (ek.pL / ek.qMag) : 1e9;
      if (pOverQ > pOverQ_upper) passep = false;
    }
    if (useLeadThetaByRegion && !passLeadRegionAndTheta) passep = false;

    if (!passep) {
      passepp = false;
      return;
    }

    if (useRecoilPMin && ek.pR < pRecoil_lower) passepp = false;
  }
};

struct Axis {
  std::string name;
  std::vector<double> edges;
};

struct HistTask {
  std::string name;
  SelectionType selection;
  std::function<bool(const EventRow&)> pass;
  std::vector<Axis*> selectorAxes;
  std::vector<std::function<double(const EventRow&)>> selectorValues;
  std::function<double(const EventRow&)> value;
  int nbins = 0;
  double xmin = 0.0;
  double xmax = 0.0;
  std::vector<double> valueEdges;
};

struct HistPack {
  std::vector<int> selectorBin;
  std::vector<double> selectorCenter;
  TH1D* nominal = nullptr;
  std::vector<TH1D*> toys;
};

struct DiffRow {
  std::string task_name;
  std::string selection;
  std::vector<int> axis_bin;
  std::vector<double> axis_center;
  int value_bin = 0;
  double value_center = 0.0;
  double value_error_low = 0.0;
  double value_error_high = 0.0;
  double count = 0.0;
  double stat_error = 0.0;
  double sys_error = 0.0;
  double sys_error_up = 0.0;
  double sys_error_down = 0.0;
};

struct IntRow {
  std::string task_name;
  std::string selection;
  std::string pattern;
  std::vector<int> axis_bin;
  std::vector<double> axis_center;
  double count = 0.0;
  double stat_error = 0.0;
  double sys_error = 0.0;
  double sys_error_up = 0.0;
  double sys_error_down = 0.0;
};

int binOf(const std::vector<double>& edges, double x) {
  if (edges.size() < 2) return -1;
  if (!(x >= edges.front() && x < edges.back())) return -1;
  for (size_t i = 0; i + 1 < edges.size(); i++) {
    if (x >= edges[i] && x < edges[i + 1]) return (int)i;
  }
  if (x == edges.back()) return (int)edges.size() - 2;
  return -1;
}

double mean(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double s = 0.0;
  for (double x : v) s += x;
  return s / (double)v.size();
}

double stddev(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double m = mean(v);
  double s2 = 0.0;
  for (double x : v) s2 += (x - m) * (x - m);
  return std::sqrt(s2 / (double)v.size());
}

TLeaf* leafReq(TTree* t, const std::vector<std::string>& names) {
  for (const auto& n : names) {
    if (TLeaf* l = t->GetLeaf(n.c_str())) return l;
  }
  std::string msg = "Missing required branch/leaf [";
  for (size_t i = 0; i < names.size(); i++) {
    if (i) msg += ", ";
    msg += names[i];
  }
  msg += "]";
  throw std::runtime_error(msg);
}

TLeaf* leafOpt(TTree* t, const std::vector<std::string>& names) {
  for (const auto& n : names) {
    if (TLeaf* l = t->GetLeaf(n.c_str())) return l;
  }
  return nullptr;
}

double lv(TLeaf* l, double dflt = 0.0) {
  return l ? l->GetValue() : dflt;
}

bool lb(TLeaf* l, bool dflt = false) {
  if (!l) return dflt;
  return l->GetValue() != 0.0;
}

struct SrcLeaves {
  TLeaf* weight = nullptr;
  TLeaf* weight_ep = nullptr;
  TLeaf* weight_epp = nullptr;
  TLeaf* xB = nullptr;
  TLeaf* Q2 = nullptr;
  TLeaf* eP = nullptr;
  TLeaf* eTheta = nullptr;
  TLeaf* ePhi = nullptr;
  TLeaf* qP = nullptr;
  TLeaf* leadP = nullptr;
  TLeaf* leadTheta = nullptr;
  TLeaf* leadPhi = nullptr;
  TLeaf* recP = nullptr;
  TLeaf* recTheta = nullptr;
  TLeaf* recPhi = nullptr;
  TLeaf* pMiss = nullptr;
  TLeaf* kMiss = nullptr;
  TLeaf* mMiss = nullptr;
  TLeaf* E0miss = nullptr;
  TLeaf* E1miss = nullptr;
  TLeaf* E2miss = nullptr;
  TLeaf* theta_PmQ = nullptr;
  TLeaf* theta_PleadQ = nullptr;
  TLeaf* theta_PmPrec = nullptr;
  TLeaf* theta_PrecQ = nullptr;
  TLeaf* theta_PleadPrec = nullptr;
  TLeaf* pRel = nullptr;
  TLeaf* pCM = nullptr;
  TLeaf* pCMx = nullptr;
  TLeaf* pCMy = nullptr;
  TLeaf* pCMz = nullptr;
  TLeaf* nProtons = nullptr;
  TLeaf* goodLead = nullptr;
  TLeaf* nGoodLeads = nullptr;
  TLeaf* leadRegion = nullptr;
};

SrcLeaves bindLeaves(TTree* t) {
  SrcLeaves L;
  L.weight = leafOpt(t, {"weight"});
  L.weight_ep = leafOpt(t, {"weight_ep"});
  L.weight_epp = leafOpt(t, {"weight_epp"});
  L.xB = leafReq(t, {"xB"});
  L.Q2 = leafReq(t, {"Q2"});
  L.eP = leafReq(t, {"eP", "pE"});
  L.eTheta = leafReq(t, {"eTheta", "thetaE"});
  L.ePhi = leafReq(t, {"ePhi", "phiE"});
  L.qP = leafReq(t, {"qP", "q_mag"});
  L.leadP = leafReq(t, {"leadP", "pLead"});
  L.leadTheta = leafReq(t, {"leadTheta", "thetaLead"});
  L.leadPhi = leafReq(t, {"leadPhi", "phiLead"});
  L.recP = leafOpt(t, {"recP", "pRec"});
  L.recTheta = leafOpt(t, {"recTheta", "thetaRec"});
  L.recPhi = leafOpt(t, {"recPhi", "phiRec"});
  L.pMiss = leafReq(t, {"pMiss"});
  L.kMiss = leafReq(t, {"kMiss"});
  L.mMiss = leafReq(t, {"mMiss"});
  L.E0miss = leafReq(t, {"E0miss"});
  L.E1miss = leafReq(t, {"E1miss"});
  L.E2miss = leafOpt(t, {"E2miss"});
  L.theta_PmQ = leafReq(t, {"theta_PmQ", "thetaPmissQ"});
  L.theta_PleadQ = leafReq(t, {"theta_PleadQ", "thetaPLeadQ"});
  L.theta_PmPrec = leafOpt(t, {"theta_PmPrec", "thetaPRecPmiss", "theta_mMissPrec"});
  L.theta_PrecQ = leafOpt(t, {"theta_PrecQ", "thetaPRecQ"});
  L.theta_PleadPrec = leafOpt(t, {"theta_PleadPrec", "thetaPLeadPrec"});
  L.pRel = leafOpt(t, {"pRel"});
  L.pCM = leafOpt(t, {"pCM"});
  L.pCMx = leafOpt(t, {"pCMx", "pcmX"});
  L.pCMy = leafOpt(t, {"pCMy", "pcmY"});
  L.pCMz = leafOpt(t, {"pCMz", "pcmZ"});
  L.nProtons = leafOpt(t, {"nProtons"});
  L.goodLead = leafOpt(t, {"goodLead"});
  L.nGoodLeads = leafOpt(t, {"nGoodLeads"});
  L.leadRegion = leafOpt(t, {"leadRegion"});
  return L;
}

int inferRegion(double leadThetaDeg) {
  if (leadThetaDeg < 37.0) return 1;
  if (leadThetaDeg > 45.0) return 2;
  return 0;
}

bool makeEvent(const SrcLeaves& L, EventRow& e, double& wEP, double& wEPP) {
  e = EventRow();

  const double kRadToDeg = 180.0 / M_PI;

  double wb = (L.weight ? lv(L.weight, 1.0) : 1.0);
  wEP = (L.weight_ep ? lv(L.weight_ep, wb) : wb);
  wEPP = (L.weight_epp ? lv(L.weight_epp, wEP) : wEP);

  e.qSq = lv(L.Q2);
  e.x = lv(L.xB);
  e.pE = lv(L.eP);
  e.tE = lv(L.eTheta) * kRadToDeg;
  e.phiE = lv(L.ePhi) * kRadToDeg;
  e.qMag = lv(L.qP);

  e.pL = lv(L.leadP);
  e.tL = lv(L.leadTheta) * kRadToDeg;
  e.phiL = lv(L.leadPhi) * kRadToDeg;

  e.pR = lv(L.recP, 0.0);
  e.tR = lv(L.recTheta, 0.0) * kRadToDeg;
  e.phiR = lv(L.recPhi, 0.0) * kRadToDeg;

  e.pM = lv(L.pMiss);
  e.kM = lv(L.kMiss);
  e.mM = lv(L.mMiss);
  e.E0 = lv(L.E0miss);
  e.E1 = lv(L.E1miss);
  e.E1adj = e.E1;
  e.E2 = lv(L.E2miss, -9.0);

  e.thPMissQ = lv(L.theta_PmQ) * kRadToDeg;
  e.thPLeadQ = lv(L.theta_PleadQ) * kRadToDeg;
  e.thPRecPMiss = lv(L.theta_PmPrec, 0.0) * kRadToDeg;
  e.thPRecQ = lv(L.theta_PrecQ, 0.0) * kRadToDeg;
  e.thPLeadPrec = lv(L.theta_PleadPrec, 0.0) * kRadToDeg;

  e.pRel = lv(L.pRel, 0.0);
  e.pCM = lv(L.pCM, 0.0);
  e.px = lv(L.pCMx, 0.0);
  e.py = lv(L.pCMy, 0.0);
  e.pz = lv(L.pCMz, 0.0);

  if (L.leadRegion) e.leadRegion = (int)lv(L.leadRegion, 0.0);
  else e.leadRegion = inferRegion(e.tL);

  bool leadOk = true;
  if (L.goodLead) leadOk = lb(L.goodLead, true);
  else if (L.nGoodLeads) leadOk = lv(L.nGoodLeads, 0.0) > 0.5;
  else leadOk = e.pL > 0.0;

  bool hasRecoil = false;
  if (L.nProtons) hasRecoil = lv(L.nProtons, 0.0) > 1.5;
  hasRecoil = hasRecoil && (e.pR > 0.0);

  e.passep = leadOk;
  e.passepp = leadOk && hasRecoil;
  return e.passep;
}

int inferMassNumberFromInputPath(const std::string& inputPath) {
  std::string base = inputPath;
  size_t slash = base.find_last_of("/");
  if (slash != std::string::npos) base = base.substr(slash + 1);

  // First, prefer any explicit integer in the basename (e.g. c12_..., 40ca_...).
  for (size_t i = 0; i < base.size();) {
    if (std::isdigit((unsigned char)base[i])) {
      int value = 0;
      size_t j = i;
      while (j < base.size() && std::isdigit((unsigned char)base[j])) {
        value = value * 10 + (base[j] - '0');
        j++;
      }
      if (value > 1 && value <= 300) return value;
      i = j;
    } else {
      i++;
    }
  }

  // Then fall back to common target tags when no number is present.
  std::string lower = base;
  for (char& c : lower) c = (char)std::tolower((unsigned char)c);

  if (lower.find("sn") != std::string::npos) return 120;
  if (lower.find("ca") != std::string::npos) return 40;
  if (lower.find("he") != std::string::npos) return 4;
  if (lower.find("c") != std::string::npos) return 12;
  if (lower.find("d") != std::string::npos) return 2;

  return 0;
}

std::string histName(const std::string& label, const HistTask& t, const std::vector<int>& selectorBin) {
  std::string n = label + "_" + t.name;
  for (size_t i = 0; i < selectorBin.size(); i++) {
    n += "_" + t.selectorAxes[i]->name + std::to_string(selectorBin[i]);
  }
  return n;
}

void buildSelectorCombosRec(const HistTask& t, int iax,
                            std::vector<int>& curBin,
                            std::vector<double>& curCenter,
                            std::vector<std::pair<std::vector<int>, std::vector<double>>>& out) {
  if (iax == (int)t.selectorAxes.size()) {
    out.push_back({curBin, curCenter});
    return;
  }
  Axis* ax = t.selectorAxes[iax];
  for (size_t b = 0; b + 1 < ax->edges.size(); b++) {
    curBin.push_back((int)b);
    curCenter.push_back(0.5 * (ax->edges[b] + ax->edges[b + 1]));
    buildSelectorCombosRec(t, iax + 1, curBin, curCenter, out);
    curBin.pop_back();
    curCenter.pop_back();
  }
}

std::vector<HistTask> buildTasks() {
  std::vector<HistTask> tasks;

  static std::vector<double> bE_Q2 = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
  static std::vector<double> bE_pmiss = {0.4, 0.55, 0.7, 0.85, 1.0};
  static std::vector<double> bE_kmiss = {0.3, 0.45, 0.6, 0.75, 0.9};
  static std::vector<double> bE_pmiss_long = {0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.85, 0.95, 1.05, 1.2};
  static std::vector<double> bE_kmiss_long = {0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};

  static Axis Q2Axis{"Q2", bE_Q2};
  static Axis pMissAxis{"pMiss", bE_pmiss};
  static Axis kMissAxis{"kMiss", bE_kmiss};

  auto passEP = [](const EventRow& e) { return e.passep; };
  auto passEPP = [](const EventRow& e) { return e.passepp; };

  auto add1D = [&](const std::string& name, SelectionType sel,
                   std::function<bool(const EventRow&)> pass,
                   std::function<double(const EventRow&)> val,
                   int nb, double xmin, double xmax,
                   const std::vector<double>& edges) {
    HistTask t;
    t.name = name;
    t.selection = sel;
    t.pass = pass;
    t.value = val;
    t.nbins = nb;
    t.xmin = xmin;
    t.xmax = xmax;
    t.valueEdges = edges;
    tasks.push_back(t);
  };

  int ep_bins = 30;
  int epp_bins = 20;

  add1D("Q2_ep", SEL_EP, passEP, [](const EventRow& e) { return e.qSq; }, ep_bins, 1.5, 5.0, {});
  add1D("Q2_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.qSq; }, epp_bins, 1.5, 5.0, {});
  add1D("xB_ep", SEL_EP, passEP, [](const EventRow& e) { return e.x; }, ep_bins, 1.0, 2.0, {});
  add1D("xB_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.x; }, epp_bins, 1.0, 2.0, {});
  add1D("pMiss_ep", SEL_EP, passEP, [](const EventRow& e) { return e.pM; }, 0, 0, 0, bE_pmiss_long);
  add1D("pMiss_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.pM; }, 0, 0, 0, bE_pmiss_long);
  add1D("kMiss_ep", SEL_EP, passEP, [](const EventRow& e) { return e.kM; }, 0, 0, 0, bE_kmiss_long);
  add1D("kMiss_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.kM; }, 0, 0, 0, bE_kmiss_long);

  add1D("pMiss_ep_note", SEL_EP, passEP, [](const EventRow& e) { return e.pM; }, ep_bins, 0.4, 1.0, {});
  add1D("pMiss_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.pM; }, epp_bins, 0.4, 1.0, {});
  add1D("kMiss_ep_note", SEL_EP, passEP, [](const EventRow& e) { return e.kM; }, ep_bins, 0.3, 1.0, {});
  add1D("kMiss_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.kM; }, epp_bins, 0.3, 1.0, {});
  add1D("mMiss_ep_note", SEL_EP, passEP, [](const EventRow& e) { return e.mM; }, ep_bins, 0.5, 1.3, {});
  add1D("mMiss_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.mM; }, epp_bins, 0.5, 1.3, {});
  add1D("EStar_ep_note", SEL_EP, passEP, [](const EventRow& e) { return e.E0; }, ep_bins, -0.1, 0.6, {});
  add1D("EStar_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.E0; }, epp_bins, -0.1, 0.6, {});
  add1D("E1miss_adj_ep_note", SEL_EP, passEP, [](const EventRow& e) { return e.E1adj; }, ep_bins, -0.3, 0.3, {});
  add1D("E1miss_adj_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.E1adj; }, epp_bins, -0.3, 0.3, {});
  add1D("E1miss_ep_note", SEL_EP, passEP, [](const EventRow& e) { return e.E1; }, ep_bins, -0.25, 0.55, {});
  add1D("E1miss_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.E1; }, epp_bins, -0.2, 0.5, {});
  add1D("E2miss_epp_note", SEL_EPP, passEPP, [](const EventRow& e) { return e.E2; }, epp_bins, -0.25, 0.5, {});

  add1D("pElectron_ep", SEL_EP, passEP, [](const EventRow& e) { return e.pE; }, ep_bins, 3.5, 6.0, {});
  add1D("pElectron_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.pE; }, epp_bins, 3.5, 6.0, {});
  add1D("thetaElectron_ep", SEL_EP, passEP, [](const EventRow& e) { return e.tE; }, ep_bins, 10, 30, {});
  add1D("thetaElectron_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.tE; }, epp_bins, 10, 30, {});
  add1D("phiElectron_ep", SEL_EP, passEP, [](const EventRow& e) { return e.phiE; }, ep_bins * 2, -180, 180, {});
  add1D("phiElectron_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.phiE; }, epp_bins * 2, -180, 180, {});
  add1D("pLead_ep", SEL_EP, passEP, [](const EventRow& e) { return e.pL; }, ep_bins, 1, 3, {});
  add1D("pLead_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.pL; }, epp_bins, 1, 3, {});
  add1D("thetaLead_ep", SEL_EP, passEP, [](const EventRow& e) { return e.tL; }, ep_bins, 0, 125, {});
  add1D("thetaLead_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.tL; }, epp_bins, 0, 125, {});
  add1D("phiLead_ep", SEL_EP, passEP, [](const EventRow& e) { return e.phiL; }, ep_bins, -180, 180, {});
  add1D("phiLead_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.phiL; }, epp_bins, -180, 180, {});
  add1D("thetaRec_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.tR; }, epp_bins, 0, 180, {});
  add1D("phiRec_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.phiR; }, epp_bins, -180, 180, {});
  add1D("theta_pRec_pmiss_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.thPRecPMiss; }, epp_bins, 80, 180, {});
  add1D("theta_pRecq_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.thPRecQ; }, epp_bins, 0, 180, {});
  add1D("theta_pLeadPrec_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.thPLeadPrec; }, epp_bins, 0, 180, {});
  add1D("q_ep", SEL_EP, passEP, [](const EventRow& e) { return e.qMag; }, ep_bins, 1, 3.5, {});
  add1D("q_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.qMag; }, epp_bins, 1, 3.5, {});
    add1D("pOverQ_ep_note", SEL_EP, passEP,
      [](const EventRow& e) { return (e.qMag > 0.0) ? (e.pL / e.qMag) : -1.0; },
      ep_bins, 0.3, 1.2, {});
    add1D("pOverQ_epp_note", SEL_EPP, passEPP,
      [](const EventRow& e) { return (e.qMag > 0.0) ? (e.pL / e.qMag) : -1.0; },
      epp_bins, 0.3, 1.2, {});
  add1D("theta_pmiss_ep", SEL_EP, passEP, [](const EventRow& e) { return e.thPMissQ; }, ep_bins, 100, 165, {});
  add1D("theta_pmiss_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.thPMissQ; }, epp_bins, 100, 165, {});
  add1D("theta_pLeadq_ep", SEL_EP, passEP, [](const EventRow& e) { return e.thPLeadQ; }, ep_bins, 0, 45, {});
  add1D("theta_pLeadq_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.thPLeadQ; }, epp_bins, 0, 45, {});
  add1D("pRel_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.pRel; }, epp_bins, 0.15, 1.0, {});
  add1D("pRec_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.pR; }, epp_bins, 0.3, 1.0, {});
  add1D("pcmx_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.px; }, epp_bins, -0.75, 0.75, {});
  add1D("pcmy_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.py; }, epp_bins, -0.75, 0.75, {});
  add1D("pcmz_epp", SEL_EPP, passEPP, [](const EventRow& e) { return e.pz; }, epp_bins, -0.3, 1.1, {});
  add1D("pcm_epp", SEL_EPP, passEPP,
        [](const EventRow& e) { return std::sqrt(e.px * e.px + e.py * e.py + e.pz * e.pz); },
        epp_bins, 0.0, 1.0, {});

  auto addSelTask = [&](const std::string& name, SelectionType sel,
                        std::function<bool(const EventRow&)> pass,
                        std::vector<Axis*> axes,
                        std::vector<std::function<double(const EventRow&)>> axisValues,
                        std::function<double(const EventRow&)> val,
                        int nb, double xmin, double xmax,
                        const std::vector<double>& edges) {
    HistTask t;
    t.name = name;
    t.selection = sel;
    t.pass = pass;
    t.selectorAxes = axes;
    t.selectorValues = axisValues;
    t.value = val;
    t.nbins = nb;
    t.xmin = xmin;
    t.xmax = xmax;
    t.valueEdges = edges;
    tasks.push_back(t);
  };

  addSelTask("Q2_ep_SRC_Q2", SEL_EP, passEP, {&Q2Axis},
             {[](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.qSq; }, 50, 1.5, 5.0, {});
  addSelTask("Q2_epp_SRC_Q2", SEL_EPP, passEPP, {&Q2Axis},
             {[](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.qSq; }, 50, 1.5, 5.0, {});

  addSelTask("Q2_ep_SRC_pmiss", SEL_EP, passEP, {&pMissAxis},
             {[](const EventRow& e) { return e.pM; }},
             [](const EventRow& e) { return e.qSq; }, 0, 0, 0, bE_Q2);
  addSelTask("Q2_epp_SRC_pmiss", SEL_EPP, passEPP, {&pMissAxis},
             {[](const EventRow& e) { return e.pM; }},
             [](const EventRow& e) { return e.qSq; }, 0, 0, 0, bE_Q2);
  addSelTask("Q2_ep_SRC_kmiss", SEL_EP, passEP, {&kMissAxis},
             {[](const EventRow& e) { return e.kM; }},
             [](const EventRow& e) { return e.qSq; }, 0, 0, 0, bE_Q2);
  addSelTask("Q2_epp_SRC_kmiss", SEL_EPP, passEPP, {&kMissAxis},
             {[](const EventRow& e) { return e.kM; }},
             [](const EventRow& e) { return e.qSq; }, 0, 0, 0, bE_Q2);

  addSelTask("E0miss_ep_SRC_pmiss", SEL_EP, passEP, {&pMissAxis},
             {[](const EventRow& e) { return e.pM; }},
             [](const EventRow& e) { return e.E0; }, 50, -0.1, 0.6, {});
  addSelTask("E1miss_ep_SRC_pmiss", SEL_EP, passEP, {&pMissAxis},
             {[](const EventRow& e) { return e.pM; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.25, 0.55, {});
  addSelTask("E1miss_epp_SRC_pmiss", SEL_EPP, passEPP, {&pMissAxis},
             {[](const EventRow& e) { return e.pM; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.2, 0.5, {});
  addSelTask("E2miss_epp_SRC_pmiss", SEL_EPP, passEPP, {&pMissAxis},
             {[](const EventRow& e) { return e.pM; }},
             [](const EventRow& e) { return e.E2; }, 15, -0.25, 0.5, {});

  addSelTask("E0miss_ep_SRC_kmiss", SEL_EP, passEP, {&kMissAxis},
             {[](const EventRow& e) { return e.kM; }},
             [](const EventRow& e) { return e.E0; }, 30, -0.1, 0.85, {});
  addSelTask("E1miss_ep_SRC_kmiss", SEL_EP, passEP, {&kMissAxis},
             {[](const EventRow& e) { return e.kM; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.25, 0.55, {});
  addSelTask("E1miss_epp_SRC_kmiss", SEL_EPP, passEPP, {&kMissAxis},
             {[](const EventRow& e) { return e.kM; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.2, 0.5, {});
  addSelTask("E2miss_epp_SRC_kmiss", SEL_EPP, passEPP, {&kMissAxis},
             {[](const EventRow& e) { return e.kM; }},
             [](const EventRow& e) { return e.E2; }, 15, -0.25, 0.5, {});

  addSelTask("E1miss_ep_SRC_pmiss_Q2", SEL_EP, passEP, {&pMissAxis, &Q2Axis},
             {[](const EventRow& e) { return e.pM; }, [](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.25, 0.55, {});
  addSelTask("E1miss_epp_SRC_pmiss_Q2", SEL_EPP, passEPP, {&pMissAxis, &Q2Axis},
             {[](const EventRow& e) { return e.pM; }, [](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.2, 0.5, {});
  addSelTask("E2miss_epp_SRC_pmiss_Q2", SEL_EPP, passEPP, {&pMissAxis, &Q2Axis},
             {[](const EventRow& e) { return e.pM; }, [](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.E2; }, 15, -0.25, 0.5, {});

  addSelTask("E1miss_ep_SRC_kmiss_Q2", SEL_EP, passEP, {&kMissAxis, &Q2Axis},
             {[](const EventRow& e) { return e.kM; }, [](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.25, 0.55, {});
  addSelTask("E1miss_epp_SRC_kmiss_Q2", SEL_EPP, passEPP, {&kMissAxis, &Q2Axis},
             {[](const EventRow& e) { return e.kM; }, [](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.E1; }, 15, -0.2, 0.5, {});
  addSelTask("E2miss_epp_SRC_kmiss_Q2", SEL_EPP, passEPP, {&kMissAxis, &Q2Axis},
             {[](const EventRow& e) { return e.kM; }, [](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.E2; }, 15, -0.25, 0.5, {});

  addSelTask("pcmx_epp_SRC_Q2", SEL_EPP, passEPP, {&Q2Axis},
             {[](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.px; }, 15, -0.75, 0.75, {});
  addSelTask("pcmy_epp_SRC_Q2", SEL_EPP, passEPP, {&Q2Axis},
             {[](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.py; }, 15, -0.75, 0.75, {});
  addSelTask("pcmz_epp_SRC_Q2", SEL_EPP, passEPP, {&Q2Axis},
             {[](const EventRow& e) { return e.qSq; }},
             [](const EventRow& e) { return e.pz; }, 15, -0.5, 1.0, {});

  return tasks;
}

void bookTask(const HistTask& t, int nToys,
              std::vector<HistPack>& packs,
              TDirectory* nomDir,
              const std::vector<TDirectory*>& toyDirs) {
  std::vector<std::pair<std::vector<int>, std::vector<double>>> combos;
  std::vector<int> cb;
  std::vector<double> cc;
  buildSelectorCombosRec(t, 0, cb, cc, combos);
  if (combos.empty()) combos.push_back({{}, {}});

  for (const auto& c : combos) {
    HistPack p;
    p.selectorBin = c.first;
    p.selectorCenter = c.second;

    std::string nNom = histName("nominal", t, p.selectorBin);
    if (!t.valueEdges.empty()) {
      p.nominal = new TH1D(nNom.c_str(), nNom.c_str(),
                           (int)t.valueEdges.size() - 1, &t.valueEdges[0]);
    } else {
      p.nominal = new TH1D(nNom.c_str(), nNom.c_str(), t.nbins, t.xmin, t.xmax);
    }
    p.nominal->Sumw2();
    nomDir->Add(p.nominal);

    p.toys.resize((size_t)nToys, nullptr);
    for (int i = 0; i < nToys; i++) {
      char label[32];
      snprintf(label, sizeof(label), "toy_%03d", i);
      std::string nToy = histName(label, t, p.selectorBin);
      if (!t.valueEdges.empty()) {
        p.toys[(size_t)i] = new TH1D(nToy.c_str(), nToy.c_str(),
                                     (int)t.valueEdges.size() - 1, &t.valueEdges[0]);
      } else {
        p.toys[(size_t)i] = new TH1D(nToy.c_str(), nToy.c_str(), t.nbins, t.xmin, t.xmax);
      }
      p.toys[(size_t)i]->Sumw2();
      toyDirs[(size_t)i]->Add(p.toys[(size_t)i]);
    }

    packs.push_back(p);
  }
}

bool selectorMatch(const HistTask& t, const EventRow& e, const HistPack& p) {
  if (t.selectorAxes.empty()) return true;
  for (size_t i = 0; i < t.selectorAxes.size(); i++) {
    double x = t.selectorValues[i](e);
    int b = binOf(t.selectorAxes[i]->edges, x);
    if (b < 0 || b != p.selectorBin[i]) return false;
  }
  return true;
}

std::string selName(SelectionType s) {
  return (s == SEL_EP) ? "ep" : "epp";
}

void addRatioRows(const HistTask& numTask,
                  const std::vector<HistPack>& numPacks,
                  const HistTask& denTask,
                  const std::vector<HistPack>& denPacks,
                  int nToys,
                  std::vector<DiffRow>& outRows) {
  std::map<std::string, const HistPack*> denMap;
  for (const auto& dp : denPacks) {
    std::string k;
    for (size_t i = 0; i < dp.selectorBin.size(); i++) {
      if (i) k += "_";
      k += std::to_string(dp.selectorBin[i]);
    }
    denMap[k] = &dp;
  }

  for (const auto& np : numPacks) {
    std::string k;
    for (size_t i = 0; i < np.selectorBin.size(); i++) {
      if (i) k += "_";
      k += std::to_string(np.selectorBin[i]);
    }
    auto it = denMap.find(k);
    if (it == denMap.end()) continue;
    const HistPack* dp = it->second;

    int nBins = np.nominal->GetNbinsX();
    for (int ib = 1; ib <= nBins; ib++) {
      double n = np.nominal->GetBinContent(ib);
      double en = np.nominal->GetBinError(ib);
      double d = dp->nominal->GetBinContent(ib);
      double ed = dp->nominal->GetBinError(ib);
      if (d == 0.0) continue;

      DiffRow r;
      r.task_name = numTask.name + "_over_" + denTask.name;
      r.selection = "ratio";
      r.axis_bin = np.selectorBin;
      r.axis_center = np.selectorCenter;
      r.value_bin = ib - 1;
      r.value_center = np.nominal->GetXaxis()->GetBinCenter(ib);
      r.value_error_low = 0.5 * np.nominal->GetXaxis()->GetBinWidth(ib);
      r.value_error_high = 0.5 * np.nominal->GetXaxis()->GetBinWidth(ib);
      r.count = n / d;
      r.stat_error = std::abs(r.count) * std::sqrt((n > 0.0 ? (en / n) * (en / n) : 0.0) +
                                                   (d > 0.0 ? (ed / d) * (ed / d) : 0.0));

      if (nToys > 0) {
        std::vector<double> toyRatio;
        toyRatio.reserve((size_t)nToys);
        for (int itoy = 0; itoy < nToys; itoy++) {
          double nt = np.toys[(size_t)itoy]->GetBinContent(ib);
          double dt = dp->toys[(size_t)itoy]->GetBinContent(ib);
          if (dt == 0.0) continue;
          toyRatio.push_back(nt / dt);
        }
        double se = stddev(toyRatio);
        r.sys_error = se;
        r.sys_error_up = se;
        r.sys_error_down = se;
      }
      outRows.push_back(r);
    }
  }
}

void writeTablesAndGraphs(const std::vector<HistTask>& tasks,
                         const std::vector<std::vector<HistPack>>& allPacks,
                         int nToys,
                         TFile* outFile,
                         std::vector<DiffRow>& outDiffRows,
                         std::vector<IntRow>& outIntRows) {
  // Differential rows from hist bins.
  for (size_t it = 0; it < tasks.size(); it++) {
    const HistTask& t = tasks[it];
    const auto& packs = allPacks[it];
    for (const auto& p : packs) {
      int nBins = p.nominal->GetNbinsX();
      for (int ib = 1; ib <= nBins; ib++) {
        DiffRow r;
        r.task_name = t.name;
        r.selection = selName(t.selection);
        r.axis_bin = p.selectorBin;
        r.axis_center = p.selectorCenter;
        r.value_bin = ib - 1;
        r.value_center = p.nominal->GetXaxis()->GetBinCenter(ib);
        r.value_error_low = 0.5 * p.nominal->GetXaxis()->GetBinWidth(ib);
        r.value_error_high = 0.5 * p.nominal->GetXaxis()->GetBinWidth(ib);
        r.count = p.nominal->GetBinContent(ib);
        r.stat_error = p.nominal->GetBinError(ib);

        if (nToys > 0) {
          std::vector<double> toyVals;
          toyVals.reserve((size_t)nToys);
          for (int itoy = 0; itoy < nToys; itoy++) {
            toyVals.push_back(p.toys[(size_t)itoy]->GetBinContent(ib));
          }
          double se = stddev(toyVals);
          r.sys_error = se;
          r.sys_error_up = se;
          r.sys_error_down = se;
        }
        outDiffRows.push_back(r);
      }
    }
  }

  // Ratio rows for tasks expected by plotting code.
  auto findTask = [&](const std::string& name) -> int {
    for (size_t i = 0; i < tasks.size(); i++) {
      if (tasks[i].name == name) return (int)i;
    }
    return -1;
  };

  std::vector<std::pair<std::string, std::string>> ratioSpecs = {
      {"pMiss_epp", "pMiss_ep"},
      {"kMiss_epp", "kMiss_ep"},
      {"Q2_epp_SRC_pmiss", "Q2_ep_SRC_pmiss"},
      {"Q2_epp_SRC_kmiss", "Q2_ep_SRC_kmiss"},
  };

  for (const auto& sp : ratioSpecs) {
    int in = findTask(sp.first);
    int id = findTask(sp.second);
    if (in < 0 || id < 0) continue;
    addRatioRows(tasks[(size_t)in], allPacks[(size_t)in],
                 tasks[(size_t)id], allPacks[(size_t)id],
                 nToys, outDiffRows);
  }

  // Integrated rows: collapse:none, per-axis collapse, collapse:all.
  for (size_t it = 0; it < tasks.size(); it++) {
    const HistTask& t = tasks[it];
    const auto& packs = allPacks[it];
    int nAxes = (int)t.selectorAxes.size();

    std::vector<std::pair<std::vector<int>, std::string>> patterns;
    patterns.push_back({{}, "collapse:none"});
    for (int ia = 0; ia < nAxes; ia++) patterns.push_back({{ia}, "collapse:" + t.selectorAxes[(size_t)ia]->name});
    if (nAxes > 1) {
      std::vector<int> all;
      for (int ia = 0; ia < nAxes; ia++) all.push_back(ia);
      patterns.push_back({all, "collapse:all"});
    }

    for (const auto& pat : patterns) {
      std::map<std::string, IntRow> groups;
      std::map<std::string, std::vector<double>> groupToyVals;

      for (const auto& p : packs) {
        std::vector<int> keepBins;
        std::vector<double> keepCenters;
        for (int ia = 0; ia < nAxes; ia++) {
          bool collapsed = false;
          for (int ca : pat.first) {
            if (ca == ia) {
              collapsed = true;
              break;
            }
          }
          if (!collapsed) {
            keepBins.push_back(p.selectorBin[(size_t)ia]);
            keepCenters.push_back(p.selectorCenter[(size_t)ia]);
          }
        }

        std::string gk;
        for (size_t i = 0; i < keepBins.size(); i++) {
          if (i) gk += "_";
          gk += std::to_string(keepBins[i]);
        }

        double statErr = 0.0;
        double nomInt = p.nominal->IntegralAndError(1, p.nominal->GetNbinsX(), statErr);

        if (!groups.count(gk)) {
          IntRow row;
          row.task_name = t.name;
          row.selection = selName(t.selection);
          row.pattern = pat.second;
          row.axis_bin = keepBins;
          row.axis_center = keepCenters;
          groups[gk] = row;
          groupToyVals[gk] = std::vector<double>((size_t)nToys, 0.0);
        }

        groups[gk].count += nomInt;
        groups[gk].stat_error = std::sqrt(groups[gk].stat_error * groups[gk].stat_error + statErr * statErr);

        for (int itoy = 0; itoy < nToys; itoy++) {
          double tint = p.toys[(size_t)itoy]->Integral(1, p.toys[(size_t)itoy]->GetNbinsX());
          groupToyVals[gk][(size_t)itoy] += tint;
        }
      }

      for (auto& kv : groups) {
        if (nToys > 0) {
          double se = stddev(groupToyVals[kv.first]);
          kv.second.sys_error = se;
          kv.second.sys_error_up = se;
          kv.second.sys_error_down = se;
        }
        outIntRows.push_back(kv.second);
      }
    }
  }

  // Write tables.
  outFile->cd();

  TTree* diffTree = new TTree("diffTable", "Differential bin table");
  std::string d_task, d_sel;
  std::vector<int> d_axisBin;
  std::vector<double> d_axisCenter;
  int d_valueBin;
  double d_valueCenter, d_valueErrLow, d_valueErrHigh;
  double d_count, d_statErr, d_sysErr, d_sysErrUp, d_sysErrDown;

  diffTree->Branch("task_name", &d_task);
  diffTree->Branch("selection", &d_sel);
  diffTree->Branch("axis_bin", &d_axisBin);
  diffTree->Branch("axis_center", &d_axisCenter);
  diffTree->Branch("value_bin", &d_valueBin, "value_bin/I");
  diffTree->Branch("value_center", &d_valueCenter, "value_center/D");
  diffTree->Branch("value_error_low", &d_valueErrLow, "value_error_low/D");
  diffTree->Branch("value_error_high", &d_valueErrHigh, "value_error_high/D");
  diffTree->Branch("count", &d_count, "count/D");
  diffTree->Branch("stat_error", &d_statErr, "stat_error/D");
  diffTree->Branch("sys_error", &d_sysErr, "sys_error/D");
  diffTree->Branch("sys_error_up", &d_sysErrUp, "sys_error_up/D");
  diffTree->Branch("sys_error_down", &d_sysErrDown, "sys_error_down/D");

  for (const auto& r : outDiffRows) {
    d_task = r.task_name;
    d_sel = r.selection;
    d_axisBin = r.axis_bin;
    d_axisCenter = r.axis_center;
    d_valueBin = r.value_bin;
    d_valueCenter = r.value_center;
    d_valueErrLow = r.value_error_low;
    d_valueErrHigh = r.value_error_high;
    d_count = r.count;
    d_statErr = r.stat_error;
    d_sysErr = r.sys_error;
    d_sysErrUp = r.sys_error_up;
    d_sysErrDown = r.sys_error_down;
    diffTree->Fill();
  }

  TTree* intTree = new TTree("integratedTable", "Integrated (collapsed-axis) table");
  std::string i_task, i_sel, i_pattern;
  std::vector<int> i_axisBin;
  std::vector<double> i_axisCenter;
  double i_count, i_statErr, i_sysErr, i_sysErrUp, i_sysErrDown;

  intTree->Branch("task_name", &i_task);
  intTree->Branch("selection", &i_sel);
  intTree->Branch("pattern", &i_pattern);
  intTree->Branch("axis_bin", &i_axisBin);
  intTree->Branch("axis_center", &i_axisCenter);
  intTree->Branch("count", &i_count, "count/D");
  intTree->Branch("stat_error", &i_statErr, "stat_error/D");
  intTree->Branch("sys_error", &i_sysErr, "sys_error/D");
  intTree->Branch("sys_error_up", &i_sysErrUp, "sys_error_up/D");
  intTree->Branch("sys_error_down", &i_sysErrDown, "sys_error_down/D");

  for (const auto& r : outIntRows) {
    i_task = r.task_name;
    i_sel = r.selection;
    i_pattern = r.pattern;
    i_axisBin = r.axis_bin;
    i_axisCenter = r.axis_center;
    i_count = r.count;
    i_statErr = r.stat_error;
    i_sysErr = r.sys_error;
    i_sysErrUp = r.sys_error_up;
    i_sysErrDown = r.sys_error_down;
    intTree->Fill();
  }

  intTree->Write();
  diffTree->Write();

  // Build graph objects from table rows.
  TDirectory* graphsDir = outFile->mkdir("graphs");
  graphsDir->cd();

  std::map<std::string, std::vector<DiffRow>> dgroups;
  for (const auto& r : outDiffRows) {
    std::string suffix;
    if (r.axis_bin.empty()) suffix = "none";
    else {
      for (size_t i = 0; i < r.axis_bin.size(); i++) {
        if (i) suffix += "_";
        suffix += std::to_string(r.axis_bin[i]);
      }
    }
    std::string gname = r.task_name + "|" + r.selection + "|" + suffix;
    dgroups[gname].push_back(r);
  }

  for (auto& kv : dgroups) {
    auto& rows = kv.second;
    std::sort(rows.begin(), rows.end(), [](const DiffRow& a, const DiffRow& b) {
      return a.value_bin < b.value_bin;
    });
    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    g->SetName(kv.first.c_str());
    for (const auto& r : rows) {
      double eu = std::sqrt(r.stat_error * r.stat_error + r.sys_error_up * r.sys_error_up);
      double ed = std::sqrt(r.stat_error * r.stat_error + r.sys_error_down * r.sys_error_down);
      int n = g->GetN();
      g->SetPoint(n, r.value_center, r.count);
      g->SetPointError(n, r.value_error_low, r.value_error_high, ed, eu);
    }
    g->Write();
  }

  std::map<std::string, std::vector<IntRow>> igroups;
  for (const auto& r : outIntRows) {
    std::string gname = r.task_name + "|" + r.selection + "|" + r.pattern;
    igroups[gname].push_back(r);
  }

  for (auto& kv : igroups) {
    auto& rows = kv.second;
    size_t naxes = rows.empty() ? 0 : rows[0].axis_bin.size();
    if (naxes > 1) continue;
    std::sort(rows.begin(), rows.end(), [](const IntRow& a, const IntRow& b) {
      int ba = a.axis_bin.empty() ? 0 : a.axis_bin[0];
      int bb = b.axis_bin.empty() ? 0 : b.axis_bin[0];
      return ba < bb;
    });

    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    g->SetName(kv.first.c_str());
    for (const auto& r : rows) {
      double x = r.axis_center.empty() ? 0.0 : r.axis_center[0];
      double eu = std::sqrt(r.stat_error * r.stat_error + r.sys_error_up * r.sys_error_up);
      double ed = std::sqrt(r.stat_error * r.stat_error + r.sys_error_down * r.sys_error_down);
      int n = g->GetN();
      g->SetPoint(n, x, r.count);
      g->SetPointError(n, 0.0, 0.0, ed, eu);
    }
    g->Write();
  }
}

struct DerivedSpec {
  std::string quantityName;
  std::string selection;
  std::string sourceTask;
  std::string axisName;
  int nAxisBins;
  int method; // 0=stddev, 1=mean, 2=fit sigma
  double fitMin;
  double fitMax;
};

double G(double x, double N, double mu, double sigma) {
  double z = (x - mu) / sigma;
  return (N / (sigma * std::sqrt(2.0 * M_PI))) * std::exp(-0.5 * z * z);
}

std::pair<double, double> fitSigmaAndError(TH1D* h, double min, double max) {
  if (!h || h->GetEntries() == 0) return {0.0, 0.0};
  if (h->GetEntries() < 20) return {h->GetStdDev(), h->GetStdDevError()};

  TF1* gFit = new TF1("srcSkimGausFit",
                      [](double* x, double* p) { return G(x[0], p[0], p[1], p[2]); },
                      min, max, 3);
  gFit->SetParameter(0, h->GetMaximum() / G(0, 1, 0, 0.1));
  gFit->SetParameter(1, 0.5 * (max + min));
  gFit->SetParLimits(1, min, max);
  gFit->SetParameter(2, 0.25 * (max - min));
  gFit->SetParLimits(2, 0.0, max - min);

  TFitResultPtr fr = h->Fit(gFit, "SrBeqn", "", min, max);
  double sigma = h->GetStdDev();
  double sigmaErr = h->GetStdDevError();
  if (fr.Get() != nullptr) {
    sigma = fr->Parameter(2);
    sigmaErr = fr->ParError(2);
  }
  delete gFit;
  if (!std::isfinite(sigmaErr) || sigmaErr < 0.0) sigmaErr = 0.0;
  return {sigma, sigmaErr};
}

std::vector<DerivedSpec> buildDerivedSpecs() {
  std::vector<DerivedSpec> out;
  out.push_back({"sigma_pcmx", "epp", "pcmx_epp_SRC_Q2", "", 1, 2, -0.2, 0.22});
  out.push_back({"sigma_pcmy", "epp", "pcmy_epp_SRC_Q2", "", 1, 2, -0.2, 0.20});

  struct EmissTask { std::string suffix; std::string selection; std::string task; };
  std::vector<EmissTask> emissTasks = {
      {"E1miss_ep_pmiss", "ep", "E1miss_ep_SRC_pmiss_Q2"},
      {"E1miss_epp_pmiss", "epp", "E1miss_epp_SRC_pmiss_Q2"},
      {"E2miss_epp_pmiss", "epp", "E2miss_epp_SRC_pmiss_Q2"},
      {"E1miss_ep_kmiss", "ep", "E1miss_ep_SRC_kmiss_Q2"},
      {"E1miss_epp_kmiss", "epp", "E1miss_epp_SRC_kmiss_Q2"},
      {"E2miss_epp_kmiss", "epp", "E2miss_epp_SRC_kmiss_Q2"},
  };

  for (const auto& e : emissTasks) {
    std::string axis = (e.suffix.find("pmiss") != std::string::npos) ? "pMiss" : "kMiss";
    out.push_back({"stddev_" + e.suffix, e.selection, e.task, axis, 4, 0, 0.0, 0.0});
    out.push_back({"mean_" + e.suffix, e.selection, e.task, axis, 4, 1, 0.0, 0.0});
  }
  return out;
}

std::string derivedHistName(const std::string& label,
                            const std::string& task,
                            const std::string& axisName,
                            int axisBin,
                            int q2Bin) {
  std::string n = label + "_" + task;
  if (!axisName.empty()) n += "_" + axisName + std::to_string(axisBin);
  n += "_Q2" + std::to_string(q2Bin);
  return n;
}

TH1D* getHist(TDirectory* histsDir, const std::string& label, const std::string& name) {
  TDirectory* d = histsDir->GetDirectory(label.c_str());
  if (!d) return nullptr;
  return dynamic_cast<TH1D*>(d->Get(name.c_str()));
}

void writeDerivedGraphs(TFile* outFile, int nToys, bool statOnlyMode) {
  TDirectory* histsDir = outFile->GetDirectory("hists");
  if (!histsDir) return;
  TDirectory* graphsDir = outFile->GetDirectory("graphs");
  if (!graphsDir) graphsDir = outFile->mkdir("graphs");

  std::vector<double> q2Edges = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
  int nQ2 = (int)q2Edges.size() - 1;

  std::vector<double> q2x;
  q2x.reserve((size_t)nQ2);
  for (int j = 0; j < nQ2; j++) {
    std::string q2Hist = derivedHistName("nominal", "Q2_epp_SRC_Q2", "", 0, j);
    TH1D* h = getHist(histsDir, "nominal", q2Hist);
    if (h && h->GetEntries() > 0) q2x.push_back(h->GetMean());
    else q2x.push_back(0.5 * (q2Edges[(size_t)j] + q2Edges[(size_t)j + 1]));
  }

  std::vector<DerivedSpec> specs = buildDerivedSpecs();
  for (const auto& sp : specs) {
    for (int ab = 0; ab < sp.nAxisBins; ab++) {
      TGraphErrors* g = new TGraphErrors();
      std::string suffix = sp.axisName.empty() ? "none" : std::to_string(ab);
      std::string gname = sp.quantityName + "|" + sp.selection + "|" + suffix;
      g->SetName(gname.c_str());

      for (int jq2 = 0; jq2 < nQ2; jq2++) {
        std::string hNomName = derivedHistName("nominal", sp.sourceTask, sp.axisName, ab, jq2);
        TH1D* hNom = getHist(histsDir, "nominal", hNomName);
        if (!hNom || hNom->GetEntries() == 0) continue;

        auto valErrFromHist = [&](TH1D* h) -> std::pair<double, double> {
          if (sp.method == 0) {
            double e = h->GetStdDevError();
            if (!std::isfinite(e) || e < 0.0) e = 0.0;
            return {h->GetStdDev(), e};
          }
          if (sp.method == 1) {
            double e = h->GetMeanError();
            if (!std::isfinite(e) || e < 0.0) e = 0.0;
            return {h->GetMean(), e};
          }
          return fitSigmaAndError(h, sp.fitMin, sp.fitMax);
        };

        double y = 0.0;
        double ey = 0.0;
        if (statOnlyMode || nToys <= 0) {
          auto ve = valErrFromHist(hNom);
          y = ve.first;
          ey = ve.second;
        } else {
          std::vector<double> toyVals;
          toyVals.reserve((size_t)nToys);
          for (int itoy = 0; itoy < nToys; itoy++) {
            char label[32];
            snprintf(label, sizeof(label), "toy_%03d", itoy);
            std::string hToyName = derivedHistName(label, sp.sourceTask, sp.axisName, ab, jq2);
            TH1D* hToy = getHist(histsDir, label, hToyName);
            if (!hToy || hToy->GetEntries() == 0) continue;
            toyVals.push_back(valErrFromHist(hToy).first);
          }
          if (toyVals.empty()) {
            auto ve = valErrFromHist(hNom);
            y = ve.first;
            ey = ve.second;
          } else {
            y = mean(toyVals);
            ey = stddev(toyVals);
          }
        }

        int n = g->GetN();
        g->SetPoint(n, q2x[(size_t)jq2], y);
        g->SetPointError(n, 0.0, ey);
      }

      if (g->GetN() > 0) {
        graphsDir->cd();
        g->Write("", TObject::kOverwrite);
      } else {
        delete g;
      }
    }
  }
}

void runPipeline(const char* inputFile,
                 const char* outputFile,
                 bool withSystematics,
                 int nToys,
                 unsigned int toySeed,
                 const char* treeName,
                 const char* leadMode,
                 const CutConfig& cutConfig,
                 int massNumber = 0) {
  if (!inputFile || !outputFile) {
    throw std::runtime_error("Input/output file path is null.");
  }

  AnalysisOptions opt;
  std::string lm = leadMode ? std::string(leadMode) : "fd";
  if (lm == "fd" || lm == "FD") opt.leadMode = LEAD_FD;
  else if (lm == "cd" || lm == "CD") opt.leadMode = LEAD_CD;
  else if (lm == "both" || lm == "BOTH") opt.leadMode = LEAD_BOTH;
  else throw std::runtime_error("leadMode must be 'fd', 'cd', or 'both'.");
  int toys = withSystematics ? nToys : 0;
  if (toys < 0) toys = 0;

  int A = massNumber;
  if (A <= 1) A = inferMassNumberFromInputPath(inputFile ? std::string(inputFile) : std::string(""));
  if (A <= 1) {
    A = 12;
    std::cerr << "[src_skim_to_plot_common] WARNING: Could not infer mass number A from input file name; "
              << "falling back to A=12. Pass massNumber explicitly to runPipeline/build macro to override." << std::endl;
  }
  const double mN_GeV = 0.939565;
  const double e1AdjCoeff = ((double)A - 2.0) / (((double)A - 1.0) * 2.0 * mN_GeV);
  std::cout << "[src_skim_to_plot_common] Using A=" << A
            << " for E1miss_adj with coefficient ((A-2)/(A-1))/(2mN)=" << e1AdjCoeff
            << " GeV^-1" << std::endl;

  TFile* in = TFile::Open(inputFile, "READ");
  if (!in || in->IsZombie()) {
    throw std::runtime_error(std::string("Could not open input file: ") + inputFile);
  }

  std::string tname = treeName ? std::string(treeName) : std::string("srcTree");
  TTree* tree = dynamic_cast<TTree*>(in->Get(tname.c_str()));
  if (!tree) {
    in->Close();
    throw std::runtime_error("Could not find tree '" + tname + "' in input file.");
  }

  SrcLeaves leaves = bindLeaves(tree);
  std::vector<HistTask> tasks = buildTasks();

  TFile* out = TFile::Open(outputFile, "RECREATE");
  if (!out || out->IsZombie()) {
    in->Close();
    throw std::runtime_error(std::string("Could not create output file: ") + outputFile);
  }

  TH1::AddDirectory(kFALSE);
  out->cd();
  TDirectory* histsDir = out->mkdir("hists");
  TDirectory* nomDir = histsDir->mkdir("nominal");
  std::vector<TDirectory*> toyDirs((size_t)toys, nullptr);
  for (int i = 0; i < toys; i++) {
    char label[32];
    snprintf(label, sizeof(label), "toy_%03d", i);
    toyDirs[(size_t)i] = histsDir->mkdir(label);
  }

  std::vector<std::vector<HistPack>> allPacks(tasks.size());
  for (size_t i = 0; i < tasks.size(); i++) {
    bookTask(tasks[i], toys, allPacks[i], nomDir, toyDirs);
  }

  CutVariation nominal = CutVariation::Nominal(cutConfig);
  bool varyAngularCut = (opt.leadMode != LEAD_BOTH);
  std::vector<CutVariation> toyCuts((size_t)toys);
  for (int i = 0; i < toys; i++) {
    toyCuts[(size_t)i] = CutVariation::Randomized(cutConfig, toySeed, i, varyAngularCut);
  }

  Long64_t nEntries = tree->GetEntries();
  Long64_t nUsed = 0;
  const Long64_t progressStep = (nEntries > 0) ? std::max<Long64_t>(1, nEntries / 50) : 1;
  auto t0 = std::chrono::steady_clock::now();
  std::cout << "[src_skim_to_plot_common] Starting event loop: 0 / " << nEntries << std::endl;

  for (Long64_t i = 0; i < nEntries; i++) {
    if ((i % progressStep) == 0 || i == nEntries - 1) {
      auto now = std::chrono::steady_clock::now();
      double elapsed = std::chrono::duration<double>(now - t0).count();
      double frac = (nEntries > 0) ? ((double)i / (double)nEntries) : 1.0;
      double eta = (frac > 1e-12) ? (elapsed * (1.0 / frac - 1.0)) : 0.0;
      std::cout << "[src_skim_to_plot_common] Iteration " << i << " / " << nEntries
                << " | elapsed " << elapsed << " s"
                << " | eta " << eta << " s"
                << "\r" << std::flush;
    }

    tree->GetEntry(i);

    EventRow ev;
    double wEP = 1.0;
    double wEPP = 1.0;
    if (!makeEvent(leaves, ev, wEP, wEPP)) continue;
    ev.E1adj = ev.E1 - e1AdjCoeff * ev.pM * ev.pM;

    bool passEPNom = false;
    bool passEPPNom = false;
    nominal.apply(ev, passEPNom, passEPPNom, opt);
    EventRow evNom = ev;
    evNom.passep = passEPNom;
    evNom.passepp = passEPPNom;

    for (size_t it = 0; it < tasks.size(); it++) {
      const HistTask& t = tasks[it];
      bool pass = t.pass(evNom);
      if (!pass) continue;
      double w = (t.selection == SEL_EP) ? wEP : wEPP;
      double v = t.value(evNom);
      for (auto& p : allPacks[it]) {
        if (!selectorMatch(t, evNom, p)) continue;
        p.nominal->Fill(v, w);
      }
    }

    for (int itoy = 0; itoy < toys; itoy++) {
      bool passEPToy = false;
      bool passEPPToy = false;
      toyCuts[(size_t)itoy].apply(ev, passEPToy, passEPPToy, opt);
      EventRow evToy = ev;
      evToy.passep = passEPToy;
      evToy.passepp = passEPPToy;

      for (size_t it = 0; it < tasks.size(); it++) {
        const HistTask& t = tasks[it];
        bool pass = t.pass(evToy);
        if (!pass) continue;
        double w = (t.selection == SEL_EP) ? wEP : wEPP;
        double v = t.value(evToy);
        for (auto& p : allPacks[it]) {
          if (!selectorMatch(t, evToy, p)) continue;
          p.toys[(size_t)itoy]->Fill(v, w);
        }
      }
    }

    nUsed++;
  }

  {
    auto t1 = std::chrono::steady_clock::now();
    double total = std::chrono::duration<double>(t1 - t0).count();
    std::cout << std::endl
              << "[src_skim_to_plot_common] Finished event loop: "
              << nEntries << " / " << nEntries
              << " in " << total << " s" << std::endl;
  }

  out->cd();
  nomDir->cd();
  for (const auto& packs : allPacks) {
    for (const auto& p : packs) p.nominal->Write("", TObject::kOverwrite);
  }
  for (int itoy = 0; itoy < toys; itoy++) {
    toyDirs[(size_t)itoy]->cd();
    for (const auto& packs : allPacks) {
      for (const auto& p : packs) p.toys[(size_t)itoy]->Write("", TObject::kOverwrite);
    }
  }

  std::vector<DiffRow> diffRows;
  std::vector<IntRow> intRows;
  writeTablesAndGraphs(tasks, allPacks, toys, out, diffRows, intRows);
  writeDerivedGraphs(out, toys, !withSystematics);

  out->cd();
  out->Write("", TObject::kOverwrite);
  out->Close();
  in->Close();

  std::cout << "[src_skim_to_plot_common] Processed " << nEntries
            << " entries, used " << nUsed
            << ". Wrote " << outputFile << std::endl;
}

}  // namespace src_skim_plot
