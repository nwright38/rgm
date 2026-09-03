#ifndef OVERLAY_SRC_COMMON_H
#define OVERLAY_SRC_COMMON_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TMath.h>
#include <TString.h>

namespace overlay_src {

struct DetectorCombination {
  const char *key;
  const char *label;
  const char *cut;
};

struct PlotVar {
  const char *name;
  const char *expr;
  const char *xTitle;
  int bins;
  double xmin;
  double xmax;
};

inline double finiteMax(double x, double fallback = 0.) {
  return std::isfinite(x) ? x : fallback;
}

inline void styleTopologyData(TH1D *hist, int idx) {
  static const int kColors[8] = {
      kBlack, kRed + 1, kBlue + 1, kGreen + 2,
      kOrange + 7, kMagenta + 1, kCyan + 2, kGray + 2};
  static const int kMarkers[8] = {20, 21, 22, 23, 24, 25, 26, 32};
  const int safeIdx = (idx >= 0 && idx < 8) ? idx : 0;
  hist->SetStats(false);
  hist->SetLineColor(kColors[safeIdx]);
  hist->SetMarkerColor(kColors[safeIdx]);
  hist->SetMarkerStyle(kMarkers[safeIdx]);
  hist->SetMarkerSize(0.80);
  hist->SetLineWidth(2);
}

inline TString trimCopy(const std::string &s) {
  const std::size_t start = s.find_first_not_of(" \t\n\r");
  if (start == std::string::npos) return "";
  const std::size_t end = s.find_last_not_of(" \t\n\r");
  return TString(s.substr(start, end - start + 1));
}

inline std::vector<TString> splitCsvToTStrings(const char *csv) {
  std::vector<TString> out;
  if (!csv || csv[0] == '\0') return out;
  std::stringstream ss(csv);
  std::string item;
  while (std::getline(ss, item, ',')) {
    TString trimmed = trimCopy(item);
    if (trimmed.Length() > 0) out.push_back(trimmed);
  }
  return out;
}

inline TString buildSelection(const char *eppCut, const char *baseCut,
                              const char *comboCut, const char *pCMyTailCut,
                              const char *q2Cut = nullptr,
                              const char *extraCut = nullptr) {
  std::vector<TString> terms;
  if (eppCut && eppCut[0] != '\0') terms.emplace_back(eppCut);
  if (baseCut && baseCut[0] != '\0') terms.emplace_back(baseCut);
  if (comboCut && comboCut[0] != '\0') terms.emplace_back(comboCut);
  if (pCMyTailCut && pCMyTailCut[0] != '\0') terms.emplace_back(pCMyTailCut);
  if (q2Cut && q2Cut[0] != '\0') terms.emplace_back(q2Cut);
  if (extraCut && extraCut[0] != '\0') terms.emplace_back(extraCut);

  if (terms.empty()) return "1";

  TString sel = "(";
  for (std::size_t i = 0; i < terms.size(); ++i) {
    if (i > 0) sel += " && ";
    sel += "(";
    sel += terms[i];
    sel += ")";
  }
  sel += ")";
  return sel;
}

inline TString applyWeightToSelection(const char *weightExpression,
                                      const TString &selection) {
  if (!weightExpression || weightExpression[0] == '\0') {
    return selection;
  }
  return Form("(%s) * (%s)", weightExpression, selection.Data());
}

inline int chooseAutoRebinFactor(int nBins, double representativeEntries,
                                 double targetAvgCountsPerBin = 8.0,
                                 int minFinalBins = 5,
                                 int maxRebinFactor = 2) {
  if (nBins <= 0 || !std::isfinite(representativeEntries) ||
      representativeEntries <= 0. || maxRebinFactor < 1) {
    return 1;
  }

  int bestMeetingTarget = 1;
  for (int factor = 2; factor <= nBins; ++factor) {
    if (factor > maxRebinFactor) break;
    if (nBins % factor != 0) continue;
    const int finalBins = nBins / factor;
    if (finalBins < minFinalBins) continue;
    const double avg = representativeEntries / static_cast<double>(finalBins);
    if (avg >= targetAvgCountsPerBin) {
      bestMeetingTarget = factor;
    }
  }
  return bestMeetingTarget;
}

inline std::vector<DetectorCombination> detectorCombinations(bool includeFdFd) {
  std::vector<DetectorCombination> all = {
      {"lead_fd_rec_fd", "Lead FD Rec FD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_fd_rec_cd", "Lead FD Rec CD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()>45.)"},
      {"lead_cd_rec_fd", "Lead CD Rec FD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_cd_rec_cd", "Lead CD Rec CD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()>45.)"}};

  if (includeFdFd) return all;
  return std::vector<DetectorCombination>(all.begin() + 1, all.end());
}

inline std::vector<PlotVar> defaultVariables() {
  return {
      {"xB", "xB", "x_{B}", 40, 1., 2.},
      {"Q2", "Q2", "Q^{2} [GeV^{2}]", 40, 1.0, 5.0},
      {"omega", "omega", "#omega [GeV]", 40, 0.5, 3.0},
      {"eP", "eP", "p_{e} [GeV]", 40, 3.0, 6.0},
      {"eTheta", "eTheta*180./TMath::Pi()", "#theta_{e} [deg]", 40, 9.0, 35.0},
      {"pLead", "leadP", "p_{lead} [GeV/c]", 40, 0.8, 3.0},
      {"recoilP", "recP", "p_{recoil} [GeV/c]", 40, 0.2, 1.2},
      {"pMiss", "pMiss", "p_{miss} [GeV/c]", 40, 0.3, 1.1},
      {"kMiss", "kMiss", "k_{miss} [GeV/c]", 40, 0.2, 1.2},
      {"mMiss", "mMiss", "m_{miss} [GeV/c^{2}]", 50, 0., 2.},
      {"EMiss", "EMiss", "E_{miss} [GeV]", 50, 0., 1.2},
      {"E0miss", "E0miss", "E_{0,miss} [GeV]", 40, 0.0, 0.5},
      {"E1miss", "E1miss", "E_{1,miss} [GeV]", 50, -0.3, 0.2},
      {"E2miss", "E2miss", "E_{2,miss} [GeV]", 50, -0.5, 0.5},
      {"leadTheta", "leadTheta*180./TMath::Pi()", "#theta_{lead} [deg]", 36,
       0.0, 110.0},
      {"recTheta", "recTheta*180./TMath::Pi()", "#theta_{rec} [deg]", 36, 0.0,
       180.0},
      {"pMissTheta", "pMissTheta*180./TMath::Pi()", "#theta_{pmiss} [deg]", 36,
       0.0, 180.0},
      {"theta_PmPrec", "theta_PmPrec*180./TMath::Pi()",
       "#theta_{p_{miss},p_{rec}} [deg]", 36, 0.0, 180.0},
    //   {"pMissMinusRecoilP", "pMiss-recP", "|p_{miss}| - |p_{recoil}| [GeV/c]",
    //    40, -0.8, 0.8},
    //   {"piMinusTheta_PmPrec", "(TMath::Pi()-theta_PmPrec)*180./TMath::Pi()",
    //    "#pi - #theta_{p_{miss},p_{rec}} [deg]", 36, 0.0, 180.0},
      {"theta_PleadPrec", "theta_PleadPrec*180./TMath::Pi()",
       "#theta_{p_{lead},p_{rec}} [deg]", 36, 0.0, 180.0},
      {"theta_PrecQ", "theta_PrecQ*180./TMath::Pi()", "#theta_{p_{rec},q} [deg]",
       36, 0.0, 180.0},
    //   {"pRelTheta", "pRelTheta*180./TMath::Pi()", "#theta_{p_{rel}} [deg]", 36,
    //    0.0, 180.0},
      {"pRel", "pRel", "p_{rel} [GeV/c]", 40, 0.0, 1.2},
      {"pRelx", "pRelx", "p_{rel,x} [GeV/c]", 40, -0.5, 0.5},
      {"pRely", "pRely", "p_{rel,y} [GeV/c]", 40, -0.5, 0.5},
      {"pRelz", "pRelz", "p_{rel,z} [GeV/c]", 40, 0., 1.},
    //   {"pLeadPlusRecOver2", "pLeadPlusRecOver2", "|p_{lead}+p_{rec}|/2 [GeV/c]",
    //    40, 0.0, 1.6},
    //   {"pLeadMinusRec", "pLeadMinusRec", "|p_{lead}-p_{rec}| [GeV/c]", 40, 0.0,
    //    2.4},
      {"pCM", "pCM", "p_{CM} [GeV/c]", 40, 0.0, 1.2},
      {"pCMx", "pCMx", "p_{CM,x} [GeV/c]", 40, -0.8, 0.8},
      {"pCMy", "pCMy", "p_{CM,y} [GeV/c]", 40, -0.8, 0.8},
      {"pCMz", "pCMz", "p_{CM,z} [GeV/c]", 40, -0.8, 0.8},
      {"pcmx_lab", "pcmx_lab", "p_{CM,x}^{lab} [GeV/c]", 40, -0.8, 0.8},
      {"pcmy_lab", "pcmy_lab", "p_{CM,y}^{lab} [GeV/c]", 40, -0.8, 0.8},
      {"pcmz_lab", "pcmz_lab", "p_{CM,z}^{lab} [GeV/c]", 40, -0.8, 0.8},
      {"alpha_1", "alpha_1", "#alpha_{1}", 40, 1.0, 1.8},
      {"alpha_2", "alpha_2", "#alpha_{2}", 40, 0.2, 1.5},
      {"alpha_CM", "alpha_CM", "#alpha_{CM}", 40, 1.5, 3.},
      {"alpha_rel", "alpha_rel", "#alpha_{rel}", 40, 0.4, 1.5},
      {"p1_plus", "p1_plus", "p_{1}^{+} [GeV]", 40, -0.5, 2.5},
      {"p2_plus", "p2_plus", "p_{2}^{+} [GeV]", 40, -0.5, 2.5},
      {"p1_perp_x", "p1_perp_x", "p_{1,T,x} [GeV/c]", 40, -0.8, 0.8},
      {"p1_perp_y", "p1_perp_y", "p_{1,T,y} [GeV/c]", 40, -0.8, 0.8},
      {"p1_perp_mag", "p1_perp_mag", "|p_{1,T}| [GeV/c]", 40, 0.0, 1.2},
      {"p2_perp_x", "p2_perp_x", "p_{2,T,x} [GeV/c]", 40, -0.8, 0.8},
      {"p2_perp_y", "p2_perp_y", "p_{2,T,y} [GeV/c]", 40, -0.8, 0.8},
      {"p2_perp_mag", "p2_perp_mag", "|p_{2,T}| [GeV/c]", 40, 0.0, 1.2},
      {"pCM_perp_x", "pCM_perp_x", "p_{CM,T,x} [GeV/c]", 40, -0.8, 0.8},
      {"pCM_perp_y", "pCM_perp_y", "p_{CM,T,y} [GeV/c]", 40, -0.8, 0.8},
      {"pCM_perp_mag", "pCM_perp_mag", "|p_{CM,T}| [GeV/c]", 40, 0.0, 1.2},
      {"prel_perp_x", "prel_perp_x", "p_{rel,T,x} [GeV/c]", 40, -0.8, 0.8},
      {"prel_perp_y", "prel_perp_y", "p_{rel,T,y} [GeV/c]", 40, -0.8, 0.8},
      {"prel_perp_mag", "prel_perp_mag", "|p_{rel,T}| [GeV/c]", 40, 0.0, 1.2},
      {"k_lc", "k", "k [GeV/c]", 40, 0.0, 1.2},
      {"k2_lc", "k2", "k^{2} [GeV^{2}]", 40, 0.0, 1.5},
      {"k_z", "k_z", "k_{z} [GeV/c]", 40, -0.8, 0.8},
      {"m_bar", "m_bar", "#bar{m} [GeV/c^{2}]", 20, 0.8, 1.1},
    //   {"phiCM", "TMath::ATan2(pCMy,pCMx)*TMath::RadToDeg()", "#phi_{CM} [deg]", 36,
    //    -180.0, 180.0},
    //   {"pCMTheta", "TMath::ACos(pCMz/pCM)*180./TMath::Pi()", "#theta_{p_{CM}} [deg]",
    //    25, 0.0, 180.0},
    //   {"pCMThetaLab",
    //    "TMath::ACos(pcmz_lab/TMath::Sqrt(pcmx_lab*pcmx_lab + "
    //    "pcmy_lab*pcmy_lab + pcmz_lab*pcmz_lab))*180./TMath::Pi()",
    //    "#theta_{p_{CM}}^{lab} [deg]", 36, 0.0, 180.0},
      {"theta_PCMlabPrel",
       "((pRel > 0.) && "
       "(TMath::Sqrt(pcmx_lab*pcmx_lab + pcmy_lab*pcmy_lab + "
       "pcmz_lab*pcmz_lab) > 0.)) ? "
       "(TMath::ACos(TMath::Max(-1.0, TMath::Min(1.0, "
       "(pcmx_lab*(pRel*TMath::Sin(pRelTheta)*TMath::Cos(pRelPhi)) + "
       " pcmy_lab*(pRel*TMath::Sin(pRelTheta)*TMath::Sin(pRelPhi)) + "
       " pcmz_lab*(pRel*TMath::Cos(pRelTheta))) / "
       "(pRel*TMath::Sqrt(pcmx_lab*pcmx_lab + pcmy_lab*pcmy_lab + "
       "pcmz_lab*pcmz_lab))))))*180./TMath::Pi() : -999.",
       "#theta_{p_{CM}^{lab},p_{rel}} [deg]", 36, 0.0, 180.0},
    //   {"phiTrento", "phiTrento", "#phi_{Trento} [rad]", 36, 0.0, TMath::Pi()},
      {"theta_PleadQ", "theta_PleadQ*180./TMath::Pi()", "#theta_{p_{lead},q} [deg]",
       36, 0.0, 60.0},
      {"theta_PmPlead", "theta_PmPlead*180./TMath::Pi()",
       "#theta_{p_{miss},p_{lead}} [deg]", 36, 0.0, 180.0},
      {"theta_PmQ", "theta_PmQ*180./TMath::Pi()", "#theta_{p_{miss},q} [deg]", 36,
       90.0, 180.0},
       {"lead_q_ratio", "leadP/qP", "|\vec{p}_{lead}|/|\vec{q}|", 40, 0.5, 1.2},
      {"qP", "qP", "|q| [GeV/c]", 40, 1., 4.5},
      {"qTheta", "qTheta*180./TMath::Pi()", "#theta_{q} [deg]", 36, 20.0, 65.0}};
}

}  // namespace overlay_src

#endif
