// overlayBasicSrcFDCD.C
//
// Draw data/simulation overlays for basic SRC observables in the four
// lead/recoil detector combinations:
//   lead FD + rec FD, lead FD + rec CD, lead CD + rec FD, lead CD + rec CD.
//
// Usage:
//   root -l -b -q 'myPlots/scratch/overlayBasicSrcFDCD.C()'
//   root -l -b -q 'myPlots/scratch/overlayBasicSrcFDCD.C("data.root","sim.root")'
//   root -l -b -q 'myPlots/scratch/overlayBasicSrcFDCD.C("~/data/RGM_DATA/c12_src_skim.root","~/data/RGM_DATA/c12_sim_skim.root","srcTree","overlay_pCMy_gt03.pdf",true,"pCM>0","(weight_epp*(weight_epp<200))","pCMy>0.3")'
//   root -l -b -q 'myPlots/scratch/overlayBasicSrcFDCD.C("~/data/RGM_DATA/c12_src_skim.root","~/data/RGM_DATA/c12_sim_skim.root","srcTree","overlay_pCMy_lt03.pdf",true,"pCM>0","(weight_epp*(weight_epp<200))","pCMy<0.3")'
//
// Notes:
//   - e'pp-only selection is enforced by eppCut (default: pCM > 0).
//   - baseCut is an additional optional selection layered on top.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

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

struct SummaryStats {
  bool valid = false;
  double nEffData = std::numeric_limits<double>::quiet_NaN();
  double nEffSim = std::numeric_limits<double>::quiet_NaN();
  double meanData = std::numeric_limits<double>::quiet_NaN();
  double meanSim = std::numeric_limits<double>::quiet_NaN();
  double meanErrData = std::numeric_limits<double>::quiet_NaN();
  double meanErrSim = std::numeric_limits<double>::quiet_NaN();
  double medData = std::numeric_limits<double>::quiet_NaN();
  double medSim = std::numeric_limits<double>::quiet_NaN();
  double medErrData = std::numeric_limits<double>::quiet_NaN();
  double medErrSim = std::numeric_limits<double>::quiet_NaN();
  double rmsData = std::numeric_limits<double>::quiet_NaN();
  double rmsSim = std::numeric_limits<double>::quiet_NaN();
};

void styleData(TH1D *hist) {
  hist->SetStats(false);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.75);
  hist->SetLineWidth(2);
}

void styleSim(TH1D *hist) {
  hist->SetStats(false);
  hist->SetLineColor(kRed + 1);
  hist->SetMarkerColor(kRed + 1);
  hist->SetMarkerStyle(24);
  hist->SetMarkerSize(0.75);
  hist->SetLineWidth(2);
}

void styleThird(TH1D *hist) {
  hist->SetStats(false);
  hist->SetLineColor(kBlue + 1);
  hist->SetMarkerColor(kBlue + 1);
  hist->SetMarkerStyle(25);
  hist->SetMarkerSize(0.75);
  hist->SetLineWidth(2);
}

void styleTopologyData(TH1D *hist, int idx) {
  static const int kColors[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
  static const int kMarkers[4] = {20, 21, 22, 23};
  const int safeIdx = (idx >= 0 && idx < 4) ? idx : 0;
  hist->SetStats(false);
  hist->SetLineColor(kColors[safeIdx]);
  hist->SetMarkerColor(kColors[safeIdx]);
  hist->SetMarkerStyle(kMarkers[safeIdx]);
  hist->SetMarkerSize(0.80);
  hist->SetLineWidth(2);
}

double finiteMax(double x, double fallback = 0.) {
  return std::isfinite(x) ? x : fallback;
}

}  // namespace

void overlayBasicSrcFDCD(
    const char *dataFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *simFileName = "~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "overlay_basic_src_fdcd_byDet.pdf",
    bool normalizeToUnity = true,
  const char *eppCut = "pCM > 0",
  const char *baseCut = "pCM > 0 && pMiss < 1. && recP < 1. && recP > .4",
  const char *simWeightExpression = "(weight_epp)",
  const char *pCMyTailCut = "",
  bool overlayByVariableDataOnly = false,
  const char *thirdFileName = "",
  const char *thirdWeightExpression = "(weight_epp)",
  const char *thirdLegendLabel = "Sim + FSI") {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile *dataFile = TFile::Open(dataFileName, "READ");
  TFile *simFile = nullptr;
  TFile *thirdFile = nullptr;
  const bool useThirdFile = thirdFileName && thirdFileName[0] != '\0';
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Could not open data file "
              << dataFileName << '\n';
    return;
  }
  if (!overlayByVariableDataOnly) {
    simFile = TFile::Open(simFileName, "READ");
    if (!simFile || simFile->IsZombie()) {
      std::cerr << "Could not open simulation file "
                << simFileName << '\n';
      dataFile->Close();
      return;
    }
    if (useThirdFile) {
      thirdFile = TFile::Open(thirdFileName, "READ");
      if (!thirdFile || thirdFile->IsZombie()) {
        std::cerr << "Could not open third file "
                  << thirdFileName << '\n';
        dataFile->Close();
        simFile->Close();
        return;
      }
    }
  }

  TTree *dataTree = dynamic_cast<TTree *>(dataFile->Get(treeName));
  TTree *simTree = nullptr;
  TTree *thirdTree = nullptr;
  if (!dataTree) {
    std::cerr << "Could not find tree \"" << treeName
              << "\" in data file.\n";
    dataFile->Close();
    if (simFile) simFile->Close();
    if (thirdFile) thirdFile->Close();
    return;
  }
  if (!overlayByVariableDataOnly) {
    simTree = dynamic_cast<TTree *>(simFile->Get(treeName));
    if (!simTree) {
      std::cerr << "Could not find tree \"" << treeName
                << "\" in simulation file.\n";
      dataFile->Close();
      simFile->Close();
      if (thirdFile) thirdFile->Close();
      return;
    }
    if (useThirdFile) {
      thirdTree = dynamic_cast<TTree *>(thirdFile->Get(treeName));
      if (!thirdTree) {
        std::cerr << "Could not find tree \"" << treeName
                  << "\" in third file.\n";
        dataFile->Close();
        simFile->Close();
        thirdFile->Close();
        return;
      }
    }
  }

  const DetectorCombination combinations[] = {
      {"lead_fd_rec_fd", "Lead FD Rec FD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_fd_rec_cd", "Lead FD Rec CD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()>45.)"},
      {"lead_cd_rec_fd", "Lead CD Rec FD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_cd_rec_cd", "Lead CD Rec CD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()>45.)"}};

  const std::vector<PlotVar> variables = {
      {"xB", "xB", "x_{B}", 40, 1., 2.},
      {"Q2", "Q2", "Q^{2} [GeV^{2}]", 40, 1.0, 5.0},
      {"pLead", "leadP", "p_{lead} [GeV/c]", 40, 0.8, 3.0},
      {"recoilP", "recP", "p_{recoil} [GeV/c]", 40, 0.2, 1.2},
      {"pMiss", "pMiss", "p_{miss} [GeV/c]", 40, 0.3, 1.1},
      {"kMiss", "kMiss", "k_{miss} [GeV/c]", 40, 0.2, 1.2},
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
      {"pMissMinusRecoilP", "pMiss-recP", "|p_{miss}| - |p_{recoil}| [GeV/c]",
       40, -0.8, 0.8},
      {"piMinusTheta_PmPrec",
       "(TMath::Pi()-theta_PmPrec)*180./TMath::Pi()",
       "#pi - #theta_{p_{miss},p_{rec}} [deg]", 36, 0.0, 180.0},
      {"theta_PleadPrec", "theta_PleadPrec*180./TMath::Pi()",
       "#theta_{p_{lead},p_{rec}} [deg]", 36, 0.0, 180.0},
      {"theta_PrecQ", "theta_PrecQ*180./TMath::Pi()",
       "#theta_{p_{rec},q} [deg]", 36, 0.0, 180.0},
      {"pRelTheta", "pRelTheta*180./TMath::Pi()", "#theta_{p_{rel}} [deg]",
       36, 0.0, 180.0},
      {"pRel", "pRel", "p_{rel} [GeV/c]", 40, 0.0, 1.2},
      {"pRelx", "pRelx", "p_{rel,x} [GeV/c]", 40, -0.8, 0.8},
      {"pRely", "pRely", "p_{rel,y} [GeV/c]", 40, -0.8, 0.8},
      {"pRelz", "pRelz", "p_{rel,z} [GeV/c]", 40, -0.8, 0.8},
      {"pLeadPlusRecOver2", "pLeadPlusRecOver2",
       "|p_{lead}+p_{rec}|/2 [GeV/c]", 40, 0.0, 1.6},
      {"pLeadMinusRec", "pLeadMinusRec",
       "|p_{lead}-p_{rec}| [GeV/c]", 40, 0.0, 2.4},
      {"pCM", "pCM", "p_{CM} [GeV/c]", 40, 0.0, 1.2},
      {"pCMx", "pCMx", "p_{CM,x} [GeV/c]", 40, -0.8, 0.8},
      {"pCMy", "pCMy", "p_{CM,y} [GeV/c]", 40, -0.8, 0.8},
      {"pCMz", "pCMz", "p_{CM,z} [GeV/c]", 40, -0.8, 0.8},
      {"phiCM", "TMath::ATan2(pCMy,pCMx)*TMath::RadToDeg()",
       "#phi_{CM} [deg]", 36, -180.0, 180.0},
      {"phiFold",
       "TMath::Min(TMath::Abs(TMath::ATan2(pCMy,pCMx)*TMath::RadToDeg()),"
       "180.0-TMath::Abs(TMath::ATan2(pCMy,pCMx)*TMath::RadToDeg()))",
       "#phi_{fold} [deg]", 25, 0.0, 90.0},
      {"pCMTheta", "TMath::ACos(pCMz/pCM)*180./TMath::Pi()",
       "#theta_{p_{CM}} [deg]", 25, 0.0, 180.0},
      {"pCMThetaLab",
       "TMath::ACos(pcmz_lab/TMath::Sqrt(pcmx_lab*pcmx_lab + pcmy_lab*pcmy_lab + pcmz_lab*pcmz_lab))*180./TMath::Pi()",
       "#theta_{p_{CM}}^{lab} [deg]", 36, 0.0, 180.0},
      {"theta_PCMlabPrel",
       "((pRel > 0.) && "
       "(TMath::Sqrt(pcmx_lab*pcmx_lab + pcmy_lab*pcmy_lab + pcmz_lab*pcmz_lab) > 0.)) ? "
       "(TMath::ACos(TMath::Max(-1.0, TMath::Min(1.0, "
       "(pcmx_lab*(pRel*TMath::Sin(pRelTheta)*TMath::Cos(pRelPhi)) + "
       " pcmy_lab*(pRel*TMath::Sin(pRelTheta)*TMath::Sin(pRelPhi)) + "
       " pcmz_lab*(pRel*TMath::Cos(pRelTheta))) / "
       "(pRel*TMath::Sqrt(pcmx_lab*pcmx_lab + pcmy_lab*pcmy_lab + pcmz_lab*pcmz_lab))))))"
       "*180./TMath::Pi() : -999.",
       "#theta_{p_{CM}^{lab},p_{rel}} [deg]", 36, 0.0, 180.0},
      {"phiTrento", "phiTrento", "#phi_{Trento} [rad]",
       36, 0.0, TMath::Pi()},
      {"theta_PleadQ", "theta_PleadQ*180./TMath::Pi()",
       "#theta_{p_{lead},q} [deg]", 36, 0.0, 180.0},
      {"theta_PmPlead", "theta_PmPlead*180./TMath::Pi()",
       "#theta_{p_{miss},p_{lead}} [deg]", 36, 0.0, 180.0},
      {"theta_PmQ", "theta_PmQ*180./TMath::Pi()", "#theta_{p_{miss},q} [deg]",
       36, 90, 180.0},
      {"qP", "qP", "|q| [GeV/c]", 40, 1., 4.5},
      {"qTheta", "qTheta*180./TMath::Pi()", "#theta_{q} [deg]", 36, 0.0,
       60.0}};

  TCanvas *canvas = new TCanvas("c_overlay_basic_src_fdcd",
                                 "SRC data/sim overlays by detector topology",
                                 1300, 950);
  canvas->SetLeftMargin(0.11);
  canvas->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  canvas->Print(pdfName + "[", "pdf");

  SummaryStats recPStats[4];
  SummaryStats pMissStats[4];
  SummaryStats pMissMinusRecPStats[4];
  SummaryStats pRelStats[4];
  SummaryStats pCMxStats[4];
  SummaryStats pCMyStats[4];
  SummaryStats pCMzStats[4];

  for (const PlotVar &v : variables) {
    canvas->Clear();

    if (overlayByVariableDataOnly) {
      canvas->cd(1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);

      std::vector<TH1D *> hDataByCombo;
      std::vector<double> rawDataIntegralsByCombo;
      std::vector<int> overlayIndices;
      hDataByCombo.reserve(4);
      rawDataIntegralsByCombo.reserve(4);
      overlayIndices.reserve(4);
      double yMax = 0.;
      double totalRawIntegral = 0.;

      // In data-only overlay mode, suppress FD+FD as requested.
      for (int idx = 0; idx < 4; ++idx) {
        if (idx == 0) continue;
        overlayIndices.push_back(idx);
      }

      for (const int idx : overlayIndices) {
        const DetectorCombination &combo = combinations[idx];
        TH1D *hData =
            new TH1D(Form("h_data_overlay_%s_%s", v.name, combo.key), "",
                     v.bins, v.xmin, v.xmax);
        hData->Sumw2();
        styleTopologyData(hData, idx);

        const bool usePCMyTailCut = pCMyTailCut && pCMyTailCut[0] != '\0';
        const TString boolSel =
            usePCMyTailCut
                ? Form("((%s) && (%s) && (%s) && (%s))", eppCut, baseCut,
                       combo.cut, pCMyTailCut)
                : Form("((%s) && (%s) && (%s))", eppCut, baseCut, combo.cut);

        dataTree->Project(hData->GetName(), v.expr, boolSel.Data());
        const double rawDataIntegral =
            hData->Integral(0, hData->GetNbinsX() + 1);
        totalRawIntegral += rawDataIntegral;
        rawDataIntegralsByCombo.push_back(rawDataIntegral);

        if (normalizeToUnity) {
          const double dataInt = hData->Integral();
          if (dataInt > 0.) hData->Scale(1. / dataInt);
        }

        yMax = std::max(yMax, finiteMax(hData->GetMaximum()));
        hDataByCombo.push_back(hData);
      }

      if (hDataByCombo.empty()) {
        continue;
      }

      if (yMax <= 0.) yMax = 1.;
      hDataByCombo.front()->SetMaximum(1.3 * yMax);
      hDataByCombo.front()->GetXaxis()->SetTitle(v.xTitle);
      hDataByCombo.front()->GetYaxis()->SetTitle(normalizeToUnity
                                                     ? "Normalized counts"
                                                     : "Counts");
      hDataByCombo.front()->GetXaxis()->SetTitleSize(0.050);
      hDataByCombo.front()->GetYaxis()->SetTitleSize(0.050);
      hDataByCombo.front()->GetXaxis()->SetLabelSize(0.043);
      hDataByCombo.front()->GetYaxis()->SetLabelSize(0.043);
      hDataByCombo.front()->GetYaxis()->SetTitleOffset(1.20);

      hDataByCombo.front()->Draw("E");
      for (std::size_t idx = 1; idx < hDataByCombo.size(); ++idx) {
        hDataByCombo[idx]->Draw("E SAME");
      }

      TLegend *legend = new TLegend(0.60, 0.7, 0.90, 0.89);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.025);
      for (std::size_t drawIdx = 0; drawIdx < hDataByCombo.size(); ++drawIdx) {
        const int comboIdx = overlayIndices[drawIdx];
        const DetectorCombination &combo = combinations[comboIdx];
        const double rawDataIntegral = rawDataIntegralsByCombo[drawIdx];
        legend->AddEntry(hDataByCombo[drawIdx],
                         Form("%s (N=%.0f)", combo.label, rawDataIntegral),
                         "lep");
      }
      legend->Draw();

      TLatex pageTitle;
      pageTitle.SetNDC();
      pageTitle.SetTextSize(0.032);
      pageTitle.SetTextFont(62);
      pageTitle.DrawLatex(0.13, 0.94,
                          Form("%s: data overlays by detector topology", v.name));

      if (totalRawIntegral <= 0.) {
        TLatex empty;
        empty.SetNDC();
        empty.SetTextSize(0.042);
        empty.SetTextColor(kGray + 2);
        empty.DrawLatex(0.20, 0.52, "No entries after cuts");
      }

      canvas->Print(pdfName, "pdf");
      continue;
    }

    canvas->Divide(2, 2);

    for (int idx = 0; idx < 4; ++idx) {
      const DetectorCombination &combo = combinations[idx];
      canvas->cd(idx + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);

      TH1D *hData = new TH1D(Form("h_data_%s_%s", v.name, combo.key), "", v.bins,
                             v.xmin, v.xmax);
      TH1D *hSim = new TH1D(Form("h_sim_%s_%s", v.name, combo.key), "", v.bins,
                            v.xmin, v.xmax);
      TH1D *hThird = useThirdFile
             ? new TH1D(Form("h_third_%s_%s", v.name, combo.key), "",
                v.bins, v.xmin, v.xmax)
             : nullptr;
      hData->Sumw2();
      hSim->Sumw2();
      if (hThird) hThird->Sumw2();
      styleData(hData);
      styleSim(hSim);
      if (hThird) styleThird(hThird);

      const bool usePCMyTailCut = pCMyTailCut && pCMyTailCut[0] != '\0';
      const TString boolSel = usePCMyTailCut
                    ? Form("((%s) && (%s) && (%s) && (%s))",
                      eppCut, baseCut, combo.cut,
                      pCMyTailCut)
                    : Form("((%s) && (%s) && (%s))", eppCut,
                      baseCut, combo.cut);
      const TString dataSel = boolSel;
      const TString simSel =
          Form("(%s) * (%s)", simWeightExpression, boolSel.Data());
        const TString thirdSel =
          Form("(%s) * (%s)", thirdWeightExpression, boolSel.Data());

      const Long64_t nDataProjected =
          dataTree->Project(hData->GetName(), v.expr, dataSel.Data());
      const Long64_t nSimProjected =
          simTree->Project(hSim->GetName(), v.expr, simSel.Data());
        const Long64_t nThirdProjected =
          hThird ? thirdTree->Project(hThird->GetName(), v.expr, thirdSel.Data())
             : 0;

      const double rawDataIntegral = hData->Integral(0, hData->GetNbinsX() + 1);
      const double rawSimIntegral = hSim->Integral(0, hSim->GetNbinsX() + 1);
        const double rawThirdIntegral =
          hThird ? hThird->Integral(0, hThird->GetNbinsX() + 1) : 0.;

      auto fillSummaryStats = [&](SummaryStats &stats) {
        stats.valid = true;
        const double nEffData = hData->GetEffectiveEntries();
        const double nEffSim = hSim->GetEffectiveEntries();
        const double sigmaData = hData->GetStdDev();
        const double sigmaSim = hSim->GetStdDev();
        const double medianScale = 1.2533141373155;  // approx. normal-theory SE factor

        stats.nEffData = nEffData;
        stats.nEffSim = nEffSim;
        stats.meanData = hData->GetMean();
        stats.meanSim = hSim->GetMean();
        stats.meanErrData =
          (std::isfinite(sigmaData) && nEffData > 0.)
            ? (sigmaData / std::sqrt(nEffData))
            : std::numeric_limits<double>::quiet_NaN();
        stats.meanErrSim =
          (std::isfinite(sigmaSim) && nEffSim > 0.)
            ? (sigmaSim / std::sqrt(nEffSim))
            : std::numeric_limits<double>::quiet_NaN();
        double q = 0.5;
        double medData = std::numeric_limits<double>::quiet_NaN();
        double medSim = std::numeric_limits<double>::quiet_NaN();
        hData->GetQuantiles(1, &medData, &q);
        hSim->GetQuantiles(1, &medSim, &q);
        stats.medData = medData;
        stats.medSim = medSim;
        stats.medErrData =
          (std::isfinite(sigmaData) && nEffData > 0.)
            ? (medianScale * sigmaData / std::sqrt(nEffData))
            : std::numeric_limits<double>::quiet_NaN();
        stats.medErrSim =
          (std::isfinite(sigmaSim) && nEffSim > 0.)
            ? (medianScale * sigmaSim / std::sqrt(nEffSim))
            : std::numeric_limits<double>::quiet_NaN();
        stats.rmsData = hData->GetRMS();
        stats.rmsSim = hSim->GetRMS();
      };
      if (!std::strcmp(v.name, "recoilP")) {
        fillSummaryStats(recPStats[idx]);
      } else if (!std::strcmp(v.name, "pMiss")) {
        fillSummaryStats(pMissStats[idx]);
      } else if (!std::strcmp(v.name, "pMissMinusRecoilP")) {
        fillSummaryStats(pMissMinusRecPStats[idx]);
      } else if (!std::strcmp(v.name, "pRel")) {
        fillSummaryStats(pRelStats[idx]);
      } else if (!std::strcmp(v.name, "pCMx")) {
        fillSummaryStats(pCMxStats[idx]);
      } else if (!std::strcmp(v.name, "pCMy")) {
        fillSummaryStats(pCMyStats[idx]);
      } else if (!std::strcmp(v.name, "pCMz")) {
        fillSummaryStats(pCMzStats[idx]);
      }

      if (!std::strcmp(v.name, "xB") || !std::strcmp(v.name, "phiTrento")) {
        const Long64_t nDataSelected = dataTree->GetEntries(boolSel.Data());
        const Long64_t nSimSelected = simTree->GetEntries(boolSel.Data());
        std::cout << "" << combo.key
                  << " boolSel=" << boolSel
                  << " | data selected=" << nDataSelected
                  << " | sim selected=" << nSimSelected
                  << " | data projected=" << nDataProjected
                  << " | sim projected=" << nSimProjected
                  << " | third projected=" << nThirdProjected
                  << " | hData integral(all bins)=" << rawDataIntegral
                  << " | hSim integral(all bins)=" << rawSimIntegral
                  << " | hThird integral(all bins)=" << rawThirdIntegral
                  << '\n';
      }


      if (normalizeToUnity) {
        const double dataInt = hData->Integral();
        const double simInt = hSim->Integral();
        const double thirdInt = hThird ? hThird->Integral() : 0.;
        if (dataInt > 0.) hData->Scale(1. / dataInt);
        if (simInt > 0.) hSim->Scale(1. / simInt);
        if (hThird && thirdInt > 0.) hThird->Scale(1. / thirdInt);
      }

      double yMax = std::max(finiteMax(hData->GetMaximum()),
                             finiteMax(hSim->GetMaximum()));
      if (hThird) yMax = std::max(yMax, finiteMax(hThird->GetMaximum()));
      if (yMax <= 0.) yMax = 1.;
      hData->SetMaximum(1.25 * yMax);

      hData->GetXaxis()->SetTitle(v.xTitle);
      hData->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts"
                                                   : "Counts (sim weighted)");
      hData->GetXaxis()->SetTitleSize(0.050);
      hData->GetYaxis()->SetTitleSize(0.050);
      hData->GetXaxis()->SetLabelSize(0.043);
      hData->GetYaxis()->SetLabelSize(0.043);
      hData->GetYaxis()->SetTitleOffset(1.20);

    //   hData->Rebin();
    //   hSim->Rebin();

      hData->Draw("E");
      hSim->Draw("E SAME");
  if (hThird) hThird->Draw("E SAME");

      TLegend *legend = new TLegend(0.52, 0.68, 0.90, 0.89);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.035);
      if (normalizeToUnity) {
        legend->AddEntry(hData, Form("Data c12 (N=%.0f)", rawDataIntegral), "lep");
        legend->AddEntry(hSim, Form("Sim c12  (sumW=%.1f)", rawSimIntegral),
      //  legend->AddEntry(hSim, Form("Sim weighted (sumW=%.1f)", rawSimIntegral),
                         "lep");
           if (hThird) {
          legend->AddEntry(hThird,
                  Form("%s (sumW=%.1f)", thirdLegendLabel,
                    rawThirdIntegral),
                  "lep");
           }
      } else {
        legend->AddEntry(hData, Form("Data (N=%.0f)", rawDataIntegral), "lep");
        legend->AddEntry(hSim, Form("Sim weighted (sumW=%.1f)", rawSimIntegral),
                         "lep");
           if (hThird) {
          legend->AddEntry(hThird,
                  Form("%s (sumW=%.1f)", thirdLegendLabel,
                    rawThirdIntegral),
                  "lep");
           }
      }
      legend->Draw();

      TLatex label;
      label.SetNDC();
      label.SetTextSize(0.040);
      label.SetTextFont(42);
      label.DrawLatex(0.14, 0.92, combo.label);

      if (rawDataIntegral <= 0. && rawSimIntegral <= 0. && rawThirdIntegral <= 0.) {
        TLatex empty;
        empty.SetNDC();
        empty.SetTextSize(0.042);
        empty.SetTextColor(kGray + 2);
        empty.DrawLatex(0.20, 0.52, "No entries after cuts");
      }
    }

    canvas->cd(0);
    TLatex pageTitle;
    pageTitle.SetNDC();
    pageTitle.SetTextSize(0.028);
    pageTitle.SetTextFont(62);
    pageTitle.DrawLatex(0.10, 0.985,
                        Form("%s: data vs simulation in detector topologies",
                             v.name));

    canvas->Print(pdfName, "pdf");
  }

  canvas->Print(pdfName + "]", "pdf");

  std::cout << "Wrote " << outputPdfName << '\n';
  std::cout << "eppCut = " << eppCut << '\n';
  std::cout << "baseCut = " << baseCut << '\n';
  std::cout << "simWeightExpression = "
            << simWeightExpression << '\n';
  std::cout << "thirdFileName = "
            << (useThirdFile ? thirdFileName : "<none>") << '\n';
  std::cout << "thirdWeightExpression = "
            << (useThirdFile ? thirdWeightExpression : "<none>") << '\n';
  std::cout << "pCMyTailCut = "
            << ((pCMyTailCut && pCMyTailCut[0] != '\0') ? pCMyTailCut
                                                         : "<none>")
            << '\n';

  if (!overlayByVariableDataOnly) {
    TString summaryFileName(outputPdfName);
    if (summaryFileName.EndsWith(".pdf")) {
      summaryFileName.ReplaceAll(".pdf", "_summary_stats.txt");
    } else {
      summaryFileName += "_summary_stats.txt";
    }

    std::ofstream summaryOut(summaryFileName.Data());
    if (!summaryOut) {
      std::cerr << "Could not open summary output file "
                << summaryFileName << '\n';
      dataFile->Close();
      if (simFile) simFile->Close();
      return;
    }

    summaryOut << "Per-topology means/medians (data vs sim weighted)\n";
    summaryOut << "Columns per variable: "
                  "mean_data +/- err / mean_sim +/- err / median_data +/- err / median_sim +/- err\n";
    summaryOut << "Medians are histogram quantiles from the plotted bins.\n";
    summaryOut << "Mean(sim) errors use sigma/sqrt(Neff), "
                  "Neff=(sumW^2/sumW2) from weighted histogram entries.\n";
    summaryOut << "Added RMS outputs for pCMx, pCMy, pCMz (data and sim).\n";
    summaryOut << std::fixed << std::setprecision(4);
    for (int idx = 0; idx < 4; ++idx) {
      const DetectorCombination &combo = combinations[idx];
      summaryOut << "" << combo.key
                 << " | Neff(data/sim): " << recPStats[idx].nEffData << " / "
                 << recPStats[idx].nEffSim
                 << " | <pRecoil>: " << recPStats[idx].meanData << " +/- "
                 << recPStats[idx].meanErrData << " / "
                 << recPStats[idx].meanSim << " +/- "
                 << recPStats[idx].meanErrSim << " / " << recPStats[idx].medData
                 << " +/- " << recPStats[idx].medErrData << " / "
                 << recPStats[idx].medSim << " +/- " << recPStats[idx].medErrSim
                 << " | <pMiss>: " << pMissStats[idx].meanData << " +/- "
                 << pMissStats[idx].meanErrData << " / "
                 << pMissStats[idx].meanSim << " +/- "
                 << pMissStats[idx].meanErrSim << " / " << pMissStats[idx].medData
                 << " +/- " << pMissStats[idx].medErrData << " / "
                 << pMissStats[idx].medSim << " +/- " << pMissStats[idx].medErrSim
                 << " | <|pMiss|-|pRecoil|>: "
                 << pMissMinusRecPStats[idx].meanData << " +/- "
                 << pMissMinusRecPStats[idx].meanErrData << " / "
                 << pMissMinusRecPStats[idx].meanSim << " +/- "
                 << pMissMinusRecPStats[idx].meanErrSim << " / "
                 << pMissMinusRecPStats[idx].medData << " +/- "
                 << pMissMinusRecPStats[idx].medErrData << " / "
                 << pMissMinusRecPStats[idx].medSim << " +/- "
                 << pMissMinusRecPStats[idx].medErrSim
                 << " | <|pRel|>: " << pRelStats[idx].meanData << " +/- "
                 << pRelStats[idx].meanErrData << " / "
                 << pRelStats[idx].meanSim << " +/- "
                 << pRelStats[idx].meanErrSim << " / " << pRelStats[idx].medData
                 << " +/- " << pRelStats[idx].medErrData << " / "
                 << pRelStats[idx].medSim << " +/- " << pRelStats[idx].medErrSim
                 << " | <pCMx>: " << pCMxStats[idx].meanData << " +/- "
                 << pCMxStats[idx].meanErrData << " / "
                 << pCMxStats[idx].meanSim << " +/- "
                 << pCMxStats[idx].meanErrSim << " / " << pCMxStats[idx].medData
                 << " +/- " << pCMxStats[idx].medErrData << " / "
                 << pCMxStats[idx].medSim << " +/- " << pCMxStats[idx].medErrSim
                 << " | <pCMy>: " << pCMyStats[idx].meanData << " +/- "
                 << pCMyStats[idx].meanErrData << " / "
                 << pCMyStats[idx].meanSim << " +/- "
                 << pCMyStats[idx].meanErrSim << " / " << pCMyStats[idx].medData
                 << " +/- " << pCMyStats[idx].medErrData << " / "
                 << pCMyStats[idx].medSim << " +/- " << pCMyStats[idx].medErrSim
                 << " | <pCMz>: " << pCMzStats[idx].meanData << " +/- "
                 << pCMzStats[idx].meanErrData << " / "
                 << pCMzStats[idx].meanSim << " +/- "
                 << pCMzStats[idx].meanErrSim << " / " << pCMzStats[idx].medData
                 << " +/- " << pCMzStats[idx].medErrData << " / "
                 << pCMzStats[idx].medSim << " +/- " << pCMzStats[idx].medErrSim
                 << " | RMS(pCMx,pCMy,pCMz) data: "
                 << pCMxStats[idx].rmsData << ", "
                 << pCMyStats[idx].rmsData << ", "
                 << pCMzStats[idx].rmsData
                 << " | sim: "
                 << pCMxStats[idx].rmsSim << ", "
                 << pCMyStats[idx].rmsSim << ", "
                 << pCMzStats[idx].rmsSim
                 << '\n';
    }
    summaryOut.close();
    std::cout << "Wrote summary stats to "
              << summaryFileName << '\n';
  } else {
    std::cout << "Data-only overlay mode enabled: skipped data/sim summary file.\n";
  }

  dataFile->Close();
  if (simFile) simFile->Close();
  if (thirdFile) thirdFile->Close();
}
