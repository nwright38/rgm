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
#include <iostream>
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

double finiteMax(double x, double fallback = 0.) {
  return std::isfinite(x) ? x : fallback;
}

}  // namespace

void overlayBasicSrcFDCD(
    const char *dataFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *simFileName = "~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "overlay_basic_src_fdcd.pdf",
    bool normalizeToUnity = true,
  const char *eppCut = "pCM > 0",
  const char *baseCut = "pCM > 0 && pMiss < 1.",
  const char *simWeightExpression = "(weight_epp*(weight_epp<200))",
  const char *pCMyTailCut = "") {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile *dataFile = TFile::Open(dataFileName, "READ");
  TFile *simFile = TFile::Open(simFileName, "READ");
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "[overlayBasicSrcFDCD] Could not open data file "
              << dataFileName << '\n';
    if (simFile) simFile->Close();
    return;
  }
  if (!simFile || simFile->IsZombie()) {
    std::cerr << "[overlayBasicSrcFDCD] Could not open simulation file "
              << simFileName << '\n';
    dataFile->Close();
    return;
  }

  TTree *dataTree = dynamic_cast<TTree *>(dataFile->Get(treeName));
  TTree *simTree = dynamic_cast<TTree *>(simFile->Get(treeName));
  if (!dataTree) {
    std::cerr << "[overlayBasicSrcFDCD] Could not find tree \"" << treeName
              << "\" in data file.\n";
    dataFile->Close();
    simFile->Close();
    return;
  }
  if (!simTree) {
    std::cerr << "[overlayBasicSrcFDCD] Could not find tree \"" << treeName
              << "\" in simulation file.\n";
    dataFile->Close();
    simFile->Close();
    return;
  }

  const DetectorCombination combinations[] = {
      {"lead_fd_rec_fd", "lead FD + rec FD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_fd_rec_cd", "lead FD + rec CD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()>45.)"},
      {"lead_cd_rec_fd", "lead CD + rec FD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_cd_rec_cd", "lead CD + rec CD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()>45.)"}};

  const std::vector<PlotVar> variables = {
      {"xB", "xB", "x_{B}", 40, 1., 2.},
      {"Q2", "Q2", "Q^{2} [GeV^{2}]", 40, 1.0, 5.0},
      {"pLead", "leadP", "p_{lead} [GeV/c]", 40, 0.8, 3.0},
      {"recoilP", "recP", "p_{recoil} [GeV/c]", 40, 0.2, 1.2},
      {"pMiss", "pMiss", "p_{miss} [GeV/c]", 40, 0.3, 1.1},
      {"kMiss", "kMiss", "k_{miss} [GeV/c]", 40, 0.2, 1.2},
      {"leadTheta", "leadTheta*180./TMath::Pi()", "#theta_{lead} [deg]", 36,
       0.0, 180.0},
      {"recTheta", "recTheta*180./TMath::Pi()", "#theta_{rec} [deg]", 36, 0.0,
       180.0},
      {"pMissTheta", "pMissTheta*180./TMath::Pi()", "#theta_{pmiss} [deg]", 36,
       0.0, 180.0},
      {"theta_PmPrec", "theta_PmPrec*180./TMath::Pi()",
       "#theta_{p_{miss},p_{rec}} [deg]", 36, 0.0, 180.0},
      {"theta_PleadPrec", "theta_PleadPrec*180./TMath::Pi()",
       "#theta_{p_{lead},p_{rec}} [deg]", 36, 0.0, 180.0},
      {"theta_PrecQ", "theta_PrecQ*180./TMath::Pi()",
       "#theta_{p_{rec},q} [deg]", 36, 0.0, 180.0},
      {"pRelTheta", "pRelTheta*180./TMath::Pi()", "#theta_{p_{rel}} [deg]",
       36, 0.0, 180.0},
      {"pRel", "pRel", "p_{rel} [GeV/c]", 40, 0.0, 1.2},
      {"pCM", "pCM", "p_{CM} [GeV/c]", 40, 0.0, 1.2},
      {"pCMx", "pCMx", "p_{CM,x} [GeV/c]", 40, -0.8, 0.8},
      {"pCMy", "pCMy", "p_{CM,y} [GeV/c]", 40, -0.8, 0.8},
      {"pCMz", "pCMz", "p_{CM,z} [GeV/c]", 40, -0.8, 0.8},
      {"pCMTheta", "TMath::ACos(pCMz/pCM)*180./TMath::Pi()",
       "#theta_{p_{CM}} [deg]", 36, 0.0, 180.0},
      {"pCMThetaLab",
       "TMath::ACos(pcmz_lab/TMath::Sqrt(pcmx_lab*pcmx_lab + pcmy_lab*pcmy_lab + pcmz_lab*pcmz_lab))*180./TMath::Pi()",
       "#theta_{p_{CM}}^{lab} [deg]", 36, 0.0, 180.0},
      {"phiTrento", "phiTrento", "#phi_{Trento} [rad]",
       36, 0.0, TMath::Pi()},
      {"theta_PleadQ", "theta_PleadQ*180./TMath::Pi()",
       "#theta_{p_{lead},q} [deg]", 36, 0.0, 180.0},
      {"theta_PmPlead", "theta_PmPlead*180./TMath::Pi()",
       "#theta_{p_{miss},p_{lead}} [deg]", 36, 0.0, 180.0},
      {"theta_PmQ", "theta_PmQ*180./TMath::Pi()", "#theta_{p_{miss},q} [deg]",
       36, 0.0, 180.0},
      {"qP", "qP", "|q| [GeV/c]", 40, 0.0, 4.5},
      {"qTheta", "qTheta*180./TMath::Pi()", "#theta_{q} [deg]", 36, 0.0,
       60.0}};

  TCanvas *canvas = new TCanvas("c_overlay_basic_src_fdcd",
                                 "SRC data/sim overlays by detector topology",
                                 1300, 950);
  canvas->SetLeftMargin(0.11);
  canvas->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  canvas->Print(pdfName + "[", "pdf");

  for (const PlotVar &v : variables) {
    canvas->Clear();
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
      hData->Sumw2();
      hSim->Sumw2();
      styleData(hData);
      styleSim(hSim);

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

      const Long64_t nDataProjected =
          dataTree->Project(hData->GetName(), v.expr, dataSel.Data());
      const Long64_t nSimProjected =
          simTree->Project(hSim->GetName(), v.expr, simSel.Data());

      const double rawDataIntegral = hData->Integral(0, hData->GetNbinsX() + 1);
      const double rawSimIntegral = hSim->Integral(0, hSim->GetNbinsX() + 1);

      if (!std::strcmp(v.name, "xB") || !std::strcmp(v.name, "phiTrento")) {
        const Long64_t nDataSelected = dataTree->GetEntries(boolSel.Data());
        const Long64_t nSimSelected = simTree->GetEntries(boolSel.Data());
        std::cout << "[overlayBasicSrcFDCD] " << combo.key
                  << " boolSel=" << boolSel
                  << " | data selected=" << nDataSelected
                  << " | sim selected=" << nSimSelected
                  << " | data projected=" << nDataProjected
                  << " | sim projected=" << nSimProjected
                  << " | hData integral(all bins)=" << rawDataIntegral
                  << " | hSim integral(all bins)=" << rawSimIntegral
                  << '\n';
      }


      if (normalizeToUnity) {
        const double dataInt = hData->Integral();
        const double simInt = hSim->Integral();
        if (dataInt > 0.) hData->Scale(1. / dataInt);
        if (simInt > 0.) hSim->Scale(1. / simInt);
      }

      double yMax = std::max(finiteMax(hData->GetMaximum()),
                             finiteMax(hSim->GetMaximum()));
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

      hData->Draw("E");
      hSim->Draw("E SAME");

      TLegend *legend = new TLegend(0.52, 0.68, 0.90, 0.89);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.035);
      if (normalizeToUnity) {
        legend->AddEntry(hData, Form("Data (N=%.0f)", rawDataIntegral), "lep");
        legend->AddEntry(hSim, Form("Sim weighted (sumW=%.1f)", rawSimIntegral),
                         "lep");
      } else {
        legend->AddEntry(hData, Form("Data (N=%.0f)", rawDataIntegral), "lep");
        legend->AddEntry(hSim, Form("Sim weighted (sumW=%.1f)", rawSimIntegral),
                         "lep");
      }
      legend->Draw();

      TLatex label;
      label.SetNDC();
      label.SetTextSize(0.040);
      label.SetTextFont(42);
      label.DrawLatex(0.14, 0.92, combo.label);

      if (rawDataIntegral <= 0. && rawSimIntegral <= 0.) {
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

  std::cout << "[overlayBasicSrcFDCD] Wrote " << outputPdfName << '\n';
  std::cout << "[overlayBasicSrcFDCD] eppCut = " << eppCut << '\n';
  std::cout << "[overlayBasicSrcFDCD] baseCut = " << baseCut << '\n';
  std::cout << "[overlayBasicSrcFDCD] simWeightExpression = "
            << simWeightExpression << '\n';
  std::cout << "[overlayBasicSrcFDCD] pCMyTailCut = "
            << ((pCMyTailCut && pCMyTailCut[0] != '\0') ? pCMyTailCut
                                                         : "<none>")
            << '\n';

  dataFile->Close();
  simFile->Close();
}
