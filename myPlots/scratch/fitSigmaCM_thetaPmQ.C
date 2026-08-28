#include "TCanvas.h"
#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

void fitSigmaCM_thetaPmQ(
    const char *pdfName = "fitSigmaCM_thetaPmQ.pdf",
    const char *dataFileName = "~/data/RGM_DATA/he4_src_skim_100MeV.root",
    const char *simFileName = "~/data/RGM_DATA/he4_sim_skim_100MeV.root",
    const char *treeName = "srcTree",
    const char *baseCutExpr = "pCM > 0 && pMiss < 1. && recP < 1. &&recTheta*180./TMath::Pi()>0. && leadTheta*180./TMath::Pi()< 37",
  const char *simWeightExpr = "weight_epp") {

  const int nSamples = 2;
  const char *sampleNames[nSamples] = {"Data", "Simulation"};
  const char *inputFileNames[nSamples] = {dataFileName, simFileName};

  TFile *inputFiles[nSamples] = {nullptr, nullptr};
  TTree *inputTrees[nSamples] = {nullptr, nullptr};
  for (int sample = 0; sample < nSamples; ++sample) {
    inputFiles[sample] = TFile::Open(inputFileNames[sample], "READ");
    if (!inputFiles[sample] || inputFiles[sample]->IsZombie()) {
      std::cerr << "[fitSigmaCM_thetaPmQ] Could not open "
                << inputFileNames[sample] << "\n";
      if (inputFiles[sample]) inputFiles[sample]->Close();
      for (int j = 0; j < sample; ++j) {
        if (inputFiles[j]) inputFiles[j]->Close();
      }
      return;
    }
    inputTrees[sample] = dynamic_cast<TTree *>(inputFiles[sample]->Get(treeName));
    if (!inputTrees[sample]) {
      std::cerr << "[fitSigmaCM_thetaPmQ] Missing tree " << treeName
                << " in " << inputFileNames[sample] << "\n";
      for (int j = 0; j <= sample; ++j) {
        if (inputFiles[j]) inputFiles[j]->Close();
      }
      return;
    }
  }

  const int nComponents = 2;
  const char *components[nComponents] = {"pCMx", "pCMy"};
  const int colors[nComponents] = {kBlue + 1, kRed + 1};

  const int nBins = 4;
  const double binEdges[nBins + 1] = {110, 130, 140,150, 180};

  const char *binExprDeg = "theta_PmQ*180./TMath::Pi()";
  const char *binVariableTitle = "#theta_{p_{miss},q} [deg]";
  const char *sigmaVariableTitle = "#sigma_{CM} [GeV/c]";

  TH1D *componentHists[nSamples][nComponents][nBins];
  TF1 *componentFits[nSamples][nComponents][nBins];
  TGraphErrors *sigmaGraphs[nSamples][nComponents];

  for (int sample = 0; sample < nSamples; ++sample) {
    for (int component = 0; component < nComponents; ++component) {
      sigmaGraphs[sample][component] = new TGraphErrors(nBins);
      sigmaGraphs[sample][component]->SetName(
          Form("sigma_%s_%s_vs_thetaPmQ", components[component], sampleNames[sample]));
      sigmaGraphs[sample][component]->SetTitle(
          Form("%s: #sigma_{%s}", sampleNames[sample], components[component]));
      sigmaGraphs[sample][component]->SetMarkerStyle(20 + component + 2 * sample);
      sigmaGraphs[sample][component]->SetMarkerColor(colors[component]);
      sigmaGraphs[sample][component]->SetLineColor(colors[component]);
      sigmaGraphs[sample][component]->SetLineWidth(2);

      for (int i = 0; i < nBins; ++i) {
        const double low = binEdges[i];
        const double high = binEdges[i + 1];
        componentHists[sample][component][i] = new TH1D(
            Form("%s_%s_thetaPmQ_bin_%d", sampleNames[sample], components[component], i),
            Form("%s %s: %.0f < %s < %.0f;%s [GeV/c];Counts",
                 sampleNames[sample], components[component], low, binVariableTitle, high,
                 components[component]),
            60, -0.8, 0.8);
        componentHists[sample][component][i]->Sumw2();

        componentFits[sample][component][i] = new TF1(
            Form("%s_%s_fit_thetaPmQ_bin_%d", sampleNames[sample], components[component], i),
            "gaus", -0.15, 0.15);
        componentFits[sample][component][i]->SetLineColor(colors[component]);
        componentFits[sample][component][i]->SetLineWidth(2);

        const TString binCutExpr =
            Form("(%s > %.12g) && (%s <= %.12g)", binExprDeg, low, binExprDeg, high);
        const TString boolCutExpr =
          Form("(%s) && (%s)", baseCutExpr, binCutExpr.Data());
        const TString selectionExpr =
          (sample == 0)
            ? boolCutExpr
            : Form("(%s) * (%s)", simWeightExpr, boolCutExpr.Data());

        inputTrees[sample]->Draw(
            Form("%s>>%s", components[component],
                 componentHists[sample][component][i]->GetName()),
          selectionExpr.Data(), "goff");

        const double entries = componentHists[sample][component][i]->GetEntries();
        const double sigmaGuess = componentHists[sample][component][i]->GetRMS();
        const double meanGuess = componentHists[sample][component][i]->GetMean();
        if (entries > 5 && std::isfinite(sigmaGuess) && sigmaGuess > 0.) {
          componentFits[sample][component][i]->SetParameters(
              componentHists[sample][component][i]->GetMaximum(), meanGuess, sigmaGuess);
          componentHists[sample][component][i]->Fit(
              componentFits[sample][component][i], "RQ0");
        }

        double sigma = TMath::QuietNaN();
        double sigmaErr = TMath::QuietNaN();
        if (entries > 5 && componentFits[sample][component][i]->GetNDF() > 0) {
          sigma = std::abs(componentFits[sample][component][i]->GetParameter(2));
          sigmaErr = componentFits[sample][component][i]->GetParError(2);
        }

        const double center = 0.5 * (low + high);
        const double halfWidth = 0.5 * (high - low);
        sigmaGraphs[sample][component]->SetPoint(i, center, sigma);
        sigmaGraphs[sample][component]->SetPointError(i, halfWidth, sigmaErr);
      }

       // sigmaGraphs[sample][component]->GetYaxis()->SetRangeUser(.1,.35);

    }
  }

  TCanvas *fitCanvas = new TCanvas("fitCanvas_thetaPmQ", "sigmaCM fits by thetaPmQ", 1200, 850);
  bool firstPdfPage = true;
  for (int sample = 0; sample < nSamples; ++sample) {
    for (int component = 0; component < nComponents; ++component) {
      fitCanvas->Clear();
      fitCanvas->Divide(3, 3);
      for (int i = 0; i < nBins; ++i) {
        fitCanvas->cd(i + 1);
        componentHists[sample][component][i]->SetLineColor(kBlack);
        componentHists[sample][component][i]->SetMarkerStyle(20);
        componentHists[sample][component][i]->SetMarkerSize(0.7);
        componentHists[sample][component][i]->Draw("E");
        if (componentFits[sample][component][i]->GetNDF() > 0) {
          componentFits[sample][component][i]->Draw("same");
        }

        TLatex txt;
        txt.SetNDC();
        txt.SetTextSize(0.05);
        if (componentFits[sample][component][i]->GetNDF() > 0) {
          txt.DrawLatex(0.14, 0.86,
                        Form("#sigma = %.4f #pm %.4f",
                             std::abs(componentFits[sample][component][i]->GetParameter(2)),
                             componentFits[sample][component][i]->GetParError(2)));
        } else {
          txt.DrawLatex(0.14, 0.86, "Fit unavailable");
        }
      }
      fitCanvas->Print(firstPdfPage ? Form("%s(", pdfName) : pdfName, "pdf");
      firstPdfPage = false;
    }
  }

  TCanvas *graphCanvas = new TCanvas("graphCanvas_thetaPmQ", "sigmaCM vs thetaPmQ", 900, 700);
  for (int sample = 0; sample < nSamples; ++sample) {
    graphCanvas->Clear();
    sigmaGraphs[sample][0]->SetTitle(
        Form("%s;%s;%s", sampleNames[sample], binVariableTitle, sigmaVariableTitle));
   // sigmaGraphs[sample][0]->SetMinimum(0.0);
    sigmaGraphs[sample][0]->Draw("AP");
    sigmaGraphs[sample][1]->Draw("P same");

    TLegend *leg = new TLegend(0.65, 0.74, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(sigmaGraphs[sample][0], "pCMx", "lp");
    leg->AddEntry(sigmaGraphs[sample][1], "pCMy", "lp");
    leg->Draw();

    graphCanvas->Print(pdfName, "pdf");
  }

  graphCanvas->Clear();
  sigmaGraphs[0][0]->SetTitle(
      Form("#sigma_{CM,x}: Data vs Simulation;%s;%s", binVariableTitle, sigmaVariableTitle));
  sigmaGraphs[0][0]->SetMarkerStyle(20);
  sigmaGraphs[0][0]->SetMarkerColor(kBlack);
  sigmaGraphs[0][0]->SetLineColor(kBlack);
  sigmaGraphs[1][0]->SetMarkerStyle(24);
  sigmaGraphs[1][0]->SetMarkerColor(kBlue + 1);
  sigmaGraphs[1][0]->SetLineColor(kBlue + 1);
  //sigmaGraphs[0][0]->SetMinimum(0.0);
  sigmaGraphs[0][0]->Draw("AP");
  sigmaGraphs[1][0]->Draw("P same");

  TLegend *cmpLeg = new TLegend(0.62, 0.74, 0.88, 0.88);
  cmpLeg->SetBorderSize(0);
  cmpLeg->SetFillStyle(0);
  cmpLeg->AddEntry(sigmaGraphs[0][0], "Data pCMx", "lp");
  cmpLeg->AddEntry(sigmaGraphs[1][0], "Simulation pCMx", "lp");
  cmpLeg->Draw();

  graphCanvas->Print(pdfName, "pdf");

  graphCanvas->Clear();
  sigmaGraphs[0][1]->SetTitle(
      Form("#sigma_{CM,y}: Data vs Simulation;%s;%s", binVariableTitle, sigmaVariableTitle));
  sigmaGraphs[0][1]->SetMarkerStyle(20);
  sigmaGraphs[0][1]->SetMarkerColor(kBlack);
  sigmaGraphs[0][1]->SetLineColor(kBlack);
  sigmaGraphs[1][1]->SetMarkerStyle(24);
  sigmaGraphs[1][1]->SetMarkerColor(kBlue + 1);
  sigmaGraphs[1][1]->SetLineColor(kBlue + 1);
  sigmaGraphs[0][1]->Draw("AP");
  sigmaGraphs[1][1]->Draw("P same");

  TLegend *cmpLegY = new TLegend(0.62, 0.74, 0.88, 0.88);
  cmpLegY->SetBorderSize(0);
  cmpLegY->SetFillStyle(0);
  cmpLegY->AddEntry(sigmaGraphs[0][1], "Data pCMy", "lp");
  cmpLegY->AddEntry(sigmaGraphs[1][1], "Simulation pCMy", "lp");
  cmpLegY->Draw();

  graphCanvas->Print(pdfName, "pdf");

  graphCanvas->Clear();
  sigmaGraphs[0][0]->SetTitle(
      Form("#sigma_{CM}: Data/Simulation combined;%s;%s", binVariableTitle,
           sigmaVariableTitle));

  sigmaGraphs[0][0]->SetMarkerStyle(20);
  sigmaGraphs[0][0]->SetMarkerColor(kBlack);
  sigmaGraphs[0][0]->SetLineColor(kBlack);

  sigmaGraphs[1][0]->SetMarkerStyle(24);
  sigmaGraphs[1][0]->SetMarkerColor(kBlue + 1);
  sigmaGraphs[1][0]->SetLineColor(kBlue + 1);

  sigmaGraphs[0][1]->SetMarkerStyle(21);
  sigmaGraphs[0][1]->SetMarkerColor(kRed + 1);
  sigmaGraphs[0][1]->SetLineColor(kRed + 1);

  sigmaGraphs[1][1]->SetMarkerStyle(25);
  sigmaGraphs[1][1]->SetMarkerColor(kMagenta + 2);
  sigmaGraphs[1][1]->SetLineColor(kMagenta + 2);

  sigmaGraphs[0][0]->Draw("AP");
  sigmaGraphs[1][0]->Draw("P same");
  sigmaGraphs[0][1]->Draw("P same");
  sigmaGraphs[1][1]->Draw("P same");

  TLegend *comboLeg = new TLegend(0.58, 0.68, 0.88, 0.88);
  comboLeg->SetBorderSize(0);
  comboLeg->SetFillStyle(0);
  comboLeg->AddEntry(sigmaGraphs[0][0], "Data pCMx", "lp");
  comboLeg->AddEntry(sigmaGraphs[1][0], "Simulation pCMx", "lp");
  comboLeg->AddEntry(sigmaGraphs[0][1], "Data pCMy", "lp");
  comboLeg->AddEntry(sigmaGraphs[1][1], "Simulation pCMy", "lp");
  comboLeg->Draw();

  graphCanvas->Print(pdfName, "pdf");

  // Extra page: ratio sigma_CMx/sigma_CMy for data and simulation.
  TGraphErrors *ratioGraphs[nSamples] = {new TGraphErrors(nBins),
                                         new TGraphErrors(nBins)};
  for (int sample = 0; sample < nSamples; ++sample) {
    ratioGraphs[sample]->SetName(Form("ratio_sigmaCMx_over_sigmaCMy_%s", sampleNames[sample]));
    ratioGraphs[sample]->SetMarkerStyle(sample == 0 ? 20 : 24);
    ratioGraphs[sample]->SetMarkerColor(sample == 0 ? kBlack : (kBlue + 1));
    ratioGraphs[sample]->SetLineColor(sample == 0 ? kBlack : (kBlue + 1));
    ratioGraphs[sample]->SetLineWidth(2);

    for (int i = 0; i < nBins; ++i) {
      const double x = sigmaGraphs[sample][0]->GetPointX(i);
      const double xErr = sigmaGraphs[sample][0]->GetErrorX(i);
      const double sigX = sigmaGraphs[sample][0]->GetPointY(i);
      const double sigY = sigmaGraphs[sample][1]->GetPointY(i);
      const double sigXErr = sigmaGraphs[sample][0]->GetErrorY(i);
      const double sigYErr = sigmaGraphs[sample][1]->GetErrorY(i);

      double ratio = TMath::QuietNaN();
      double ratioErr = TMath::QuietNaN();
      if (std::isfinite(sigX) && std::isfinite(sigY) && std::isfinite(sigXErr) &&
          std::isfinite(sigYErr) && sigX > 0. && sigY > 0.) {
        ratio = sigX / sigY;
        const double relX = sigXErr / sigX;
        const double relY = sigYErr / sigY;
        ratioErr = ratio * std::sqrt(relX * relX + relY * relY);
      }

      ratioGraphs[sample]->SetPoint(i, x, ratio);
      ratioGraphs[sample]->SetPointError(i, xErr, ratioErr);
    }
  }

  graphCanvas->Clear();
  ratioGraphs[0]->SetTitle(Form("#sigma_{CM,x}/#sigma_{CM,y}: Data vs Simulation;%s;#sigma_{CM,x}/#sigma_{CM,y}",
                                binVariableTitle));


  ratioGraphs[0]->GetYaxis()->SetRangeUser(0.5, 1.2);
  ratioGraphs[0]->Draw("AP");
  ratioGraphs[1]->Draw("P same");
  TLegend *ratioLeg = new TLegend(0.62, 0.74, 0.88, 0.88);
  ratioLeg->SetBorderSize(0);
  ratioLeg->SetFillStyle(0);
  ratioLeg->AddEntry(ratioGraphs[0], "Data", "lp");
  ratioLeg->AddEntry(ratioGraphs[1], "Simulation", "lp");
  ratioLeg->Draw();

  graphCanvas->Print(Form("%s)", pdfName), "pdf");

  std::cout << "[fitSigmaCM_thetaPmQ] Wrote " << pdfName << "\n";

  for (int sample = 0; sample < nSamples; ++sample) {
    if (inputFiles[sample]) inputFiles[sample]->Close();
  }
}