#include "OverlaySrcCommon.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

void styleDataHist(TH1D *hist) {
  hist->SetStats(false);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.75);
  hist->SetLineWidth(2);
}

void styleSimHist(TH1D *hist) {
  hist->SetStats(false);
  hist->SetLineColor(kRed + 1);
  hist->SetMarkerColor(kRed + 1);
  hist->SetMarkerStyle(24);
  hist->SetMarkerSize(0.75);
  hist->SetLineWidth(2);
}

std::vector<TPad *> makeFivePanelPads(TCanvas *canvas) {
  canvas->Clear();

  const double xLeft = 0.06;
  const double xMidLeft = 0.365;
  const double xMidRight = 0.635;
  const double xRight = 0.94;

  const double yTop = 0.94;
  const double yMid = 0.51;
  const double yBottom = 0.08;

  std::vector<TPad *> pads;
  pads.reserve(5);

  pads.push_back(new TPad("pad1", "pad1", xLeft, yMid, xMidLeft, yTop));
  pads.push_back(new TPad("pad2", "pad2", xMidLeft, yMid, xMidRight, yTop));
  pads.push_back(new TPad("pad3", "pad3", xMidRight, yMid, xRight, yTop));
  pads.push_back(new TPad("pad4", "pad4", 0.19, yBottom, 0.50, yMid));
  pads.push_back(new TPad("pad5", "pad5", 0.50, yBottom, 0.81, yMid));

  for (TPad *pad : pads) {
    pad->SetLeftMargin(0.13);
    pad->SetBottomMargin(0.12);
    pad->SetRightMargin(0.04);
    pad->SetTopMargin(0.08);
    pad->Draw();
  }

  return pads;
}

}  // namespace

void overlay_data_sim_q2_by_detector(
    const char *dataFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *simFileName = "~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "pdf/c12_data_sim_q2_by_detector.pdf",
    bool normalizeToUnity = true,
    const char *eppCut = "pCM > 0",
    const char *baseCut = "pCM > 0 && pMiss < 1. && recP < 1.",
    const char *dataWeightExpression = "",
    const char *simWeightExpression = "(weight_epp)*(weight_epp<200)",
    const char *dataLabel = "Data",
    const char *simLabel = "Sim",
    const char *pCMyTailCut = "",
    const char *q2EdgesCsv = "1.5,1.80,2.10,2.40,3.00,5.0",
    bool autoRebin = true,
    int maxRebinFactor = 2,
    bool includeFdFd = false) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile *dataFile = TFile::Open(dataFileName, "READ");
  TFile *simFile = TFile::Open(simFileName, "READ");
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Could not open data file " << dataFileName << "\n";
    if (simFile) simFile->Close();
    return;
  }
  if (!simFile || simFile->IsZombie()) {
    std::cerr << "Could not open sim file " << simFileName << "\n";
    dataFile->Close();
    return;
  }

  TTree *dataTree = dynamic_cast<TTree *>(dataFile->Get(treeName));
  TTree *simTree = dynamic_cast<TTree *>(simFile->Get(treeName));
  if (!dataTree || !simTree) {
    std::cerr << "Could not find tree \"" << treeName << "\" in both files\n";
    dataFile->Close();
    simFile->Close();
    return;
  }

  std::vector<TString> q2EdgeStrings = overlay_src::splitCsvToTStrings(q2EdgesCsv);
  std::vector<double> q2Edges;
  q2Edges.reserve(q2EdgeStrings.size());
  for (const auto &edge : q2EdgeStrings) q2Edges.push_back(edge.Atof());
  if (q2Edges.size() < 2) {
    std::cerr << "Need at least two Q2 edges.\n";
    dataFile->Close();
    simFile->Close();
    return;
  }

  const std::vector<overlay_src::DetectorCombination> combos =
      overlay_src::detectorCombinations(includeFdFd);
  const std::vector<overlay_src::PlotVar> vars = overlay_src::defaultVariables();

  TCanvas *c = new TCanvas("c_overlay_data_sim_q2_by_detector",
                           "Data/Sim Q2 panels by detector", 1300, 950);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  const std::size_t nQ2Bins = q2Edges.size() - 1;
  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  for (const auto &v : vars) {
    for (const auto &combo : combos) {
      double representativeEntries = 0.;
      int representativeCount = 0;
      for (std::size_t iq2 = 0; iq2 < nQ2Bins; ++iq2) {
        const TString q2Cut = Form("(Q2>=%g) && (Q2<%g)", q2Edges[iq2], q2Edges[iq2 + 1]);
        const TString sel = overlay_src::buildSelection(
            eppCut, baseCut, combo.cut, pCMyTailCut, q2Cut.Data());
        const Long64_t nDataSel = dataTree->GetEntries(sel.Data());
        const Long64_t nSimSel = simTree->GetEntries(sel.Data());
        const Long64_t nSel = std::max(nDataSel, nSimSel);
        if (nSel > 0) {
          representativeEntries += static_cast<double>(nSel);
          ++representativeCount;
        }
      }
      if (representativeCount > 0) {
        representativeEntries /= static_cast<double>(representativeCount);
      }
      const int rebinFactor = autoRebin
                                  ? overlay_src::chooseAutoRebinFactor(
                                        v.bins, representativeEntries, 8.0, 5,
                                        maxRebinFactor)
                                  : 1;

      for (std::size_t firstBin = 0; firstBin < nQ2Bins; firstBin += 5) {
        std::vector<TPad *> pads = makeFivePanelPads(c);

        for (std::size_t panel = 0; panel < 5; ++panel) {
          pads[panel]->cd();

          const std::size_t iq2 = firstBin + panel;
          if (iq2 >= nQ2Bins) {
            TLatex omitted;
            omitted.SetNDC();
            omitted.SetTextSize(0.040);
            omitted.SetTextColor(kGray + 1);
            omitted.DrawLatex(0.18, 0.52, "Unused panel");
            continue;
          }

          TH1D *hData = new TH1D(Form("h_data_%s_%s_bin%zu", v.name, combo.key, iq2), "",
                                 v.bins, v.xmin, v.xmax);
          TH1D *hSim = new TH1D(Form("h_sim_%s_%s_bin%zu", v.name, combo.key, iq2), "",
                                v.bins, v.xmin, v.xmax);
          hData->Sumw2();
          hSim->Sumw2();
          styleDataHist(hData);
          styleSimHist(hSim);

          const TString q2Cut = Form("(Q2>=%g) && (Q2<%g)", q2Edges[iq2], q2Edges[iq2 + 1]);
          const TString boolSel = overlay_src::buildSelection(
              eppCut, baseCut, combo.cut, pCMyTailCut, q2Cut.Data());
          const TString dataSel = overlay_src::applyWeightToSelection(dataWeightExpression, boolSel);
          const TString simSel = overlay_src::applyWeightToSelection(simWeightExpression, boolSel);

          dataTree->Project(hData->GetName(), v.expr, dataSel.Data());
          simTree->Project(hSim->GetName(), v.expr, simSel.Data());

          const double rawDataIntegral = hData->Integral(0, hData->GetNbinsX() + 1);
          const double rawSimIntegral = hSim->Integral(0, hSim->GetNbinsX() + 1);

          if (rebinFactor > 1) {
            hData->Rebin(rebinFactor);
            hSim->Rebin(rebinFactor);
          }

          if (normalizeToUnity) {
            const double dataInRange = hData->Integral();
            const double simInRange = hSim->Integral();
            if (dataInRange > 0.) hData->Scale(1. / dataInRange);
            if (simInRange > 0.) hSim->Scale(1. / simInRange);
          }

          double yMax = std::max(overlay_src::finiteMax(hData->GetMaximum()),
                                 overlay_src::finiteMax(hSim->GetMaximum()));
          if (yMax <= 0.) yMax = 1.;

          hData->SetMinimum(0.);
          hData->SetMaximum(1.25 * yMax);
          hData->GetXaxis()->SetTitle(v.xTitle);
          hData->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts" : "Counts");
          hData->GetXaxis()->SetTitleSize(0.050);
          hData->GetYaxis()->SetTitleSize(0.050);
          hData->GetXaxis()->SetLabelSize(0.043);
          hData->GetYaxis()->SetLabelSize(0.043);
          hData->GetYaxis()->SetTitleOffset(1.20);

          hData->Draw("E");
          hSim->Draw("E SAME");

          TLegend *leg = new TLegend(0.53, 0.64, 0.90, 0.89);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.028);
          leg->AddEntry(hData, Form("%s (sumW=%.1f)", dataLabel, rawDataIntegral), "lep");
          leg->AddEntry(hSim, Form("%s (sumW=%.1f)", simLabel, rawSimIntegral), "lep");
          leg->Draw();

          TLatex lab;
          lab.SetNDC();
          lab.SetTextSize(0.038);
          lab.SetTextFont(42);
          lab.DrawLatex(0.14, 0.92,
                        Form("%.2f < Q^{2} < %.2f", q2Edges[iq2], q2Edges[iq2 + 1]));
          if (rebinFactor > 1) {
            lab.SetTextSize(0.030);
            lab.DrawLatex(0.14, 0.86, Form("auto rebin x%d", rebinFactor));
          }
          if (rawDataIntegral <= 0. && rawSimIntegral <= 0.) {
            TLatex empty;
            empty.SetNDC();
            empty.SetTextSize(0.038);
            empty.SetTextColor(kGray + 2);
            empty.DrawLatex(0.20, 0.52, "No entries after cuts");
          }
        }

        c->cd(0);
        TLatex title;
        title.SetNDC();
        title.SetTextSize(0.028);
        title.SetTextFont(62);
        title.DrawLatex(0.10, 0.985,
                        Form("%s: %s", v.name, combo.label));

        c->Print(pdfName, "pdf");
      }
    }
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << pdfName << "\n";

  dataFile->Close();
  simFile->Close();
}