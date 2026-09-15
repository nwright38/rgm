#include "OverlaySrcCommon.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

void overlay_default_multi_cut_columns(
    const char *fileNamesCsv =
        "~/data/RGM_DATA/c12_src_skim.root,~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "pdf/c12_data_sim_overlay_default_multi_cut_columns.pdf",
    bool normalizeToUnity = true,
    const char *eppCut = "pCM > 0",
    const char *baseCut = "pMiss < 1. && recP < 1.",
    const char *weightsCsv = "(weight_epp),(weight_epp)*(weight_epp < 150.)",
    const char *labelsCsv = "C12 Data,C12 PWIA,C12 FSI",
    const char *pCMyTailCut = "",
    bool includeFdFd = false,
    Long64_t maxEvents = -1,
    Long64_t firstEvent = 0,
    const char *epCut = "goodLead",
    const char *epBaseCut = "pMiss < 1.",
    const char *epWeightsCsv = "(weight_ep),(weight_ep)",
    const char *extraCut = "xB > 1.3",
    const char *extraCutLabel = "xB > 1.3",
    bool omitLeadCdRecFd = true) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const std::vector<TString> fileNames = overlay_src::splitCsvToTStrings(fileNamesCsv);
  if (fileNames.size() < 2) {
    std::cerr << "Need at least two files in fileNamesCsv.\n";
    return;
  }

  std::vector<TString> weights = overlay_src::splitCsvToTStrings(weightsCsv);
  std::vector<TString> epWeights = overlay_src::splitCsvToTStrings(epWeightsCsv);
  std::vector<TString> labels = overlay_src::splitCsvToTStrings(labelsCsv);

  while (weights.size() < fileNames.size()) weights.push_back("(weight_epp)");
  while (epWeights.size() < fileNames.size()) epWeights.push_back("(weight_ep)");
  while (labels.size() < fileNames.size()) {
    labels.push_back(Form("Sample %zu", labels.size() + 1));
  }

  std::vector<TFile *> files;
  std::vector<TTree *> trees;
  files.reserve(fileNames.size());
  trees.reserve(fileNames.size());

  for (std::size_t i = 0; i < fileNames.size(); ++i) {
    TFile *f = TFile::Open(fileNames[i].Data(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "Could not open file " << fileNames[i] << "\n";
      for (TFile *ff : files) ff->Close();
      return;
    }
    TTree *t = dynamic_cast<TTree *>(f->Get(treeName));
    if (!t) {
      std::cerr << "Could not find tree \"" << treeName << "\" in " << fileNames[i]
                << "\n";
      f->Close();
      for (TFile *ff : files) ff->Close();
      return;
    }
    files.push_back(f);
    trees.push_back(t);
  }

  std::vector<overlay_src::DetectorCombination> eppCombos =
      overlay_src::detectorCombinations(includeFdFd);
  if (omitLeadCdRecFd) {
    eppCombos.erase(
        std::remove_if(eppCombos.begin(), eppCombos.end(),
                       [](const overlay_src::DetectorCombination &combo) {
                         return std::string(combo.key) == "lead_cd_rec_fd";
                       }),
        eppCombos.end());
  }
  if (eppCombos.empty()) {
    std::cerr << "No detector combinations selected.\n";
    for (TFile *f : files) f->Close();
    return;
  }

  const std::vector<overlay_src::DetectorCombination> epCombos = {
      {"lead_fd", "Lead FD", "(leadTheta*180./TMath::Pi()<37.)"},
      {"lead_cd", "Lead CD", "(leadTheta*180./TMath::Pi()>45.)"}};

  const std::vector<overlay_src::PlotVar> vars = overlay_src::defaultVariables();
  const std::vector<std::string> epCompatibleVars = {
      "xB",          "Q2",           "omega",         "eP",      "eTheta",
      "pLead",       "pMiss",        "kMiss",         "mMiss",   "EMiss",
      "E0miss",      "E1miss",       "leadTheta",     "pMissTheta",
      "theta_PmQ",   "theta_PmPlead", "theta_PleadQ", "alpha_1",
      "alpha_q",     "alpha_pLead",  "p1_plus",       "p1_lc_z",
      "p1_perp_mag", "lead_q_ratio", "qP",            "qTheta"};

  const auto isEpCompatible = [&epCompatibleVars](const char *name) {
    return std::find(epCompatibleVars.begin(), epCompatibleVars.end(),
                     std::string(name)) != epCompatibleVars.end();
  };

  const int nRows = 2;
  const int canvasHeight = std::max(880, 320 * nRows);
  TCanvas *c = new TCanvas("c_overlay_default_multi_cut_columns",
                           "Default multi overlay with extra-cut columns", 1450,
                           canvasHeight);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  const auto drawCutColumns = [&](const overlay_src::PlotVar &v,
                                  const std::vector<overlay_src::DetectorCombination> &pageCombos,
                                  const char *selCut,
                                  const char *selBaseCut,
                                  const char *pageTag,
                                  const std::vector<TString> &pageWeights) {
    c->Clear();
    c->Divide(2, 2);

    for (std::size_t ic = 0; ic < pageCombos.size() && ic < 2; ++ic) {
      const auto &combo = pageCombos[ic];

      std::vector<TH1D *> hNominal;
      std::vector<TH1D *> hExtra;
      std::vector<double> rawNominal;
      std::vector<double> rawExtra;
      hNominal.reserve(fileNames.size());
      hExtra.reserve(fileNames.size());
      rawNominal.reserve(fileNames.size());
      rawExtra.reserve(fileNames.size());

      double yMax = 0.;

      for (std::size_t is = 0; is < trees.size(); ++is) {
        TH1D *hLeft = new TH1D(Form("h_nom_%s_%s_s%zu", v.name, combo.key, is), "",
                               v.bins, v.xmin, v.xmax);
        TH1D *hRight = new TH1D(Form("h_extra_%s_%s_s%zu", v.name, combo.key, is), "",
                                v.bins, v.xmin, v.xmax);
        hLeft->Sumw2();
        hRight->Sumw2();
        overlay_src::styleTopologyData(hLeft, static_cast<int>(is));
        overlay_src::styleTopologyData(hRight, static_cast<int>(is));

        const TString nominalSel =
          overlay_src::buildSelection(selCut, selBaseCut, combo.cut, pCMyTailCut);
        const TString extraSel = overlay_src::buildSelection(
          selCut, selBaseCut, combo.cut, pCMyTailCut, nullptr, extraCut);

        const TString nominalWeighted =
          overlay_src::applyWeightToSelection(pageWeights[is].Data(), nominalSel);
        const TString extraWeighted =
          overlay_src::applyWeightToSelection(pageWeights[is].Data(), extraSel);

        if (maxEvents > 0) {
          trees[is]->Project(hLeft->GetName(), v.expr, nominalWeighted.Data(), "",
                             maxEvents, firstEvent);
          trees[is]->Project(hRight->GetName(), v.expr, extraWeighted.Data(), "",
                             maxEvents, firstEvent);
        } else {
          trees[is]->Project(hLeft->GetName(), v.expr, nominalWeighted.Data());
          trees[is]->Project(hRight->GetName(), v.expr, extraWeighted.Data());
        }

        rawNominal.push_back(hLeft->Integral(0, hLeft->GetNbinsX() + 1));
        rawExtra.push_back(hRight->Integral(0, hRight->GetNbinsX() + 1));

        if (normalizeToUnity) {
          const double nominalInRange = hLeft->Integral();
          const double extraInRange = hRight->Integral();
          if (nominalInRange > 0.) hLeft->Scale(1. / nominalInRange);
          if (extraInRange > 0.) hRight->Scale(1. / extraInRange);
        }

        yMax = std::max(yMax, overlay_src::finiteMax(hLeft->GetMaximum()));
        yMax = std::max(yMax, overlay_src::finiteMax(hRight->GetMaximum()));

        hNominal.push_back(hLeft);
        hExtra.push_back(hRight);
      }

      if (hNominal.empty()) continue;
      if (yMax <= 0.) yMax = 1.;

      const int leftPad = 1 + 2 * static_cast<int>(ic);
      const int rightPad = leftPad + 1;

      c->cd(leftPad);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.10);

      hNominal.front()->SetMinimum(0.);
      hNominal.front()->SetMaximum(1.25 * yMax);
      hNominal.front()->SetTitle(Form("%s (%s) | Nominal", v.name, pageTag));
      hNominal.front()->GetXaxis()->SetTitle(v.xTitle);
      hNominal.front()->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts"
                                                              : "Counts");
      hNominal.front()->GetXaxis()->SetTitleOffset(1.05);
      hNominal.front()->GetXaxis()->SetTitleSize(0.050);
      hNominal.front()->GetYaxis()->SetTitleSize(0.050);
      hNominal.front()->GetXaxis()->SetLabelSize(0.043);
      hNominal.front()->GetYaxis()->SetLabelSize(0.043);
      hNominal.front()->GetYaxis()->SetTitleOffset(1.20);
      hNominal.front()->Draw("E");
      for (std::size_t is = 1; is < hNominal.size(); ++is) hNominal[is]->Draw("E SAME");

        TLegend *legLeft = new TLegend(0.50, 0.64, 0.90, 0.89);
      legLeft->SetBorderSize(0);
      legLeft->SetFillStyle(0);
        legLeft->SetTextSize(0.030);
      for (std::size_t is = 0; is < hNominal.size(); ++is) {
        legLeft->AddEntry(hNominal[is],
                          Form("%s (sumW=%.1f)", labels[is].Data(), rawNominal[is]),
                          "lep");
      }
      legLeft->Draw();

        TLatex labLeft;
        labLeft.SetNDC();
        labLeft.SetTextSize(0.030);
        labLeft.SetTextFont(42);
        labLeft.DrawLatex(0.14, 0.92, combo.label);

      c->cd(rightPad);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.10);

      hExtra.front()->SetMinimum(0.);
      hExtra.front()->SetMaximum(1.25 * yMax);
      hExtra.front()->SetTitle(
          Form("%s (%s) | Nominal + %s", v.name, pageTag,
               (extraCutLabel && extraCutLabel[0] != '\0') ? extraCutLabel : extraCut));
      hExtra.front()->GetXaxis()->SetTitle(v.xTitle);
      hExtra.front()->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts"
                                                            : "Counts");
      hExtra.front()->GetXaxis()->SetTitleOffset(1.05);
      hExtra.front()->GetXaxis()->SetTitleSize(0.050);
      hExtra.front()->GetYaxis()->SetTitleSize(0.050);
      hExtra.front()->GetXaxis()->SetLabelSize(0.043);
      hExtra.front()->GetYaxis()->SetLabelSize(0.043);
      hExtra.front()->GetYaxis()->SetTitleOffset(1.20);
      hExtra.front()->Draw("E");
      for (std::size_t is = 1; is < hExtra.size(); ++is) hExtra[is]->Draw("E SAME");

      TLegend *legRight = new TLegend(0.50, 0.64, 0.90, 0.89);
      legRight->SetBorderSize(0);
      legRight->SetFillStyle(0);
      legRight->SetTextSize(0.030);
      for (std::size_t is = 0; is < hExtra.size(); ++is) {
        legRight->AddEntry(hExtra[is],
                           Form("%s (sumW=%.1f)", labels[is].Data(), rawExtra[is]),
                           "lep");
      }
      legRight->Draw();

      TLatex labRight;
      labRight.SetNDC();
      labRight.SetTextSize(0.030);
      labRight.SetTextFont(42);
      labRight.DrawLatex(0.14, 0.92, combo.label);
    }

    for (std::size_t ip = pageCombos.size(); ip < 2; ++ip) {
      const int leftPad = 1 + 2 * static_cast<int>(ip);
      const int rightPad = leftPad + 1;

      c->cd(leftPad);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      TLatex omitted;
      omitted.SetNDC();
      omitted.SetTextSize(0.040);
      omitted.SetTextColor(kGray + 1);
      omitted.DrawLatex(0.22, 0.50, "Unused panel");

      c->cd(rightPad);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      TLatex omittedRight;
      omittedRight.SetNDC();
      omittedRight.SetTextSize(0.040);
      omittedRight.SetTextColor(kGray + 1);
      omittedRight.DrawLatex(0.22, 0.50, "Unused panel");
    }

    c->Print(pdfName, "pdf");
  };

  for (const auto &v : vars) {
    if (isEpCompatible(v.name)) {
      drawCutColumns(v, epCombos, epCut, epBaseCut, "e'p", epWeights);
      drawCutColumns(v, eppCombos, eppCut, baseCut, "e'pp", weights);
    } else {
      drawCutColumns(v, eppCombos, eppCut, baseCut, "e'pp", weights);
    }
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";

  for (TFile *f : files) f->Close();
}