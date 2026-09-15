#include "OverlaySrcCommon.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

void overlay_default_multi(
    const char *fileNamesCsv =
    "~/data/RGM_DATA/c12_src_skim_onlyPlead.root,~/data/RGM_DATA/c12_sim_skim_onlyPlead.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "pdf/c12_data_sim_overlay_default_multi_ep_test.pdf",
    bool normalizeToUnity = true,
    const char *eppCut = "pCM > 0",
    const char *baseCut = "Q2 < 5.",
    const char *weightsCsv = "(weight_epp),(weight_epp)*(weight_epp < 150.)",
    const char *labelsCsv = "C12 Data,C12 PWIA,C12 FSI",
    const char *pCMyTailCut = "",
    bool includeFdFd = false,
    Long64_t maxEvents = -1,
    Long64_t firstEvent = 0,
    const char *epCut = "goodLead",
    const char *epBaseCut = "Q2 < 5.",
    const char *epWeightsCsv = "(weight_ep),(weight_ep)*(weight_ep < 200.)") {

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

  const std::vector<overlay_src::DetectorCombination> combos =
      overlay_src::detectorCombinations(includeFdFd);
  const std::vector<overlay_src::DetectorCombination> epCombos = {
      {"lead_fd", "Lead FD", "(leadTheta*180./TMath::Pi()<37.)"},
      {"lead_cd", "Lead CD", "(leadTheta*180./TMath::Pi()>45.)"}};
  const std::vector<overlay_src::PlotVar> vars = overlay_src::defaultVariables();

  const std::vector<std::string> epCompatibleVars = {
      "xB",          "Q2",          "omega",       "eP",      "eTheta",
      "pLead",       "pMiss",       "kMiss",       "mMiss",   "EMiss",
      "E0miss",      "E1miss",      "leadTheta",   "pMissTheta",
      "theta_PmQ",   "theta_PmPlead","theta_PleadQ","alpha_1",
      "alpha_q",     "alpha_pLead", "p1_plus",     "p1_lc_z",
      "p1_perp_mag", "lead_q_ratio",
      "qP",          "qTheta"};

  const auto isEpCompatible = [&epCompatibleVars](const char *name) {
    return std::find(epCompatibleVars.begin(), epCompatibleVars.end(),
                     std::string(name)) != epCompatibleVars.end();
  };

  TCanvas *c = new TCanvas("c_overlay_default_multi", "Default multi overlay", 1300, 950);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  const auto drawPage = [&](const overlay_src::PlotVar &v,
                            const std::vector<overlay_src::DetectorCombination> &pageCombos,
                            const char *selCut,
                            const char *selBaseCut,
                            const char *pageTag,
                            const std::vector<TString> &pageWeights) {
    c->Clear();
    c->Divide(2, 2);

    for (std::size_t ic = 0; ic < pageCombos.size() && ic < 4; ++ic) {
      const auto &combo = pageCombos[ic];
      c->cd(static_cast<int>(ic) + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.10);

      std::vector<TH1D *> hs;
      std::vector<double> rawIntegrals;
      hs.reserve(fileNames.size());
      rawIntegrals.reserve(fileNames.size());
      double yMax = 0.;

      for (std::size_t is = 0; is < trees.size(); ++is) {
        TH1D *h = new TH1D(Form("h_%s_%s_%s_s%zu", v.name, pageTag, combo.key, is), "",
                           v.bins, v.xmin, v.xmax);
        h->Sumw2();
        overlay_src::styleTopologyData(h, static_cast<int>(is));

        const TString boolSel =
            overlay_src::buildSelection(selCut, selBaseCut, combo.cut, pCMyTailCut);
        const TString weightedSel =
          overlay_src::applyWeightToSelection(pageWeights[is].Data(), boolSel);

        if (maxEvents > 0) {
          trees[is]->Project(h->GetName(), v.expr, weightedSel.Data(), "", maxEvents,
                             firstEvent);
        } else {
          trees[is]->Project(h->GetName(), v.expr, weightedSel.Data());
        }
        const double rawInt = h->Integral(0, h->GetNbinsX() + 1);
        rawIntegrals.push_back(rawInt);

        if (normalizeToUnity) {
          const double intInRange = h->Integral();
          if (intInRange > 0.) h->Scale(1. / intInRange);
        }

        yMax = std::max(yMax, overlay_src::finiteMax(h->GetMaximum()));
        hs.push_back(h);
      }

      if (hs.empty()) continue;

      if (yMax <= 0.) yMax = 1.;
      hs.front()->SetMinimum(0.);
      hs.front()->SetMaximum(1.25 * yMax);
      hs.front()->GetXaxis()->SetTitle(v.xTitle);
      hs.front()->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts" : "Counts");
      hs.front()->GetXaxis()->SetTitleSize(0.050);
      hs.front()->GetYaxis()->SetTitleSize(0.050);
      hs.front()->GetXaxis()->SetLabelSize(0.043);
      hs.front()->GetYaxis()->SetLabelSize(0.043);
      hs.front()->GetYaxis()->SetTitleOffset(1.20);

      hs.front()->Draw("E");
      for (std::size_t is = 1; is < hs.size(); ++is) hs[is]->Draw("E SAME");

      TLegend *leg = new TLegend(0.50, 0.64, 0.90, 0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.030);
      for (std::size_t is = 0; is < hs.size(); ++is) {
        leg->AddEntry(hs[is], Form("%s (sumW=%.1f)", labels[is].Data(), rawIntegrals[is]), "lep");
      }
      leg->Draw();

      TLatex lab;
      lab.SetNDC();
      lab.SetTextSize(0.040);
      lab.SetTextFont(42);
      lab.DrawLatex(0.14, 0.92, combo.label);
    }

    // for (std::size_t ip = pageCombos.size(); ip < 4; ++ip) {
    //   c->cd(static_cast<int>(ip) + 1);
    //   gPad->SetLeftMargin(0.13);
    //   gPad->SetBottomMargin(0.12);
    //   gPad->SetTopMargin(0.10);
    //   TLatex omitted;
    //   omitted.SetNDC();
    //   omitted.SetTextSize(0.040);
    //   omitted.SetTextColor(kGray + 1);
    //   omitted.DrawLatex(0.18, 0.52, "Unused panel");
    // }

    c->cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.030);
    title.SetTextFont(62);
    title.DrawLatex(0.4, 0.965, Form("%s (%s)", v.name, pageTag));

    c->Print(pdfName, "pdf");
  };

  for (const auto &v : vars) {
    if (isEpCompatible(v.name)) {
      drawPage(v, epCombos, epCut, epBaseCut, "e'p", epWeights);
      drawPage(v, combos, eppCut, baseCut, "e'pp", weights);
    } else {
      drawPage(v, combos, eppCut, baseCut, "e'pp", weights);
    }
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";

  for (TFile *f : files) f->Close();
}
