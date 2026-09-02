#include "OverlaySrcCommon.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

void overlay_data_by_detector(
    const char *fileName = "~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "pdf/overlay_data_by_detector_c12_sim.pdf",
    bool normalizeToUnity = true,
    const char *eppCut = "pCM > 0",
    const char *baseCut = "pCM > 0 && pMiss < 1. && recP < 1.",
    const char *dataWeightExpression = "(weight_epp)*(weight_epp < 200)",
    const char *pCMyTailCut = "",
    bool includeFdFd = false) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile *file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Could not open file " << fileName << "\n";
    return;
  }
  TTree *tree = dynamic_cast<TTree *>(file->Get(treeName));
  if (!tree) {
    std::cerr << "Could not find tree \"" << treeName << "\"\n";
    file->Close();
    return;
  }

  const std::vector<overlay_src::DetectorCombination> combos =
      overlay_src::detectorCombinations(includeFdFd);
  const std::vector<overlay_src::PlotVar> vars = overlay_src::defaultVariables();

  TCanvas *c = new TCanvas("c_overlay_data_by_detector", "Data by detector", 1200, 900);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  for (const auto &v : vars) {
    c->Clear();
    c->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);

    std::vector<TH1D *> hs;
    std::vector<double> rawIntegrals;
    hs.reserve(combos.size());
    rawIntegrals.reserve(combos.size());
    double yMax = 0.;

    for (std::size_t ic = 0; ic < combos.size(); ++ic) {
      const auto &combo = combos[ic];
      TH1D *h = new TH1D(Form("h_data_%s_%s", v.name, combo.key), "",
                         v.bins, v.xmin, v.xmax);
      h->Sumw2();
      overlay_src::styleTopologyData(h, static_cast<int>(ic));

      const TString boolSel =
          overlay_src::buildSelection(eppCut, baseCut, combo.cut, pCMyTailCut);
      const TString weightedSel =
          overlay_src::applyWeightToSelection(dataWeightExpression, boolSel);

      tree->Project(h->GetName(), v.expr, weightedSel.Data());
      const double rawInt = h->Integral(0, h->GetNbinsX() + 1);
      rawIntegrals.push_back(rawInt);

      if (normalizeToUnity) {
        const double inRange = h->Integral();
        if (inRange > 0.) h->Scale(1. / inRange);
      }
      yMax = std::max(yMax, overlay_src::finiteMax(h->GetMaximum()));
      hs.push_back(h);
    }

    if (hs.empty()) continue;

    if (yMax <= 0.) yMax = 1.;
    hs.front()->SetMinimum(0.);
    hs.front()->SetMaximum(1.30 * yMax);
    hs.front()->GetXaxis()->SetTitle(v.xTitle);
    hs.front()->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts" : "Counts");
    hs.front()->GetXaxis()->SetTitleSize(0.050);
    hs.front()->GetYaxis()->SetTitleSize(0.050);
    hs.front()->GetXaxis()->SetLabelSize(0.043);
    hs.front()->GetYaxis()->SetLabelSize(0.043);
    hs.front()->GetYaxis()->SetTitleOffset(1.20);

    hs.front()->Draw("E");
    for (std::size_t ic = 1; ic < hs.size(); ++ic) hs[ic]->Draw("E SAME");

    TLegend *leg = new TLegend(0.60, 0.70, 0.90, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.026);
    for (std::size_t ic = 0; ic < hs.size(); ++ic) {
      leg->AddEntry(hs[ic], Form("%s (sumW=%.1f)", combos[ic].label, rawIntegrals[ic]), "lep");
    }
    leg->Draw();

    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.032);
    title.SetTextFont(62);
    title.DrawLatex(0.13, 0.94, Form("%s: data overlays by detector", v.name));

    c->Print(pdfName, "pdf");
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";
  file->Close();
}
