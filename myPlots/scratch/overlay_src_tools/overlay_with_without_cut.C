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

void overlay_with_without_cut(
    const char *fileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "overlay_with_without_cut.pdf",
    const char *eppCut = "pCM > 0",
    const char *baseCut = "pCM > 0 && pMiss < 1. && recP < 1.",
    const char *extraCut = "pCMz > .4",
    const char *extraCutLabel = "pCMz > .4",
    const char *dataWeightExpression = "(weight_epp)",
    const char *pCMyTailCut = "",
    bool includeFdFd = false,
    bool normalizeToUnity = false) {

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

  TCanvas *c = new TCanvas("c_overlay_with_without_cut", "With vs without cut", 1300, 950);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  for (const auto &v : vars) {
    c->Clear();
    c->Divide(2, 2);

    for (std::size_t ic = 0; ic < combos.size() && ic < 4; ++ic) {
      const auto &combo = combos[ic];
      c->cd(static_cast<int>(ic) + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);

      TH1D *hNoCut = new TH1D(Form("h_nocut_%s_%s", v.name, combo.key), "",
                              v.bins, v.xmin, v.xmax);
      TH1D *hWithCut = new TH1D(Form("h_withcut_%s_%s", v.name, combo.key), "",
                                v.bins, v.xmin, v.xmax);
      hNoCut->Sumw2();
      hWithCut->Sumw2();
      overlay_src::styleTopologyData(hNoCut, 0);
      overlay_src::styleTopologyData(hWithCut, 1);

      const TString noCutSel =
          overlay_src::buildSelection(eppCut, baseCut, combo.cut, pCMyTailCut);
      const TString withCutSel =
          overlay_src::buildSelection(eppCut, baseCut, combo.cut, pCMyTailCut,
                                      nullptr, extraCut);

      const TString noCutWeighted =
          overlay_src::applyWeightToSelection(dataWeightExpression, noCutSel);
      const TString withCutWeighted =
          overlay_src::applyWeightToSelection(dataWeightExpression, withCutSel);

      tree->Project(hNoCut->GetName(), v.expr, noCutWeighted.Data());
      tree->Project(hWithCut->GetName(), v.expr, withCutWeighted.Data());

      const double rawNoCut = hNoCut->Integral(0, hNoCut->GetNbinsX() + 1);
      const double rawWithCut = hWithCut->Integral(0, hWithCut->GetNbinsX() + 1);

      if (normalizeToUnity) {
        const double noCutInRange = hNoCut->Integral();
        const double withCutInRange = hWithCut->Integral();
        if (noCutInRange > 0.) hNoCut->Scale(1. / noCutInRange);
        if (withCutInRange > 0.) hWithCut->Scale(1. / withCutInRange);
      }

      double yMax = std::max(overlay_src::finiteMax(hNoCut->GetMaximum()),
                             overlay_src::finiteMax(hWithCut->GetMaximum()));
      if (yMax <= 0.) yMax = 1.;
      hNoCut->SetMinimum(0.);
      hNoCut->SetMaximum(1.25 * yMax);
      hNoCut->GetXaxis()->SetTitle(v.xTitle);
      hNoCut->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts" : "Counts");
      hNoCut->GetXaxis()->SetTitleSize(0.050);
      hNoCut->GetYaxis()->SetTitleSize(0.050);
      hNoCut->GetXaxis()->SetLabelSize(0.043);
      hNoCut->GetYaxis()->SetLabelSize(0.043);
      hNoCut->GetYaxis()->SetTitleOffset(1.20);

      hNoCut->Draw("E");
      hWithCut->Draw("E SAME");

      TLegend *leg = new TLegend(0.50, 0.68, 0.90, 0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.028);
      leg->AddEntry(hNoCut, Form("No extra cut (sumW=%.1f)", rawNoCut), "lep");
      leg->AddEntry(hWithCut,
                    Form("%s (sumW=%.1f)",
                         (extraCutLabel && extraCutLabel[0] != '\0') ? extraCutLabel
                                                                       : "With extra cut",
                         rawWithCut),
                    "lep");
      leg->Draw();

      TLatex lab;
      lab.SetNDC();
      lab.SetTextSize(0.040);
      lab.SetTextFont(42);
      lab.DrawLatex(0.14, 0.92, combo.label);
    }

    for (std::size_t ip = combos.size(); ip < 4; ++ip) {
      c->cd(static_cast<int>(ip) + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      TLatex omitted;
      omitted.SetNDC();
      omitted.SetTextSize(0.040);
      omitted.SetTextColor(kGray + 1);
      omitted.DrawLatex(0.18, 0.52, includeFdFd ? "Unused panel" : "FD+FD omitted");
    }

    c->cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.028);
    title.SetTextFont(62);
    title.DrawLatex(0.10, 0.985, Form("%s: with-vs-without cut", v.name));

    c->Print(pdfName, "pdf");
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";
  file->Close();
}
