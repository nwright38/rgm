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

void overlay_q2_by_detector(
    const char *fileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "pdf/overlay_q2_by_detector_c12_data.pdf",
    bool normalizeToUnity = true,
    const char *eppCut = "pCM > 0",
    const char *baseCut = "pCM > 0 && pMiss < 1. && recP < 1.",
    const char *dataWeightExpression = "(weight_epp)",
    const char *pCMyTailCut = "",
    const char *q2EdgesCsv = "1.5,1.80,2.10,2.40,3.00,5.0",
    bool autoRebin = true,
    int maxRebinFactor = 2,
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

  std::vector<TString> q2EdgeStrings = overlay_src::splitCsvToTStrings(q2EdgesCsv);
  std::vector<double> q2Edges;
  q2Edges.reserve(q2EdgeStrings.size());
  for (const auto &edge : q2EdgeStrings) q2Edges.push_back(edge.Atof());
  if (q2Edges.size() < 2) {
    std::cerr << "Need at least two Q2 edges.\n";
    file->Close();
    return;
  }

  const std::vector<overlay_src::DetectorCombination> combos =
      overlay_src::detectorCombinations(includeFdFd);
  const std::vector<overlay_src::PlotVar> vars = overlay_src::defaultVariables();

  TCanvas *c = new TCanvas("c_overlay_q2_by_detector", "Q2 by detector", 1300, 950);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  for (const auto &v : vars) {
    c->Clear();
    c->Divide(2, 2);

    double representativeEntries = 0.;
    int representativeCount = 0;
    for (const auto &combo : combos) {
      for (std::size_t iq2 = 0; iq2 + 1 < q2Edges.size(); ++iq2) {
        const TString q2Cut = Form("(Q2>=%g) && (Q2<%g)", q2Edges[iq2], q2Edges[iq2 + 1]);
        const TString sel =
            overlay_src::buildSelection(eppCut, baseCut, combo.cut, pCMyTailCut, q2Cut.Data());
        const Long64_t nSel = tree->GetEntries(sel.Data());
        if (nSel > 0) {
          representativeEntries += static_cast<double>(nSel);
          ++representativeCount;
        }
      }
    }
    if (representativeCount > 0) {
      representativeEntries /= static_cast<double>(representativeCount);
    }
    const int rebinFactor = autoRebin
                                ? overlay_src::chooseAutoRebinFactor(v.bins, representativeEntries,
                                                                     8.0, 5, maxRebinFactor)
                                : 1;

    for (std::size_t ic = 0; ic < combos.size() && ic < 4; ++ic) {
      const auto &combo = combos[ic];
      c->cd(static_cast<int>(ic) + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);

      std::vector<TH1D *> hs;
      std::vector<double> rawIntegrals;
      hs.reserve(q2Edges.size() - 1);
      rawIntegrals.reserve(q2Edges.size() - 1);
      double yMax = 0.;

      for (std::size_t iq2 = 0; iq2 + 1 < q2Edges.size(); ++iq2) {
        TH1D *h = new TH1D(Form("h_q2_%s_%s_bin%zu", v.name, combo.key, iq2), "",
                           v.bins, v.xmin, v.xmax);
        h->Sumw2();
        overlay_src::styleTopologyData(h, static_cast<int>(iq2));

        const TString q2Cut = Form("(Q2>=%g) && (Q2<%g)", q2Edges[iq2], q2Edges[iq2 + 1]);
        const TString sel =
            overlay_src::buildSelection(eppCut, baseCut, combo.cut, pCMyTailCut, q2Cut.Data());
        const TString weightedSel =
            overlay_src::applyWeightToSelection(dataWeightExpression, sel);

        tree->Project(h->GetName(), v.expr, weightedSel.Data());
        const double rawInt = h->Integral(0, h->GetNbinsX() + 1);
        rawIntegrals.push_back(rawInt);

        if (rebinFactor > 1) h->Rebin(rebinFactor);
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
      hs.front()->SetMaximum(1.25 * yMax);
      hs.front()->GetXaxis()->SetTitle(v.xTitle);
      hs.front()->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts" : "Counts");
      hs.front()->GetXaxis()->SetTitleSize(0.050);
      hs.front()->GetYaxis()->SetTitleSize(0.050);
      hs.front()->GetXaxis()->SetLabelSize(0.043);
      hs.front()->GetYaxis()->SetLabelSize(0.043);
      hs.front()->GetYaxis()->SetTitleOffset(1.20);

      hs.front()->Draw("E");
      for (std::size_t iq2 = 1; iq2 < hs.size(); ++iq2) hs[iq2]->Draw("E SAME");

      TLegend *leg = new TLegend(0.56, 0.52, 0.90, 0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.022);
      for (std::size_t iq2 = 0; iq2 < hs.size(); ++iq2) {
        leg->AddEntry(hs[iq2],
                      Form("%.2f<Q^{2}<%.2f (sumW=%.1f)", q2Edges[iq2], q2Edges[iq2 + 1],
                           rawIntegrals[iq2]),
                      "lep");
      }
      leg->Draw();

      TLatex lab;
      lab.SetNDC();
      lab.SetTextSize(0.040);
      lab.SetTextFont(42);
      lab.DrawLatex(0.14, 0.92, combo.label);
      if (rebinFactor > 1) {
        lab.SetTextSize(0.030);
        lab.DrawLatex(0.14, 0.86, Form("auto rebin x%d", rebinFactor));
      }
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
    title.DrawLatex(0.10, 0.985, Form("%s: Q2-bin overlays by detector", v.name));

    c->Print(pdfName, "pdf");
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";
  file->Close();
}
