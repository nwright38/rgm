// h2d_src_fdcd.C
//
// Draw TH2D scatter plots of SRC variable pairs split across the four
// lead/recoil detector topology panels (same 2x2 grouping as overlay_default_multi).
// Operates on a single input file.
//
// Usage (data):
//   root -l -b -q 'h2d_src_fdcd.C()'
//   root -l -b -q 'h2d_src_fdcd.C("~/data/RGM_DATA/c12_src_skim.root")'
//
// Usage (sim with weight):
//   root -l -b -q 'h2d_src_fdcd.C("~/data/RGM_DATA/c12_sim_skim.root","srcTree","pdf/h2d_sim.pdf","pCM > 0","pCM > 0 && pMiss < 1. && recP < 1. && recP > .4","weight_epp*(weight_epp<200)")'

#include "OverlaySrcCommon.h"

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

struct Plot2D {
  const char *name;      // unique key (used in histogram name)
  const char *yExpr;     // TTree expression for Y axis
  const char *xExpr;     // TTree expression for X axis
  const char *yTitle;
  const char *xTitle;
  int    yBins; double yMin; double yMax;
  int    xBins; double xMin; double xMax;
};

std::vector<Plot2D> default2DPlots() {
  return {
    {"E2miss_omega",
     "E2miss", "omega",
     "E_{2,miss} [GeV]", "#omega [GeV]",
     50, -0.5,  0.5,   40, 0.5,  3.0},
    {"pMiss_recP",
     "pMiss", "recP",
     "p_{miss} [GeV/c]", "p_{recoil} [GeV/c]",
     40, 0.3,  1.1,   40, 0.2,  1.2},
    {"pMiss_pCM",
     "pMiss", "pCM",
     "p_{miss} [GeV/c]", "p_{CM} [GeV/c]",
     40, 0.3,  1.1,   40, 0.0,  1.2},
    {"pCM_pRel",
     "pCM", "pRel",
     "p_{CM} [GeV/c]", "p_{rel} [GeV/c]",
     40, 0.0,  1.2,   40, 0.0,  1.2},
    {"pMiss_pRel",
     "pMiss", "pRel",
     "p_{miss} [GeV/c]", "p_{rel} [GeV/c]",
     40, 0.3,  1.1,   40, 0.0,  1.2},
    {"theta_PleadQ_leadOverQ",
     "theta_PleadQ*180./TMath::Pi()", "leadP/qP",
     "#theta_{p_{lead},q} [deg]", "|p_{lead}| / |q|",
     36, 0.0, 60.0,   40, 0.5,  1.2},
    {"xB_mMiss",
     "xB", "mMiss",
     "x_{B}", "m_{miss} [GeV/c^{2}]",
     40, 1.0,  2.0,   50, 0.0,  2.0},
    {"theta_p_pMiss",
     "pMissTheta*180./TMath::Pi()", "pMiss",
     "#theta_{p_{miss}} [deg]", "p_{miss} [GeV/c]",
     36, 0.0, 180.0,   40, 0.3,  1.1},
    {"theta_p_pLead",
     "leadTheta*180./TMath::Pi()", "leadP",
     "#theta_{p_{lead}} [deg]", "p_{lead} [GeV/c]",
     36, 0.0, 110.0,   40, 0.8,  3.0},
    {"theta_p_pRecoil",
     "recTheta*180./TMath::Pi()", "recP",
     "#theta_{p_{recoil}} [deg]", "p_{recoil} [GeV/c]",
     36, 0.0, 180.0,   40, 0.2,  1.2},
  };
}

}  // namespace

void h2d_src_fdcd(
    const char *fileName  = "~/data/RGM_DATA/he4_src_skim_100MeV.root",
    const char *treeName  = "srcTree",
    const char *outputPdfName = "pdf/h2d_src_fdcd_he4_data.pdf",
    const char *eppCut    = "pCM > 0",
    const char *baseCut   = "pCM > 0 && pMiss < 1. && recP < 1.",
    const char *weightExpression = "(weight_epp)*(weight_epp < 200)",      // empty = unweighted (data)
    const char *pCMyTailCut = "",
    bool includeFdFd = true) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[h2d_src_fdcd] Could not open " << fileName << "\n";
    return;
  }
  TTree *tree = dynamic_cast<TTree *>(f->Get(treeName));
  if (!tree) {
    std::cerr << "[h2d_src_fdcd] Could not find tree \"" << treeName << "\"\n";
    f->Close();
    return;
  }

  const std::vector<overlay_src::DetectorCombination> combos =
      overlay_src::detectorCombinations(includeFdFd);
  const std::vector<Plot2D> plots = default2DPlots();

  TCanvas *c = new TCanvas("c_h2d_src_fdcd", "SRC 2D distributions by topology",
                            1300, 950);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  for (const auto &p : plots) {
    c->Clear();
    c->Divide(2, 2);

    for (std::size_t ic = 0; ic < combos.size() && ic < 4; ++ic) {
      const auto &combo = combos[ic];
      c->cd(static_cast<int>(ic) + 1);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.12);
      gPad->SetRightMargin(0.14);

      TH2D *h = new TH2D(Form("h2d_%s_%s", p.name, combo.key), "",
                         p.xBins, p.xMin, p.xMax,
                         p.yBins, p.yMin, p.yMax);
      h->Sumw2();
      h->SetStats(false);

      const TString boolSel =
          overlay_src::buildSelection(eppCut, baseCut, combo.cut, pCMyTailCut);
      const TString drawSel =
          overlay_src::applyWeightToSelection(weightExpression, boolSel);

      // TTree::Project with "Y:X" fills a TH2 as (xaxis=X, yaxis=Y).
      const TString projExpr = Form("%s:%s", p.yExpr, p.xExpr);
      tree->Project(h->GetName(), projExpr.Data(), drawSel.Data());

      h->GetXaxis()->SetTitle(p.xTitle);
      h->GetYaxis()->SetTitle(p.yTitle);
      h->GetXaxis()->SetTitleSize(0.050);
      h->GetYaxis()->SetTitleSize(0.050);
      h->GetXaxis()->SetLabelSize(0.043);
      h->GetYaxis()->SetLabelSize(0.043);
      h->GetYaxis()->SetTitleOffset(1.30);

      h->Draw("COLZ");

      TLatex lab;
      lab.SetNDC();
      lab.SetTextSize(0.038);
      lab.SetTextFont(42);
      lab.DrawLatex(0.15, 0.92, combo.label);

      TLatex cntLabel;
      cntLabel.SetNDC();
      cntLabel.SetTextSize(0.030);
      cntLabel.SetTextColor(kWhite);
      const double total = h->Integral(0, h->GetNbinsX() + 1, 0, h->GetNbinsY() + 1);
      cntLabel.DrawLatex(0.15, 0.85, Form("N=%.0f", total));
    }

    // blank unused panels when FD-FD is excluded (3 combos)
    for (std::size_t ip = combos.size(); ip < 4; ++ip) {
      c->cd(static_cast<int>(ip) + 1);
      TLatex omitted;
      omitted.SetNDC();
      omitted.SetTextSize(0.040);
      omitted.SetTextColor(kGray + 1);
      omitted.DrawLatex(0.18, 0.52, "Unused panel");
    }

    c->cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.027);
    title.SetTextFont(62);
    title.DrawLatex(0.10, 0.985,
                    Form("%s vs %s : by detector topology", p.yExpr, p.xExpr));

    c->Print(pdfName, "pdf");
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "[h2d_src_fdcd] Wrote " << outputPdfName << "\n";

  f->Close();
}
