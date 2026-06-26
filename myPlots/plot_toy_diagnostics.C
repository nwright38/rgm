// plot_toy_diagnostics.C
//
// Additive toy-ensemble diagnostics (does not replace existing workflows).
//
// Primary page implemented:
//   - Distribution of toy values for (e,e'pp)/(e,e'p) in a selected
//     pMiss bin and Q2 value bin, with nominal central value overlay.
//
// Usage:
//   root -l -b -q 'plot_toy_diagnostics.C("Data_He_hists_6GeV.root")'
//   root -l -b -q 'plot_toy_diagnostics.C("Data_He_hists_6GeV.root", "pdf/Q2_pmiss/toy_diagnostics.pdf", 3, 6)'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

using std::string;
using std::vector;

namespace {

const vector<double> kQ2Edges = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};

vector<string> toyLabels(TDirectory* histsDir) {
  vector<string> labels;
  if (!histsDir) return labels;

  TList* keys = histsDir->GetListOfKeys();
  if (!keys) return labels;

  TIter next(keys);
  while (TObject* obj = next()) {
    TKey* key = dynamic_cast<TKey*>(obj);
    if (!key) continue;
    const string name = key->GetName();
    const string cls = key->GetClassName();
    if (cls == "TDirectoryFile" && name.rfind("toy_", 0) == 0) {
      labels.push_back(name);
    }
  }
  std::sort(labels.begin(), labels.end());
  return labels;
}

TH1D* getHist(TDirectory* histsDir,
              const string& label,
              const string& task,
              const string& suffix) {
  if (!histsDir) return nullptr;
  TDirectory* d = histsDir->GetDirectory(label.c_str());
  if (!d) return nullptr;
  const string hname = label + "_" + task + suffix;
  return dynamic_cast<TH1D*>(d->Get(hname.c_str()));
}

void drawRatioToyDistributionPage(TCanvas* c,
                                  TDirectory* histsDir,
                                  const string& outPdf,
                                  int pMissBin,
                                  int q2ValueBin) {
  const string numeratorTask = "Q2_epp_SRC_pmiss";
  const string denominatorTask = "Q2_ep_SRC_pmiss";
  const string suffix = string("_pMiss") + std::to_string(pMissBin);

  // Q2 value bin is 0-based from analysis code, TH1 bins are 1-based.
  const int rootBin = q2ValueBin + 1;

  TH1D* hNumNom = getHist(histsDir, "nominal", numeratorTask, suffix);
  TH1D* hDenNom = getHist(histsDir, "nominal", denominatorTask, suffix);
  if (!hNumNom || !hDenNom) {
    c->Clear();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.DrawLatex(0.08, 0.60, "Missing nominal histograms for ratio toy diagnostics.");
    t.DrawLatex(0.08, 0.52, Form("Expected: nominal_%s%s and nominal_%s%s",
                                 numeratorTask.c_str(), suffix.c_str(),
                                 denominatorTask.c_str(), suffix.c_str()));
    c->Print(outPdf.c_str(), "pdf");
    return;
  }

  const double nomNum = hNumNom->GetBinContent(rootBin);
  const double nomDen = hDenNom->GetBinContent(rootBin);
  const bool nomValid = (nomDen > 0.0);
  const double nomRatio = nomValid ? (nomNum / nomDen) : 0.0;

  vector<string> labels = toyLabels(histsDir);
  vector<double> toyVals;
  toyVals.reserve(labels.size());

  for (const auto& lab : labels) {
    TH1D* hNum = getHist(histsDir, lab, numeratorTask, suffix);
    TH1D* hDen = getHist(histsDir, lab, denominatorTask, suffix);
    if (!hNum || !hDen) continue;
    const double num = hNum->GetBinContent(rootBin);
    const double den = hDen->GetBinContent(rootBin);
    if (den <= 0.0) continue;
    toyVals.push_back(num / den);
  }

  c->Clear();
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.06);

  if (toyVals.empty()) {
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.DrawLatex(0.08, 0.60, "No valid toy ratio values found.");
    c->Print(outPdf.c_str(), "pdf");
    return;
  }

  double vmin = *std::min_element(toyVals.begin(), toyVals.end());
  double vmax = *std::max_element(toyVals.begin(), toyVals.end());
  const double span = std::max(1e-6, vmax - vmin);
  vmin -= 0.15 * span;
  vmax += 0.15 * span;
  if (nomValid) {
    vmin = std::min(vmin, nomRatio - 0.10 * span);
    vmax = std::max(vmax, nomRatio + 0.10 * span);
  }

  TH1D* hToy = new TH1D("hToyRatio", "", 40, vmin, vmax);
  for (double v : toyVals) hToy->Fill(v);

  hToy->SetLineColor(kAzure + 2);
  hToy->SetFillColorAlpha(kAzure + 1, 0.30);
  hToy->SetLineWidth(2);
  hToy->GetXaxis()->SetTitle("Toy value of (e,e^{#prime}pp)/(e,e^{#prime}p)");
  hToy->GetYaxis()->SetTitle("Number of toys");
  hToy->GetXaxis()->SetTitleSize(0.045);
  hToy->GetYaxis()->SetTitleSize(0.045);
  hToy->GetXaxis()->SetLabelSize(0.040);
  hToy->GetYaxis()->SetLabelSize(0.040);
  hToy->Draw("HIST");

  const double toyMean = hToy->GetMean();
  const double toyRms = hToy->GetStdDev();

  TLine* lToyMean = new TLine(toyMean, 0.0, toyMean, hToy->GetMaximum() * 1.05);
  lToyMean->SetLineColor(kBlue + 2);
  lToyMean->SetLineStyle(2);
  lToyMean->SetLineWidth(3);
  lToyMean->Draw("SAME");

  TLine* lNom = nullptr;
  if (nomValid) {
    lNom = new TLine(nomRatio, 0.0, nomRatio, hToy->GetMaximum() * 1.05);
    lNom->SetLineColor(kRed + 1);
    lNom->SetLineStyle(2);
    lNom->SetLineWidth(3);
    lNom->Draw("SAME");
  }

  TLegend* leg = new TLegend(0.58, 0.70, 0.92, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->AddEntry(hToy, "Toy distribution", "lf");
  leg->AddEntry(lToyMean, "Toy mean", "l");
  if (lNom) leg->AddEntry(lNom, "Nominal central value", "l");
  leg->Draw();

  TLatex txt;
  txt.SetNDC();
  txt.SetTextSize(0.032);
  txt.DrawLatex(0.12, 0.92,
                Form("Toy diagnostics: ratio in pMiss bin %d, Q2 bin %d", pMissBin, q2ValueBin));
  if (q2ValueBin >= 0 && q2ValueBin < (int)kQ2Edges.size() - 1) {
    txt.DrawLatex(0.12, 0.87,
                  Form("Q2 range: %.2f < Q^{2} < %.2f", kQ2Edges[q2ValueBin], kQ2Edges[q2ValueBin + 1]));
  }

  txt.DrawLatex(0.12, 0.81, Form("N toys used: %d", (int)toyVals.size()));
  txt.DrawLatex(0.12, 0.76, Form("Toy mean: %.5f", toyMean));
  txt.DrawLatex(0.12, 0.71, Form("Toy RMS: %.5f", toyRms));
  if (nomValid) {
    const double pull = (toyRms > 0.0) ? ((nomRatio - toyMean) / toyRms) : 0.0;
    txt.DrawLatex(0.12, 0.66, Form("Nominal: %.5f", nomRatio));
    txt.DrawLatex(0.12, 0.61, Form("(Nominal - ToyMean)/ToyRMS: %.3f", pull));
  } else {
    txt.DrawLatex(0.12, 0.66, "Nominal ratio undefined (denominator <= 0)");
  }

  c->Print(outPdf.c_str(), "pdf");
}

}  // namespace

void plot_toy_diagnostics(const char* inFileName,
                          const char* outPdfName = "pdf/Q2_pmiss/toy_diagnostics.pdf",
                          int pMissBin = 3,
                          int q2ValueBin = 6) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* f = TFile::Open(inFileName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open input file: " << inFileName << std::endl;
    return;
  }

  TDirectory* histsDir = f->GetDirectory("hists");
  if (!histsDir) {
    std::cerr << "Missing hists directory in: " << inFileName << std::endl;
    f->Close();
    return;
  }

  TCanvas* c = new TCanvas("c_toy", "Toy diagnostics", 1100, 800);
  string outPdf(outPdfName);

  c->Print((outPdf + "[").c_str());
  drawRatioToyDistributionPage(c, histsDir, outPdf, pMissBin, q2ValueBin);
  c->Print((outPdf + "]").c_str());

  std::cout << "Wrote " << outPdf << std::endl;
  f->Close();
}
