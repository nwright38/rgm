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
#include <TBox.h>
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

struct ToySummary {
  double mean = 0.0;
  double stddev = 0.0;
  double p16 = 0.0;
  double p50 = 0.0;
  double p84 = 0.0;
};

void meanStddev(const vector<double>& v, double& mean, double& stddev) {
  mean = 0.0;
  stddev = 0.0;
  if (v.empty()) return;
  for (double x : v) mean += x;
  mean /= (double)v.size();
  double s2 = 0.0;
  for (double x : v) s2 += (x - mean) * (x - mean);
  stddev = std::sqrt(s2 / (double)v.size());
}

double percentile(vector<double> vals, double p) {
  if (vals.empty()) return 0.0;
  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;
  std::sort(vals.begin(), vals.end());
  const double x = p * (vals.size() - 1);
  const int i0 = static_cast<int>(std::floor(x));
  const int i1 = static_cast<int>(std::ceil(x));
  if (i0 == i1) return vals[i0];
  const double t = x - i0;
  return (1.0 - t) * vals[i0] + t * vals[i1];
}

ToySummary summarizeToyValues(const vector<double>& vals) {
  ToySummary s;
  meanStddev(vals, s.mean, s.stddev);
  s.p16 = percentile(vals, 0.16);
  s.p50 = percentile(vals, 0.50);
  s.p84 = percentile(vals, 0.84);
  return s;
}

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

void drawToySummary(TLegend* leg,
                    TH1D* hToy,
                    const ToySummary& s,
                    bool drawNominal,
                    double nominalValue) {
  const double ymax = hToy->GetMaximum();

  TBox* percentileBand = new TBox(s.p16, 0.0, s.p84, ymax);
  percentileBand->SetFillColorAlpha(kAzure + 1, 0.16);
  percentileBand->SetLineColorAlpha(kAzure + 1, 0.0);
  percentileBand->Draw("SAME");

  TBox* stddevBand = new TBox(s.mean - s.stddev, 0.0, s.mean + s.stddev, ymax);
  stddevBand->SetFillColorAlpha(kOrange + 1, 0.20);
  stddevBand->SetLineColorAlpha(kOrange + 1, 0.0);
  stddevBand->Draw("SAME");

  hToy->Draw("HIST SAME");

  TLine* lMean = new TLine(s.mean, 0.0, s.mean, ymax);
  lMean->SetLineColor(kOrange + 8);
  lMean->SetLineStyle(1);
  lMean->SetLineWidth(3);
  lMean->Draw("SAME");

  TLine* lMedian = new TLine(s.p50, 0.0, s.p50, ymax);
  lMedian->SetLineColor(kBlue + 2);
  lMedian->SetLineStyle(2);
  lMedian->SetLineWidth(3);
  lMedian->Draw("SAME");

  TLine* lNom = nullptr;
  if (drawNominal) {
    lNom = new TLine(nominalValue, 0.0, nominalValue, ymax);
    lNom->SetLineColor(kRed + 1);
    lNom->SetLineStyle(1);
    lNom->SetLineWidth(3);
    lNom->Draw("SAME");
  }

  leg->AddEntry(hToy, "Toy distribution", "lf");
  leg->AddEntry(lMean, Form("Mean %.5g", s.mean), "l");
  leg->AddEntry(stddevBand, Form("Mean #pm std dev %.5g", s.stddev), "f");
  leg->AddEntry(lMedian, Form("Median %.5g", s.p50), "l");
  leg->AddEntry(percentileBand, Form("p16-p84 %.5g-%.5g", s.p16, s.p84), "f");
  if (lNom) leg->AddEntry(lNom, Form("Nominal %.5g", nominalValue), "l");
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
  hToy->GetXaxis()->SetTitle("Value of (e,e'pp)/(e,e'p)");
  hToy->GetYaxis()->SetTitle("Counts");
  hToy->GetXaxis()->SetTitleSize(0.045);
  hToy->GetYaxis()->SetTitleSize(0.045);
  hToy->GetXaxis()->SetLabelSize(0.040);
  hToy->GetYaxis()->SetLabelSize(0.040);
  hToy->Draw("HIST");

  const ToySummary s = summarizeToyValues(toyVals);

  TLegend* leg = new TLegend(0.56, 0.66, 0.92, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.027);
  drawToySummary(leg, hToy, s, nomValid, nomRatio);
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

  c->Print(outPdf.c_str(), "pdf");
}

void drawYieldToyDistributionPage(TCanvas* c,
                                  TDirectory* histsDir,
                                  const string& outPdf,
                                  const string& taskName,
                                  const string& pageLabel,
                                  int pMissBin,
                                  int q2ValueBin) {
  const string suffix = string("_pMiss") + std::to_string(pMissBin);
  const int rootBin = q2ValueBin + 1;

  TH1D* hNom = getHist(histsDir, "nominal", taskName, suffix);
  if (!hNom) {
    c->Clear();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.DrawLatex(0.08, 0.60, Form("Missing nominal histogram: %s", taskName.c_str()));
    c->Print(outPdf.c_str(), "pdf");
    return;
  }

  const double nomVal = hNom->GetBinContent(rootBin);

  vector<string> labels = toyLabels(histsDir);
  vector<double> toyVals;
  toyVals.reserve(labels.size());

  for (const auto& lab : labels) {
    TH1D* hToy = getHist(histsDir, lab, taskName, suffix);
    if (!hToy) continue;
    toyVals.push_back(hToy->GetBinContent(rootBin));
  }

  c->Clear();
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.06);

  if (toyVals.empty()) {
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.DrawLatex(0.08, 0.60, "No valid toy values found.");
    c->Print(outPdf.c_str(), "pdf");
    return;
  }

  double vmin = *std::min_element(toyVals.begin(), toyVals.end());
  double vmax = *std::max_element(toyVals.begin(), toyVals.end());
  const double span = std::max(1e-6, vmax - vmin);
  vmin -= 0.15 * span;
  vmax += 0.15 * span;
  vmin = std::min(vmin, nomVal - 0.10 * span);
  vmax = std::max(vmax, nomVal + 0.10 * span);

  TH1D* hToyDist = new TH1D(Form("hToy_%s", taskName.c_str()), "", 40, vmin, vmax);
  for (double v : toyVals) hToyDist->Fill(v);

  hToyDist->SetLineColor(kGreen + 3);
  hToyDist->SetFillColorAlpha(kGreen + 2, 0.30);
  hToyDist->SetLineWidth(2);
  hToyDist->GetXaxis()->SetTitle(Form("Value of %s yield", pageLabel.c_str()));
  hToyDist->GetYaxis()->SetTitle("Counts");
  hToyDist->GetXaxis()->SetTitleSize(0.045);
  hToyDist->GetYaxis()->SetTitleSize(0.045);
  hToyDist->GetXaxis()->SetLabelSize(0.040);
  hToyDist->GetYaxis()->SetLabelSize(0.040);
  hToyDist->Draw("HIST");

  const ToySummary s = summarizeToyValues(toyVals);

  TLegend* leg = new TLegend(0.15, 0.6, 0.5, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  drawToySummary(leg, hToyDist, s, true, nomVal);
  leg->Draw();

  TLatex txt;
  txt.SetNDC();
  txt.SetTextSize(0.032);
  txt.DrawLatex(0.12, 0.92,
                Form("Toy diagnostics: %s in pMiss bin %d, Q2 bin %d",
                     pageLabel.c_str(), pMissBin, q2ValueBin));
  if (q2ValueBin >= 0 && q2ValueBin < (int)kQ2Edges.size() - 1) {
    txt.DrawLatex(0.12, 0.87,
                  Form("Q2 range: %.2f < Q^{2} < %.2f", kQ2Edges[q2ValueBin], kQ2Edges[q2ValueBin + 1]));
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
  drawYieldToyDistributionPage(c, histsDir, outPdf,
                               "Q2_ep_SRC_pmiss", "e,e'p", pMissBin, q2ValueBin);
  drawYieldToyDistributionPage(c, histsDir, outPdf,
                               "Q2_epp_SRC_pmiss", "e,e'pp", pMissBin, q2ValueBin);
  drawRatioToyDistributionPage(c, histsDir, outPdf, pMissBin, q2ValueBin);
  c->Print((outPdf + "]").c_str());

  std::cout << "Wrote " << outPdf << std::endl;
  f->Close();
}
