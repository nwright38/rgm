// plot_extraction_diagnostics.C
//
// Additive diagnostics package for extraction sanity checks.
// Does NOT replace any existing plotting scripts.
//
// Produces a multipage PDF with:
//  1) E_miss distributions in pMiss bins (2x2 pages), annotated with mean/sigma.
//  2) E_miss distributions split by pMiss and Q2 (3x3 pages), annotated with mean/sigma.
//  3) sigma_CM fit diagnostics: per-Q2-bin pcm histograms with Gaussian fit overlay,
//     and fit metrics (mu, sigma, chi2/ndf) shown on each panel.
//
// Usage:
//   root -l -b -q 'plot_extraction_diagnostics.C("Data_He_hists_6GeV.root")'
//   root -l -b -q 'plot_extraction_diagnostics.C("Data_He_hists_6GeV.root", "pdf/Q2_pmiss/extraction_diagnostics.pdf")'

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace {

const vector<double> kQ2Edges = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};

void styleHist(TH1D* h) {
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
  h->SetMarkerColor(kBlack);
}

TH1D* getNominalHist(TDirectory* nominalDir, const string& histName) {
  if (!nominalDir) return nullptr;
  TH1D* h = (TH1D*)nominalDir->Get(histName.c_str());
  return h;
}

void drawStatsBox(TH1D* h, double xNdc, double yNdc, double size = 0.030) {
  if (!h) return;
  TLatex t;
  t.SetNDC();
  t.SetTextSize(size);

  const double mean = h->GetMean();
  const double sigma = h->GetStdDev();
  const double meanErr = h->GetMeanError();
  const double sigmaErr = h->GetStdDevError();

  t.DrawLatex(xNdc, yNdc,
              Form("#mu = %.4f #pm %.4f", mean, meanErr));
  t.DrawLatex(xNdc, yNdc - 0.055,
              Form("#sigma = %.4f #pm %.4f", sigma, sigmaErr));
}

void drawEmissGuideLines(TH1D* h) {
  if (!h) return;

  const double mean = h->GetMean();
  const double sigma = h->GetStdDev();
  const double yTop = h->GetMaximum();

  TLine* lMean = new TLine(mean, 0.0, mean, yTop * 1.05);
  lMean->SetLineColor(kMagenta + 1);
  lMean->SetLineStyle(3);
  lMean->SetLineWidth(2);
  lMean->Draw("SAME");

  TLine* lMinus = new TLine(mean - sigma, 0.0, mean - sigma, yTop * 1.05);
  lMinus->SetLineColor(kAzure + 1);
  lMinus->SetLineStyle(3);
  lMinus->SetLineWidth(2);
  lMinus->Draw("SAME");

  TLine* lPlus = new TLine(mean + sigma, 0.0, mean + sigma, yTop * 1.05);
  lPlus->SetLineColor(kAzure + 1);
  lPlus->SetLineStyle(3);
  lPlus->SetLineWidth(2);
  lPlus->Draw("SAME");
}

void drawPanelLabel(const string& label, double xNdc = 0.15, double yNdc = 0.88,
                    double size = 0.038) {
  TLatex t;
  t.SetNDC();
  t.SetTextSize(size);
  t.DrawLatex(xNdc, yNdc, label.c_str());
}

void drawPageTitle(TCanvas* c, const string& title) {
  c->cd(0);
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.035);
  t.SetTextFont(62);
  t.DrawLatex(0.06, 0.98, title.c_str());
}

void drawEmissPMissPage(TCanvas* c,
                        TDirectory* nominalDir,
                        const string& histPrefix,
                        const string& outPdf,
                        const string& pageTitle,
                        const string& xTitle) {
  c->Clear();
  c->Divide(2, 2);

  for (int p = 0; p < 4; p++) {
    c->cd(p + 1);
    gPad->SetTicks(1, 1);

    const string hname = histPrefix + Form("_pMiss%d", p);
    TH1D* h = getNominalHist(nominalDir, hname);
    if (!h) {
      drawPanelLabel("Missing: " + hname, 0.08, 0.5, 0.030);
      continue;
    }

    styleHist(h);
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle("Counts");
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->Draw("E");
    drawEmissGuideLines(h);

    const string pmissLabel = Form("pMiss bin %d", p);
    drawPanelLabel(pmissLabel, 0.16, 0.88, 0.040);
    drawStatsBox(h, 0.16, 0.80, 0.033);
  }

  drawPageTitle(c, pageTitle);
  c->Print(outPdf.c_str(), "pdf");
}

void drawEmissQ2DeepDivePages(TCanvas* c,
                              TDirectory* nominalDir,
                              const string& histPrefix,
                              const string& outPdf,
                              const string& pageBaseTitle,
                              const string& xTitle) {
  for (int p = 0; p < 4; p++) {
    c->Clear();
    c->Divide(3, 3);

    for (int q = 0; q < 7; q++) {
      c->cd(q + 1);
      gPad->SetTicks(1, 1);

      const string hname = histPrefix + Form("_pMiss%d_Q2%d", p, q);
      TH1D* h = getNominalHist(nominalDir, hname);
      if (!h) {
        drawPanelLabel("Missing: " + hname, 0.08, 0.5, 0.028);
        continue;
      }

      styleHist(h);
      h->GetXaxis()->SetTitle(xTitle.c_str());
      h->GetYaxis()->SetTitle("Counts");
      h->GetXaxis()->SetTitleSize(0.050);
      h->GetYaxis()->SetTitleSize(0.050);
      h->GetXaxis()->SetLabelSize(0.040);
      h->GetYaxis()->SetLabelSize(0.040);
      h->Draw("E");
      drawEmissGuideLines(h);

      const string q2Label = Form("%.2f < Q^{2} < %.2f", kQ2Edges[q], kQ2Edges[q + 1]);
      drawPanelLabel(q2Label, 0.12, 0.88, 0.030);
      drawStatsBox(h, 0.12, 0.79, 0.027);
    }

    c->cd(8);
    gPad->Clear();
    drawPanelLabel("Per-panel: nominal histogram", 0.08, 0.70, 0.035);
    drawPanelLabel("numbers are direct TH1 mean/stddev", 0.08, 0.58, 0.032);
    drawPanelLabel("(not fit-derived in this section)", 0.08, 0.46, 0.032);

    c->cd(9);
    gPad->Clear();
    drawPanelLabel(Form("pMiss bin %d", p), 0.20, 0.65, 0.055);

    drawPageTitle(c, pageBaseTitle + Form(" (pMiss bin %d)", p));
    c->Print(outPdf.c_str(), "pdf");
  }
}

void drawSigmaFitPages(TCanvas* c,
                       TDirectory* nominalDir,
                       const string& histPrefix,
                       const string& outPdf,
                       const string& pageTitle,
                       const string& xTitle,
                       double fitMin,
                       double fitMax) {
  c->Clear();
  c->Divide(3, 3);

  for (int q = 0; q < 7; q++) {
    c->cd(q + 1);
    gPad->SetTicks(1, 1);

    const string hname = histPrefix + Form("_Q2%d", q);
    TH1D* h = getNominalHist(nominalDir, hname);
    if (!h) {
      drawPanelLabel("Missing: " + hname, 0.08, 0.5, 0.028);
      continue;
    }

    styleHist(h);
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle("Counts");
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->Draw("E");

    TF1* gfit = nullptr;
    if (h->GetEntries() > 20) {
      gfit = new TF1(Form("gfit_%s_%d", histPrefix.c_str(), q), "gaus", fitMin, fitMax);
      gfit->SetLineColor(kRed + 1);
      gfit->SetLineWidth(2);
      h->Fit(gfit, "QRS");
      gfit->Draw("SAME");
    }

    const string q2Label = Form("%.2f < Q^{2} < %.2f", kQ2Edges[q], kQ2Edges[q + 1]);
    drawPanelLabel(q2Label, 0.11, 0.88, 0.030);

    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.026);

    if (gfit) {
      const double mu = gfit->GetParameter(1);
      const double muErr = gfit->GetParError(1);
      const double sigma = gfit->GetParameter(2);
      const double sigmaErr = gfit->GetParError(2);
      const double chi2 = gfit->GetChisquare();
      const double ndf = gfit->GetNDF();
      const double chi2ndf = (ndf > 0.0) ? (chi2 / ndf) : 0.0;

      t.DrawLatex(0.11, 0.79, Form("fit: #mu = %.4f #pm %.4f", mu, muErr));
      t.DrawLatex(0.11, 0.72, Form("fit: #sigma = %.4f #pm %.4f", sigma, sigmaErr));
      t.DrawLatex(0.11, 0.65, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      t.DrawLatex(0.11, 0.58, Form("fit range: [%.2f, %.2f]", fitMin, fitMax));
    } else {
      t.DrawLatex(0.11, 0.72, "fit skipped (low statistics)");
    }
  }

  c->cd(8);
  gPad->Clear();
  drawPanelLabel("Gaussian-fit diagnostics", 0.10, 0.72, 0.040);
  drawPanelLabel("Same fit window as extraction", 0.10, 0.58, 0.034);
  drawPanelLabel(Form("[%.2f, %.2f]", fitMin, fitMax), 0.10, 0.46, 0.036);

  c->cd(9);
  gPad->Clear();
  drawPanelLabel("Red curve = Gaussian fit", 0.10, 0.66, 0.036);
  drawPanelLabel("Black points = nominal histogram", 0.10, 0.54, 0.034);

  drawPageTitle(c, pageTitle);
  c->Print(outPdf.c_str(), "pdf");
}

void drawPcmXYOverlayByQ2(TCanvas* c,
                          TDirectory* nominalDir,
                          const string& outPdf,
                          bool normalizeToArea,
                          const string& pageTitle) {
  c->Clear();
  c->Divide(3, 3);

  for (int q = 0; q < 7; q++) {
    c->cd(q + 1);
    gPad->SetTicks(1, 1);

    const string hxName = Form("nominal_pcmx_epp_SRC_Q2_Q2%d", q);
    const string hyName = Form("nominal_pcmy_epp_SRC_Q2_Q2%d", q);
    TH1D* hxSrc = getNominalHist(nominalDir, hxName);
    TH1D* hySrc = getNominalHist(nominalDir, hyName);

    if (!hxSrc || !hySrc) {
      drawPanelLabel("Missing pcmx/pcmy hist(s)", 0.08, 0.5, 0.028);
      continue;
    }

    TH1D* hx = (TH1D*)hxSrc->Clone(Form("hx_overlay_q2_%d", q));
    TH1D* hy = (TH1D*)hySrc->Clone(Form("hy_overlay_q2_%d", q));
    hx->SetDirectory(nullptr);
    hy->SetDirectory(nullptr);

    if (normalizeToArea) {
      const double ix = hx->Integral();
      const double iy = hy->Integral();
      if (ix > 0.0) hx->Scale(1.0 / ix);
      if (iy > 0.0) hy->Scale(1.0 / iy);
    }

    hx->SetLineColor(kBlue + 1);
    hx->SetMarkerColor(kBlue + 1);
    hx->SetLineWidth(2);
    hx->SetMarkerStyle(20);
    hx->SetMarkerSize(0.5);

    hy->SetLineColor(kRed + 1);
    hy->SetMarkerColor(kRed + 1);
    hy->SetLineWidth(2);
    hy->SetMarkerStyle(21);
    hy->SetMarkerSize(0.5);

    const double ymax = std::max(hx->GetMaximum(), hy->GetMaximum());
    hx->SetMaximum(ymax * 1.25);

    hx->GetXaxis()->SetTitle("p_{C.M.} component [GeV]");
    hx->GetYaxis()->SetTitle(normalizeToArea ? "A.U." : "Counts");
    hx->GetXaxis()->SetTitleSize(0.050);
    hx->GetYaxis()->SetTitleSize(0.050);
    hx->GetXaxis()->SetLabelSize(0.040);
    hx->GetYaxis()->SetLabelSize(0.040);

    hx->Draw("HIST E");
    hy->Draw("HIST E SAME");

    drawPanelLabel(Form("%.2f < Q^{2} < %.2f", kQ2Edges[q], kQ2Edges[q + 1]),
                   0.10, 0.88, 0.030);

    TLegend* leg = new TLegend(0.58, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->AddEntry(hx, "pcmx", "lep");
    leg->AddEntry(hy, "pcmy", "lep");
    leg->Draw();
  }

  c->cd(8);
  gPad->Clear();
  drawPanelLabel(normalizeToArea ? "Overlaid with unit-area normalization"
                                 : "Overlaid raw counts",
                 0.08, 0.65, 0.036);
  drawPanelLabel("Data nominal histograms in each Q2 bin", 0.08, 0.52, 0.033);

  c->cd(9);
  gPad->Clear();
  drawPanelLabel("Blue: pcmx", 0.12, 0.66, 0.040);
  drawPanelLabel("Red: pcmy", 0.12, 0.54, 0.040);

  drawPageTitle(c, pageTitle);
  c->Print(outPdf.c_str(), "pdf");
}

void drawSingleComponentAllQ2Overlay(TCanvas* c,
                                     TDirectory* nominalDir,
                                     const string& outPdf,
                                     const string& histPrefix,
                                     const string& xTitle,
                                     const string& pageTitle,
                                     bool normalizeToArea = true) {
  c->Clear();
  c->Divide(1, 1);
  c->cd(1);
  gPad->SetTicks(1, 1);

  vector<TH1D*> hs;
  hs.reserve(7);
  double ymax = 0.0;

  for (int q = 0; q < 7; q++) {
    const string hname = histPrefix + Form("_Q2%d", q);
    TH1D* hsrc = getNominalHist(nominalDir, hname);
    if (!hsrc) continue;

    TH1D* h = (TH1D*)hsrc->Clone(Form("overlay_%s_q2_%d", histPrefix.c_str(), q));
    h->SetDirectory(nullptr);

    if (normalizeToArea) {
      const double integ = h->Integral();
      if (integ > 0.0) h->Scale(1.0 / integ);
    }

    int color = kBlue + q;
    if (q >= 4) color = kRed + (q - 4);
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(20 + (q % 4));
    h->SetMarkerSize(0.5);

    ymax = std::max(ymax, h->GetMaximum());
    hs.push_back(h);
  }

  if (hs.empty()) {
    drawPanelLabel("No histograms found for all-Q2 overlay.", 0.12, 0.55, 0.040);
    drawPageTitle(c, pageTitle);
    c->Print(outPdf.c_str(), "pdf");
    return;
  }

  hs[0]->SetMaximum(ymax * 1.25);
  hs[0]->GetXaxis()->SetTitle(xTitle.c_str());
  hs[0]->GetYaxis()->SetTitle(normalizeToArea ? "A.U." : "Counts");
  hs[0]->GetXaxis()->SetTitleSize(0.045);
  hs[0]->GetYaxis()->SetTitleSize(0.045);
  hs[0]->GetXaxis()->SetLabelSize(0.038);
  hs[0]->GetYaxis()->SetLabelSize(0.038);
  hs[0]->Draw("HIST E");

  for (size_t i = 1; i < hs.size(); i++) {
    hs[i]->Draw("HIST E SAME");
  }

  TLegend* leg = new TLegend(0.58, 0.55, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  for (int q = 0; q < (int)hs.size(); q++) {
    leg->AddEntry(hs[q], Form("%.2f < Q^{2} < %.2f", kQ2Edges[q], kQ2Edges[q + 1]), "lep");
  }
  leg->Draw();

  drawPageTitle(c, pageTitle);
  c->Print(outPdf.c_str(), "pdf");
}

}  // namespace

void plot_extraction_diagnostics(const char* inFileName,
                                 const char* outPdfName = "pdf/Q2_pmiss/extraction_diagnostics.pdf") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* f = TFile::Open(inFileName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open input file: " << inFileName << std::endl;
    return;
  }

  TDirectory* nominalDir = f->GetDirectory("hists/nominal");
  if (!nominalDir) {
    std::cerr << "Missing directory hists/nominal in: " << inFileName << std::endl;
    f->Close();
    return;
  }

  TCanvas* c = new TCanvas("c_diag", "Extraction diagnostics", 1400, 1000);

  string outPdf(outPdfName);
  c->Print((outPdf + "[").c_str());

  // Tier 1: pMiss-binned E_miss overview pages.
  drawEmissPMissPage(
      c, nominalDir,
      "nominal_E1miss_ep_SRC_pmiss",
      outPdf,
      "E1miss (e'p), pMiss-binned nominal distributions",
      "E_{1,miss} [GeV]");

  drawEmissPMissPage(
      c, nominalDir,
      "nominal_E1miss_epp_SRC_pmiss",
      outPdf,
      "E1miss (e'pp), pMiss-binned nominal distributions",
      "E_{1,miss} [GeV]");

  drawEmissPMissPage(
      c, nominalDir,
      "nominal_E2miss_epp_SRC_pmiss",
      outPdf,
      "E2miss (e'pp), pMiss-binned nominal distributions",
      "E_{2,miss} [GeV]");

  // Tier 2: Q2-split deep-dive pages (one page per pMiss bin).
  drawEmissQ2DeepDivePages(
      c, nominalDir,
      "nominal_E1miss_ep_SRC_pmiss_Q2",
      outPdf,
      "E1miss (e'p), Q2-split diagnostics",
      "E_{1,miss} [GeV]");

  drawEmissQ2DeepDivePages(
      c, nominalDir,
      "nominal_E1miss_epp_SRC_pmiss_Q2",
      outPdf,
      "E1miss (e'pp), Q2-split diagnostics",
      "E_{1,miss} [GeV]");

  drawEmissQ2DeepDivePages(
      c, nominalDir,
      "nominal_E2miss_epp_SRC_pmiss_Q2",
      outPdf,
      "E2miss (e'pp), Q2-split diagnostics",
      "E_{2,miss} [GeV]");

  // sigma_CM fit diagnostics (same fit windows used in ExtractFitQuantities.C).
  drawSigmaFitPages(
      c, nominalDir,
      "nominal_pcmx_epp_SRC_Q2",
      outPdf,
      "sigma_{x,C.M.} extraction diagnostics: per-Q2 Gaussian fits",
      "p_{C.M.,x} [GeV]",
      -0.2, 0.22);

  drawSigmaFitPages(
      c, nominalDir,
      "nominal_pcmy_epp_SRC_Q2",
      outPdf,
      "sigma_{y,C.M.} extraction diagnostics: per-Q2 Gaussian fits",
      "p_{C.M.,y} [GeV]",
      -0.2, 0.2);

    // pcmx vs pcmy overlays in each Q2 bin (raw and normalized).
    drawPcmXYOverlayByQ2(
      c, nominalDir, outPdf, false,
      "pcmx/pcmy overlays by Q2 bin (raw counts)");

    drawPcmXYOverlayByQ2(
      c, nominalDir, outPdf, true,
      "pcmx/pcmy overlays by Q2 bin (unit-area normalized)");

  // One-page-per-component overlays with all Q2 bins on one axis.
  drawSingleComponentAllQ2Overlay(
      c, nominalDir, outPdf,
      "nominal_pcmx_epp_SRC_Q2",
      "p_{C.M.,x} [GeV]",
      "pcmx all-Q2 overlay (single panel, unit-area normalized)",
      true);

  drawSingleComponentAllQ2Overlay(
      c, nominalDir, outPdf,
      "nominal_pcmy_epp_SRC_Q2",
      "p_{C.M.,y} [GeV]",
      "pcmy all-Q2 overlay (single panel, unit-area normalized)",
      true);

  c->Print((outPdf + "]").c_str());
  cout << "Wrote " << outPdf << endl;

  f->Close();
}
