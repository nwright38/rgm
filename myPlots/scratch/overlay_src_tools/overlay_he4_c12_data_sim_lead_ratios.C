#include "OverlaySrcCommon.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

struct SampleConfig {
  TString path;
  TString weight;
};

struct HistSummary {
  TH1D *raw;
  double rawIntegral;
};

void styleNucleusHist(TH1D *hist, int nucleusIndex, bool isSim) {
  static const int kColors[2] = {kAzure + 1, kOrange + 7};
  const int safeIndex = (nucleusIndex == 0) ? 0 : 1;

  hist->SetStats(false);
  hist->SetLineColor(kColors[safeIndex]);
  hist->SetMarkerColor(kColors[safeIndex]);
  hist->SetLineWidth(2);
  hist->SetMarkerSize(0.80);
  hist->SetMarkerStyle(isSim ? 24 + safeIndex : 20 + safeIndex);
  hist->SetLineStyle(isSim ? 2 : 1);
}

std::vector<overlay_src::DetectorCombination> leadCombinations() {
  return {{"lead_fd", "Lead FD", "(leadTheta*180./TMath::Pi()<37.)"},
          {"lead_cd", "Lead CD", "(leadTheta*180./TMath::Pi()>45.)"}};
}

overlay_src::PlotVar pmissVariable() {
  return {"pMiss", "pMiss", "p_{miss} [GeV/c]", 20, 0.3, 1.1};
}

HistSummary projectHistogram(TTree *tree, const overlay_src::PlotVar &var,
                             const char *histName, const char *selectionCut,
                             const char *baseCut, const char *comboCut,
                             const char *pCMyTailCut, const char *weightExpression) {
  TH1D *hist = new TH1D(histName, "", var.bins, var.xmin, var.xmax);
  hist->Sumw2();

  const TString boolSel =
      overlay_src::buildSelection(selectionCut, baseCut, comboCut, pCMyTailCut);
  const TString weightedSel =
      overlay_src::applyWeightToSelection(weightExpression, boolSel);
  tree->Project(hist->GetName(), var.expr, weightedSel.Data());

  return {hist, hist->Integral(0, hist->GetNbinsX() + 1)};
}

TH1D *cloneForDisplay(const TH1D *source, const char *name, bool normalizeToUnity) {
  TH1D *clone = dynamic_cast<TH1D *>(source->Clone(name));
  clone->SetDirectory(nullptr);
  if (normalizeToUnity) {
    const double integral = clone->Integral();
    if (integral > 0.) clone->Scale(1. / integral);
  }
  return clone;
}

TH1D *buildUnityNormalizedRatio(const TH1D *numeratorSource,
                                const TH1D *denominatorSource,
                                const char *name) {
  TH1D *numerator = cloneForDisplay(numeratorSource, Form("%s_num", name), true);
  TH1D *denominator = cloneForDisplay(denominatorSource, Form("%s_den", name), true);
  TH1D *ratio = dynamic_cast<TH1D *>(numerator->Clone(name));
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICESM");
  ratio->Divide(numerator, denominator);
  delete numerator;
  delete denominator;
  return ratio;
}

double histogramMax(const std::vector<TH1D *> &hists) {
  double yMax = 0.;
  for (TH1D *hist : hists) {
    if (!hist) continue;
    yMax = std::max(yMax, overlay_src::finiteMax(hist->GetMaximum()));
  }
  return (yMax > 0.) ? yMax : 1.;
}

void histogramRange(const std::vector<TH1D *> &hists, double &yMin, double &yMax) {
  yMin = 1e9;
  yMax = -1e9;
  bool found = false;

  for (TH1D *hist : hists) {
    if (!hist) continue;
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
      const double value = hist->GetBinContent(bin);
      if (!std::isfinite(value)) continue;
      yMin = std::min(yMin, value);
      yMax = std::max(yMax, value);
      found = true;
    }
  }

  if (!found) {
    yMin = 0.5;
    yMax = 1.5;
    return;
  }

  if (yMin == yMax) {
    yMin -= 0.2;
    yMax += 0.2;
  } else {
    const double padding = 0.15 * (yMax - yMin);
    yMin -= padding;
    yMax += padding;
  }

  if (yMin > 1.) yMin = 1.;
}

void setPadStyle() {
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.10);
  gPad->SetRightMargin(0.04);
}

}  // namespace

void overlay_he4_c12_data_sim_lead_ratios(
    const char *heDataFileName = "~/data/RGM_DATA/he4_src_skim_100MeV.root",
    const char *c12DataFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *heSimFileName = "~/data/RGM_DATA/he4_sim_skim_allD.root",
    const char *c12SimFileName = "~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
  const char *outputPdfName = "pdf/he4_c12_pmiss_data_sim_lead_ratios.pdf",
    bool normalizeOverlayPage = true,
    const char *leadCut = "goodLead",
    const char *baseCut = "pMiss < 1.",
    const char *dataWeightsCsv = "(weight_ep),(weight_ep)",
    const char *simWeightsCsv = "(weight_ep),(weight_ep)",
  const char *pCMyTailCut = "") {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  std::vector<TString> dataWeights = overlay_src::splitCsvToTStrings(dataWeightsCsv);
  std::vector<TString> simWeights = overlay_src::splitCsvToTStrings(simWeightsCsv);
  while (dataWeights.size() < 2) dataWeights.push_back("(weight_ep)");
  while (simWeights.size() < 2) simWeights.push_back("(weight_ep)");

  const std::vector<SampleConfig> samples = {
      {heDataFileName, dataWeights[0]},
      {c12DataFileName, dataWeights[1]},
      {heSimFileName, simWeights[0]},
      {c12SimFileName, simWeights[1]}};

  std::vector<TFile *> files;
  std::vector<TTree *> trees;
  files.reserve(samples.size());
  trees.reserve(samples.size());

  for (const auto &sample : samples) {
    TFile *file = TFile::Open(sample.path.Data(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Could not open file " << sample.path << "\n";
      for (TFile *openFile : files) openFile->Close();
      return;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get(treeName));
    if (!tree) {
      std::cerr << "Could not find tree \"" << treeName << "\" in " << sample.path
                << "\n";
      file->Close();
      for (TFile *openFile : files) openFile->Close();
      return;
    }

    files.push_back(file);
    trees.push_back(tree);
  }

  const std::vector<overlay_src::DetectorCombination> combos = leadCombinations();
//   const overlay_src::PlotVar var = pmissVariable();
  const overlay_src::PlotVar var = {"pMiss", "pMiss", "p_{miss} [GeV/c]", 10, 0.25, 1.};

  TCanvas *canvas = new TCanvas("c_overlay_he4_c12_data_sim_lead_ratios",
                                "He4/C12 lead overlays and ratios", 1350, 950);
  canvas->SetLeftMargin(0.11);
  canvas->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  canvas->Print(pdfName + "[", "pdf");

    struct ComboProducts {
    TH1D *heDataOverlay;
    TH1D *cDataOverlay;
    TH1D *heSimOverlay;
    TH1D *cSimOverlay;
    TH1D *dataRatio;
    TH1D *simRatio;
    TH1D *doubleRatio;
    double heDataIntegral;
    double cDataIntegral;
    double heSimIntegral;
    double cSimIntegral;
    } comboProducts[2];

    for (std::size_t ic = 0; ic < combos.size(); ++ic) {
    const auto &combo = combos[ic];

    const HistSummary heDataRaw = projectHistogram(
      trees[0], var, Form("h_%s_%s_he_data_raw", var.name, combo.key), leadCut,
      baseCut, combo.cut, pCMyTailCut, samples[0].weight.Data());
    const HistSummary cDataRaw = projectHistogram(
      trees[1], var, Form("h_%s_%s_c_data_raw", var.name, combo.key), leadCut,
      baseCut, combo.cut, pCMyTailCut, samples[1].weight.Data());
    const HistSummary heSimRaw = projectHistogram(
      trees[2], var, Form("h_%s_%s_he_sim_raw", var.name, combo.key), leadCut,
      baseCut, combo.cut, pCMyTailCut, samples[2].weight.Data());
    const HistSummary cSimRaw = projectHistogram(
      trees[3], var, Form("h_%s_%s_c_sim_raw", var.name, combo.key), leadCut,
      baseCut, combo.cut, pCMyTailCut, samples[3].weight.Data());

    comboProducts[ic].heDataOverlay = cloneForDisplay(
      heDataRaw.raw, Form("h_%s_%s_he_data_overlay", var.name, combo.key),
      normalizeOverlayPage);
    comboProducts[ic].cDataOverlay = cloneForDisplay(
      cDataRaw.raw, Form("h_%s_%s_c_data_overlay", var.name, combo.key),
      normalizeOverlayPage);
    comboProducts[ic].heSimOverlay = cloneForDisplay(
      heSimRaw.raw, Form("h_%s_%s_he_sim_overlay", var.name, combo.key),
      normalizeOverlayPage);
    comboProducts[ic].cSimOverlay = cloneForDisplay(
      cSimRaw.raw, Form("h_%s_%s_c_sim_overlay", var.name, combo.key),
      normalizeOverlayPage);
    comboProducts[ic].dataRatio = buildUnityNormalizedRatio(
      heDataRaw.raw, cDataRaw.raw, Form("h_%s_%s_data_ratio", var.name, combo.key));
    comboProducts[ic].simRatio = buildUnityNormalizedRatio(
      heSimRaw.raw, cSimRaw.raw, Form("h_%s_%s_sim_ratio", var.name, combo.key));
    comboProducts[ic].doubleRatio = dynamic_cast<TH1D *>(
      comboProducts[ic].dataRatio->Clone(Form("h_%s_%s_double_ratio", var.name, combo.key)));
    comboProducts[ic].doubleRatio->SetDirectory(nullptr);
    comboProducts[ic].doubleRatio->Reset("ICESM");
    comboProducts[ic].doubleRatio->Divide(comboProducts[ic].dataRatio,
                                          comboProducts[ic].simRatio);

    comboProducts[ic].heDataIntegral = heDataRaw.rawIntegral;
    comboProducts[ic].cDataIntegral = cDataRaw.rawIntegral;
    comboProducts[ic].heSimIntegral = heSimRaw.rawIntegral;
    comboProducts[ic].cSimIntegral = cSimRaw.rawIntegral;

    styleNucleusHist(comboProducts[ic].heDataOverlay, 0, false);
    styleNucleusHist(comboProducts[ic].cDataOverlay, 1, false);
    styleNucleusHist(comboProducts[ic].heSimOverlay, 0, true);
    styleNucleusHist(comboProducts[ic].cSimOverlay, 1, true);

    styleNucleusHist(comboProducts[ic].dataRatio, 0, false);
    comboProducts[ic].dataRatio->SetLineColor(kBlack);
    comboProducts[ic].dataRatio->SetMarkerColor(kBlack);
    comboProducts[ic].dataRatio->SetMarkerStyle(20);
    comboProducts[ic].simRatio->SetLineColor(kRed + 1);
    comboProducts[ic].simRatio->SetMarkerColor(kRed + 1);
    comboProducts[ic].simRatio->SetMarkerStyle(24);
    comboProducts[ic].simRatio->SetLineStyle(2);
    comboProducts[ic].simRatio->SetLineWidth(2);
    comboProducts[ic].simRatio->SetStats(false);
    comboProducts[ic].doubleRatio->SetStats(false);
    comboProducts[ic].doubleRatio->SetLineWidth(2);
    comboProducts[ic].doubleRatio->SetLineColor(kBlue + 2);
    comboProducts[ic].doubleRatio->SetMarkerColor(kBlue + 2);
    comboProducts[ic].doubleRatio->SetMarkerStyle(21);

    delete heDataRaw.raw;
    delete cDataRaw.raw;
    delete heSimRaw.raw;
    delete cSimRaw.raw;
    }

    canvas->Clear();
    canvas->Divide(2, 2);

    for (std::size_t ic = 0; ic < combos.size(); ++ic) {
    const auto &combo = combos[ic];

    canvas->cd(static_cast<int>(ic) + 1);
    setPadStyle();
    const double yMaxData = histogramMax(
      {comboProducts[ic].heDataOverlay, comboProducts[ic].cDataOverlay});
    comboProducts[ic].heDataOverlay->SetMinimum(0.);
    comboProducts[ic].heDataOverlay->SetMaximum(1.25 * yMaxData);
    comboProducts[ic].heDataOverlay->GetXaxis()->SetTitle(var.xTitle);
    comboProducts[ic].heDataOverlay->GetYaxis()->SetTitle(
      normalizeOverlayPage ? "Normalized counts" : "Counts");
    comboProducts[ic].heDataOverlay->Draw("E");
    comboProducts[ic].cDataOverlay->Draw("E SAME");

    TLegend *legend = new TLegend(0.48, 0.67, 0.90, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.028);
    legend->AddEntry(comboProducts[ic].heDataOverlay,
             Form("He4 data (sumW=%.1f)", comboProducts[ic].heDataIntegral),
             "lep");
    legend->AddEntry(comboProducts[ic].cDataOverlay,
             Form("C12 data (sumW=%.1f)", comboProducts[ic].cDataIntegral),
             "lep");
    legend->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.038);
    label.DrawLatex(0.14, 0.92, Form("%s data", combo.label));

    canvas->cd(static_cast<int>(ic) + 3);
    setPadStyle();
    const double yMaxSim = histogramMax(
      {comboProducts[ic].heSimOverlay, comboProducts[ic].cSimOverlay});
    comboProducts[ic].heSimOverlay->SetMinimum(0.);
    comboProducts[ic].heSimOverlay->SetMaximum(1.25 * yMaxSim);
    comboProducts[ic].heSimOverlay->GetXaxis()->SetTitle(var.xTitle);
    comboProducts[ic].heSimOverlay->GetYaxis()->SetTitle(
      normalizeOverlayPage ? "Normalized counts" : "Counts");
    comboProducts[ic].heSimOverlay->Draw("E");
    comboProducts[ic].cSimOverlay->Draw("E SAME");

    TLegend *legendSim = new TLegend(0.48, 0.67, 0.90, 0.89);
    legendSim->SetBorderSize(0);
    legendSim->SetFillStyle(0);
    legendSim->SetTextSize(0.028);
    legendSim->AddEntry(comboProducts[ic].heSimOverlay,
              Form("He4 sim (sumW=%.1f)", comboProducts[ic].heSimIntegral),
              "lep");
    legendSim->AddEntry(comboProducts[ic].cSimOverlay,
              Form("C12 sim (sumW=%.1f)", comboProducts[ic].cSimIntegral),
              "lep");
    legendSim->Draw();

    TLatex simLabel;
    simLabel.SetNDC();
    simLabel.SetTextSize(0.038);
    simLabel.DrawLatex(0.14, 0.92, Form("%s sim", combo.label));
    }

    canvas->cd(0);
    TLatex pageTitle;
    pageTitle.SetNDC();
    pageTitle.SetTextFont(62);
    pageTitle.SetTextSize(0.030);
    pageTitle.DrawLatex(0.18, 0.965, "pMiss: He4/C12 overlays by lead detector");
    canvas->Print(pdfName, "pdf");

    canvas->Clear();
    canvas->Divide(2, 2);

    for (std::size_t ic = 0; ic < combos.size(); ++ic) {
    const auto &combo = combos[ic];

    canvas->cd(static_cast<int>(ic) + 1);
    setPadStyle();
    double yMin = 0.;
    double yMax = 0.;
    histogramRange({comboProducts[ic].dataRatio}, yMin, yMax);
    comboProducts[ic].dataRatio->SetMinimum(yMin);
    comboProducts[ic].dataRatio->SetMaximum(yMax);
    comboProducts[ic].dataRatio->GetXaxis()->SetTitle(var.xTitle);
    comboProducts[ic].dataRatio->GetYaxis()->SetTitle("He4/C12 (each norm. to 1)");
    comboProducts[ic].dataRatio->Draw("E");

    TLine *unity = new TLine(var.xmin, 1.0, var.xmax, 1.0);
    unity->SetLineStyle(7);
    unity->SetLineColor(kGray + 2);
    unity->Draw();
    comboProducts[ic].dataRatio->Draw("E SAME");

    TLegend *legend = new TLegend(0.52, 0.74, 0.90, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.028);
    legend->AddEntry(comboProducts[ic].dataRatio, "Data He4/C12", "lep");
    legend->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.038);
    label.DrawLatex(0.14, 0.92, Form("%s data ratio", combo.label));

    canvas->cd(static_cast<int>(ic) + 3);
    setPadStyle();
    histogramRange({comboProducts[ic].simRatio}, yMin, yMax);
    comboProducts[ic].simRatio->SetMinimum(yMin);
    comboProducts[ic].simRatio->SetMaximum(yMax);
    comboProducts[ic].simRatio->GetXaxis()->SetTitle(var.xTitle);
    comboProducts[ic].simRatio->GetYaxis()->SetTitle("He4/C12 (each norm. to 1)");
    comboProducts[ic].simRatio->Draw("E");

    TLine *unitySim = new TLine(var.xmin, 1.0, var.xmax, 1.0);
    unitySim->SetLineStyle(7);
    unitySim->SetLineColor(kGray + 2);
    unitySim->Draw();
    comboProducts[ic].simRatio->Draw("E SAME");

    TLegend *legendSim = new TLegend(0.52, 0.74, 0.90, 0.89);
    legendSim->SetBorderSize(0);
    legendSim->SetFillStyle(0);
    legendSim->SetTextSize(0.028);
    legendSim->AddEntry(comboProducts[ic].simRatio, "Sim He4/C12", "lep");
    legendSim->Draw();

    TLatex simLabel;
    simLabel.SetNDC();
    simLabel.SetTextSize(0.038);
    simLabel.DrawLatex(0.14, 0.92, Form("%s sim ratio", combo.label));
    }

    canvas->cd(0);
    pageTitle.DrawLatex(0.14, 0.965, "pMiss: normalized He4/C12 ratios by lead detector");
    canvas->Print(pdfName, "pdf");

    for (std::size_t ic = 0; ic < combos.size(); ++ic) {
    const auto &combo = combos[ic];
    canvas->Clear();
    canvas->cd(1);
    setPadStyle();

    double yMin = 0.;
    double yMax = 0.;
    histogramRange({comboProducts[ic].dataRatio, comboProducts[ic].simRatio}, yMin, yMax);
    comboProducts[ic].dataRatio->SetMinimum(yMin);
    comboProducts[ic].dataRatio->SetMaximum(yMax);
    comboProducts[ic].dataRatio->GetXaxis()->SetTitle(var.xTitle);
    comboProducts[ic].dataRatio->GetYaxis()->SetTitle("He4/C12 (each norm. to 1)");
    comboProducts[ic].dataRatio->Draw("E");

    TLine *unity = new TLine(var.xmin, 1.0, var.xmax, 1.0);
    unity->SetLineStyle(7);
    unity->SetLineColor(kGray + 2);
    unity->Draw();
    comboProducts[ic].dataRatio->Draw("E SAME");
    comboProducts[ic].simRatio->Draw("E SAME");

    TLegend *legend = new TLegend(0.58, 0.72, 0.90, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.028);
    legend->AddEntry(comboProducts[ic].dataRatio, "Data He4/C12", "lep");
    legend->AddEntry(comboProducts[ic].simRatio, "Sim He4/C12", "lep");
    legend->Draw();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(62);
    title.SetTextSize(0.032);
    title.DrawLatex(0.18, 0.95,
            Form("pMiss: data/sim overlay of He4/C12 ratio (%s)", combo.label));

    canvas->Print(pdfName, "pdf");
    }

      canvas->Clear();
      canvas->cd(1);
      setPadStyle();

      double yMin = 0.;
      double yMax = 0.;
      histogramRange({comboProducts[0].doubleRatio, comboProducts[1].doubleRatio}, yMin, yMax);
      comboProducts[0].doubleRatio->SetMinimum(yMin);
      comboProducts[0].doubleRatio->SetMaximum(yMax);
      comboProducts[0].doubleRatio->GetXaxis()->SetTitle(var.xTitle);
      comboProducts[0].doubleRatio->GetYaxis()->SetTitle("(Data He4/C12) / (Sim He4/C12)");
      comboProducts[0].doubleRatio->GetYaxis()->SetRangeUser(.2,1.8);

      comboProducts[0].doubleRatio->SetLineColor(kBlue + 2);
      comboProducts[0].doubleRatio->SetMarkerColor(kBlue + 2);
      comboProducts[0].doubleRatio->SetMarkerStyle(21);
      comboProducts[1].doubleRatio->SetLineColor(kGreen + 2);
      comboProducts[1].doubleRatio->SetMarkerColor(kGreen + 2);
      comboProducts[1].doubleRatio->SetMarkerStyle(22);

      comboProducts[0].doubleRatio->Draw("E");
      TLine *unity = new TLine(var.xmin, 1.0, var.xmax, 1.0);
      unity->SetLineStyle(7);
      unity->SetLineColor(kGray + 2);
      unity->Draw();
      comboProducts[0].doubleRatio->Draw("E SAME");
      comboProducts[1].doubleRatio->Draw("E SAME");

      TLegend *legend = new TLegend(0.52, 0.73, 0.90, 0.89);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.030);
      legend->AddEntry(comboProducts[0].doubleRatio, "Double ratio (Lead FD)", "lep");
      legend->AddEntry(comboProducts[1].doubleRatio, "Double ratio (Lead CD)", "lep");
      legend->Draw();

      TLatex doubleTitle;
      doubleTitle.SetNDC();
      doubleTitle.SetTextFont(62);
      doubleTitle.SetTextSize(0.032);
      doubleTitle.DrawLatex(0.18, 0.95,
                  "pMiss: (Data He4/C12) / (Sim He4/C12), FD and CD overlay");
      canvas->Print(pdfName, "pdf");

    for (std::size_t ic = 0; ic < combos.size(); ++ic) {
    delete comboProducts[ic].heDataOverlay;
    delete comboProducts[ic].cDataOverlay;
    delete comboProducts[ic].heSimOverlay;
    delete comboProducts[ic].cSimOverlay;
    delete comboProducts[ic].dataRatio;
    delete comboProducts[ic].simRatio;
    delete comboProducts[ic].doubleRatio;
  }

  canvas->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";

  for (TFile *file : files) file->Close();
}