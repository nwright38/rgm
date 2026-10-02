#include "OverlaySrcCommon.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
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

namespace {

struct FileCutPair {
  std::size_t fileIdx;
  std::size_t cutIdx;
};

bool parseIndex(const std::string &s, std::size_t *outIdx) {
  if (!outIdx) return false;
  if (s.empty()) return false;

  char *endPtr = nullptr;
  const long v = std::strtol(s.c_str(), &endPtr, 10);
  if (!endPtr || *endPtr != '\0' || v < 0) return false;
  if (static_cast<unsigned long>(v) >
      static_cast<unsigned long>(std::numeric_limits<std::size_t>::max())) {
    return false;
  }

  *outIdx = static_cast<std::size_t>(v);
  return true;
}

std::vector<FileCutPair> buildPairs(std::size_t nFiles, std::size_t nCuts,
                                    const char *fileCutPairsCsv) {
  std::vector<FileCutPair> pairs;

  if (fileCutPairsCsv && fileCutPairsCsv[0] != '\0') {
    const std::vector<TString> rawPairs =
        overlay_src::splitCsvToTStrings(fileCutPairsCsv);
    pairs.reserve(rawPairs.size());

    for (const auto &raw : rawPairs) {
      const std::string token = raw.Data();
      const std::size_t colon = token.find(':');
      if (colon == std::string::npos) {
        std::cerr << "Ignoring invalid pair token '" << token
                  << "' (expected fileIndex:cutIndex).\n";
        continue;
      }

      const TString leftTrim = overlay_src::trimCopy(token.substr(0, colon));
      const TString rightTrim = overlay_src::trimCopy(token.substr(colon + 1));

      std::size_t fileIdx = 0;
      std::size_t cutIdx = 0;
      const bool fileOk = parseIndex(leftTrim.Data(), &fileIdx);
      const bool cutOk = parseIndex(rightTrim.Data(), &cutIdx);
      if (!fileOk || !cutOk) {
        std::cerr << "Ignoring invalid pair token '" << token
                  << "' (non-integer index).\n";
        continue;
      }

      if (fileIdx >= nFiles || cutIdx >= nCuts) {
        std::cerr << "Ignoring out-of-range pair token '" << token << "' for "
                  << nFiles << " file(s) and " << nCuts << " cut(s).\n";
        continue;
      }

      pairs.push_back({fileIdx, cutIdx});
    }
  }

  if (!pairs.empty()) return pairs;

  pairs.reserve(nFiles * nCuts);
  for (std::size_t iFile = 0; iFile < nFiles; ++iFile) {
    for (std::size_t iCut = 0; iCut < nCuts; ++iCut) {
      pairs.push_back({iFile, iCut});
    }
  }
  return pairs;
}

TString combineCuts(const char *a, const char *b) {
  std::vector<TString> terms;
  if (a && a[0] != '\0') terms.emplace_back(a);
  if (b && b[0] != '\0') terms.emplace_back(b);

  if (terms.empty()) return "";

  TString out = "(";
  for (std::size_t i = 0; i < terms.size(); ++i) {
    if (i > 0) out += " && ";
    out += "(";
    out += terms[i];
    out += ")";
  }
  out += ")";
  return out;
}

void styleOverlaySeries(TH1D *hist, int idx) {
  static const int kColors[] = {
      kBlack, kRed + 1, kBlue + 1,  kGreen + 2,  kOrange + 7, kMagenta + 1,
      kCyan + 2, kGray + 2, kViolet + 1, kTeal + 3,  kAzure + 2, kPink + 7};
  static const int kMarkers[] = {20, 21, 22, 23, 24, 25, 26, 32, 33, 34, 29, 30};

  const int nStyles = static_cast<int>(sizeof(kColors) / sizeof(kColors[0]));
  const int safe = (idx >= 0) ? (idx % nStyles) : 0;

  hist->SetStats(false);
  hist->SetLineColor(kColors[safe]);
  hist->SetMarkerColor(kColors[safe]);
  hist->SetMarkerStyle(kMarkers[safe]);
  hist->SetMarkerSize(0.80);
  hist->SetLineWidth(2);
}

}  // namespace

void overlay_singlefile_multicut(
    const char *fileNamesCsv =
        "~/data/RGM_DATA/c12_src_skim.root,~/data/RGM_DATA/c12_sim_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "pdf/c12_overlay_multicut_test.pdf",
    const char *fileLabelsCsv = "C12 Data,C12 Sim",
    const char *eppWeightsCsv = "(weight_epp),(weight_epp)*(weight_epp < 150.)",
    const char *epWeightsCsv = "(weight_ep),(weight_ep)*(weight_ep < 200.)",
    const char *cutsCsv = "mMiss > 0.65, leadP/qP < 0.95",
    const char *cutLabelsCsv = "mMiss > 0.65, p/q < 0.95",
    const char *fileCutPairsCsv = "",
    const char *epCut = "Q2 > 1.5 && leadP > 1. && xB > 1.2 && mMiss < 1.1 && kMiss > 0.3",
    const char *eppExtraCut = "pCM > 0 && recP > .35",
    const char *baseCut =
        "Q2 > 1.5 && leadP > 1. && xB > 1.2 && mMiss < 1.1 && kMiss > 0.3",
    const char *pCMyTailCut = "",
    bool includeFdFd = false,
    bool normalizeToUnity = true,
    bool printSelectionDiagnostics = true,
    Long64_t maxEvents = -1,
    Long64_t firstEvent = 0) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kTRUE);

  const std::vector<TString> fileNames = overlay_src::splitCsvToTStrings(fileNamesCsv);
  if (fileNames.empty()) {
    std::cerr << "fileNamesCsv is empty. Provide at least one file.\n";
    return;
  }

  std::vector<TString> fileLabels = overlay_src::splitCsvToTStrings(fileLabelsCsv);
  std::vector<TString> eppWeights = overlay_src::splitCsvToTStrings(eppWeightsCsv);
  std::vector<TString> epWeights = overlay_src::splitCsvToTStrings(epWeightsCsv);
  std::vector<TString> cuts = overlay_src::splitCsvToTStrings(cutsCsv);
  std::vector<TString> cutLabels = overlay_src::splitCsvToTStrings(cutLabelsCsv);

  if (cuts.empty()) {
    std::cerr << "cutsCsv is empty. Provide at least one cut expression.\n";
    return;
  }

  while (fileLabels.size() < fileNames.size()) {
    fileLabels.push_back(Form("File %zu", fileLabels.size() + 1));
  }
  while (eppWeights.size() < fileNames.size()) {
    eppWeights.push_back("(weight_epp)");
  }
  while (epWeights.size() < fileNames.size()) {
    epWeights.push_back("(weight_ep)");
  }
  while (cutLabels.size() < cuts.size()) {
    cutLabels.push_back(Form("Cut %zu", cutLabels.size() + 1));
  }

  const std::vector<FileCutPair> pairs =
      buildPairs(fileNames.size(), cuts.size(), fileCutPairsCsv);
  if (pairs.empty()) {
    std::cerr << "No valid file/cut pairs to draw.\n";
    return;
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
      "xB",          "Q2",           "omega",       "eP",      "eTheta",
      "pLead",       "pMiss",        "kMiss",       "mMiss",   "EMiss",
      "E0miss",      "E1miss",       "leadTheta",   "pMissTheta",
      "theta_PmQ",   "theta_PmPlead", "theta_PleadQ", "alpha_1",
      "alpha_q",     "alpha_pLead",  "p1_plus",     "p1_lc_z",
      "p1_perp_mag", "lead_q_ratio", "qP",          "qTheta"};

  const auto isEpCompatible = [&epCompatibleVars](const char *name) {
    return std::find(epCompatibleVars.begin(), epCompatibleVars.end(),
                     std::string(name)) != epCompatibleVars.end();
  };

  const TString eppCutCombined = combineCuts(epCut, eppExtraCut);

  TCanvas *c =
      new TCanvas("c_overlay_singlefile_multicut", "Multi-file multi-cut overlay", 1300,
                  950);
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.11);

  TString pdfName(outputPdfName);
  c->Print(pdfName + "[", "pdf");

  const auto drawPage = [&](const overlay_src::PlotVar &v,
                            const std::vector<overlay_src::DetectorCombination> &pageCombos,
                            const char *channelCut,
                            const std::vector<TString> &pageWeights,
                            const char *pageTag,
                            bool printDiagnosticsForThisPage) {
    c->Clear();
    c->Divide(2, 2);

    for (std::size_t ic = 0; ic < pageCombos.size() && ic < 4; ++ic) {
      const auto &combo = pageCombos[ic];
      c->cd(static_cast<int>(ic) + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.10);

      std::vector<TH1D *> hists;
      std::vector<double> rawIntegrals;
      hists.reserve(pairs.size());
      rawIntegrals.reserve(pairs.size());

      for (std::size_t iPair = 0; iPair < pairs.size(); ++iPair) {
        const std::size_t iFile = pairs[iPair].fileIdx;
        const std::size_t iCut = pairs[iPair].cutIdx;

        TH1D *h = new TH1D(
            Form("h_%s_%s_%s_f%zu_c%zu", v.name, pageTag, combo.key, iFile, iCut),
            "", v.bins, v.xmin, v.xmax);
        h->Sumw2();
        styleOverlaySeries(h, static_cast<int>(iPair));

        const TString fullBoolSel =
            overlay_src::buildSelection(channelCut, baseCut, combo.cut, pCMyTailCut,
                                        nullptr, cuts[iCut].Data());
        const TString weightedSel =
            overlay_src::applyWeightToSelection(pageWeights[iFile].Data(), fullBoolSel);

        if (maxEvents > 0) {
          trees[iFile]->Project(h->GetName(), v.expr, weightedSel.Data(), "", maxEvents,
                                firstEvent);
        } else {
          trees[iFile]->Project(h->GetName(), v.expr, weightedSel.Data());
        }

        const double rawInt = h->Integral(0, h->GetNbinsX() + 1);
        rawIntegrals.push_back(rawInt);

        if (printDiagnosticsForThisPage && ic == 0) {
          const Long64_t unweightedEntries = trees[iFile]->GetEntries(fullBoolSel.Data());
          std::cout << "Diag[" << pageTag << "] file=" << fileLabels[iFile]
                    << " cut=" << cutLabels[iCut] << " entries=" << unweightedEntries
                    << " sumW=" << rawInt << " weightExpr=" << pageWeights[iFile]
                    << "\n";
        }

        if (normalizeToUnity) {
          const double inRange = h->Integral();
          if (inRange > 0.) h->Scale(1. / inRange);
        }

        hists.push_back(h);
      }

      if (hists.empty()) continue;

      double panelYMax = 0.;
      for (TH1D *h : hists) {
        panelYMax = std::max(panelYMax, overlay_src::finiteMax(h->GetMaximum()));
      }
      if (panelYMax <= 0.) panelYMax = 1.;

      hists.front()->SetMinimum(0.);
      hists.front()->SetMaximum(1.25 * panelYMax);
      hists.front()->GetXaxis()->SetTitle(v.xTitle);
      hists.front()->GetYaxis()->SetTitle(normalizeToUnity ? "Normalized counts"
                                                           : "Counts");
      hists.front()->GetXaxis()->SetTitleSize(0.050);
      hists.front()->GetYaxis()->SetTitleSize(0.050);
      hists.front()->GetXaxis()->SetLabelSize(0.043);
      hists.front()->GetYaxis()->SetLabelSize(0.043);
      hists.front()->GetYaxis()->SetTitleOffset(1.20);

      hists.front()->Draw("E");
      for (std::size_t i = 1; i < hists.size(); ++i) {
        hists[i]->Draw("E SAME");
      }

      TLegend *leg = new TLegend(0.36, 0.52, 0.90, 0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.024);
      for (std::size_t iPair = 0; iPair < hists.size(); ++iPair) {
        const std::size_t iFile = pairs[iPair].fileIdx;
        const std::size_t iCut = pairs[iPair].cutIdx;
        leg->AddEntry(
            hists[iPair],
            Form("%s | %s (sumW=%.1f)", fileLabels[iFile].Data(),
                 cutLabels[iCut].Data(), rawIntegrals[iPair]),
            "lep");
      }
      leg->Draw();

      TLatex panelLabel;
      panelLabel.SetNDC();
      panelLabel.SetTextSize(0.040);
      panelLabel.SetTextFont(42);
      panelLabel.DrawLatex(0.14, 0.92, combo.label);
    }

    for (std::size_t ip = pageCombos.size(); ip < 4; ++ip) {
      c->cd(static_cast<int>(ip) + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.10);
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
    title.DrawLatex(0.10, 0.985, Form("%s: multi-cut overlay (%s)", v.name, pageTag));

    c->Print(pdfName, "pdf");
  };

  for (std::size_t iv = 0; iv < vars.size(); ++iv) {
    const auto &v = vars[iv];
    const bool printDiagThisVar = printSelectionDiagnostics && iv == 0;

    if (isEpCompatible(v.name)) {
      drawPage(v, epCombos, epCut, epWeights, "e'p", printDiagThisVar);
      drawPage(v, combos, eppCutCombined.Data(), eppWeights, "e'pp",
               printDiagThisVar);
    } else {
      drawPage(v, combos, eppCutCombined.Data(), eppWeights, "e'pp",
               printDiagThisVar);
    }
  }

  c->Print(pdfName + "]", "pdf");
  std::cout << "Wrote " << outputPdfName << "\n";

  for (TFile *f : files) f->Close();
}
