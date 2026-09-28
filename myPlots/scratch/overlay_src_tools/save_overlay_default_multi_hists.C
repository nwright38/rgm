#include "OverlaySrcCommon.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

void save_overlay_default_multi_hists(
    const char *fileNamesCsv =
    "~/data/RGM_DATA/c12_src_skim.root,~/data/RGM_DATA/he4_src_skim_100MeV.root",
    const char *treeName = "srcTree",
    const char *outputRootName =
      "scratch/overlay_src_tools/c12_he4_src_skim_hists.root",
  const char *runTag = "",
  bool appendToFile = false,
    bool normalizeToUnity = true,
    const char *eppCut = "pCM > 0",
    const char *baseCut = "Q2 < 5. && leadP > 1.",
    const char *weightsCsv = "(weight_epp),(weight_epp)*(weight_epp < 150.)",
    const char *labelsCsv = "C12 Data,He4 Data",
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

  const char *outMode = appendToFile ? "UPDATE" : "RECREATE";
  TFile *outFile = TFile::Open(outputRootName, outMode);
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Could not create output file " << outputRootName << "\n";
    for (TFile *f : files) f->Close();
    return;
  }

  const TString runTagStr =
      (runTag && runTag[0] != '\0') ? TString(runTag) : TString("");
  TDirectory *rootDir = outFile;
  if (runTagStr.Length() > 0) {
    TDirectory *existingRunDir = outFile->GetDirectory(runTagStr.Data());
    rootDir = existingRunDir ? existingRunDir : outFile->mkdir(runTagStr.Data());
  }

  TDirectory *metaDir = rootDir->GetDirectory("meta");
  if (!metaDir) metaDir = rootDir->mkdir("meta");
  metaDir->cd();
  for (std::size_t i = 0; i < fileNames.size(); ++i) {
    TH1D *meta = new TH1D(Form("sample_%zu", i), labels[i].Data(), 1, 0.0, 1.0);
    meta->SetBinContent(1, static_cast<double>(i));
    meta->GetXaxis()->SetTitle(fileNames[i]);
    meta->Write();
  }

  const std::vector<overlay_src::DetectorCombination> eppDrawnCombos =
      overlay_src::detectorCombinations(includeFdFd);
  const std::vector<overlay_src::DetectorCombination> epLeadCombos = {
      {"lead_fd", "Lead FD", "(leadTheta*180./TMath::Pi()<37.)"},
      {"lead_cd", "Lead CD", "(leadTheta*180./TMath::Pi()>45.)"}};
  const std::vector<overlay_src::DetectorCombination> leadOnlyCombos = {
      {"lead_fd_anyrec", "Lead FD (all recoil)", "(leadTheta*180./TMath::Pi()<37.)"},
      {"lead_cd_anyrec", "Lead CD (all recoil)", "(leadTheta*180./TMath::Pi()>45.)"}};
  const std::vector<overlay_src::DetectorCombination> inclusiveCombo = {
      {"inclusive", "Inclusive", "1"}};

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

  const std::vector<overlay_src::PlotVar> vars = overlay_src::defaultVariables();

  const auto writeComboSet =
      [&](const overlay_src::PlotVar &v,
          const std::vector<overlay_src::DetectorCombination> &comboSet,
          const char *selectionTag,
          const char *groupTag,
          const char *selCut,
          const char *selBaseCut,
          const std::vector<TString> &selWeights) {
        TDirectory *dirVar = rootDir->GetDirectory(v.name);
        if (!dirVar) dirVar = rootDir->mkdir(v.name);

        TDirectory *dirSel = dirVar->GetDirectory(selectionTag);
        if (!dirSel) dirSel = dirVar->mkdir(selectionTag);

        TDirectory *dirGroup = dirSel->GetDirectory(groupTag);
        if (!dirGroup) dirGroup = dirSel->mkdir(groupTag);

        for (const auto &combo : comboSet) {
          TDirectory *dirTopo = dirGroup->GetDirectory(combo.key);
          if (!dirTopo) dirTopo = dirGroup->mkdir(combo.key);
          dirTopo->cd();

          const TString boolSel =
              overlay_src::buildSelection(selCut, selBaseCut, combo.cut, pCMyTailCut);

          for (std::size_t is = 0; is < trees.size(); ++is) {
            const TString weightedSel =
                overlay_src::applyWeightToSelection(selWeights[is].Data(), boolSel);

            TH1D *h = new TH1D(
                Form("h_%s_%s_%s_s%zu", v.name, selectionTag, combo.key, is),
                Form("%s ;%s;Counts", combo.label, v.xTitle),
                v.bins, v.xmin, v.xmax);
            h->Sumw2();

            if (maxEvents > 0) {
              trees[is]->Project(h->GetName(), v.expr, weightedSel.Data(), "", maxEvents,
                                 firstEvent);
            } else {
              trees[is]->Project(h->GetName(), v.expr, weightedSel.Data());
            }

            if (normalizeToUnity) {
              const double intInRange = h->Integral();
              if (intInRange > 0.) h->Scale(1. / intInRange);
            }

            h->SetDirectory(dirTopo);
            h->Write();
            delete h;
          }
        }
      };

  for (const auto &v : vars) {
    if (isEpCompatible(v.name)) {
      writeComboSet(v, epLeadCombos, "ep", "drawn", epCut, epBaseCut, epWeights);
      writeComboSet(v, inclusiveCombo, "ep", "inclusive", epCut, epBaseCut, epWeights);
    }

    writeComboSet(v, eppDrawnCombos, "epp", "drawn", eppCut, baseCut, weights);
    writeComboSet(v, leadOnlyCombos, "epp", "lead_only", eppCut, baseCut, weights);
    writeComboSet(v, inclusiveCombo, "epp", "inclusive", eppCut, baseCut, weights);
  }

  outFile->Write();
  outFile->Close();

  if (runTagStr.Length() > 0) {
    std::cout << "Wrote " << outputRootName << " (runTag=" << runTagStr << ")\n";
  } else {
    std::cout << "Wrote " << outputRootName << "\n";
  }

  for (TFile *f : files) f->Close();
}
