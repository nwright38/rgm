// STEP 1 OF 2: create the ROOT file used by plotCMByMiss.C.
//
// Default workflow:
//   root -l -b -q 'myPlots/scratch/fillCMByMiss.C()'
//     Input:  ~/data/RGM_DATA/c12_src_skim.root
//     Output: cm_by_miss.root
//
// Then make the PDF with:
//   root -l -b -q 'myPlots/scratch/plotCMByMiss.C()'
//
// Custom input and output ROOT file:
//   root -l -b -q 'myPlots/scratch/fillCMByMiss.C("input.root","my_hists.root")'
//   root -l -b -q 'myPlots/scratch/plotCMByMiss.C("my_hists.root","my_plots.pdf")'
//
// Change the CONFIGURATION block below to select a detector topology.  All
// enabled selections are ANDed with baseCut, so pCM > 0 is always retained.

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TTree.h>
#include <TTreeFormula.h>

namespace {

// --------------------------- CONFIGURATION ---------------------------
// Leave all four false for no detector-topology restriction.
// Set at most one lead option and at most one recoil option to true.
const bool requireLeadFD = false;    // leadTheta < 37 degrees
const bool requireLeadCD = false;    // leadTheta > 45 degrees
const bool requireRecoilFD = false;  // recTheta  < 37 degrees
const bool requireRecoilCD = false;  // recTheta  > 45 degrees

const char *baseCut = "pCM > 0 && weight_epp < 300 && pMiss < 1. && recP < 1.";
const char *weightExpression = "weight_epp";

std::vector<double> pmissBins() {
  return {0.2, 0.4, 0.50, 0.60, 0.70, 0.85, 1.};
}

std::vector<double> kmissBins() {
  return {0.30, 0.45, 0.60, 0.75, 0.90, 1.};
}
// ---------------------------------------------------------------------

std::string topologyCut() {
  std::vector<std::string> cuts;
  if (requireLeadFD)
    cuts.push_back("leadTheta*180./TMath::Pi() < 37.");
  if (requireLeadCD)
    cuts.push_back("leadTheta*180./TMath::Pi() > 45.");
  if (requireRecoilFD)
    cuts.push_back("recTheta*180./TMath::Pi() < 37.");
  if (requireRecoilCD)
    cuts.push_back("recTheta*180./TMath::Pi() > 45.");
  if (cuts.empty()) return "1";

  std::string result;
  for (const std::string &cut : cuts) {
    if (!result.empty()) result += " && ";
    result += "(" + cut + ")";
  }
  return result;
}

struct HistPair {
  TH1D *x;
  TH1D *y;
};

bool hasBranch(TTree *tree, const char *name) {
  return tree->GetBranch(name) || tree->GetLeaf(name);
}

int findBin(double value, const std::vector<double> &edges) {
  for (int i = 0; i + 1 < static_cast<int>(edges.size()); ++i)
    if (value >= edges[i] &&
        (value < edges[i + 1] ||
         (i + 2 == static_cast<int>(edges.size()) && value == edges[i + 1])))
      return i;
  return -1;
}

std::vector<HistPair> makeHistograms(const char *variable,
                                     const std::vector<double> &edges) {
  std::vector<HistPair> result;
  for (int i = 0; i + 1 < static_cast<int>(edges.size()); ++i) {
    const TString label =
        Form("%.3g #leq %s < %.3g GeV/c", edges[i], variable, edges[i + 1]);
    TH1D *hx = new TH1D(Form("h_pCMx_%s_bin%02d", variable, i),
                        Form("%s;p_{CM,x} [GeV/c];Weighted counts", label.Data()),
                        40, -0.75, 0.75);
    TH1D *hy = new TH1D(Form("h_pCMy_%s_bin%02d", variable, i),
                        Form("%s;p_{CM,y} [GeV/c];Weighted counts", label.Data()),
                        40, -0.75, 0.75);
    hx->Sumw2();
    hy->Sumw2();
    result.push_back({hx, hy});
  }
  return result;
}

}  // namespace

void fillCMByMiss(
    const char *inputFileName = "~/data/RGM_DATA/c12_sim_skim.root",
    const char *outputFileName = "cm_by_miss.root",
    const char *treeName = "srcTree") {

  if ((requireLeadFD && requireLeadCD) ||
      (requireRecoilFD && requireRecoilCD)) {
    std::cerr << "[fillCMByMiss] Choose only one detector region (FD or CD) "
                 "for each proton.\n";
    return;
  }
  const std::string detectorCut = topologyCut();

  TFile *input = TFile::Open(inputFileName, "READ");
  TTree *tree = input && !input->IsZombie()
                    ? dynamic_cast<TTree *>(input->Get(treeName))
                    : nullptr;
  if (!tree) {
    std::cerr << "[fillCMByMiss] Cannot read tree " << treeName << " from "
              << inputFileName << '\n';
    if (input) input->Close();
    return;
  }
  for (const char *branch : {"pCMx", "pCMy", "pMiss", "kMiss"}) {
    if (!hasBranch(tree, branch)) {
      std::cerr << "[fillCMByMiss] Missing required branch " << branch << '\n';
      input->Close();
      return;
    }
  }

  TTreeFormula topology("topology_formula", detectorCut.c_str(), tree);
  TTreeFormula common("base_formula", baseCut, tree);
  TTreeFormula weight("weight_formula", weightExpression, tree);
  TTreeFormula pCMxValue("pCMx_formula", "pCMx", tree);
  TTreeFormula pCMyValue("pCMy_formula", "pCMy", tree);
  TTreeFormula pMissValue("pMiss_formula", "pMiss", tree);
  TTreeFormula kMissValue("kMiss_formula", "kMiss", tree);
  if (topology.GetNdim() <= 0 || common.GetNdim() <= 0 ||
      weight.GetNdim() <= 0 || pCMxValue.GetNdim() <= 0 ||
      pCMyValue.GetNdim() <= 0 || pMissValue.GetNdim() <= 0 ||
      kMissValue.GetNdim() <= 0) {
    std::cerr << "[fillCMByMiss] Invalid cut or weight expression.\n";
    input->Close();
    return;
  }

  const std::vector<double> pEdges = pmissBins();
  const std::vector<double> kEdges = kmissBins();
  std::vector<HistPair> pHists = makeHistograms("pmiss", pEdges);
  std::vector<HistPair> kHists = makeHistograms("kmiss", kEdges);

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    const Long64_t local = tree->LoadTree(entry);
    if (local < 0) break;
    tree->GetEntry(entry);
    topology.GetNdata();
    common.GetNdata();
    weight.GetNdata();
    if (topology.EvalInstance() == 0. || common.EvalInstance() == 0.) continue;
    const double w = weight.EvalInstance();
    if (!std::isfinite(w) || w == 0.) continue;

    const double pCMx = pCMxValue.EvalInstance();
    const double pCMy = pCMyValue.EvalInstance();
    const double pMiss = pMissValue.EvalInstance();
    const double kMiss = kMissValue.EvalInstance();
    const int pBin = findBin(pMiss, pEdges);
    const int kBin = findBin(kMiss, kEdges);
    if (pBin >= 0) {
      pHists[pBin].x->Fill(pCMx, w);
      pHists[pBin].y->Fill(pCMy, w);
    }
    if (kBin >= 0) {
      kHists[kBin].x->Fill(pCMx, w);
      kHists[kBin].y->Fill(pCMy, w);
    }
  }

  TFile *output = TFile::Open(outputFileName, "RECREATE");
  if (!output || output->IsZombie()) {
    std::cerr << "[fillCMByMiss] Cannot create " << outputFileName << '\n';
    input->Close();
    return;
  }
  TNamed("topology_cut", detectorCut.c_str()).Write();
  TNamed("base_cut", baseCut).Write();
  TNamed("weight_expression", weightExpression).Write();

  char variable[16] = "";
  int bin = 0;
  double low = 0., high = 0.;
  TTree binInfo("bin_info", "Binning used for the pCM component histograms");
  binInfo.Branch("variable", variable, "variable/C");
  binInfo.Branch("bin", &bin, "bin/I");
  binInfo.Branch("low", &low, "low/D");
  binInfo.Branch("high", &high, "high/D");
  for (const auto &spec : {std::make_pair("pmiss", &pEdges),
                           std::make_pair("kmiss", &kEdges)}) {
    snprintf(variable, sizeof(variable), "%s", spec.first);
    for (bin = 0; bin + 1 < static_cast<int>(spec.second->size()); ++bin) {
      low = (*spec.second)[bin];
      high = (*spec.second)[bin + 1];
      binInfo.Fill();
    }
  }
  binInfo.Write();
  for (const HistPair &pair : pHists) {
    pair.x->Write();
    pair.y->Write();
  }
  for (const HistPair &pair : kHists) {
    pair.x->Write();
    pair.y->Write();
  }
  output->Close();
  input->Close();
  std::cout << "[fillCMByMiss] Wrote " << outputFileName << '\n';
}
