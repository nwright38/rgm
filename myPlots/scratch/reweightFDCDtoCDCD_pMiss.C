// reweightFDCDtoCDCD_pMiss.C
//
// Source -> target reweighting in pMiss only:
//   source: lead FD + rec CD
//   target: lead CD + rec CD
//
// Event-by-event source weight:
//   w = h_cd_cd->GetBinContent(h_cd_cd->FindBin(pMiss)) /
//       h_fd_cd->GetBinContent(h_fd_cd->FindBin(pMiss))
//
// Output:
//   <input>_reweight_pmiss.root with added branch
//     w_fdcd_to_cdcd_pmiss
//
// Usage:
//   root -l -b -q 'myPlots/scratch/reweightFDCDtoCDCD_pMiss.C("in.root")'
//   root -l -b -q 'myPlots/scratch/reweightFDCDtoCDCD_pMiss.C("in.root","srcTree")'

#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

namespace {

bool isFinite(double x) { return std::isfinite(x); }

double safeLeafValue(TLeaf *leaf) {
  if (!leaf) return 1.0;
  const double v = leaf->GetValue();
  if (!isFinite(v)) return 0.0;
  return v;
}

TString makeOutputName(const char *inputFileName) {
  TString outName(inputFileName);
  if (outName.EndsWith(".root")) {
    outName.ReplaceAll(".root", "_reweight_pmiss.root");
  } else {
    outName += "_reweight_pmiss.root";
  }
  return outName;
}

double toDeg(double thetaRaw, bool thetaIsDeg) {
  return thetaIsDeg ? thetaRaw : thetaRaw * TMath::RadToDeg();
}

bool isLeadFDRecCD(double leadThetaDeg, double recThetaDeg) {
  return (leadThetaDeg < 37.0 && recThetaDeg > 45.0);
}

bool isLeadCDRecCD(double leadThetaDeg, double recThetaDeg) {
  return (leadThetaDeg > 45.0 && recThetaDeg > 45.0);
}

}  // namespace

void reweightFDCDtoCDCD_pMiss(const char *inputFileName = "~/data/RGM_DATA/c12_src_skim.root",
                              const char *treeName = "srcTree",
                              const char *baseWeightBranch = "",
                              int nPMissBins = 60,
                              double pMissMin = 0.3,
                              double pMissMax = 1.0) {
  TFile *inFile = TFile::Open(inputFileName, "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "Could not open input file: " << inputFileName << "\n";
    return;
  }

  TTree *inTree = dynamic_cast<TTree *>(inFile->Get(treeName));
  if (!inTree) {
    std::cerr << "Could not find tree \"" << treeName << "\" in "
              << inputFileName << "\n";
    inFile->Close();
    return;
  }

  Float_t pMiss = -9.f;
  Float_t recTheta = -9.f;
  Float_t leadTheta = -9.f;

  inTree->SetBranchAddress("pMiss", &pMiss);
  inTree->SetBranchAddress("recTheta", &recTheta);
  inTree->SetBranchAddress("leadTheta", &leadTheta);

  TLeaf *baseWeightLeaf = nullptr;
  if (baseWeightBranch && baseWeightBranch[0] != '\0') {
    baseWeightLeaf = inTree->GetLeaf(baseWeightBranch);
    if (!baseWeightLeaf) {
      std::cerr << "Requested base weight branch not found: "
                << baseWeightBranch << "\n";
      inFile->Close();
      return;
    }
  }

  const Long64_t nEntries = inTree->GetEntries();

  double maxLeadThetaRaw = -1.0;
  double maxRecThetaRaw = -1.0;
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    if (leadTheta > maxLeadThetaRaw) maxLeadThetaRaw = leadTheta;
    if (recTheta > maxRecThetaRaw) maxRecThetaRaw = recTheta;
  }

  const bool leadThetaIsDeg = (maxLeadThetaRaw > 3.3);
  const bool recThetaIsDeg = (maxRecThetaRaw > 3.3);

  TH1D *h_fd_cd = new TH1D("h_fd_cd", "FD/CD pMiss;|p_{miss}| [GeV/c];Weighted counts",
                           nPMissBins, pMissMin, pMissMax);
  TH1D *h_cd_cd = new TH1D("h_cd_cd", "CD/CD pMiss;|p_{miss}| [GeV/c];Weighted counts",
                           nPMissBins, pMissMin, pMissMax);
  h_fd_cd->Sumw2();
  h_cd_cd->Sumw2();

  Long64_t nFDCD = 0;
  Long64_t nCDCD = 0;

  // Pass 1: only pMiss projections by detector scenario.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);
    const double baseW = safeLeafValue(baseWeightLeaf);

    if (!isFinite(pMiss) || !isFinite(leadThetaDeg) || !isFinite(recThetaDeg)) continue;
    if (!(baseW > 0.0)) continue;
    if (!(pMiss >= pMissMin && pMiss <= pMissMax)) continue;

    if (isLeadFDRecCD(leadThetaDeg, recThetaDeg)) {
      h_fd_cd->Fill(pMiss, baseW);
      ++nFDCD;
    } else if (isLeadCDRecCD(leadThetaDeg, recThetaDeg)) {
      h_cd_cd->Fill(pMiss, baseW);
      ++nCDCD;
    }
  }

  TString outName = makeOutputName(inputFileName);
  TFile *outFile = TFile::Open(outName, "RECREATE");
  TTree *outTree = inTree->CloneTree(0);

  Float_t w_fdcd_to_cdcd_pmiss = 1.0f;
  outTree->Branch("w_fdcd_to_cdcd_pmiss", &w_fdcd_to_cdcd_pmiss,
                  "w_fdcd_to_cdcd_pmiss/F");

  Long64_t nWeightedFDCD = 0;
  Long64_t nZeroDenominator = 0;

  // Pass 2: write event-by-event branch from pMiss-bin ratio.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    w_fdcd_to_cdcd_pmiss = 1.0f;

    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);

    if (isFinite(pMiss) && isFinite(leadThetaDeg) && isFinite(recThetaDeg) &&
        pMiss >= pMissMin && pMiss <= pMissMax &&
        isLeadFDRecCD(leadThetaDeg, recThetaDeg)) {
      const int binCDCD = h_cd_cd->FindBin(pMiss);
      const int binFDCD = h_fd_cd->FindBin(pMiss);
      const double num = h_cd_cd->GetBinContent(binCDCD);
      const double den = h_fd_cd->GetBinContent(binFDCD);

      if (den > 0.0) {
        w_fdcd_to_cdcd_pmiss = static_cast<Float_t>(num / den);
      } else {
        w_fdcd_to_cdcd_pmiss = 0.0f;
        ++nZeroDenominator;
      }
      ++nWeightedFDCD;
    }

    outTree->Fill();
  }

  outTree->Write();
  h_fd_cd->Write();
  h_cd_cd->Write();

  std::cout << "Input entries:             " << nEntries << "\n";
  std::cout << "Base weight branch:        "
            << ((baseWeightLeaf) ? baseWeightBranch : "(none, weight=1)") << "\n";
  std::cout << "FD/CD entries in pMiss:    " << nFDCD << "\n";
  std::cout << "CD/CD entries in pMiss:    " << nCDCD << "\n";
  std::cout << "FD/CD weighted events:     " << nWeightedFDCD << "\n";
  std::cout << "Zero denominator events:   " << nZeroDenominator << "\n";
  std::cout << "Wrote output file: " << outName << "\n";

  outFile->Close();
  inFile->Close();
}
