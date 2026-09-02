// reweightFDCDtoCDCD_pMiss.C
//
// Source -> target reweighting in (pMiss, thetaMiss):
//   source: lead FD + rec CD
//   target: lead CD + rec CD
//
// Event-by-event source weight:
//   w = h_cd_cd->GetBinContent(h_cd_cd->FindBin(pMiss, thetaMiss)) /
//       h_fd_cd->GetBinContent(h_fd_cd->FindBin(pMiss, thetaMiss))
//
// Output:
//   <input>_reweight_pmiss.root with added branch
//     w_fdcd2cd
//     w_fdcd2cd_has_tgt
//
// Usage:
//   root -l -b -q 'myPlots/scratch/reweightFDCDtoCDCD_pMiss.C("in.root")'
//   root -l -b -q 'myPlots/scratch/reweightFDCDtoCDCD_pMiss.C("in.root","srcTree")'

#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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
    outName.ReplaceAll(".root", "_reweight_pmisstheta.root");
  } else {
    outName += "_reweight_pmisstheta.root";
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
                              int nPMissBins = 20,
                              double pMissMin = 0.5,
                              double pMissMax = 1.0,
                              const char *thetaMissBranch = "pMissTheta",
                              int nThetaMissBins = 10,
                              double thetaMissMinDeg = 60.0,
                              double thetaMissMaxDeg = 110.0) {
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
  Float_t thetaMiss = -9.f;
  Float_t recTheta = -9.f;
  Float_t leadTheta = -9.f;

  inTree->SetBranchAddress("pMiss", &pMiss);
  inTree->SetBranchAddress("recTheta", &recTheta);
  inTree->SetBranchAddress("leadTheta", &leadTheta);

  const char *thetaBranchToUse = thetaMissBranch;
  if (!inTree->GetLeaf(thetaBranchToUse)) {
    const char *candidates[] = {"pMissTheta", "thetaMiss", "theta_miss", "theta_PmQ", nullptr};
    thetaBranchToUse = nullptr;
    for (int i = 0; candidates[i] != nullptr; ++i) {
      if (inTree->GetLeaf(candidates[i])) {
        thetaBranchToUse = candidates[i];
        break;
      }
    }
  }

  if (!thetaBranchToUse) {
    std::cerr << "Could not find theta-miss branch. Tried requested branch: "
              << thetaMissBranch << " and fallback candidates." << "\n";
    inFile->Close();
    return;
  }

  inTree->SetBranchAddress(thetaBranchToUse, &thetaMiss);

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
  double maxThetaMissRaw = -1.0;
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    if (leadTheta > maxLeadThetaRaw) maxLeadThetaRaw = leadTheta;
    if (recTheta > maxRecThetaRaw) maxRecThetaRaw = recTheta;
    if (thetaMiss > maxThetaMissRaw) maxThetaMissRaw = thetaMiss;
  }

  const bool leadThetaIsDeg = (maxLeadThetaRaw > 3.3);
  const bool recThetaIsDeg = (maxRecThetaRaw > 3.3);
  const bool thetaMissIsDeg = (maxThetaMissRaw > 3.3);

  TH2D *h_fd_cd = new TH2D("h_fd_cd", "FD/CD;|p_{miss}| [GeV/c];#theta_{miss} [deg]",
                           nPMissBins, pMissMin, pMissMax,
                           nThetaMissBins, thetaMissMinDeg, thetaMissMaxDeg);
  TH2D *h_cd_cd = new TH2D("h_cd_cd", "CD/CD;|p_{miss}| [GeV/c];#theta_{miss} [deg]",
                           nPMissBins, pMissMin, pMissMax,
                           nThetaMissBins, thetaMissMinDeg, thetaMissMaxDeg);
  h_fd_cd->Sumw2();
  h_cd_cd->Sumw2();

  Long64_t nFDCD = 0;
  Long64_t nCDCD = 0;

  // Pass 1: only (pMiss, thetaMiss) projections by detector scenario.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);
    const double thetaMissDeg = toDeg(thetaMiss, thetaMissIsDeg);
    const double baseW = safeLeafValue(baseWeightLeaf);

    if (!isFinite(pMiss) || !isFinite(thetaMissDeg) || !isFinite(leadThetaDeg) || !isFinite(recThetaDeg)) continue;
    if (!(baseW > 0.0)) continue;
    if (!(pMiss >= pMissMin && pMiss <= pMissMax)) continue;
    if (!(thetaMissDeg >= thetaMissMinDeg && thetaMissDeg <= thetaMissMaxDeg)) continue;

    if (isLeadFDRecCD(leadThetaDeg, recThetaDeg)) {
      h_fd_cd->Fill(pMiss, thetaMissDeg, baseW);
      ++nFDCD;
    } else if (isLeadCDRecCD(leadThetaDeg, recThetaDeg)) {
      h_cd_cd->Fill(pMiss, thetaMissDeg, baseW);
      ++nCDCD;
    }
  }

  TString outName = makeOutputName(inputFileName);
  TFile *outFile = TFile::Open(outName, "RECREATE");
  TTree *outTree = inTree->CloneTree(0);

  Float_t w_fdcd2cd = 1.0f;
  UChar_t w_fdcd2cd_has_tgt = 0;
  outTree->Branch("w_fdcd2cd", &w_fdcd2cd, "w_fdcd2cd/F");
  outTree->Branch("w_fdcd2cd_has_tgt", &w_fdcd2cd_has_tgt, "w_fdcd2cd_has_tgt/b");

  Long64_t nWeightedFDCD = 0;
  Long64_t nZeroDenominator = 0;
  Long64_t nZeroNumerator = 0;

  // Pass 2: write event-by-event branch from (pMiss, thetaMiss)-bin ratio.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    w_fdcd2cd = 1.0f;
    w_fdcd2cd_has_tgt = 0;

    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);
    const double thetaMissDeg = toDeg(thetaMiss, thetaMissIsDeg);

    if (isFinite(pMiss) && isFinite(thetaMissDeg) && isFinite(leadThetaDeg) && isFinite(recThetaDeg) &&
        pMiss >= pMissMin && pMiss <= pMissMax &&
        thetaMissDeg >= thetaMissMinDeg && thetaMissDeg <= thetaMissMaxDeg &&
        isLeadFDRecCD(leadThetaDeg, recThetaDeg)) {
      const int binCDCD = h_cd_cd->FindBin(pMiss, thetaMissDeg);
      const int binFDCD = h_fd_cd->FindBin(pMiss, thetaMissDeg);
      const double num = h_cd_cd->GetBinContent(binCDCD);
      const double den = h_fd_cd->GetBinContent(binFDCD);
      w_fdcd2cd_has_tgt = (num > 0.0) ? 1 : 0;

      if (den > 0.0) {
        w_fdcd2cd = static_cast<Float_t>(num / den);
        if (num <= 0.0) ++nZeroNumerator;
      } else {
        w_fdcd2cd = 0.0f;
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
  std::cout << "theta-miss branch:         " << thetaBranchToUse << "\n";
  std::cout << "Base weight branch:        "
            << ((baseWeightLeaf) ? baseWeightBranch : "(none, weight=1)") << "\n";
  std::cout << "FD/CD entries in 2D:       " << nFDCD << "\n";
  std::cout << "CD/CD entries in 2D:       " << nCDCD << "\n";
  std::cout << "FD/CD weighted events:     " << nWeightedFDCD << "\n";
  std::cout << "Zero numerator events:     " << nZeroNumerator
            << " (CD/CD bin empty)" << "\n";
  std::cout << "Zero denominator events:   " << nZeroDenominator << "\n";
  std::cout << "Wrote output file: " << outName << "\n";

  outFile->Close();
  inFile->Close();
}
