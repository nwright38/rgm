// reweightFDCDtoCDCD_4D.C
//
// Source -> target reweighting:
//   source: lead FD + rec CD
//   target: lead CD + rec CD
//
// Weight model (reduced dimensionality):
//   3D in (pMiss, recP, recThetaLab)
//
// Validity used in all loops:
//   pMiss in [0.3, 1.0]
//   recP  in [0.3, 1.0]
//   recTheta > 0, leadTheta > 0, pCM > 0
//
// Output:
//   <input>_reweight.root with added branches
//     w_fdcd_to_cdcd_4d
//     w_fdcd_to_cdcd_4d_valid
//
// Notes:
//   - Branch name is kept for compatibility, but the model is 3D.
//   - Unmapped/invalid source events get weight 0 and valid=0.
//
// Usage:
//   root -l -b -q 'myPlots/scratch/reweightFDCDtoCDCD_4D.C("in.root")'
//   root -l -b -q 'myPlots/scratch/reweightFDCDtoCDCD_4D.C("in.root","srcTree")'

#include <cmath>
#include <iostream>
#include <limits>

#include <TFile.h>
#include <TH1D.h>
#include <THnSparse.h>
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
    outName.ReplaceAll(".root", "_reweight.root");
  } else {
    outName += "_reweight.root";
  }
  return outName;
}

bool isLeadFDRecCD(double leadThetaDeg, double recThetaDeg) {
  return (leadThetaDeg < 37.0 && recThetaDeg > 45.0);
}

bool isLeadCDRecCD(double leadThetaDeg, double recThetaDeg) {
  return (leadThetaDeg > 45.0 && recThetaDeg > 45.0);
}

bool passKinematicValidity(float pCM, float pMiss, float recP,
                           float recThetaDeg, float leadThetaDeg) {
  if (!(pCM > 0.0f)) return false;
  if (!(pMiss >= 0.3f && pMiss <= 1.0f)) return false;
  if (!(recP >= 0.3f && recP <= 1.0f)) return false;
  if (!(recThetaDeg > 0.0f)) return false;
  if (!(leadThetaDeg > 0.0f)) return false;
  return true;
}

double toDeg(double thetaRaw, bool thetaIsDeg) {
  return thetaIsDeg ? thetaRaw : thetaRaw * TMath::RadToDeg();
}

void fill3(double pMiss, double recP, double recThetaDeg,
           TH1D *hPMiss, TH1D *hRecP, TH1D *hRecTheta,
           double w) {
  hPMiss->Fill(pMiss, w);
  hRecP->Fill(recP, w);
  hRecTheta->Fill(recThetaDeg, w);
}

double chi2NdfShape(const TH1D *after, const TH1D *target, int &ndfOut) {
  ndfOut = 0;
  if (!after || !target) return std::numeric_limits<double>::quiet_NaN();

  TH1D *a = static_cast<TH1D *>(after->Clone("h_tmp_after_shape"));
  TH1D *t = static_cast<TH1D *>(target->Clone("h_tmp_target_shape"));
  const double ia = a->Integral();
  const double it = t->Integral();
  if (ia > 0.0 && it > 0.0) a->Scale(it / ia);

  double chi2 = 0.0;
  for (int b = 1; b <= a->GetNbinsX(); ++b) {
    const double ya = a->GetBinContent(b);
    const double yt = t->GetBinContent(b);
    const double ea = a->GetBinError(b);
    const double et = t->GetBinError(b);
    const double var = ea * ea + et * et;
    if (var <= 0.0) continue;
    const double d = ya - yt;
    chi2 += d * d / var;
    ++ndfOut;
  }

  delete a;
  delete t;
  if (ndfOut <= 1) return std::numeric_limits<double>::quiet_NaN();
  // one effective constraint from normalization
  --ndfOut;
  if (ndfOut <= 0) return std::numeric_limits<double>::quiet_NaN();
  return chi2 / ndfOut;
}

}  // namespace

void reweightFDCDtoCDCD_4D(const char *inputFileName = "~/data/RGM_DATA/c12_src_skim.root",
                           const char *treeName = "srcTree",
                           const char *baseWeightBranch = "",
                           double maxWeight = 20.0) {
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
  Float_t recP = -9.f;
  Float_t theta_PleadQ = -9.f;
  Float_t recTheta = -9.f;
  Float_t leadTheta = -9.f;
  Float_t pCM = -9.f;

  inTree->SetBranchAddress("pMiss", &pMiss);
  inTree->SetBranchAddress("recP", &recP);
  inTree->SetBranchAddress("theta_PleadQ", &theta_PleadQ);
  inTree->SetBranchAddress("recTheta", &recTheta);
  inTree->SetBranchAddress("leadTheta", &leadTheta);
  inTree->SetBranchAddress("pCM", &pCM);

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

  // Unit verification before anything else.
  double maxLeadThetaRaw = -1.0;
  double maxRecThetaRaw = -1.0;
  double maxThetaPleadQRaw = -1.0;
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    if (leadTheta > maxLeadThetaRaw) maxLeadThetaRaw = leadTheta;
    if (recTheta > maxRecThetaRaw) maxRecThetaRaw = recTheta;
    if (theta_PleadQ > maxThetaPleadQRaw) maxThetaPleadQRaw = theta_PleadQ;
  }

  const bool leadThetaIsDeg = (maxLeadThetaRaw > 3.3);
  const bool recThetaIsDeg = (maxRecThetaRaw > 3.3);
  const bool thetaPleadQIsDeg = (maxThetaPleadQRaw > 3.3);

  const int nDim = 3;
  const int nBins[nDim] = {6, 6, 6};
  double xMin[nDim] = {
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  double xMax[nDim] = {
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};

  Long64_t nSrc = 0;
  Long64_t nTgt = 0;
  double srcSumW = 0.0;
  double tgtSumW = 0.0;
  double minThetaPleadQDeg = std::numeric_limits<double>::infinity();
  double maxThetaPleadQDeg = -std::numeric_limits<double>::infinity();

  // Pass 1a: discover phase-space bounds from FD/CD and CD/CD only.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);

    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);
    const double thetaPleadQDeg = toDeg(theta_PleadQ, thetaPleadQIsDeg);
    const double baseW = safeLeafValue(baseWeightLeaf);

    if (!passKinematicValidity(pCM, pMiss, recP, recThetaDeg, leadThetaDeg)) continue;
    if (thetaPleadQDeg < minThetaPleadQDeg) minThetaPleadQDeg = thetaPleadQDeg;
    if (thetaPleadQDeg > maxThetaPleadQDeg) maxThetaPleadQDeg = thetaPleadQDeg;
    if (!(baseW > 0.0)) continue;

    const bool src = isLeadFDRecCD(leadThetaDeg, recThetaDeg);
    const bool tgt = isLeadCDRecCD(leadThetaDeg, recThetaDeg);
    if (!(src || tgt)) continue;

    const double x[nDim] = {pMiss, recP, recThetaDeg};
    for (int d = 0; d < nDim; ++d) {
      if (x[d] < xMin[d]) xMin[d] = x[d];
      if (x[d] > xMax[d]) xMax[d] = x[d];
    }
  }

  for (int d = 0; d < nDim; ++d) {
    if (!isFinite(xMin[d]) || !isFinite(xMax[d])) {
      std::cerr << "Invalid 3D range discovered on axis " << d << "\n";
      inFile->Close();
      return;
    }
    if (xMax[d] <= xMin[d]) {
      xMin[d] -= 0.5;
      xMax[d] += 0.5;
    }
    const double span = xMax[d] - xMin[d];
    xMin[d] -= 1e-6 * span;
    xMax[d] += 1e-6 * span;
  }

  THnSparseD *hSrc = new THnSparseD(
      "hSrc",
      "FD/CD source;|p_{miss}|;|p_{rec}|;#theta_{rec}^{lab} [deg]",
      nDim, nBins, xMin, xMax);
  THnSparseD *hTgt = new THnSparseD(
      "hTgt",
      "CD/CD target;|p_{miss}|;|p_{rec}|;#theta_{rec}^{lab} [deg]",
      nDim, nBins, xMin, xMax);
  hSrc->Sumw2();
  hTgt->Sumw2();

  // Pass 1b: fill source and target 3D histograms.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);

    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);
    const double baseW = safeLeafValue(baseWeightLeaf);

    if (!passKinematicValidity(pCM, pMiss, recP, recThetaDeg, leadThetaDeg)) continue;
    if (!(baseW > 0.0)) continue;

    const double x[nDim] = {pMiss, recP, recThetaDeg};
    if (isLeadFDRecCD(leadThetaDeg, recThetaDeg)) {
      hSrc->Fill(x, baseW);
      ++nSrc;
      srcSumW += baseW;
    } else if (isLeadCDRecCD(leadThetaDeg, recThetaDeg)) {
      hTgt->Fill(x, baseW);
      ++nTgt;
      tgtSumW += baseW;
    }
  }

  if (!(nSrc > 0 && nTgt > 0 && srcSumW > 0.0 && tgtSumW > 0.0)) {
    std::cerr << "Insufficient source/target statistics to build weights."
              << " nSrc=" << nSrc << ", nTgt=" << nTgt
              << ", srcSumW=" << srcSumW << ", tgtSumW=" << tgtSumW << "\n";
    inFile->Close();
    return;
  }

  const double srcSum = srcSumW;
  const double tgtSum = tgtSumW;
  long long nCells = 1;
  for (int d = 0; d < nDim; ++d) nCells *= nBins[d];
  const double alpha = 1.0;
  const double normFactor = (srcSum + alpha * nCells) / (tgtSum + alpha * nCells);

  TString outName = makeOutputName(inputFileName);
  TFile *outFile = TFile::Open(outName, "RECREATE");
  TTree *outTree = inTree->CloneTree(0);

  Float_t w_fdcd_to_cdcd_4d = 1.0f;
  UChar_t w_fdcd_to_cdcd_4d_valid = 0;
  outTree->Branch("w_fdcd_to_cdcd_4d", &w_fdcd_to_cdcd_4d,
                  "w_fdcd_to_cdcd_4d/F");
  outTree->Branch("w_fdcd_to_cdcd_4d_valid", &w_fdcd_to_cdcd_4d_valid,
                  "w_fdcd_to_cdcd_4d_valid/b");

  Long64_t nWeighted = 0;
  Long64_t nNoMap = 0;
  Long64_t nClipped = 0;
  Long64_t nSupportZeroTgt = 0;
  double sw = 0.0;
  double sw2 = 0.0;

  TH1D *hPMissSrcBefore = new TH1D("hPMiss_src_before", "FD/CD pMiss before;|p_{miss}| [GeV/c];Weighted counts", 60, 0.3, 1.0);
  TH1D *hPMissSrcAfter = new TH1D("hPMiss_src_after", "FD/CD pMiss after;|p_{miss}| [GeV/c];Weighted counts", 60, 0.3, 1.0);
  TH1D *hPMissTgt = new TH1D("hPMiss_tgt", "CD/CD pMiss target;|p_{miss}| [GeV/c];Weighted counts", 60, 0.3, 1.0);

  TH1D *hRecPSrcBefore = new TH1D("hRecP_src_before", "FD/CD recP before;|p_{rec}| [GeV/c];Weighted counts", 60, 0.3, 1.0);
  TH1D *hRecPSrcAfter = new TH1D("hRecP_src_after", "FD/CD recP after;|p_{rec}| [GeV/c];Weighted counts", 60, 0.3, 1.0);
  TH1D *hRecPTgt = new TH1D("hRecP_tgt", "CD/CD recP target;|p_{rec}| [GeV/c];Weighted counts", 60, 0.3, 1.0);

  TH1D *hRecThetaSrcBefore = new TH1D("hRecTheta_src_before", "FD/CD recTheta before;#theta_{rec}^{lab} [deg];Weighted counts", 60, xMin[2], xMax[2]);
  TH1D *hRecThetaSrcAfter = new TH1D("hRecTheta_src_after", "FD/CD recTheta after;#theta_{rec}^{lab} [deg];Weighted counts", 60, xMin[2], xMax[2]);
  TH1D *hRecThetaTgt = new TH1D("hRecTheta_tgt", "CD/CD recTheta target;#theta_{rec}^{lab} [deg];Weighted counts", 60, xMin[2], xMax[2]);

  // Pass 2: apply source->target 3D ratio event-by-event.
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    w_fdcd_to_cdcd_4d = 1.0f;
    w_fdcd_to_cdcd_4d_valid = 0;

    const double leadThetaDeg = toDeg(leadTheta, leadThetaIsDeg);
    const double recThetaDeg = toDeg(recTheta, recThetaIsDeg);
    const double baseW = safeLeafValue(baseWeightLeaf);

    if (passKinematicValidity(pCM, pMiss, recP, recThetaDeg, leadThetaDeg)) {
      if (isLeadFDRecCD(leadThetaDeg, recThetaDeg)) {
        fill3(pMiss, recP, recThetaDeg,
              hPMissSrcBefore, hRecPSrcBefore, hRecThetaSrcBefore,
              baseW);

        const double x[nDim] = {pMiss, recP, recThetaDeg};

        Int_t idx[nDim];
        bool inRange = true;
        for (int d = 0; d < nDim; ++d) {
          idx[d] = hSrc->GetAxis(d)->FindBin(x[d]);
          if (idx[d] < 1 || idx[d] > nBins[d]) {
            inRange = false;
            break;
          }
        }

        double srcBin = 0.0;
        double tgtBin = 0.0;
        if (inRange) {
          srcBin = hSrc->GetBinContent(idx);
          tgtBin = hTgt->GetBinContent(idx);
        }

        if (inRange && srcBin > 0.0 && baseW > 0.0) {
          double rw = ((tgtBin + alpha) / (srcBin + alpha)) * normFactor;
          if (rw > maxWeight) {
            rw = maxWeight;
            ++nClipped;
          }
          if (isFinite(rw) && rw >= 0.0) {
            w_fdcd_to_cdcd_4d = static_cast<Float_t>(rw);
            w_fdcd_to_cdcd_4d_valid = (tgtBin > 0.0) ? 1 : 0;
            if (tgtBin <= 0.0) ++nSupportZeroTgt;
            ++nWeighted;
            sw += rw;
            sw2 += rw * rw;
          } else {
            w_fdcd_to_cdcd_4d = 0.0f;
            w_fdcd_to_cdcd_4d_valid = 0;
            ++nNoMap;
          }
        } else {
          // A7: unmapped events must not remain in the weighted sample.
          w_fdcd_to_cdcd_4d = 0.0f;
          w_fdcd_to_cdcd_4d_valid = 0;
          ++nNoMap;
        }

        fill3(pMiss, recP, recThetaDeg,
              hPMissSrcAfter, hRecPSrcAfter, hRecThetaSrcAfter,
              baseW * w_fdcd_to_cdcd_4d);
      } else if (isLeadCDRecCD(leadThetaDeg, recThetaDeg) && baseW > 0.0) {
        fill3(pMiss, recP, recThetaDeg,
              hPMissTgt, hRecPTgt, hRecThetaTgt,
              baseW);
      }
    }

    outTree->Fill();
  }

  int ndfPMiss = 0;
  int ndfRecP = 0;
  int ndfRecTheta = 0;
  const double chi2PMiss = chi2NdfShape(hPMissSrcAfter, hPMissTgt, ndfPMiss);
  const double chi2RecP = chi2NdfShape(hRecPSrcAfter, hRecPTgt, ndfRecP);
  const double chi2RecTheta = chi2NdfShape(hRecThetaSrcAfter, hRecThetaTgt, ndfRecTheta);

  const double nEff = (sw2 > 0.0) ? (sw * sw / sw2) : 0.0;

  outTree->Write();
  hSrc->Write();
  hTgt->Write();
  hPMissSrcBefore->Write();
  hPMissSrcAfter->Write();
  hPMissTgt->Write();
  hRecPSrcBefore->Write();
  hRecPSrcAfter->Write();
  hRecPTgt->Write();
  hRecThetaSrcBefore->Write();
  hRecThetaSrcAfter->Write();
  hRecThetaTgt->Write();

  std::cout << "Theta unit detection:"
            << " leadTheta=" << (leadThetaIsDeg ? "deg" : "rad")
            << " recTheta=" << (recThetaIsDeg ? "deg" : "rad")
            << " theta_PleadQ=" << (thetaPleadQIsDeg ? "deg" : "rad")
            << "\n";
  std::cout << "Raw max values:"
            << " leadTheta=" << maxLeadThetaRaw
            << " recTheta=" << maxRecThetaRaw
            << " theta_PleadQ=" << maxThetaPleadQRaw << "\n";
  std::cout << "theta_PleadQ range (deg): [" << minThetaPleadQDeg
            << ", " << maxThetaPleadQDeg << "]\n";

  std::cout << "Input entries:             " << nEntries << "\n";
  std::cout << "Base weight branch:        "
            << ((baseWeightLeaf) ? baseWeightBranch : "(none, weight=1)") << "\n";
  std::cout << "Source events (FD/CD):     " << nSrc << "\n";
  std::cout << "Target events (CD/CD):     " << nTgt << "\n";
  std::cout << "FD/CD events weighted:     " << nWeighted << "\n";
  std::cout << "FD/CD unmapped (w=0):      " << nNoMap << "\n";
  std::cout << "FD/CD support tgtBin==0:   " << nSupportZeroTgt
            << " (" << (nSrc > 0 ? 100.0 * nSupportZeroTgt / nSrc : 0.0) << "%)\n";
  std::cout << "Weights clipped at maxW:   " << nClipped
            << " (maxWeight=" << maxWeight << ")\n";
  std::cout << "alpha * nCells:            " << (alpha * nCells)
            << "  srcSum=" << srcSum << "  tgtSum=" << tgtSum << "\n";
  std::cout << "N_eff = " << nEff << " of " << nSrc
            << " (" << (nSrc > 0 ? 100.0 * nEff / nSrc : 0.0) << "%)\n";
  if (nSrc > 0 && nEff < 0.3 * nSrc) {
    std::cout << "WARNING: N_eff is below 30% of source sample;"
              << " reweighted comparison has limited statistical power.\n";
  }

  std::cout << "3D range used:"
            << " pMiss[" << xMin[0] << ", " << xMax[0] << "]"
            << " recP[" << xMin[1] << ", " << xMax[1] << "]"
            << " recThetaLab[" << xMin[2] << ", " << xMax[2] << "]\n";

  std::cout << "Closure chi2/ndf (after vs target, shape-normalized):\n";
  std::cout << "  pMiss:    " << chi2PMiss << " (ndf=" << ndfPMiss << ")\n";
  std::cout << "  recP:     " << chi2RecP << " (ndf=" << ndfRecP << ")\n";
  std::cout << "  recTheta: " << chi2RecTheta << " (ndf=" << ndfRecTheta << ")\n";

  std::cout << "Wrote output file: " << outName << "\n";

  outFile->Close();
  inFile->Close();
}
