#include "SigmaCMResultIO.h"

#include "SigmaCMConfig.h"

#include <TFile.h>
#include <TTree.h>

#include <memory>
#include <stdexcept>

namespace sigmacm {

void writeResultsTree(const std::string& path, const std::vector<Result>& results,
                      const std::string& treeName) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) throw std::runtime_error("Could not create output ROOT file '" + path + "'");
  TTree tree(treeName.c_str(), "sigma_CM extraction results");

  double sigmaX, sigmaY, sigmaZ, sigmaXErrLow, sigmaXErrHigh, sigmaYErrLow, sigmaYErrHigh;
  double sigmaZErrLow, sigmaZErrHigh, scale, scaleErr, rhoScaleSigmaX, rhoScaleSigmaY, rhoScaleSigmaZ;
  double chi2, chi2X, chi2Y, chi2Z, effectiveMCEntries;
  int ndf, q2BinIndex, leadMode, sharedScale, statOnly, converged;
  long long nEventsData;
  unsigned long long seed;
  double xBLower, q2Lower, q2Upper, q2BinLow, q2BinHigh, mMissLower, mMissUpper, kMissLower;
  double pLeadLower, thetaFDUpper;
  double thetaCDLower, pRecLower;
  double cutRangeXY, fitZMin, fitZMax;
  int binsX, binsY, binsZ, fdLeadRegionValue, cdLeadRegionValue, requirePcmLtPrel;
  std::string configJson, auxWeightBranch, status;

  tree.Branch("sigmaX", &sigmaX);
  tree.Branch("sigmaY", &sigmaY);
  tree.Branch("sigmaZ", &sigmaZ);
  tree.Branch("sigmaXErrLow", &sigmaXErrLow);
  tree.Branch("sigmaXErrHigh", &sigmaXErrHigh);
  tree.Branch("sigmaYErrLow", &sigmaYErrLow);
  tree.Branch("sigmaYErrHigh", &sigmaYErrHigh);
  tree.Branch("sigmaZErrLow", &sigmaZErrLow);
  tree.Branch("sigmaZErrHigh", &sigmaZErrHigh);
  tree.Branch("scale", &scale);
  tree.Branch("scaleErr", &scaleErr);
  tree.Branch("rhoScaleSigmaX", &rhoScaleSigmaX);
  tree.Branch("rhoScaleSigmaY", &rhoScaleSigmaY);
  tree.Branch("rhoScaleSigmaZ", &rhoScaleSigmaZ);
  tree.Branch("chi2", &chi2);
  tree.Branch("ndf", &ndf);
  tree.Branch("chi2X", &chi2X);
  tree.Branch("chi2Y", &chi2Y);
  tree.Branch("chi2Z", &chi2Z);
  tree.Branch("nEventsData", &nEventsData);
  tree.Branch("effectiveMCEntries", &effectiveMCEntries);
  tree.Branch("converged", &converged);
  tree.Branch("status", &status);
  tree.Branch("configJson", &configJson);
  tree.Branch("auxWeightBranch", &auxWeightBranch);
  tree.Branch("seed", &seed);
  tree.Branch("leadMode", &leadMode);
  tree.Branch("q2BinIndex", &q2BinIndex);
  tree.Branch("sharedScale", &sharedScale);
  tree.Branch("statOnly", &statOnly);
  tree.Branch("xBLower", &xBLower);
  tree.Branch("q2Lower", &q2Lower);
  tree.Branch("q2Upper", &q2Upper);
  tree.Branch("q2BinLow", &q2BinLow);
  tree.Branch("q2BinHigh", &q2BinHigh);
  tree.Branch("mMissLower", &mMissLower);
  tree.Branch("mMissUpper", &mMissUpper);
  tree.Branch("kMissLower", &kMissLower);
  tree.Branch("pLeadLower", &pLeadLower);
  tree.Branch("thetaFDUpper", &thetaFDUpper);
  tree.Branch("thetaCDLower", &thetaCDLower);
  tree.Branch("pRecLower", &pRecLower);
  tree.Branch("requirePcmLtPrel", &requirePcmLtPrel);
  tree.Branch("cutRangeXY", &cutRangeXY);
  tree.Branch("fitZMin", &fitZMin);
  tree.Branch("fitZMax", &fitZMax);
  tree.Branch("binsX", &binsX);
  tree.Branch("binsY", &binsY);
  tree.Branch("binsZ", &binsZ);
  tree.Branch("fdLeadRegionValue", &fdLeadRegionValue);
  tree.Branch("cdLeadRegionValue", &cdLeadRegionValue);

  TTree profile("profile", "profile chi2 scan");
  int resultIndex, axis;
  double profileSigma, profileChi2, profileSigmaX, profileSigmaY, profileSigmaZ, profileScale;
  profile.Branch("resultIndex", &resultIndex);
  profile.Branch("axis", &axis);
  profile.Branch("sigma", &profileSigma);
  profile.Branch("chi2", &profileChi2);
  profile.Branch("sigmaX", &profileSigmaX);
  profile.Branch("sigmaY", &profileSigmaY);
  profile.Branch("sigmaZ", &profileSigmaZ);
  profile.Branch("scale", &profileScale);

  for (size_t i = 0; i < results.size(); ++i) {
    const auto& r = results[i];
    sigmaX = r.sigmaX; sigmaY = r.sigmaY; sigmaZ = r.sigmaZ;
    sigmaXErrLow = r.sigmaXErrLow; sigmaXErrHigh = r.sigmaXErrHigh;
    sigmaYErrLow = r.sigmaYErrLow; sigmaYErrHigh = r.sigmaYErrHigh;
    sigmaZErrLow = r.sigmaZErrLow; sigmaZErrHigh = r.sigmaZErrHigh;
    scale = r.scale; scaleErr = r.scaleErr;
    rhoScaleSigmaX = r.rhoScaleSigmaX; rhoScaleSigmaY = r.rhoScaleSigmaY; rhoScaleSigmaZ = r.rhoScaleSigmaZ;
    chi2 = r.chi2; ndf = r.ndf;
    chi2X = r.projectionChi2[0]; chi2Y = r.projectionChi2[1]; chi2Z = r.projectionChi2[2];
    nEventsData = r.nEventsData; effectiveMCEntries = r.effectiveMCEntries;
    converged = r.converged ? 1 : 0; status = r.status;
    configJson = toJson(r.config); auxWeightBranch = r.config.auxWeightBranch;
    seed = r.config.seed; leadMode = static_cast<int>(r.config.leadMode);
    q2BinIndex = r.config.integratedQ2 ? -1 : r.config.q2BinIndex;
    sharedScale = r.config.sharedScale ? 1 : 0; statOnly = r.config.statOnly ? 1 : 0;
    xBLower = r.config.xBLower; q2Lower = r.config.q2Lower; q2Upper = r.config.q2Upper;
    q2BinLow = r.config.q2Lower;
    q2BinHigh = r.config.q2Upper;
    if (!r.config.integratedQ2 &&
        r.config.q2BinIndex >= 0 &&
        r.config.q2BinIndex + 1 < static_cast<int>(r.config.q2Edges.size())) {
      q2BinLow = r.config.q2Edges[r.config.q2BinIndex];
      q2BinHigh = r.config.q2Edges[r.config.q2BinIndex + 1];
    }
    mMissLower = r.config.mMissLower; mMissUpper = r.config.mMissUpper;
    kMissLower = r.config.kMissLower; pLeadLower = r.config.pLeadLower;
    thetaFDUpper = r.config.thetaFDUpper; thetaCDLower = r.config.thetaCDLower;
    pRecLower = r.config.pRecLower; requirePcmLtPrel = r.config.requirePcmLtPrel ? 1 : 0;
    cutRangeXY = r.config.cutRangeXY; fitZMin = r.config.fitZMin; fitZMax = r.config.fitZMax;
    binsX = r.config.binsX; binsY = r.config.binsY; binsZ = r.config.binsZ;
    fdLeadRegionValue = r.config.fdLeadRegionValue; cdLeadRegionValue = r.config.cdLeadRegionValue;
    tree.Fill();

    resultIndex = static_cast<int>(i);
    axis = r.profileAxis;
    for (const auto& p : r.profile) {
      profileSigma = p.sigma; profileChi2 = p.chi2;
      profileSigmaX = p.sigmaX; profileSigmaY = p.sigmaY; profileSigmaZ = p.sigmaZ;
      profileScale = p.scale;
      profile.Fill();
    }
  }
  tree.Write();
  profile.Write();
  file->Close();
}

std::vector<Result> readResultsTree(const std::string&, const std::string&) {
  throw std::runtime_error("readResultsTree is intentionally left to Python uproot budget tooling");
}

}  // namespace sigmacm
