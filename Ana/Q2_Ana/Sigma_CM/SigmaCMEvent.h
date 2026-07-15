#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace sigmacm {

struct Event {
  int run = 0;
  long long event = 0;
  double Q2 = 0.0;
  double xB = 0.0;
  double mMiss = 0.0;
  double kMiss = 0.0;
  double pMiss = 0.0;
  double pLead = 0.0;
  double thetaLead = 0.0;
  int leadRegion = 0;
  bool hasRecoil = false;
  double pRec = 0.0;
  double thetaRec = 0.0;
  double pcmX = 0.0;
  double pcmY = 0.0;
  double pcmZ = 0.0;
  double pRel = 0.0;
  double pCM = 0.0;

  bool isMC = false;
  double genWeight = 1.0;
  double pcmXTruth = 0.0;
  double pcmYTruth = 0.0;
  double pcmZTruth = 0.0;
  std::map<std::string, double> auxWeights;

  double bootstrapWeight = 1.0;
  std::uint64_t stableIndex = 0;
};

struct GcfToyParams {
  std::string branchName;
  double Cpp0 = 0.0;
  double Cpn0 = 0.0;
  double Cnn0 = 0.0;
  double Cpn1 = 0.0;
  double sigmaCM = 0.0;
  double P00 = 0.0;
  double P01 = 0.0;
  double P10 = 0.0;
  double P11 = 0.0;
  double P20 = 0.0;
  double P21 = 0.0;
  double P30 = 0.0;
  double P31 = 0.0;
  double TN = 0.0;
  double TNN = 0.0;
};

struct Sample {
  std::vector<Event> events;
  double sigmaGen = 0.0;
  int fdLeadRegionValue = 2000;
  int cdLeadRegionValue = 4000;
  std::vector<std::string> auxWeightBranches;
  std::vector<GcfToyParams> gcfToyParams;
  std::string metadataDump;
};

}  // namespace sigmacm
