#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace sigmacm {

enum class LeadMode { FD, CD, BOTH };

struct Config {
  double xBLower = 1.2;
  double q2Lower = 1.5;
  double q2Upper = 5.0;
  double mMissLower = 0.65;
  double mMissUpper = 1.10;
  double kMissLower = 0.3;
  double pLeadLower = 1.0;
  double thetaFDUpper = 37.0;
  double thetaCDLower = 45.0;
  double pRecLower = 0.3;
  bool requirePcmLtPrel = false;

  LeadMode leadMode = LeadMode::FD;
  int fdLeadRegionValue = 2000;
  int cdLeadRegionValue = 4000;

  bool integratedQ2 = true;
  int q2BinIndex = -1;
  std::vector<double> q2Edges{1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};

  int binsX = 21;
  double histXMin = -0.55;
  double histXMax = 0.55;
  int binsY = 21;
  double histYMin = -0.55;
  double histYMax = 0.55;
  int binsZ = 24;
  double histZMin = -0.5;
  double histZMax = 1.0;

  double cutRangeXY = 0.55;
  double fitZMin = -0.5;
  double fitZMax = 1.0;

  double sigmaInit = 0.16;
  double sigmaStep = 0.003;
  double sigmaMin = 0.04;
  double sigmaMax = 0.40;

  std::string auxWeightBranch;
  bool sharedScale = true;
  bool statOnly = true;
  std::uint64_t seed = 1;

  bool useBootstrapWeights = false;
};

struct FitRange {
  double cutRangeXY = 0.55;
  double fitZMin = -0.5;
  double fitZMax = 1.0;
};

std::string leadModeName(LeadMode mode);
LeadMode parseLeadMode(const std::string& text);

std::string toJson(const Config& cfg);
Config configFromArgs(int argc, char** argv, int startIndex, std::vector<std::string>& positional);
void printCommonUsage(const char* program);

}  // namespace sigmacm
