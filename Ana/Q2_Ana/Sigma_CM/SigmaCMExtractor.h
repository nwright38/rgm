#pragma once

#include "SigmaCMConfig.h"
#include "SigmaCMEvent.h"

#include <array>
#include <limits>
#include <vector>

namespace sigmacm {

struct ProfilePoint {
  double sigma = 0.0;
  double chi2 = 0.0;
  double sigmaX = 0.0;
  double sigmaY = 0.0;
  double sigmaZ = 0.0;
  double scale = 0.0;
};

struct Result {
  Config config;
  double sigmaX = 0.0;
  double sigmaY = 0.0;
  double sigmaZ = 0.0;
  double sigmaXErrLow = 0.0;
  double sigmaXErrHigh = 0.0;
  double sigmaYErrLow = 0.0;
  double sigmaYErrHigh = 0.0;
  double sigmaZErrLow = 0.0;
  double sigmaZErrHigh = 0.0;
  double scale = 0.0;
  double scaleErr = 0.0;
  double rhoScaleSigmaX = 0.0;
  double rhoScaleSigmaY = 0.0;
  double rhoScaleSigmaZ = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
  std::array<double, 3> projectionChi2{{0.0, 0.0, 0.0}};
  long long nEventsData = 0;
  double effectiveMCEntries = 0.0;
  bool converged = false;
  std::string status;
  double closureInjectedSigma = std::numeric_limits<double>::quiet_NaN();
  int profileAxis = -1;
  std::vector<ProfilePoint> profile;
};

double gaussianRatio(double p, double sigmaTarget, double sigmaGen);
double eventWeight(const Event& event, double sigmaX, double sigmaY, double sigmaZ,
                   double sigmaGen, const std::string& auxWeightBranch);
bool passesCuts(const Event& event, const Config& cfg);
Result extract(const std::vector<Event>& dataEvents, const std::vector<Event>& mcEvents,
               double sigmaGen, const Config& cfg);
Result extractProfileScan(const std::vector<Event>& dataEvents, const std::vector<Event>& mcEvents,
                          double sigmaGen, const Config& cfg, int axis,
                          double scanMin, double scanMax, int nPoints);

}  // namespace sigmacm
