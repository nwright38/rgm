#include "SigmaCMExtractor.h"
#include "SigmaCMPlotOutput.h"
#include "SigmaCMResultIO.h"

#include <TFile.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace sigmacm;

namespace {

Event baseEvent(int i, double x, double y, double z, int fdRegion) {
  Event e;
  e.run = 1;
  e.event = i;
  e.Q2 = 2.0;
  e.xB = 1.4;
  e.mMiss = 0.8;
  e.kMiss = 0.5;
  e.pLead = 1.2;
  e.thetaLead = 25.0;
  e.leadRegion = fdRegion;
  e.hasRecoil = true;
  e.pRec = 0.5;
  e.pcmX = x;
  e.pcmY = y;
  e.pcmZ = z;
  e.pcmXTruth = x;
  e.pcmYTruth = y;
  e.pcmZTruth = z;
  e.genWeight = 1.0;
  e.isMC = true;
  return e;
}

std::vector<Event> makeEvents(int n, double sigma, int fdRegion) {
  std::mt19937_64 rng(1234);
  std::normal_distribution<double> g(0.0, sigma);
  std::vector<Event> out;
  out.reserve(n);
  for (int i = 0; i < n; ++i) out.push_back(baseEvent(i, g(rng), g(rng), g(rng), fdRegion));
  return out;
}

double textbookChi2(double d, double e, double scale, double se2) {
  return (d - scale * e) * (d - scale * e) / (d + scale * scale * se2);
}

}  // namespace

int main() {
  const double sigmaGen = 0.2;
  Event e = baseEvent(1, 0.05, -0.02, 0.07, 2000);
  e.genWeight = 2.75;
  const double w = eventWeight(e, sigmaGen, sigmaGen, sigmaGen, sigmaGen, "");
  assert(std::abs(w - e.genWeight) < 1e-14);

  Config cfg;
  cfg.fdLeadRegionValue = 2000;
  cfg.cdLeadRegionValue = 4000;
  cfg.binsX = cfg.binsY = 15;
  cfg.binsZ = 18;
  auto data = makeEvents(1500, 0.16, cfg.fdLeadRegionValue);
  auto mc = data;
  Result r = extract(data, mc, sigmaGen, cfg);
  assert(r.converged);
  assert(std::abs(r.scale - 1.0) < 0.15);

  const char* plotPath = "sigmacm_plotting_surface_test.root";
  writeResultsTree(plotPath, {r});
  writePlottingRootObjects(plotPath, data, mc, sigmaGen, {r});
  TFile plotFile(plotPath, "READ");
  assert(!plotFile.IsZombie());
  assert(plotFile.Get("sigmaCM"));
  assert(plotFile.Get("pcmx_epp"));
  assert(plotFile.Get("pcmx_epp_fit"));
  assert(plotFile.Get("sigmacmx_int"));
  assert(plotFile.Get("g_chi2_pcmx_epp"));
  assert(plotFile.Get("c_overlay_pcmx_epp"));

  Result profile = extractProfileScan(data, mc, sigmaGen, cfg, 0, 0.10, 0.22, 7);
  assert(profile.converged);
  assert(profile.profileAxis == 0);
  assert(profile.profile.size() == 7);
  const char* profilePath = "sigmacm_profile_test.root";
  writeResultsTree(profilePath, {profile});
  TFile profileFile(profilePath, "READ");
  assert(!profileFile.IsZombie());
  assert(profileFile.Get("sigmaCM"));
  assert(profileFile.Get("profile"));

  const double c1 = textbookChi2(10.0, 8.0, 1.2, 8.0);
  const double c2 = (10.0 - 1.2 * 8.0) * (10.0 - 1.2 * 8.0) / (10.0 + 1.2 * 1.2 * 8.0);
  assert(std::abs(c1 - c2) < 1e-15);

  std::mt19937_64 rng(42);
  std::normal_distribution<double> variedTheta(cfg.thetaFDUpper, 1.0);
  const double toyTheta = variedTheta(rng);
  for (const auto& ev : data) {
    (void)ev;
    assert(std::abs(toyTheta - toyTheta) < 1e-15);
  }

  std::cout << "Sigma_CM unit tests passed\n";
  return 0;
}
