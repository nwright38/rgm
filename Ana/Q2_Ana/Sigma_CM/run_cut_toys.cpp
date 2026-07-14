#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <algorithm>
#include <iostream>
#include <random>

using namespace sigmacm;

namespace {
Config varied(Config cfg, std::mt19937_64& rng) {
  auto g = [&](double mean, double width) { return std::normal_distribution<double>(mean, width)(rng); };
  cfg.xBLower = g(cfg.xBLower, 0.01);
  cfg.q2Lower = g(cfg.q2Lower, 0.01);
  cfg.mMissLower = g(cfg.mMissLower, 0.03);
  cfg.mMissUpper = g(cfg.mMissUpper, 0.03);
  cfg.kMissLower = g(cfg.kMissLower, 0.03);
  cfg.pLeadLower = g(cfg.pLeadLower, 0.03);
  cfg.thetaFDUpper = g(cfg.thetaFDUpper, 1.0);
  cfg.pRecLower = g(cfg.pRecLower, 0.045);
  return cfg;
}
}

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  int nCutToys = 100;
  int nBoot = 200;
  for (auto it = pos.begin(); it != pos.end();) {
    if (it->rfind("--n-cut-toys=", 0) == 0) { nCutToys = std::atoi(it->substr(13).c_str()); it = pos.erase(it); }
    else if (it->rfind("--n-bootstrap=", 0) == 0) { nBoot = std::atoi(it->substr(14).c_str()); it = pos.erase(it); }
    else ++it;
  }
  if (pos.size() != 3) {
    printCommonUsage(argv[0]);
    std::cerr << "Extra options: --n-cut-toys=N --n-bootstrap=N\n";
    return 2;
  }
  try {
    Sample data = loadSkim(pos[0], false);
    Sample mc = loadSkim(pos[1], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    std::mt19937_64 rng(cfg.seed);
    std::vector<Result> results;
    for (int i = 0; i < nCutToys; ++i) {
      Config toyCfg = varied(cfg, rng);
      toyCfg.seed = cfg.seed + i;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, toyCfg));
    }
    for (int i = 0; i < nBoot; ++i) {
      auto bootData = data.events;
      std::poisson_distribution<int> pois(1.0);
      for (auto& e : bootData) e.bootstrapWeight = pois(rng);
      Config bootCfg = cfg;
      bootCfg.useBootstrapWeights = true;
      bootCfg.seed = cfg.seed + 100000 + i;
      results.push_back(extract(bootData, mc.events, mc.sigmaGen, bootCfg));
    }
    writeResultsTree(pos[2], results);
  } catch (const std::exception& e) {
    std::cerr << "run_cut_toys failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
