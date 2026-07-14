#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <iostream>
#include <random>

using namespace sigmacm;

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  int nToys = 100;
  for (auto it = pos.begin(); it != pos.end();) {
    if (it->rfind("--n-toys=", 0) == 0) { nToys = std::atoi(it->substr(9).c_str()); it = pos.erase(it); }
    else ++it;
  }
  if (pos.size() != 3) {
    printCommonUsage(argv[0]);
    std::cerr << "Extra option: --n-toys=N\n";
    return 2;
  }
  try {
    Sample data = loadSkim(pos[0], false);
    Sample mc = loadSkim(pos[1], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    std::mt19937_64 rng(cfg.seed);
    std::vector<Result> results;
    std::uniform_int_distribution<size_t> branchPick(0, mc.auxWeightBranches.empty() ? 0 : mc.auxWeightBranches.size() - 1);
    auto g = [&](double mean, double width) { return std::normal_distribution<double>(mean, width)(rng); };
    for (int i = 0; i < nToys; ++i) {
      Config toyCfg = cfg;
      toyCfg.xBLower = g(cfg.xBLower, 0.01);
      toyCfg.q2Lower = g(cfg.q2Lower, 0.01);
      toyCfg.mMissLower = g(cfg.mMissLower, 0.03);
      toyCfg.mMissUpper = g(cfg.mMissUpper, 0.03);
      toyCfg.kMissLower = g(cfg.kMissLower, 0.03);
      toyCfg.pLeadLower = g(cfg.pLeadLower, 0.03);
      toyCfg.thetaFDUpper = g(cfg.thetaFDUpper, 1.0);
      toyCfg.pRecLower = g(cfg.pRecLower, 0.045);
      if (!mc.auxWeightBranches.empty()) toyCfg.auxWeightBranch = mc.auxWeightBranches[branchPick(rng)];
      toyCfg.seed = cfg.seed + i;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, toyCfg));
    }
    writeResultsTree(pos[2], results);
  } catch (const std::exception& e) {
    std::cerr << "run_combined_toys failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
