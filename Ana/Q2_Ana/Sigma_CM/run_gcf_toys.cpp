#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMSkimIO.h"

#include <iostream>

using namespace sigmacm;

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  if (pos.size() != 3) {
    printCommonUsage(argv[0]);
    return 2;
  }
  try {
    Sample data = loadSkim(pos[0], false);
    Sample mc = loadSkim(pos[1], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    std::vector<Result> results;
    for (const auto& branch : mc.auxWeightBranches) {
      Config toyCfg = cfg;
      toyCfg.auxWeightBranch = branch;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, toyCfg));
    }
    writeResultsTree(pos[2], results);
  } catch (const std::exception& e) {
    std::cerr << "run_gcf_toys failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
