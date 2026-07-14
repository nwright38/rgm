#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMLegacyOutput.h"
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
    cfg.integratedQ2 = true;
    cfg.q2BinIndex = -1;
    results.push_back(extract(data.events, mc.events, mc.sigmaGen, cfg));
    for (int i = 0; i + 1 < static_cast<int>(cfg.q2Edges.size()); ++i) {
      Config binCfg = cfg;
      binCfg.integratedQ2 = false;
      binCfg.q2BinIndex = i;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, binCfg));
    }
    writeResultsTree(pos[2], results);
    writePlottingRootObjects(pos[2], data.events, mc.events, mc.sigmaGen, results);
  } catch (const std::exception& e) {
    std::cerr << "run_nominal failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
