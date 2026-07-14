#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

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
    for (int i = 0; i + 1 < static_cast<int>(cfg.q2Edges.size()); ++i) {
      Config c = cfg;
      c.integratedQ2 = false;
      c.q2BinIndex = i;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, c));
    }
    for (double theta = 32.0; theta <= 42.0001; theta += 1.0) {
      Config c = cfg;
      c.thetaFDUpper = theta;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, c));
    }
    writeResultsTree(pos[2], results);
  } catch (const std::exception& e) {
    std::cerr << "run_diagnostics failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
