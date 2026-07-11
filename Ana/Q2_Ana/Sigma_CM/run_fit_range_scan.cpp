#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMSkimIO.h"

#include <iostream>
#include <sstream>

using namespace sigmacm;

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  std::vector<double> xyRanges{0.45, 0.50, 0.55};
  for (auto it = pos.begin(); it != pos.end();) {
    if (it->rfind("--xy-ranges=", 0) == 0) {
      xyRanges.clear();
      std::stringstream ss(it->substr(12));
      std::string item;
      while (std::getline(ss, item, ',')) xyRanges.push_back(std::atof(item.c_str()));
      it = pos.erase(it);
    } else ++it;
  }
  if (pos.size() != 3) {
    printCommonUsage(argv[0]);
    std::cerr << "Extra option: --xy-ranges=0.45,0.50,0.55\n";
    return 2;
  }
  try {
    Sample data = loadSkim(pos[0], false);
    Sample mc = loadSkim(pos[1], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    std::vector<Result> results;
    for (double rxy : xyRanges) {
      Config c = cfg;
      c.cutRangeXY = rxy;
      c.fitZMin = -rxy;
      c.fitZMax = 2.0 * rxy;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, c));
    }
    writeResultsTree(pos[2], results);
  } catch (const std::exception& e) {
    std::cerr << "run_fit_range_scan failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
