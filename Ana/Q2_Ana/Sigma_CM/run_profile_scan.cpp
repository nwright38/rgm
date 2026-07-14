#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <iostream>

using namespace sigmacm;

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  int axis = 0;
  double min = 0.08, max = 0.26;
  int n = 61;
  for (auto it = pos.begin(); it != pos.end();) {
    if (it->rfind("--axis=", 0) == 0) { axis = std::atoi(it->substr(7).c_str()); it = pos.erase(it); }
    else if (it->rfind("--scan-min=", 0) == 0) { min = std::atof(it->substr(11).c_str()); it = pos.erase(it); }
    else if (it->rfind("--scan-max=", 0) == 0) { max = std::atof(it->substr(11).c_str()); it = pos.erase(it); }
    else if (it->rfind("--n-points=", 0) == 0) { n = std::atoi(it->substr(11).c_str()); it = pos.erase(it); }
    else ++it;
  }
  if (pos.size() != 3) {
    printCommonUsage(argv[0]);
    std::cerr << "Extra options: --axis=0|1|2 --scan-min=v --scan-max=v --n-points=N\n";
    return 2;
  }
  try {
    Sample data = loadSkim(pos[0], false);
    Sample mc = loadSkim(pos[1], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    Result r = extractProfileScan(data.events, mc.events, mc.sigmaGen, cfg, axis, min, max, n);
    writeResultsTree(pos[2], {r});
  } catch (const std::exception& e) {
    std::cerr << "run_profile_scan failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
