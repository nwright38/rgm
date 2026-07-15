#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <iostream>

using namespace sigmacm;

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  for (auto it = pos.begin(); it != pos.end();) {
    if (*it == "--from-hipo") {
      std::cerr << "run_gcf_toys cannot read hipo directly because GCF toys require "
                << "w_gcf_toy_* auxiliary weight branches. Use a skim ROOT file that "
                << "contains those branches.\n";
      return 2;
    } else if (it->rfind("--n-toys=", 0) == 0) {
      std::cerr << "run_gcf_toys uses every w_gcf_toy_* branch in the MC file; ignoring "
                << *it << "\n";
      it = pos.erase(it);
    } else {
      ++it;
    }
  }
  if (pos.size() != 3) {
    printCommonUsage(argv[0]);
    return 2;
  }
  try {
    Sample data = loadSkim(pos[0], false);
    Sample mc = loadSkim(pos[1], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    if (mc.auxWeightBranches.empty()) {
      std::cout << "No w_gcf_toy_* branches found in " << pos[1]
                << "; skipping GCF toys and not writing " << pos[2] << "\n"
                << "Make the MC skim with, for example:\n"
                << "  sigmacm_make_skim A mc mc_skim.root sim.hipo --gcf-toys=100\n"
                << "If you already passed --gcf-toys, rebuild from the full repository "
                << "so SIGMACM_WITH_GCF_TOYS is enabled.\n";
      return 0;
    }
    std::vector<Result> results;
    for (const auto& branch : mc.auxWeightBranches) {
      Config toyCfg = cfg;
      toyCfg.auxWeightBranch = branch;
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, toyCfg));
    }
    writeResultsTree(pos[2], results);
    std::cout << "Wrote " << pos[2] << " with " << results.size()
              << " GCF toy result rows\n";
  } catch (const std::exception& e) {
    std::cerr << "run_gcf_toys failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
