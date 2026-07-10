#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMSkimIO.h"

#include <iostream>
#include <random>

using namespace sigmacm;

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  if (pos.size() != 2) {
    std::cerr << "Usage: " << argv[0] << " <mc.root> <out.root> [options]\n";
    return 2;
  }
  try {
    Sample mc = loadSkim(pos[0], true);
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    std::vector<Event> evenPool, oddTemplates;
    for (size_t i = 0; i < mc.events.size(); ++i) {
      if (i % 2 == 0) evenPool.push_back(mc.events[i]);
      else oddTemplates.push_back(mc.events[i]);
    }
    std::vector<Result> results;
    std::mt19937_64 rng(cfg.seed);
    const double targets[] = {0.10, 0.14, 0.18, 0.22};
    for (double target : targets) {
      auto pseudo = evenPool;
      double sumW = 0.0;
      for (const auto& e : pseudo) sumW += eventWeight(e, target, target, target, mc.sigmaGen, "");
      const double norm = sumW > 0.0 ? static_cast<double>(pseudo.size()) / sumW : 1.0;
      for (auto& e : pseudo) {
        const double lambda = std::max(0.0, norm * eventWeight(e, target, target, target, mc.sigmaGen, ""));
        e.bootstrapWeight = std::poisson_distribution<int>(lambda)(rng);
        e.isMC = false;
      }
      Config c = cfg;
      c.useBootstrapWeights = true;
      c.seed = cfg.seed + static_cast<std::uint64_t>(target * 1000.0);
      Result r = extract(pseudo, oddTemplates, mc.sigmaGen, c);
      r.status += "; closure injected sigma=" + std::to_string(target) +
                  "; split=even entries pseudo-data, odd entries templates";
      results.push_back(r);
    }
    writeResultsTree(pos[1], results);
  } catch (const std::exception& e) {
    std::cerr << "run_closure failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
