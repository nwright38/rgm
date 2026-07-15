#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace sigmacm;

namespace {
bool takeFlag(std::vector<std::string>& args, const std::string& flag) {
  for (auto it = args.begin(); it != args.end(); ++it) {
    if (*it == flag) {
      args.erase(it);
      return true;
    }
  }
  return false;
}

std::string takeOption(std::vector<std::string>& args, const std::string& key) {
  for (auto it = args.begin(); it != args.end();) {
    if (*it == key) {
      if (it + 1 == args.end()) throw std::runtime_error("Missing value after " + key);
      std::string value = *(it + 1);
      it = args.erase(it, it + 2);
      return value;
    }
    const std::string prefix = key + "=";
    if (it->rfind(prefix, 0) == 0) {
      std::string value = it->substr(prefix.size());
      args.erase(it);
      return value;
    }
    ++it;
  }
  return "";
}

void usage(const char* program) {
  std::cerr << "Usage:\n"
            << "  " << program << " mc.root out.root [options]\n"
            << "  " << program << " --from-hipo A sim.hipo out.root [options]\n"
            << "Extra options: --beam-energy=v --max-events=N\n";
}
}

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  const bool fromHipo = takeFlag(pos, "--from-hipo");
  HipoLoadOptions hipoOptions;
  const std::string beamEnergy = takeOption(pos, "--beam-energy");
  if (!beamEnergy.empty()) hipoOptions.beamEnergy = std::atof(beamEnergy.c_str());
  const std::string maxEvents = takeOption(pos, "--max-events");
  if (!maxEvents.empty()) hipoOptions.maxEvents = std::atoll(maxEvents.c_str());
  if ((!fromHipo && pos.size() != 2) || (fromHipo && pos.size() != 3)) {
    usage(argv[0]);
    return 2;
  }
  try {
    std::string outPath;
    Sample mc;
    if (fromHipo) {
      hipoOptions.nucleusA = std::atoi(pos[0].c_str());
      mc = loadHipo({pos[1]}, true, cfg, hipoOptions);
      outPath = pos[2];
    } else {
      mc = loadSkim(pos[0], true);
      outPath = pos[1];
    }
    cfg.fdLeadRegionValue = mc.fdLeadRegionValue;
    cfg.cdLeadRegionValue = mc.cdLeadRegionValue;
    std::vector<Event> evenPool, oddTemplates;
    for (size_t i = 0; i < mc.events.size(); ++i) {
      if (i % 2 == 0) evenPool.push_back(mc.events[i]);
      else oddTemplates.push_back(mc.events[i]);
    }
    if (evenPool.empty() || oddTemplates.empty()) {
      throw std::runtime_error("Closure needs at least one even and one odd MC event after loading");
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
      r.closureInjectedSigma = target;
      r.status += "; closure injected sigma=" + std::to_string(target) +
                  "; split=even entries pseudo-data, odd entries templates";
      results.push_back(r);
    }
    writeResultsTree(outPath, results);
    std::cout << "Wrote " << outPath << " with " << results.size()
              << " closure result rows\n";
  } catch (const std::exception& e) {
    std::cerr << "run_closure failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
