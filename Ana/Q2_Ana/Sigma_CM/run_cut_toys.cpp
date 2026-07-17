#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>

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
  cfg.thetaCDLower = g(cfg.thetaCDLower, 1.0);
  cfg.pRecLower = g(cfg.pRecLower, 0.045);
  return cfg;
}

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
            << "  " << program << " data.root mc.root out.root [options]\n"
            << "  " << program << " --from-hipo A data.hipo sim.hipo out.root [options]\n"
            << "Extra options: --n-cut-toys=N --n-bootstrap=N --beam-energy=v --max-events=N\n";
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
  int nCutToys = 100;
  int nBoot = 200;
  for (auto it = pos.begin(); it != pos.end();) {
    if (it->rfind("--n-cut-toys=", 0) == 0) { nCutToys = std::atoi(it->substr(13).c_str()); it = pos.erase(it); }
    else if (it->rfind("--n-bootstrap=", 0) == 0) { nBoot = std::atoi(it->substr(14).c_str()); it = pos.erase(it); }
    else ++it;
  }
  if ((!fromHipo && pos.size() != 3) || (fromHipo && pos.size() != 4)) {
    usage(argv[0]);
    return 2;
  }
  try {
    std::string outPath;
    Sample data;
    Sample mc;
    if (fromHipo) {
      hipoOptions.nucleusA = std::atoi(pos[0].c_str());
      data = loadHipo({pos[1]}, false, cfg, hipoOptions);
      mc = loadHipo({pos[2]}, true, cfg, hipoOptions);
      outPath = pos[3];
    } else {
      data = loadSkim(pos[0], false);
      mc = loadSkim(pos[1], true);
      outPath = pos[2];
    }
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
    writeResultsTree(outPath, results);
    std::cout << "Wrote " << outPath << " with " << results.size()
              << " cut/bootstrap result rows\n";
  } catch (const std::exception& e) {
    std::cerr << "run_cut_toys failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
