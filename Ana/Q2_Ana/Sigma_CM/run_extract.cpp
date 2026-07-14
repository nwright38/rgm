#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMHipoIO.h"
#include "SigmaCMLegacyOutput.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMSkimIO.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sigmacm;

namespace {

void usage(const char* program) {
  std::cerr << "Usage:\n"
            << "  " << program << " --from-skim out.root data.root mc.root [options]\n\n"
            << "Hipo mode:\n"
            << "  " << program << " A out.root data.hipo sim.hipo [options]\n\n"
            << "Common options: --lead-mode FD|CD|BOTH --aux-weight name --seed N\n"
            << "  --independent-scales --pcm-lt-prel --cut-range-xy=v --fit-z-min=v --fit-z-max=v\n"
            << "  --beam-energy=v --max-events=N\n";
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

void summarize(const std::vector<Result>& results) {
  std::cout << "\nSigma_CM summary\n"
            << "  bin      sigmaX        sigmaY        sigmaZ        chi2/ndf   status\n";
  for (const auto& r : results) {
    const std::string bin = r.config.integratedQ2 ? "int" : ("Q2_" + std::to_string(r.config.q2BinIndex));
    const double chi2ndf = r.ndf > 0 ? r.chi2 / r.ndf : 0.0;
    std::cout << "  " << std::setw(6) << bin << "  "
              << std::fixed << std::setprecision(4)
              << r.sigmaX << " +/- " << r.sigmaXErrHigh << "  "
              << r.sigmaY << " +/- " << r.sigmaYErrHigh << "  "
              << r.sigmaZ << " +/- " << r.sigmaZErrHigh << "  "
              << std::setprecision(2) << chi2ndf << "    "
              << (r.converged ? "OK" : r.status) << "\n";
  }
}

std::vector<Result> runNominal(const Sample& data, const Sample& mc, Config cfg) {
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
  return results;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  const bool fromSkim = takeFlag(pos, "--from-skim");
  HipoLoadOptions hipoOptions;
  const std::string beamEnergy = takeOption(pos, "--beam-energy");
  if (!beamEnergy.empty()) hipoOptions.beamEnergy = std::atof(beamEnergy.c_str());
  const std::string maxEvents = takeOption(pos, "--max-events");
  if (!maxEvents.empty()) hipoOptions.maxEvents = std::atoll(maxEvents.c_str());

  if (fromSkim && pos.size() != 3) {
    usage(argv[0]);
    return 2;
  }
  if (!fromSkim && pos.size() != 4) {
    usage(argv[0]);
    return 2;
  }

  try {
    std::string outPath;
    Sample data;
    Sample mc;
    if (fromSkim) {
      outPath = pos[0];
      data = loadSkim(pos[1], false);
      mc = loadSkim(pos[2], true);
    } else {
      hipoOptions.nucleusA = std::atoi(pos[0].c_str());
      outPath = pos[1];
      data = loadHipo({pos[2]}, false, cfg, hipoOptions);
      mc = loadHipo({pos[3]}, true, cfg, hipoOptions);
    }
    std::vector<Result> results = runNominal(data, mc, cfg);
    writeResultsTree(outPath, results);
    writePlottingRootObjects(outPath, data.events, mc.events, mc.sigmaGen, results);
    summarize(results);
    std::cout << "\nWrote " << outPath << "\n";
  } catch (const std::exception& e) {
    std::cerr << "sigmacm_extract failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
