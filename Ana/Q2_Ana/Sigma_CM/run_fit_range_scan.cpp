#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <utility>
#include <sstream>
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
            << "  " << program << " data.root mc.root out.root [options]\n"
            << "  " << program << " --from-hipo A data.hipo sim.hipo out.root [options]\n"
            << "Extra option: --xy-ranges=0.45,0.50,0.55 --z-ranges=lo:hi,lo:hi "
            << "--beam-energy=v --max-events=N\n";
}

bool defaultZWindow(const Config& cfg) {
  return std::abs(cfg.fitZMin - (-0.5)) < 1e-12 && std::abs(cfg.fitZMax - 1.0) < 1e-12;
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
  std::vector<double> xyRanges{0.45, 0.50, 0.55};
  std::vector<std::pair<double, double>> zRanges;
  for (auto it = pos.begin(); it != pos.end();) {
    if (it->rfind("--xy-ranges=", 0) == 0) {
      xyRanges.clear();
      std::stringstream ss(it->substr(12));
      std::string item;
      while (std::getline(ss, item, ',')) xyRanges.push_back(std::atof(item.c_str()));
      it = pos.erase(it);
    } else if (it->rfind("--z-ranges=", 0) == 0) {
      zRanges.clear();
      std::stringstream ss(it->substr(11));
      std::string item;
      while (std::getline(ss, item, ',')) {
        const auto colon = item.find(':');
        if (colon == std::string::npos) {
          throw std::runtime_error("--z-ranges entries must be lo:hi pairs");
        }
        zRanges.emplace_back(std::atof(item.substr(0, colon).c_str()),
                             std::atof(item.substr(colon + 1).c_str()));
      }
      it = pos.erase(it);
    } else ++it;
  }
  if ((!fromHipo && pos.size() != 3) || (fromHipo && pos.size() != 4)) {
    usage(argv[0]);
    return 2;
  }
  if (xyRanges.empty()) {
    std::cerr << "run_fit_range_scan failed: --xy-ranges did not contain any values\n";
    return 1;
  }
  if (!zRanges.empty() && zRanges.size() != xyRanges.size()) {
    std::cerr << "run_fit_range_scan failed: --z-ranges must contain the same number "
              << "of entries as --xy-ranges\n";
    return 1;
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
    std::vector<Result> results;
    const bool keepExplicitZWindow = zRanges.empty() && !defaultZWindow(cfg);
    for (size_t i = 0; i < xyRanges.size(); ++i) {
      const double rxy = xyRanges[i];
      Config c = cfg;
      c.cutRangeXY = rxy;
      if (!zRanges.empty()) {
        c.fitZMin = zRanges[i].first;
        c.fitZMax = zRanges[i].second;
      } else if (!keepExplicitZWindow) {
        c.fitZMin = -rxy;
        c.fitZMax = 2.0 * rxy;
      }
      results.push_back(extract(data.events, mc.events, mc.sigmaGen, c));
    }
    writeResultsTree(outPath, results);
    std::cout << "Wrote " << outPath << " with " << results.size()
              << " fit-range result rows\n";
  } catch (const std::exception& e) {
    std::cerr << "run_fit_range_scan failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
