#include "SigmaCMConfig.h"
#include "SigmaCMExtractor.h"
#include "SigmaCMResultIO.h"
#include "SigmaCMInput.h"

#include <cstdlib>
#include <iostream>
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
            << "Extra options: --axis=0|1|2 --scan-min=v --scan-max=v --n-points=N "
            << "--beam-energy=v --max-events=N\n";
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
    Result r = extractProfileScan(data.events, mc.events, mc.sigmaGen, cfg, axis, min, max, n);
    writeResultsTree(outPath, {r});
  } catch (const std::exception& e) {
    std::cerr << "run_profile_scan failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
