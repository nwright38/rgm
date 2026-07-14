#include "SigmaCMConfig.h"
#include "SigmaCMInput.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sigmacm;

namespace {

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
            << "  " << program << " A data|mc out.root input.hipo [more.hipo ...] "
            << "[--beam-energy=v] [--max-events=N] [--gcf-toys=N]\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> pos;
  Config cfg = configFromArgs(argc, argv, 1, pos);
  HipoLoadOptions hipoOptions;
  const std::string beamEnergy = takeOption(pos, "--beam-energy");
  if (!beamEnergy.empty()) hipoOptions.beamEnergy = std::atof(beamEnergy.c_str());
  const std::string maxEvents = takeOption(pos, "--max-events");
  if (!maxEvents.empty()) hipoOptions.maxEvents = std::atoll(maxEvents.c_str());
  const std::string gcfToys = takeOption(pos, "--gcf-toys");
  if (!gcfToys.empty()) hipoOptions.nGcfToys = std::atoi(gcfToys.c_str());

  if (pos.size() < 4 || (pos[1] != "data" && pos[1] != "mc")) {
    usage(argv[0]);
    return 2;
  }

  try {
    hipoOptions.nucleusA = std::atoi(pos[0].c_str());
    const bool isMC = pos[1] == "mc";
    if (!isMC) hipoOptions.nGcfToys = 0;
    const std::string outPath = pos[2];
    std::vector<std::string> inputs(pos.begin() + 3, pos.end());
    Sample sample = loadHipo(inputs, isMC, cfg, hipoOptions);
    writeSkim(outPath, sample, isMC);
    std::cout << "Wrote " << outPath << " with " << sample.events.size()
              << " cached events\n";
  } catch (const std::exception& e) {
    std::cerr << "sigmacm_make_skim failed: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
