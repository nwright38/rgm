#pragma once

#include "SigmaCMConfig.h"
#include "SigmaCMEvent.h"

#include <string>
#include <vector>

namespace sigmacm {

struct HipoLoadOptions {
  int nucleusA = 4;
  double beamEnergy = 5.98636;
  long long maxEvents = -1;
  int nGcfToys = 0;
};

Sample loadSkim(const std::string& path, bool requireMC);
void writeSkim(const std::string& path, const Sample& sample, bool writeMCBranches);
Sample loadHipo(const std::vector<std::string>& paths, bool requireMC,
                const Config& cfg, const HipoLoadOptions& options = {});
bool hipoSupportAvailable();

}  // namespace sigmacm
