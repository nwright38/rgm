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
};

Sample loadHipo(const std::vector<std::string>& paths, bool requireMC,
                const Config& cfg, const HipoLoadOptions& options = {});

bool hipoSupportAvailable();

}  // namespace sigmacm
