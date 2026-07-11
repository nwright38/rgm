#pragma once

#include "SigmaCMEvent.h"

#include <string>

namespace sigmacm {

Sample loadSkim(const std::string& path, bool requireMC);

}  // namespace sigmacm
