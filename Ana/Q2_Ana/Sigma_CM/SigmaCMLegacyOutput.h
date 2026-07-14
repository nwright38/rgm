#pragma once

#include "SigmaCMEvent.h"
#include "SigmaCMExtractor.h"

#include <string>
#include <vector>

namespace sigmacm {

void writePlottingRootObjects(const std::string& path,
                              const std::vector<Event>& dataEvents,
                              const std::vector<Event>& mcEvents,
                              double sigmaGen,
                              const std::vector<Result>& results);

void writeLegacyRootObjects(const std::string& path,
                            const std::vector<Event>& dataEvents,
                            const std::vector<Event>& mcEvents,
                            double sigmaGen,
                            const std::vector<Result>& results);

}  // namespace sigmacm
