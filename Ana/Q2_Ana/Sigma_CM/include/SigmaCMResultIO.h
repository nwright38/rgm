#pragma once

#include "SigmaCMExtractor.h"

#include <string>
#include <vector>

namespace sigmacm {

void writeResultsTree(const std::string& path, const std::vector<Result>& results,
                      const std::string& treeName = "sigmaCM");
std::vector<Result> readResultsTree(const std::string& path,
                                    const std::string& treeName = "sigmaCM");

}  // namespace sigmacm
