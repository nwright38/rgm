#include "SigmaCMConfig.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace sigmacm {

std::string leadModeName(LeadMode mode) {
  if (mode == LeadMode::CD) return "CD";
  if (mode == LeadMode::BOTH) return "BOTH";
  return "FD";
}

LeadMode parseLeadMode(const std::string& text) {
  if (text == "fd" || text == "FD") return LeadMode::FD;
  if (text == "cd" || text == "CD") return LeadMode::CD;
  if (text == "both" || text == "BOTH") return LeadMode::BOTH;
  throw std::runtime_error("Unknown leadMode '" + text + "'; expected FD, CD, or BOTH");
}

std::string toJson(const Config& c) {
  std::ostringstream os;
  os << "{"
     << "\"xBLower\":" << c.xBLower
     << ",\"q2Lower\":" << c.q2Lower
     << ",\"q2Upper\":" << c.q2Upper
     << ",\"mMissLower\":" << c.mMissLower
     << ",\"mMissUpper\":" << c.mMissUpper
     << ",\"kMissLower\":" << c.kMissLower
     << ",\"pLeadLower\":" << c.pLeadLower
     << ",\"thetaFDUpper\":" << c.thetaFDUpper
     << ",\"thetaCDLower\":" << c.thetaCDLower
     << ",\"pRecLower\":" << c.pRecLower
     << ",\"requirePcmLtPrel\":" << (c.requirePcmLtPrel ? "true" : "false")
     << ",\"leadMode\":\"" << leadModeName(c.leadMode) << "\""
     << ",\"fdLeadRegionValue\":" << c.fdLeadRegionValue
     << ",\"cdLeadRegionValue\":" << c.cdLeadRegionValue
     << ",\"integratedQ2\":" << (c.integratedQ2 ? "true" : "false")
     << ",\"q2BinIndex\":" << c.q2BinIndex
     << ",\"binsX\":" << c.binsX
     << ",\"histXMin\":" << c.histXMin
     << ",\"histXMax\":" << c.histXMax
     << ",\"binsY\":" << c.binsY
     << ",\"histYMin\":" << c.histYMin
     << ",\"histYMax\":" << c.histYMax
     << ",\"binsZ\":" << c.binsZ
     << ",\"histZMin\":" << c.histZMin
     << ",\"histZMax\":" << c.histZMax
     << ",\"cutRangeXY\":" << c.cutRangeXY
     << ",\"fitZMin\":" << c.fitZMin
     << ",\"fitZMax\":" << c.fitZMax
     << ",\"sigmaInit\":" << c.sigmaInit
     << ",\"sigmaStep\":" << c.sigmaStep
     << ",\"sigmaMin\":" << c.sigmaMin
     << ",\"sigmaMax\":" << c.sigmaMax
     << ",\"auxWeightBranch\":\"" << c.auxWeightBranch << "\""
     << ",\"sharedScale\":" << (c.sharedScale ? "true" : "false")
     << ",\"statOnly\":" << (c.statOnly ? "true" : "false")
     << ",\"seed\":" << c.seed
     << ",\"q2Edges\":[";
  for (size_t i = 0; i < c.q2Edges.size(); ++i) {
    if (i) os << ",";
    os << c.q2Edges[i];
  }
  os << "]}";
  return os.str();
}

static double valueAfterEquals(const std::string& arg) {
  const auto pos = arg.find('=');
  if (pos == std::string::npos) throw std::runtime_error("Expected --key=value for " + arg);
  return std::atof(arg.substr(pos + 1).c_str());
}

Config configFromArgs(int argc, char** argv, int startIndex, std::vector<std::string>& positional) {
  Config cfg;
  for (int i = startIndex; i < argc; ++i) {
    std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::runtime_error("Missing value after " + arg);
      return argv[++i];
    };
    if (arg == "--lead-mode") cfg.leadMode = parseLeadMode(next());
    else if (arg.rfind("--lead-mode=", 0) == 0) cfg.leadMode = parseLeadMode(arg.substr(12));
    else if (arg == "--aux-weight") cfg.auxWeightBranch = next();
    else if (arg.rfind("--aux-weight=", 0) == 0) cfg.auxWeightBranch = arg.substr(13);
    else if (arg == "--q2-bin") { cfg.integratedQ2 = false; cfg.q2BinIndex = std::atoi(next().c_str()); }
    else if (arg.rfind("--q2-bin=", 0) == 0) { cfg.integratedQ2 = false; cfg.q2BinIndex = std::atoi(arg.substr(9).c_str()); }
    else if (arg == "--independent-scales") cfg.sharedScale = false;
    else if (arg == "--full-errors") cfg.statOnly = false;
    else if (arg == "--seed") cfg.seed = std::strtoull(next().c_str(), nullptr, 10);
    else if (arg.rfind("--seed=", 0) == 0) cfg.seed = std::strtoull(arg.substr(7).c_str(), nullptr, 10);
    else if (arg.rfind("--xB-lower=", 0) == 0) cfg.xBLower = valueAfterEquals(arg);
    else if (arg.rfind("--q2-lower=", 0) == 0) cfg.q2Lower = valueAfterEquals(arg);
    else if (arg.rfind("--q2-upper=", 0) == 0) cfg.q2Upper = valueAfterEquals(arg);
    else if (arg.rfind("--mmiss-lower=", 0) == 0) cfg.mMissLower = valueAfterEquals(arg);
    else if (arg.rfind("--mmiss-upper=", 0) == 0) cfg.mMissUpper = valueAfterEquals(arg);
    else if (arg.rfind("--kmiss-lower=", 0) == 0) cfg.kMissLower = valueAfterEquals(arg);
    else if (arg.rfind("--plead-lower=", 0) == 0) cfg.pLeadLower = valueAfterEquals(arg);
    else if (arg.rfind("--theta-fd-upper=", 0) == 0) cfg.thetaFDUpper = valueAfterEquals(arg);
    else if (arg.rfind("--theta-cd-lower=", 0) == 0) cfg.thetaCDLower = valueAfterEquals(arg);
    else if (arg.rfind("--prec-lower=", 0) == 0) cfg.pRecLower = valueAfterEquals(arg);
    else if (arg == "--pcm-lt-prel") cfg.requirePcmLtPrel = true;
    else if (arg.rfind("--cut-range-xy=", 0) == 0) cfg.cutRangeXY = valueAfterEquals(arg);
    else if (arg.rfind("--fit-z-min=", 0) == 0) cfg.fitZMin = valueAfterEquals(arg);
    else if (arg.rfind("--fit-z-max=", 0) == 0) cfg.fitZMax = valueAfterEquals(arg);
    else if (arg.rfind("--sigma-init=", 0) == 0) cfg.sigmaInit = valueAfterEquals(arg);
    else if (arg.rfind("--sigma-step=", 0) == 0) cfg.sigmaStep = valueAfterEquals(arg);
    else if (arg.rfind("--sigma-min=", 0) == 0) cfg.sigmaMin = valueAfterEquals(arg);
    else if (arg.rfind("--sigma-max=", 0) == 0) cfg.sigmaMax = valueAfterEquals(arg);
    else positional.push_back(arg);
  }
  return cfg;
}

void printCommonUsage(const char* program) {
  std::cerr << "Usage: " << program << " <data.root> <mc.root> <out.root> [options]\n"
            << "Common options: --lead-mode FD|CD|BOTH --q2-bin N --aux-weight name --seed N\n"
            << "  --independent-scales --pcm-lt-prel --cut-range-xy=v --fit-z-min=v --fit-z-max=v\n"
            << "  --sigma-min=v --sigma-max=v --sigma-init=v --sigma-step=v\n";
}

}  // namespace sigmacm
