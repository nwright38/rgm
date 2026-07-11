#include "SigmaCMSkimIO.h"

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace sigmacm {
namespace {

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

bool regexNumber(const std::string& text, const std::string& key, double& value) {
  const std::regex r(key + R"(\s*[:=]\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?))",
                     std::regex::icase);
  std::smatch m;
  if (std::regex_search(text, m, r)) {
    value = std::stod(m[1].str());
    return true;
  }
  return false;
}

bool regexInt(const std::string& text, const std::string& key, int& value) {
  double d = 0.0;
  if (!regexNumber(text, key, d)) return false;
  value = static_cast<int>(d);
  return true;
}

void requireBranch(TTree* tree, const std::string& name) {
  if (!tree->GetBranch(name.c_str())) {
    throw std::runtime_error("Input skim is missing required branch '" + name + "'");
  }
}

std::vector<std::string> discoverAuxBranches(TTree* tree) {
  std::vector<std::string> out;
  TObjArray* branches = tree->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i) {
    const std::string name = branches->At(i)->GetName();
    if (name.rfind("w_gcf_toy_", 0) == 0) out.push_back(name);
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::string collectNamedMetadata(TFile& file) {
  std::ostringstream os;
  TIter next(file.GetListOfKeys());
  while (TKey* key = static_cast<TKey*>(next())) {
    std::unique_ptr<TObject> obj(key->ReadObj());
    if (auto* named = dynamic_cast<TNamed*>(obj.get())) {
      os << named->GetName() << "\n" << named->GetTitle() << "\n";
    }
  }
  return os.str();
}

bool readParamsTree(TFile& file, Sample& sample, bool& fdFound, bool& cdFound) {
  TTree* params = dynamic_cast<TTree*>(file.Get("params"));
  if (!params) return false;
  bool found = false;
  auto readDouble = [&](const char* name, double& target) {
    TLeaf* leaf = params->GetLeaf(name);
    if (!leaf) return;
    params->GetEntry(0);
    target = leaf->GetValue();
    found = true;
  };
  auto readInt = [&](const char* name, int& target) {
    TLeaf* leaf = params->GetLeaf(name);
    if (!leaf) return;
    params->GetEntry(0);
    target = static_cast<int>(leaf->GetValue());
    found = true;
  };
  readDouble("sigma_gen", sample.sigmaGen);
  auto readFd = [&](const char* name) {
    const int before = sample.fdLeadRegionValue;
    readInt(name, sample.fdLeadRegionValue);
    if (params->GetLeaf(name) || sample.fdLeadRegionValue != before) fdFound = true;
  };
  auto readCd = [&](const char* name) {
    const int before = sample.cdLeadRegionValue;
    readInt(name, sample.cdLeadRegionValue);
    if (params->GetLeaf(name) || sample.cdLeadRegionValue != before) cdFound = true;
  };
  readFd("leadRegionFD");
  readFd("lead_region_fd");
  readFd("fdLeadRegionValue");
  readCd("leadRegionCD");
  readCd("lead_region_cd");
  readCd("cdLeadRegionValue");
  return found;
}

}  // namespace

Sample loadSkim(const std::string& path, bool requireMC) {
  TH1::AddDirectory(kFALSE);
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("Could not open skim ROOT file '" + path + "'");
  TTree* tree = dynamic_cast<TTree*>(file->Get("srcTree"));
  if (!tree) tree = dynamic_cast<TTree*>(file->Get("skim"));
  if (!tree) throw std::runtime_error("Input file '" + path + "' is missing TTree 'srcTree' or 'skim'");

  const char* common[] = {"run", "event", "Q2", "xB", "mMiss", "kMiss", "pMiss", "pLead",
                          "thetaLead", "leadRegion", "hasRecoil", "pRec", "thetaRec",
                          "pcmX", "pcmY", "pcmZ", "pRel", "pCM"};
  for (const char* b : common) requireBranch(tree, b);
  if (requireMC) {
    const char* mc[] = {"isMC", "genWeight", "pcmX_truth", "pcmY_truth", "pcmZ_truth"};
    for (const char* b : mc) requireBranch(tree, b);
  }

  Sample sample;
  sample.auxWeightBranches = discoverAuxBranches(tree);
  sample.metadataDump = collectNamedMetadata(*file);
  bool fdFound = false;
  bool cdFound = false;
  readParamsTree(*file, sample, fdFound, cdFound);

  double sigmaGen = sample.sigmaGen;
  if (sigmaGen <= 0.0) regexNumber(sample.metadataDump, "sigma_gen", sigmaGen);
  sample.sigmaGen = sigmaGen;

  fdFound = regexInt(sample.metadataDump, "leadRegionFD", sample.fdLeadRegionValue) ||
            regexInt(sample.metadataDump, "lead_region_fd", sample.fdLeadRegionValue) ||
            regexInt(sample.metadataDump, "fdLeadRegionValue", sample.fdLeadRegionValue) || fdFound;
  cdFound = regexInt(sample.metadataDump, "leadRegionCD", sample.cdLeadRegionValue) ||
            regexInt(sample.metadataDump, "lead_region_cd", sample.cdLeadRegionValue) ||
            regexInt(sample.metadataDump, "cdLeadRegionValue", sample.cdLeadRegionValue) || cdFound;
  if (sample.sigmaGen <= 0.0) {
    throw std::runtime_error("Input file '" + path + "' is missing required metadata 'sigma_gen'");
  }
  if (!fdFound || !cdFound) {
    throw std::runtime_error("Input file '" + path +
                             "' is missing required lead-region enum metadata "
                             "(expected leadRegionFD/leadRegionCD or lead_region_fd/lead_region_cd)");
  }

  Event e;
  Bool_t hasRecoil = false;
  Bool_t isMC = false;
  tree->SetBranchAddress("run", &e.run);
  tree->SetBranchAddress("event", &e.event);
  tree->SetBranchAddress("Q2", &e.Q2);
  tree->SetBranchAddress("xB", &e.xB);
  tree->SetBranchAddress("mMiss", &e.mMiss);
  tree->SetBranchAddress("kMiss", &e.kMiss);
  tree->SetBranchAddress("pMiss", &e.pMiss);
  tree->SetBranchAddress("pLead", &e.pLead);
  tree->SetBranchAddress("thetaLead", &e.thetaLead);
  tree->SetBranchAddress("leadRegion", &e.leadRegion);
  tree->SetBranchAddress("hasRecoil", &hasRecoil);
  tree->SetBranchAddress("pRec", &e.pRec);
  tree->SetBranchAddress("thetaRec", &e.thetaRec);
  tree->SetBranchAddress("pcmX", &e.pcmX);
  tree->SetBranchAddress("pcmY", &e.pcmY);
  tree->SetBranchAddress("pcmZ", &e.pcmZ);
  tree->SetBranchAddress("pRel", &e.pRel);
  tree->SetBranchAddress("pCM", &e.pCM);
  if (requireMC) {
    tree->SetBranchAddress("isMC", &isMC);
    tree->SetBranchAddress("genWeight", &e.genWeight);
    tree->SetBranchAddress("pcmX_truth", &e.pcmXTruth);
    tree->SetBranchAddress("pcmY_truth", &e.pcmYTruth);
    tree->SetBranchAddress("pcmZ_truth", &e.pcmZTruth);
  }
  std::vector<double> auxValues(sample.auxWeightBranches.size(), 1.0);
  for (size_t i = 0; i < sample.auxWeightBranches.size(); ++i) {
    tree->SetBranchAddress(sample.auxWeightBranches[i].c_str(), &auxValues[i]);
  }

  const Long64_t n = tree->GetEntries();
  sample.events.reserve(static_cast<size_t>(n));
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    e.hasRecoil = hasRecoil;
    e.isMC = requireMC ? static_cast<bool>(isMC) : false;
    e.stableIndex = static_cast<std::uint64_t>(i);
    e.auxWeights.clear();
    for (size_t j = 0; j < sample.auxWeightBranches.size(); ++j) {
      e.auxWeights.emplace(sample.auxWeightBranches[j], auxValues[j]);
    }
    sample.events.push_back(e);
  }
  return sample;
}

}  // namespace sigmacm
