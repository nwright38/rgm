#include "SigmaCMInput.h"
#include "SigmaCMExtractor.h"

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <iomanip>
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

void writeSkim(const std::string& path, const Sample& sample, bool writeMCBranches) {
  TH1::AddDirectory(kFALSE);
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) throw std::runtime_error("Could not create skim ROOT file '" + path + "'");

  TTree tree("srcTree", "Sigma_CM cached event skim");
  Event e;
  Bool_t hasRecoil = false;
  Bool_t isMC = false;
  tree.Branch("run", &e.run);
  tree.Branch("event", &e.event);
  tree.Branch("Q2", &e.Q2);
  tree.Branch("xB", &e.xB);
  tree.Branch("mMiss", &e.mMiss);
  tree.Branch("kMiss", &e.kMiss);
  tree.Branch("pMiss", &e.pMiss);
  tree.Branch("pLead", &e.pLead);
  tree.Branch("thetaLead", &e.thetaLead);
  tree.Branch("leadRegion", &e.leadRegion);
  tree.Branch("hasRecoil", &hasRecoil);
  tree.Branch("pRec", &e.pRec);
  tree.Branch("thetaRec", &e.thetaRec);
  tree.Branch("pcmX", &e.pcmX);
  tree.Branch("pcmY", &e.pcmY);
  tree.Branch("pcmZ", &e.pcmZ);
  tree.Branch("pRel", &e.pRel);
  tree.Branch("pCM", &e.pCM);
  if (writeMCBranches) {
    tree.Branch("isMC", &isMC);
    tree.Branch("genWeight", &e.genWeight);
    tree.Branch("pcmX_truth", &e.pcmXTruth);
    tree.Branch("pcmY_truth", &e.pcmYTruth);
    tree.Branch("pcmZ_truth", &e.pcmZTruth);
  }
  std::vector<double> auxValues(sample.auxWeightBranches.size(), 1.0);
  for (size_t i = 0; i < sample.auxWeightBranches.size(); ++i) {
    tree.Branch(sample.auxWeightBranches[i].c_str(), &auxValues[i]);
  }

  for (const auto& event : sample.events) {
    e = event;
    hasRecoil = event.hasRecoil;
    isMC = event.isMC;
    for (size_t i = 0; i < sample.auxWeightBranches.size(); ++i) {
      auto it = event.auxWeights.find(sample.auxWeightBranches[i]);
      auxValues[i] = it == event.auxWeights.end() ? 1.0 : it->second;
    }
    tree.Fill();
  }

  TTree params("params", "Sigma_CM skim metadata");
  double sigmaGen = sample.sigmaGen;
  int fdLeadRegionValue = sample.fdLeadRegionValue;
  int cdLeadRegionValue = sample.cdLeadRegionValue;
  params.Branch("sigma_gen", &sigmaGen);
  params.Branch("fdLeadRegionValue", &fdLeadRegionValue);
  params.Branch("cdLeadRegionValue", &cdLeadRegionValue);
  params.Branch("leadRegionFD", &fdLeadRegionValue);
  params.Branch("leadRegionCD", &cdLeadRegionValue);
  params.Fill();

  const std::string metadata =
      sample.metadataDump +
      "\nsigma_gen=" + std::to_string(sample.sigmaGen) +
      "\nfdLeadRegionValue=" + std::to_string(sample.fdLeadRegionValue) +
      "\ncdLeadRegionValue=" + std::to_string(sample.cdLeadRegionValue) + "\n";
  TNamed meta("sigmacm_skim_metadata", metadata.c_str());

  tree.Write();
  params.Write();
  if (!sample.gcfToyParams.empty()) {
    TTree gcfToyParams("gcfToyParams", "GCF toy randomized parameters");
    std::string branchName;
    double Cpp0, Cpn0, Cnn0, Cpn1, sigmaCM;
    double P00, P01, P10, P11, P20, P21, P30, P31;
    double TN, TNN;
    gcfToyParams.Branch("branchName", &branchName);
    gcfToyParams.Branch("Cpp0", &Cpp0);
    gcfToyParams.Branch("Cpn0", &Cpn0);
    gcfToyParams.Branch("Cnn0", &Cnn0);
    gcfToyParams.Branch("Cpn1", &Cpn1);
    gcfToyParams.Branch("sigmaCM", &sigmaCM);
    gcfToyParams.Branch("P00", &P00);
    gcfToyParams.Branch("P01", &P01);
    gcfToyParams.Branch("P10", &P10);
    gcfToyParams.Branch("P11", &P11);
    gcfToyParams.Branch("P20", &P20);
    gcfToyParams.Branch("P21", &P21);
    gcfToyParams.Branch("P30", &P30);
    gcfToyParams.Branch("P31", &P31);
    gcfToyParams.Branch("TN", &TN);
    gcfToyParams.Branch("TNN", &TNN);
    for (const auto& p : sample.gcfToyParams) {
      branchName = p.branchName;
      Cpp0 = p.Cpp0; Cpn0 = p.Cpn0; Cnn0 = p.Cnn0; Cpn1 = p.Cpn1;
      sigmaCM = p.sigmaCM;
      P00 = p.P00; P01 = p.P01; P10 = p.P10; P11 = p.P11;
      P20 = p.P20; P21 = p.P21; P30 = p.P30; P31 = p.P31;
      TN = p.TN; TNN = p.TNN;
      gcfToyParams.Fill();
    }
    gcfToyParams.Write();
  }
  meta.Write();
  file->Close();
}

}  // namespace sigmacm

#ifdef SIGMACM_WITH_HIPO

#include "Corrections.h"
#include "HipoChain.h"
#include "clas12ana.h"
#ifdef SIGMACM_WITH_GCF_TOYS
#include "reweighter.h"
#endif

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>

namespace sigmacm {
namespace {

constexpr double kNominalGcfSigmaCM = 0.15;
constexpr double kMassD = 1.8756;
constexpr double kMassElectron = 0.000511;
constexpr double kAtomicMassUnit = 0.9314941024;

double helium4Mass() {
  return 4.00260325415 * kAtomicMassUnit - 2.0 * kMassElectron;
}

int protonNumber(int a) {
  if (a == 48) return 20;
  if (a == 120) return 50;
  return a / 2;
}

int neutronNumber(int a) {
  if (a == 48) return 20;
  if (a == 120) return 70;
  return a / 2;
}

void setElectronP4LikeSkimEp(TLorentzVector& p4, clas12::region_part_ptr electron, bool isMC) {
  GetLorentzVector_ReconVector(p4, electron);
  if (!isMC) {
    SetLorentzVector_ThetaCorrection(p4, electron);
    SetLorentzVector_MomentumCorrection(p4, electron);
  }
  if (isMC) {
    SetLorentzVector_MomentumSimulationSmear(p4, electron);
  }
}

void setProtonP4LikeSkimEp(TLorentzVector& p4, clas12::region_part_ptr proton, bool isMC) {
  GetLorentzVector_ReconVector(p4, proton);
  if (!isMC) {
    SetLorentzVector_ThetaCorrection(p4, proton);
  }
  SetLorentzVector_EnergyLossCorrection(p4, proton);
  if (!isMC) {
    SetLorentzVector_MomentumCorrection(p4, proton);
  }
  if (isMC) {
    SetLorentzVector_MomentumSimulationSmear(p4, proton);
  }
}

bool fillHipoEvent(const std::unique_ptr<clas12::clas12reader>& c12,
                   clas12ana& clasAna,
                   bool isMC,
                   double beamEnergy,
                   Event& out) {
  auto* db = TDatabasePDG::Instance();
  const double massP = db->GetParticle(2212)->Mass();
  const double massN = db->GetParticle(2112)->Mass();

  TLorentzVector beam(0.0, 0.0, beamEnergy, beamEnergy);
  TLorentzVector deuteron(0.0, 0.0, 0.0, kMassD);
  TLorentzVector electron(0.0, 0.0, 0.0, kMassElectron);
  TLorentzVector leadP4(0.0, 0.0, 0.0, massP);
  TLorentzVector recoilP4(0.0, 0.0, 0.0, massP);

  clasAna.Run(c12);
  auto electrons = clasAna.getByPid(11);
  if (electrons.size() != 1) return false;

  setElectronP4LikeSkimEp(electron, electrons[0], isMC);
  TLorentzVector q = beam - electron;
  const double q2 = -q.M2();
  const double omega = q.E();
  if (omega <= 0.0) return false;
  const double xB = q2 / (2.0 * massP * omega);

  clasAna.getLeadRecoilSRC(beam, deuteron, electron);
  auto lead = clasAna.getLeadSRC();
  auto recoil = clasAna.getRecoilSRC();
  if (lead.size() != 1) return false;

  setProtonP4LikeSkimEp(leadP4, lead[0], isMC);
  TLorentzVector miss = q + deuteron - leadP4;
  const double mmiss2 = miss.M2();
  if (mmiss2 < 0.0) return false;
  TVector3 missNeg = -miss.Vect();

  TLorentzVector missLC = leadP4 - q;
  TVector3 qUnit = q.Vect().Unit();
  const double pmm = missLC.E() - missLC.Vect().Dot(qUnit);
  const double pmp = missLC.Vect().Perp(qUnit);
  const double denom = pmm * (2.0 * massN - pmm);
  if (denom <= 0.0) return false;
  const double kmiss2 = massN * massN * ((pmp * pmp + massN * massN) / denom) - massN * massN;
  if (kmiss2 < 0.0) return false;

  out.run = c12->runconfig()->getRun();
  out.event = c12->runconfig()->getEvent();
  out.Q2 = q2;
  out.xB = xB;
  out.mMiss = std::sqrt(mmiss2);
  out.kMiss = std::sqrt(kmiss2);
  out.pMiss = miss.P();
  out.pLead = leadP4.P();
  out.thetaLead = leadP4.Theta() * 180.0 / M_PI;
  out.leadRegion = static_cast<int>(lead[0]->getRegion());
  out.isMC = isMC;
  out.genWeight = isMC ? c12->mcevent()->getWeight() : 1.0;

  if (recoil.size() == 1) {
    setProtonP4LikeSkimEp(recoilP4, recoil[0], isMC);
    TVector3 vRec = recoilP4.Vect();
    TVector3 vRel = (missNeg - vRec) * 0.5;
    TVector3 vCM = missNeg + vRec;

    if (missNeg.Mag() > 0.0 && missNeg.Cross(q.Vect()).Mag() > 0.0) {
      TVector3 vz = missNeg.Unit();
      TVector3 vy = missNeg.Cross(q.Vect()).Unit();
      TVector3 vx = vz.Cross(vy).Unit();
      out.pcmX = vCM.Dot(vx);
      out.pcmY = vCM.Dot(vy);
      out.pcmZ = vCM.Dot(vz);
    }
    out.hasRecoil = true;
    out.pRec = vRec.Mag();
    out.thetaRec = recoilP4.Theta() * 180.0 / M_PI;
    out.pRel = vRel.Mag();
    out.pCM = vCM.Mag();
  }

  if (isMC) {
    auto mcInfo = c12->mcparts();
    if (!mcInfo || mcInfo->getRows() < 3) return false;
    TVector3 vBeam(0.0, 0.0, beamEnergy);
    TVector3 vElectron(mcInfo->getPx(0), mcInfo->getPy(0), mcInfo->getPz(0));
    TVector3 vLead(mcInfo->getPx(1), mcInfo->getPy(1), mcInfo->getPz(1));
    TVector3 vRec(mcInfo->getPx(2), mcInfo->getPy(2), mcInfo->getPz(2));
    TVector3 vQ = vBeam - vElectron;
    TVector3 vMiss = vLead - vQ;
    TVector3 vCM = vMiss + vRec;
    if (vMiss.Mag() <= 0.0 || vMiss.Cross(vQ).Mag() <= 0.0) return false;
    TVector3 tz = vMiss.Unit();
    TVector3 ty = vMiss.Cross(vQ).Unit();
    TVector3 tx = tz.Cross(ty).Unit();
    out.pcmXTruth = vCM.Dot(tx);
    out.pcmYTruth = vCM.Dot(ty);
    out.pcmZTruth = vCM.Dot(tz);
  }

  return true;
}

#ifdef SIGMACM_WITH_GCF_TOYS
GcfToyParams snapshotGcfToyParams(const std::string& branchName, const reweighter& w) {
  GcfToyParams p;
  p.branchName = branchName;
  p.Cpp0 = w.get_Cpp0();
  p.Cpn0 = w.get_Cpn0();
  p.Cnn0 = w.get_Cnn0();
  p.Cpn1 = w.get_Cpn1();
  p.sigmaCM = w.get_sigma_cm_fin();
  p.P00 = w.get_P(0, 0);
  p.P01 = w.get_P(0, 1);
  p.P10 = w.get_P(1, 0);
  p.P11 = w.get_P(1, 1);
  p.P20 = w.get_P(2, 0);
  p.P21 = w.get_P(2, 1);
  p.P30 = w.get_P(3, 0);
  p.P31 = w.get_P(3, 1);
  p.TN = w.get_TN();
  p.TNN = w.get_TNN();
  return p;
}
#endif

}  // namespace

Sample loadHipo(const std::vector<std::string>& paths, bool requireMC,
                const Config& cfg, const HipoLoadOptions& options) {
  if (paths.empty()) throw std::runtime_error("No hipo input files were provided");

  clas12root::HipoChain chain;
  for (const auto& path : paths) chain.Add(path);
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto configC12 = chain.GetC12Reader();
  (void)configC12;
  auto& c12 = chain.C12ref();

  clas12ana clasAna;
  thisRand->SetSeed(static_cast<UInt_t>(cfg.seed));
  Sample sample;
  sample.sigmaGen = 0.2;
  sample.fdLeadRegionValue = static_cast<int>(clas12::FD);
  sample.cdLeadRegionValue = static_cast<int>(clas12::CD);

#ifdef SIGMACM_WITH_GCF_TOYS
  std::vector<reweighter> gcfToyWeights;
  std::unique_ptr<reweighter> nominalGcfWeight;
  if (requireMC) {
    const int z = protonNumber(options.nucleusA);
    const int n = neutronNumber(options.nucleusA);
    nominalGcfWeight = std::make_unique<reweighter>(
        options.beamEnergy, z, n, kelly, const_cast<char*>("AV18"), kNominalGcfSigmaCM);
  }
  if (requireMC && options.nGcfToys > 0) {
    gcfToyWeights.reserve(static_cast<size_t>(options.nGcfToys));
    const int z = protonNumber(options.nucleusA);
    const int n = neutronNumber(options.nucleusA);
    for (int i = 0; i < options.nGcfToys; ++i) {
      reweighter w(options.beamEnergy, z, n, kelly, const_cast<char*>("AV18"), kNominalGcfSigmaCM);
      w.randomize_Config(false, false);
      std::ostringstream name;
      name << "w_gcf_toy_" << std::setw(3) << std::setfill('0') << i;
      const std::string branchName = name.str();
      sample.auxWeightBranches.push_back(branchName);
      sample.gcfToyParams.push_back(snapshotGcfToyParams(branchName, w));
      gcfToyWeights.push_back(w);
    }
  }
#else
  if (requireMC && options.nGcfToys > 0) {
    throw std::runtime_error("This Sigma_CM build does not include GCF toy reweighting support. "
                             "Build from the full repository so the reweighter target is available.");
  }
#endif

  std::ostringstream meta;
  std::ostringstream auxList;
  for (size_t i = 0; i < sample.auxWeightBranches.size(); ++i) {
    if (i) auxList << ",";
    auxList << sample.auxWeightBranches[i];
  }
  meta << "sigma_gen=0.2\n"
       << "leadRegionFD=" << sample.fdLeadRegionValue << "\n"
       << "leadRegionCD=" << sample.cdLeadRegionValue << "\n"
       << "lead_region_fd=" << sample.fdLeadRegionValue << "\n"
       << "lead_region_cd=" << sample.cdLeadRegionValue << "\n"
       << "beam_energy=" << options.beamEnergy << "\n"
       << "nucleus_A=" << options.nucleusA << "\n"
       << "nucleus_Z=" << protonNumber(options.nucleusA) << "\n"
      << "nucleus_N=" << neutronNumber(options.nucleusA) << "\n"
      << "target_mass_4He=" << helium4Mass() << "\n"
      << "smearing_rng_seed=" << cfg.seed << "\n"
      << "momentum_corrections=skim_ep_explicit_recon_data_corrections_mc_smearing\n"
      << "gcf_nominal_potential=AV18\n"
      << "gcf_nominal_weight_applied=" << (requireMC ? 1 : 0) << "\n"
      << "gcf_nominal_sigma_cm_factor_removed=" << (requireMC ? 1 : 0) << "\n"
      << "gcf_nominal_sigma_cm=" << kNominalGcfSigmaCM << "\n"
      << "gcf_toy_weight_mode=toy_over_nominal_ratio\n"
      << "gcf_toy_transparency_randomized=0\n"
      << "gcf_toy_sigma_cm_randomized=0\n"
      << "gcf_toy_params_tree=gcfToyParams\n"
      << "gcf_toy_branch_count=" << sample.auxWeightBranches.size() << "\n"
      << "gcf_toy_branches=" << auxList.str() << "\n"
      << "n_gcf_toys=" << (requireMC ? options.nGcfToys : 0) << "\n";
  sample.metadataDump = meta.str();

  long long scanned = 0;
  while (chain.Next()) {
    if (options.maxEvents >= 0 && scanned >= options.maxEvents) break;
    Event event;
    event.stableIndex = static_cast<std::uint64_t>(scanned);
    if (fillHipoEvent(c12, clasAna, requireMC, options.beamEnergy, event)) {
#ifdef SIGMACM_WITH_GCF_TOYS
      if (requireMC && nominalGcfWeight) {
        auto* mcInfo = c12->mcparts();
        const double nominal = nominalGcfWeight->get_weight_epp(mcInfo);
        if (!std::isfinite(nominal) || nominal == 0.0) {
          ++scanned;
          continue;
        }
        const double nominalSigma =
            gaussianRatio(event.pcmXTruth, kNominalGcfSigmaCM, sample.sigmaGen) *
            gaussianRatio(event.pcmYTruth, kNominalGcfSigmaCM, sample.sigmaGen) *
            gaussianRatio(event.pcmZTruth, kNominalGcfSigmaCM, sample.sigmaGen);
        if (!std::isfinite(nominalSigma) || nominalSigma == 0.0) {
          ++scanned;
          continue;
        }
        event.genWeight *= nominal / nominalSigma;
        for (size_t i = 0; i < gcfToyWeights.size(); ++i) {
          const double toy = gcfToyWeights[i].get_weight_epp(mcInfo);
          event.auxWeights[sample.auxWeightBranches[i]] =
              std::isfinite(toy) ? toy / nominal : 1.0;
        }
      }
#endif
      sample.events.push_back(event);
    }
    ++scanned;
    if (scanned % 100000 == 0) {
      std::cout << "  scanned " << scanned << " hipo events, kept "
                << sample.events.size() << " candidate events\n";
    }
  }
  std::cout << "Loaded " << sample.events.size() << " candidate events from "
            << scanned << " hipo events\n";
  return sample;
}

bool hipoSupportAvailable() {
  return true;
}

}  // namespace sigmacm

#else

namespace sigmacm {

Sample loadHipo(const std::vector<std::string>&, bool, const Config&, const HipoLoadOptions&) {
  throw std::runtime_error("This Sigma_CM build does not include hipo support. "
                           "Build from the full repository environment or use --from-skim.");
}

bool hipoSupportAvailable() {
  return false;
}

}  // namespace sigmacm

#endif
