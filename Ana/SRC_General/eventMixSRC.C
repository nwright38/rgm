// Build an event-mixed e'pp sample from a flat SRC skim.
//
// Usage:
//   root -l -b -q \
//     'Ana/SRC_General/eventMixSRC.C+("input.root","mixed.root","FD","all")'
//
// Only entries with pCM > 0 are used.  For N such entries, the output has
// N*(N-1) entries: the electron and leading proton come from entry i, and
// the recoil proton comes from every entry j != i.  Event weights are copied
// unchanged from entry i; no normalization is applied.
//
// leadDetector and recoilDetector independently accept:
//   "all" : no detector-angle selection
//   "FD"  : theta < 37 degrees
//   "CD"  : theta > 45 degrees

#include <ROOT/RConfig.hxx>
#include <TBranch.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TObject.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr double kProtonMass = 0.9382720813; // GeV
// Nuclear mass of 12C (atomic mass minus six electron masses), in GeV.
constexpr double kCarbon12Mass = 12.0 * 0.93149410242 - 6.0 * 0.00051099895;

struct Recoil {
  Long64_t entry;
  Float_t p;
  Float_t theta;
  Float_t phi;
  Float_t beta;
  Float_t tof;
};

enum class DetectorSelection { All, FD, CD, Invalid };

DetectorSelection parseDetector(const char *value)
{
  std::string text = value ? value : "";
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (text.empty() || text == "all") return DetectorSelection::All;
  if (text == "fd" || text == "forward") return DetectorSelection::FD;
  if (text == "cd" || text == "central") return DetectorSelection::CD;
  return DetectorSelection::Invalid;
}

const char *detectorName(DetectorSelection selection)
{
  switch (selection) {
    case DetectorSelection::All: return "all";
    case DetectorSelection::FD: return "FD (theta < 37 deg)";
    case DetectorSelection::CD: return "CD (theta > 45 deg)";
    default: return "invalid";
  }
}

bool passesDetector(Float_t theta, DetectorSelection selection)
{
  const double thetaDegrees = theta * 180.0 / M_PI;
  switch (selection) {
    case DetectorSelection::All: return true;
    case DetectorSelection::FD: return thetaDegrees < 37.0;
    case DetectorSelection::CD: return thetaDegrees > 45.0;
    default: return false;
  }
}

bool requireBranch(TTree *tree, const char *name)
{
  if (tree->GetBranch(name)) return true;
  std::cerr << "eventMixSRC: required branch '" << name << "' is missing\n";
  return false;
}

TVector3 sphericalVector(double magnitude, double theta, double phi)
{
  TVector3 result;
  result.SetMagThetaPhi(magnitude, theta, phi);
  return result;
}

Float_t angleOrSentinel(const TVector3 &a, const TVector3 &b)
{
  if (a.Mag2() == 0.0 || b.Mag2() == 0.0) return -9.0F;
  return static_cast<Float_t>(a.Angle(b));
}

Float_t mixedChiFrame(const TVector3 &pMiss, const TVector3 &q,
                      const TVector3 &pRec)
{
  if (pMiss.Mag2() == 0.0) return -9.0F;
  const TVector3 y = pMiss.Cross(q);
  if (y.Mag2() == 0.0) return -9.0F;
  const TVector3 zHat = pMiss.Unit();
  const TVector3 yHat = y.Unit();
  const TVector3 xHat = zHat.Cross(yHat).Unit();
  return static_cast<Float_t>(
      std::atan2(pRec.Dot(yHat), pRec.Dot(xHat)));
}

} // namespace

int eventMixSRC(const char *inputName,
                const char *outputName = "",
                const char *leadDetector = "all",
                const char *recoilDetector = "all",
                double targetMass = kCarbon12Mass)
{
  if (!inputName || !*inputName) {
    std::cerr << "Usage: eventMixSRC(input.root, output.root"
                 "[, leadDetector, recoilDetector, targetMassGeV])\n";
    return 1;
  }

  const DetectorSelection leadSelection = parseDetector(leadDetector);
  const DetectorSelection recoilSelection = parseDetector(recoilDetector);
  if (leadSelection == DetectorSelection::Invalid ||
      recoilSelection == DetectorSelection::Invalid) {
    std::cerr << "eventMixSRC: detector selectors must be 'all', 'FD', or "
                 "'CD'\n";
    return 1;
  }

  const std::string expandedInput = gSystem->ExpandPathName(inputName);
  std::string expandedOutput;
  if (outputName && *outputName) {
    expandedOutput = gSystem->ExpandPathName(outputName);
  } else {
    expandedOutput = expandedInput;
    const std::string suffix = ".root";
    const auto pos = expandedOutput.rfind(suffix);
    if (pos == std::string::npos) expandedOutput += "_mixed.root";
    else expandedOutput.insert(pos, "_mixed");
  }

  if (expandedInput == expandedOutput) {
    std::cerr << "eventMixSRC: input and output must be different files\n";
    return 1;
  }
  if (!(targetMass > 0.0) || !std::isfinite(targetMass)) {
    std::cerr << "eventMixSRC: targetMass must be a positive value in GeV\n";
    return 1;
  }

  TFile input(expandedInput.c_str(), "READ");
  if (input.IsZombie()) {
    std::cerr << "eventMixSRC: cannot open " << expandedInput << '\n';
    return 1;
  }
  auto *source = dynamic_cast<TTree *>(input.Get("srcTree"));
  if (!source) {
    std::cerr << "eventMixSRC: no TTree named 'srcTree' in input\n";
    return 1;
  }

  const char *required[] = {
      "pCM", "pMiss", "pMissTheta", "pMissPhi",
      "qP", "qTheta", "qPhi", "leadP", "leadTheta", "leadPhi",
      "recP", "recTheta", "recPhi", "recBeta", "recToF",
      "pRel", "pRelTheta", "pRelPhi", "pCMx", "pCMy", "pCMz",
      "chi_frame", "E2miss", "theta_PleadPrec", "theta_PmPrec",
      "theta_PrecQ", "omega"};
  for (const char *name : required)
    if (!requireBranch(source, name)) return 1;

  // Cache only the small recoil record.  This avoids rereading donor entries
  // N times during the N*(N-1) mixing loop.
  Float_t selectPCM = -9.0F;
  Float_t selectLeadTheta = -9.0F;
  Recoil scratch{};
  source->SetBranchStatus("*", 0);
  for (const char *name : {"pCM", "leadTheta", "recP", "recTheta",
                           "recPhi", "recBeta", "recToF"})
    source->SetBranchStatus(name, 1);
  source->SetBranchAddress("pCM", &selectPCM);
  source->SetBranchAddress("leadTheta", &selectLeadTheta);
  source->SetBranchAddress("recP", &scratch.p);
  source->SetBranchAddress("recTheta", &scratch.theta);
  source->SetBranchAddress("recPhi", &scratch.phi);
  source->SetBranchAddress("recBeta", &scratch.beta);
  source->SetBranchAddress("recToF", &scratch.tof);

  std::vector<Recoil> recoils;
  recoils.reserve(static_cast<std::size_t>(source->GetEntries() / 10));
  std::uint64_t eppCount = 0;
  for (Long64_t entry = 0; entry < source->GetEntries(); ++entry) {
    source->GetEntry(entry);
    if (selectPCM > 0.0F) {
      ++eppCount;
      if (!passesDetector(selectLeadTheta, leadSelection) ||
          !passesDetector(scratch.theta, recoilSelection))
        continue;
      scratch.entry = entry;
      recoils.push_back(scratch);
    }
  }

  const auto n = static_cast<std::uint64_t>(recoils.size());
  std::cout << "Input entries: " << source->GetEntries() << '\n'
            << "e'pp entries (pCM > 0): " << eppCount << '\n'
            << "Lead selection: " << detectorName(leadSelection) << '\n'
            << "Recoil selection: " << detectorName(recoilSelection) << '\n'
            << "Events in selected mixing pool (N): " << n << '\n';
  if (n < 2) {
    std::cerr << "eventMixSRC: need at least two entries after all "
                 "selections; found " << n << '\n';
    return 1;
  }
  if (n > std::numeric_limits<std::uint64_t>::max() / (n - 1)) {
    std::cerr << "eventMixSRC: mixed entry count overflows uint64_t\n";
    return 1;
  }
  const std::uint64_t expected = n * (n - 1);
  std::cout << "Mixed output entries N*(N-1): " << expected << '\n';

  source->ResetBranchAddresses();
  source->SetBranchStatus("*", 1);

  Float_t pMiss, pMissTheta, pMissPhi;
  Float_t qP, qTheta, qPhi;
  Float_t leadP, leadTheta, leadPhi;
  Float_t recP, recTheta, recPhi, recBeta, recToF;
  Float_t pRel, pRelTheta, pRelPhi;
  Float_t pCM, pCMx, pCMy, pCMz, chiFrame, e2miss;
  Float_t thetaLeadRec, thetaMissRec, thetaRecQ, omega;

  source->SetBranchAddress("pMiss", &pMiss);
  source->SetBranchAddress("pMissTheta", &pMissTheta);
  source->SetBranchAddress("pMissPhi", &pMissPhi);
  source->SetBranchAddress("qP", &qP);
  source->SetBranchAddress("qTheta", &qTheta);
  source->SetBranchAddress("qPhi", &qPhi);
  source->SetBranchAddress("leadP", &leadP);
  source->SetBranchAddress("leadTheta", &leadTheta);
  source->SetBranchAddress("leadPhi", &leadPhi);
  source->SetBranchAddress("recP", &recP);
  source->SetBranchAddress("recTheta", &recTheta);
  source->SetBranchAddress("recPhi", &recPhi);
  source->SetBranchAddress("recBeta", &recBeta);
  source->SetBranchAddress("recToF", &recToF);
  source->SetBranchAddress("pRel", &pRel);
  source->SetBranchAddress("pRelTheta", &pRelTheta);
  source->SetBranchAddress("pRelPhi", &pRelPhi);
  source->SetBranchAddress("pCM", &pCM);
  source->SetBranchAddress("pCMx", &pCMx);
  source->SetBranchAddress("pCMy", &pCMy);
  source->SetBranchAddress("pCMz", &pCMz);
  source->SetBranchAddress("chi_frame", &chiFrame);
  source->SetBranchAddress("E2miss", &e2miss);
  source->SetBranchAddress("theta_PleadPrec", &thetaLeadRec);
  source->SetBranchAddress("theta_PmPrec", &thetaMissRec);
  source->SetBranchAddress("theta_PrecQ", &thetaRecQ);
  source->SetBranchAddress("omega", &omega);

  TFile output(expandedOutput.c_str(), "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "eventMixSRC: cannot create " << expandedOutput << '\n';
    return 1;
  }
  output.SetCompressionSettings(input.GetCompressionSettings());
  auto *mixed = source->CloneTree(0);
  mixed->SetName(source->GetName());
  mixed->SetTitle(source->GetTitle());
  mixed->SetAutoFlush(-100000000); // about 100 MB per cluster

  const TLorentzVector target(0.0, 0.0, 0.0, targetMass);
  for (std::size_t i = 0; i < recoils.size(); ++i) {
    source->GetEntry(recoils[i].entry);

    // These quantities are entirely determined by the primary event.
    const TVector3 miss = sphericalVector(pMiss, pMissTheta, pMissPhi);
    const TVector3 q = sphericalVector(qP, qTheta, qPhi);
    const TVector3 lead = sphericalVector(leadP, leadTheta, leadPhi);

    TLorentzVector q4;
    q4.SetPxPyPzE(q.X(), q.Y(), q.Z(), omega);
    TLorentzVector lead4;
    lead4.SetVectM(lead, kProtonMass);

    const TVector3 zHat = miss.Unit();
    const TVector3 yRaw = miss.Cross(q);
    const bool validFrame = miss.Mag2() > 0.0 && yRaw.Mag2() > 0.0;
    const TVector3 yHat = validFrame ? yRaw.Unit() : TVector3();
    const TVector3 xHat = validFrame ? zHat.Cross(yHat).Unit() : TVector3();

    for (std::size_t j = 0; j < recoils.size(); ++j) {
      if (j == i) continue;
      const Recoil &donor = recoils[j];
      recP = donor.p;
      recTheta = donor.theta;
      recPhi = donor.phi;
      recBeta = donor.beta;
      recToF = donor.tof;

      const TVector3 rec = sphericalVector(recP, recTheta, recPhi);
      const TVector3 rel = 0.5 * (miss - rec);
      const TVector3 cm = miss + rec;

      pRel = static_cast<Float_t>(rel.Mag());
      pRelTheta = static_cast<Float_t>(rel.Theta());
      pRelPhi = static_cast<Float_t>(rel.Phi());
      pCM = static_cast<Float_t>(cm.Mag());
      pCMx = validFrame ? static_cast<Float_t>(cm.Dot(xHat)) : -9.0F;
      pCMy = validFrame ? static_cast<Float_t>(cm.Dot(yHat)) : -9.0F;
      pCMz = validFrame ? static_cast<Float_t>(cm.Dot(zHat)) : -9.0F;
      chiFrame = mixedChiFrame(miss, q, rec);
      thetaLeadRec = angleOrSentinel(lead, rec);
      thetaMissRec = angleOrSentinel(miss, rec);
      thetaRecQ = angleOrSentinel(rec, q);

      TLorentzVector rec4;
      rec4.SetVectM(rec, kProtonMass);
      const TLorentzVector residual = q4 + target - lead4 - rec4;
      const double tLead = lead4.E() - lead4.M();
      const double tRec = rec4.E() - rec4.M();
      const double tResidual = residual.E() - residual.M();
      e2miss = static_cast<Float_t>(omega - tLead - tRec - tResidual);

      mixed->Fill();
    }

    if ((i + 1) % 100 == 0 || i + 1 == recoils.size()) {
      std::cout << "\rMixed primary event " << (i + 1) << '/'
                << recoils.size() << std::flush;
    }
  }
  std::cout << '\n';

  output.cd();
  mixed->Write("", TObject::kOverwrite);
  output.Close();
  input.Close();

  std::cout << "Wrote " << expected << " entries to " << expandedOutput
            << " (no normalization applied)\n";
  return 0;
}
