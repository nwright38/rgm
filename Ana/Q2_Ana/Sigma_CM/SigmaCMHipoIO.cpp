#include "SigmaCMHipoIO.h"

#include <stdexcept>

#ifdef SIGMACM_WITH_HIPO

#include "Corrections.h"
#include "HipoChain.h"
#include "clas12ana.h"

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

namespace sigmacm {
namespace {

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

bool fillEvent(const std::unique_ptr<clas12::clas12reader>& c12,
               clas12ana& clasAna,
               bool isMC,
               double beamEnergy,
               Event& out) {
  auto* db = TDatabasePDG::Instance();
  const double massP = db->GetParticle(2212)->Mass();
  const double massN = db->GetParticle(2112)->Mass();

  TLorentzVector beam(0.0, 0.0, beamEnergy, beamEnergy);
  TLorentzVector nucleus(0.0, 0.0, 0.0, helium4Mass());
  TLorentzVector deuteron(0.0, 0.0, 0.0, kMassD);
  TLorentzVector electron(0.0, 0.0, 0.0, massP);
  TLorentzVector leadP4(0.0, 0.0, 0.0, massP);
  TLorentzVector recoilP4(0.0, 0.0, 0.0, massP);

  electron.SetE(kMassElectron);

  clasAna.Run(c12);
  auto electrons = clasAna.getByPid(11);
  if (electrons.size() != 1) return false;

  GetLorentzVector_Corrected(electron, electrons[0], isMC);
  TLorentzVector q = beam - electron;
  const double q2 = -q.M2();
  const double omega = q.E();
  if (omega <= 0.0) return false;
  const double xB = q2 / (2.0 * massP * omega);

  clasAna.getLeadRecoilSRC(beam, deuteron, electron);
  auto lead = clasAna.getLeadSRC();
  auto recoil = clasAna.getRecoilSRC();
  if (lead.size() != 1) return false;

  GetLorentzVector_Corrected(leadP4, lead[0], isMC);
  TLorentzVector miss = q + deuteron - leadP4;
  const double mmiss2 = miss.M2();
  if (mmiss2 < 0.0) return false;
  const double mmiss = std::sqrt(mmiss2);
  TVector3 missNeg = -miss.Vect();
  const double pmiss = miss.P();

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
  out.mMiss = mmiss;
  out.kMiss = std::sqrt(kmiss2);
  out.pMiss = pmiss;
  out.pLead = leadP4.P();
  out.thetaLead = leadP4.Theta() * 180.0 / M_PI;
  out.leadRegion = static_cast<int>(lead[0]->getRegion());
  out.hasRecoil = false;
  out.pRec = 0.0;
  out.thetaRec = 0.0;
  out.pcmX = 0.0;
  out.pcmY = 0.0;
  out.pcmZ = 0.0;
  out.pRel = 0.0;
  out.pCM = 0.0;
  out.isMC = isMC;
  out.genWeight = isMC ? c12->mcevent()->getWeight() : 1.0;

  if (recoil.size() == 1) {
    GetLorentzVector_Corrected(recoilP4, recoil[0], isMC);
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

}  // namespace

Sample loadHipo(const std::vector<std::string>& paths, bool requireMC,
                const Config&, const HipoLoadOptions& options) {
  if (paths.empty()) throw std::runtime_error("No hipo input files were provided");

  clas12root::HipoChain chain;
  for (const auto& path : paths) chain.Add(path);
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto configC12 = chain.GetC12Reader();
  (void)configC12;
  auto& c12 = chain.C12ref();

  clas12ana clasAna;
  Sample sample;
  sample.sigmaGen = 0.2;
  sample.fdLeadRegionValue = static_cast<int>(clas12::FD);
  sample.cdLeadRegionValue = static_cast<int>(clas12::CD);

  std::ostringstream meta;
  meta << "sigma_gen=0.2\n"
       << "leadRegionFD=" << sample.fdLeadRegionValue << "\n"
       << "leadRegionCD=" << sample.cdLeadRegionValue << "\n"
       << "lead_region_fd=" << sample.fdLeadRegionValue << "\n"
       << "lead_region_cd=" << sample.cdLeadRegionValue << "\n"
       << "beam_energy=" << options.beamEnergy << "\n"
       << "nucleus_A=" << options.nucleusA << "\n"
       << "nucleus_Z=" << protonNumber(options.nucleusA) << "\n"
       << "nucleus_N=" << neutronNumber(options.nucleusA) << "\n";
  sample.metadataDump = meta.str();

  long long scanned = 0;
  while (chain.Next()) {
    if (options.maxEvents >= 0 && scanned >= options.maxEvents) break;
    Event event;
    event.stableIndex = static_cast<std::uint64_t>(scanned);
    if (fillEvent(c12, clasAna, requireMC, options.beamEnergy, event)) {
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
