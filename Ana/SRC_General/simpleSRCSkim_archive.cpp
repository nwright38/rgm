#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>

#include "clas12reader.h"
#include "clas12ana.h"
#include "HipoChain.h"
#include "Corrections.h"
#include "reweighter.h"

using namespace std;
using namespace clas12;

const double mP = 0.938;
const double mU = 0.9314941024;
const double me = 0.000511;
const double m_4He = 4.00260325415 * mU - 2 * me;
const double m_12C = 12.0 * 0.93149410242 - 6.0 * 0.00051099895;
const double kLightConeEdgeEpsilon = 1e-6;
const double kLightConeRoundoffTolerance = 1e-9;

const int MAXP = 4;  // max number of proton candidates considered per event

struct LightConeBasis {
  TVector3 xHat;
  TVector3 yHat;
  TVector3 zHat;
  bool valid = false;
};

struct LightConeKinematics {
  bool validBasis = false;
  bool pairDefined = false;
  bool qualityBit = false;
  bool k2Clipped = false;
  double mBar = -9.;
  double alpha1 = -9.;
  double alpha2 = -9.;
  double alphaCM = -9.;
  double alphaRel = -9.;
  double qMinus = -9.;
  double k = -9.;
  double k2 = -9.;
  double kZ = -9.;
  TVector3 p1Perp;
  TVector3 p2Perp;
  TVector3 pCMPerp;
  TVector3 pRelPerp;
  double p1PerpX = -9.;
  double p1PerpY = -9.;
  double p1PerpMag = -9.;
  double p2PerpX = -9.;
  double p2PerpY = -9.;
  double p2PerpMag = -9.;
  double pCMPerpX = -9.;
  double pCMPerpY = -9.;
  double pCMPerpMag = -9.;
  double pRelPerpX = -9.;
  double pRelPerpY = -9.;
  double pRelPerpMag = -9.;
};

double square(double x) { return x * x; }

bool nearlyEqual(double a, double b, double tol = 1e-10)
{
  return std::abs(a - b) <= tol;
}

bool nearlyEqualVec(const TVector3 &a, const TVector3 &b, double tol = 1e-10)
{
  return (a - b).Mag() <= tol;
}

LightConeBasis makeLightConeBasis(const TVector3 &qVec)
{
  LightConeBasis basis;
  if(qVec.Mag2() <= 0.) return basis;

  basis.zHat = qVec.Unit();

  const TVector3 beamAxis(0., 0., 1.);
  TVector3 xSeed = beamAxis - basis.zHat * beamAxis.Dot(basis.zHat);
  if(xSeed.Mag2() <= 1e-12){
    const TVector3 xLab(1., 0., 0.);
    xSeed = xLab - basis.zHat * xLab.Dot(basis.zHat);
  }
  if(xSeed.Mag2() <= 1e-12){
    const TVector3 yLab(0., 1., 0.);
    xSeed = yLab - basis.zHat * yLab.Dot(basis.zHat);
  }
  if(xSeed.Mag2() <= 1e-12) return basis;

  basis.xHat = xSeed.Unit();
  basis.yHat = basis.zHat.Cross(basis.xHat);
  if(basis.yHat.Mag2() <= 1e-12) return basis;

  basis.yHat = basis.yHat.Unit();
  basis.xHat = basis.yHat.Cross(basis.zHat).Unit();
  basis.valid = true;
  return basis;
}

/*
 * Light-cone conventions used for the added GCF variables:
 *   z || q on an event-by-event basis,
 *   p^{\pm} = p^{0} \pm p^{3},
 *   alpha = p^{-} / m_bar,
 *   m_bar = m_A / A.
 * Nucleon 1 is the struck initial-state nucleon and nucleon 2 is the recoil nucleon.
 */
LightConeKinematics computeLightConeKinematics(const TVector3 &qVec, double nu,
                                               const TVector3 &pLead,
                                               const TVector3 &pRec,
                                               double mLead,
                                               double mRec,
                                               double targetMass,
                                               int massNumber)
{
  LightConeKinematics out;
  if(targetMass <= 0. || massNumber <= 0) return out;

  out.mBar = targetMass / static_cast<double>(massNumber);
  const LightConeBasis basis = makeLightConeBasis(qVec);
  if(!basis.valid) return out;

  out.validBasis = true;
  out.qualityBit = true;

  const double qMag = qVec.Mag();
  out.qMinus = nu - qMag;

  const double eLead = std::sqrt(pLead.Mag2() + mLead * mLead);
  const double eRec = std::sqrt(pRec.Mag2() + mRec * mRec);
  const double pLeadZ = pLead.Dot(basis.zHat);
  const double pRecZ = pRec.Dot(basis.zHat);

  // alpha_1 uses p_1^- = p_lead^- - q^- to avoid assigning an on-shell energy to the struck initial nucleon.
  out.alpha1 = (eLead - pLeadZ - out.qMinus) / out.mBar;
  out.alpha2 = (eRec - pRecZ) / out.mBar;

  out.p1Perp = pLead - basis.zHat * pLeadZ;
  out.p2Perp = pRec - basis.zHat * pRecZ;
  out.p1PerpX = out.p1Perp.Dot(basis.xHat);
  out.p1PerpY = out.p1Perp.Dot(basis.yHat);
  out.p1PerpMag = out.p1Perp.Mag();
  out.p2PerpX = out.p2Perp.Dot(basis.xHat);
  out.p2PerpY = out.p2Perp.Dot(basis.yHat);
  out.p2PerpMag = out.p2Perp.Mag();

  out.alphaCM = out.alpha1 + out.alpha2;
  if(std::abs(out.alphaCM) <= kLightConeRoundoffTolerance){
    out.qualityBit = false;
    return out;
  }

  out.pairDefined = true;
  out.pCMPerp = out.p1Perp + out.p2Perp;
  out.pCMPerpX = out.pCMPerp.Dot(basis.xHat);
  out.pCMPerpY = out.pCMPerp.Dot(basis.yHat);
  out.pCMPerpMag = out.pCMPerp.Mag();

  out.alphaRel = 2. * out.alpha2 / out.alphaCM;
  out.pRelPerp = out.p2Perp - (out.alpha2 / out.alphaCM) * out.pCMPerp;
  out.pRelPerpX = out.pRelPerp.Dot(basis.xHat);
  out.pRelPerpY = out.pRelPerp.Dot(basis.yHat);
  out.pRelPerpMag = out.pRelPerp.Mag();

  const double alphaFactor = out.alphaRel * (2. - out.alphaRel);
  if(out.alphaRel <= kLightConeEdgeEpsilon || out.alphaRel >= 2. - kLightConeEdgeEpsilon){
    out.qualityBit = false;
  }
  if(alphaFactor <= 0.){
    out.qualityBit = false;
    return out;
  }

  const double pairMass = 0.5 * (mLead + mRec);
  const double rawK2 = (pairMass * pairMass + out.pRelPerp.Mag2()) / alphaFactor - pairMass * pairMass;
  if(rawK2 < 0.){
    out.k2Clipped = true;
    if(rawK2 < -1e-7) out.qualityBit = false;
  }
  out.k2 = std::max(0., rawK2);
  out.k = std::sqrt(out.k2);
  out.kZ = (out.alphaRel - 1.) * std::sqrt(pairMass * pairMass + out.k2);
  return out;
}

TVector3 makeLeadMomentumFromInitialState(const TVector3 &p1Initial,
                                         double qMinus,
                                         double mLead,
                                         double mBar)
{
  const double alpha1 = (std::sqrt(p1Initial.Mag2() + mLead * mLead) - p1Initial.Z()) / mBar;
  const double pLeadMinus = alpha1 * mBar + qMinus;
  const double pPerp2 = square(p1Initial.X()) + square(p1Initial.Y());
  const double pLeadZ = (mLead * mLead + pPerp2 - pLeadMinus * pLeadMinus) / (2. * pLeadMinus);
  return TVector3(p1Initial.X(), p1Initial.Y(), pLeadZ);
}

void runLightConeSelfTests()
{
  const double protonMass = mP;
  const double targetMass = 2. * protonMass;
  const int massNumber = 2;

  {
    const TVector3 qVec(0., 0., 1.1);
    const double nu = 1.0;
    const TVector3 p1Initial(0., 0., 0.);
    const TVector3 pLead = makeLeadMomentumFromInitialState(
        p1Initial, nu - qVec.Mag(), protonMass, targetMass / massNumber);
    const TVector3 pRec(0., 0., 0.);
    const LightConeKinematics lc = computeLightConeKinematics(
        qVec, nu, pLead, pRec, protonMass, protonMass, targetMass, massNumber);

    assert(lc.validBasis);
    assert(lc.pairDefined);
    assert(nearlyEqual(lc.alpha1, 1.0, 1e-10));
    assert(nearlyEqual(lc.alpha2, 1.0, 1e-10));
    assert(nearlyEqual(lc.alphaRel, 1.0, 1e-10));
    assert(nearlyEqual(lc.pRelPerpMag, 0.0, 1e-10));
    assert(nearlyEqual(lc.k, 0.0, 1e-10));
  }

  {
    const TVector3 qVec(0., 0., 1.1);
    const double nu = 1.0;
    const TVector3 p1Initial(0.06, -0.03, 0.04);
    const TVector3 pLead = makeLeadMomentumFromInitialState(
        p1Initial, nu - qVec.Mag(), protonMass, targetMass / massNumber);
    const TVector3 pRec(-0.04, 0.02, -0.01);
    const LightConeKinematics lc = computeLightConeKinematics(
        qVec, nu, pLead, pRec, protonMass, protonMass, targetMass, massNumber);

    assert(lc.validBasis);
    assert(lc.pairDefined);

    const TVector3 pRelPerpAlt =
        (lc.alpha1 * lc.p2Perp - lc.alpha2 * lc.p1Perp) * (1. / lc.alphaCM);
    assert(nearlyEqualVec(lc.pRelPerp, pRelPerpAlt, 1e-10));

    const double alphaRelSwapped = 2. * lc.alpha1 / lc.alphaCM;
    const TVector3 pRelPerpSwapped =
        lc.p1Perp - (lc.alpha1 / lc.alphaCM) * lc.pCMPerp;
    const double k2Swapped =
        (protonMass * protonMass + pRelPerpSwapped.Mag2()) /
            (alphaRelSwapped * (2. - alphaRelSwapped)) -
        protonMass * protonMass;
    assert(nearlyEqual(alphaRelSwapped, 2. - lc.alphaRel, 1e-10));
    assert(nearlyEqualVec(pRelPerpSwapped, -lc.pRelPerp, 1e-10));
    assert(nearlyEqual(k2Swapped, lc.k2, 1e-10));

    const TVector3 pRelLab = (p1Initial - pRec) * 0.5;
    assert(std::abs(lc.k - pRelLab.Mag()) / pRelLab.Mag() < 0.05);
  }
}

void Usage()
{
  cerr << "Usage: ./code <MC=1,Data=0> <Ebeam(GeV)> <output.root> <A> <input.hipo> [more hipo...]\n";
}

bool inFD(int status){ return (abs(status) >= 2000 && abs(status) < 4000); }
bool inCD(int status){ return (abs(status) >= 4000 && abs(status) < 8000); }

double getChi(const TVector3 &pMiss, const TVector3 &q)
{
  // yHat is normal to the (pMiss,q) plane, while phiHat_m is normal to the
  // (pMiss,beam) plane.  chi is their dihedral angle about pMiss.
  const TVector3 zLab(0., 0., 1.);
  const TVector3 yHatV = pMiss.Cross(q);
  const TVector3 phiHatMV = zLab.Cross(pMiss);
  if(yHatV.Mag2() == 0. || phiHatMV.Mag2() == 0.) return -9.;

  double cosChi = yHatV.Unit().Dot(phiHatMV.Unit());
  if(cosChi > 1.) cosChi = 1.;
  if(cosChi < -1.) cosChi = -1.;
  return acos(cosChi);
}

double getChiFrame(const TVector3 &pMiss, const TVector3 &q,
                   const TVector3 &pRec)
{
  // Signed azimuth of the recoil about pMiss, measured from the
  // (q,pMiss) plane.  These are the same axes used for pCMx and pCMy.
  if(pMiss.Mag2() == 0.) return -9.;

  const TVector3 zHat = pMiss.Unit();
  const TVector3 yHatV = pMiss.Cross(q);
  if(yHatV.Mag2() == 0.) return -9.;

  const TVector3 yHat = yHatV.Unit();
  const TVector3 xHat = zHat.Cross(yHat).Unit();
  return atan2(pRec.Dot(yHat), pRec.Dot(xHat));
}

double getPhiTrentoZeroCM(const TVector3 &pLead, const TVector3 &pMiss,
                          const TVector3 &pRec)
{
  // Infer the third nucleon from zero-CM: p3 = -(pRec + pMiss).
  if(pLead.Mag2() == 0. || pRec.Mag2() == 0. || pMiss.Mag2() == 0.) return -9.;

  const TVector3 p3 = -(pRec + pMiss);
  const TVector3 n1 = pLead.Cross(pRec);
  const TVector3 n2 = pLead.Cross(p3);
  if(n1.Mag2() == 0. || n2.Mag2() == 0.) return -9.;

  double cosPhi = n1.Unit().Dot(n2.Unit());
  if(cosPhi > 1.) cosPhi = 1.;
  if(cosPhi < -1.) cosPhi = -1.;
  return acos(cosPhi);
}

void setElectronP4(TLorentzVector &p4, clas12::region_part_ptr electron, bool isMC)
{
  GetLorentzVector_ReconVector(p4, electron);
  if(!isMC){
    SetLorentzVector_ThetaCorrection(p4, electron);
    SetLorentzVector_MomentumCorrection(p4, electron);
  }
  if(isMC){
    SetLorentzVector_MomentumSimulationSmear(p4, electron);
  }
}

void setProtonP4(TLorentzVector &p4, clas12::region_part_ptr proton, bool isMC)
{
  GetLorentzVector_ReconVector(p4, proton);
  if(!isMC){
    SetLorentzVector_ThetaCorrection(p4, proton);
  }
  SetLorentzVector_EnergyLossCorrection(p4, proton);
  if(!isMC){
    SetLorentzVector_MomentumCorrection(p4, proton);
  }
  if(isMC){
    SetLorentzVector_MomentumSimulationSmear(p4, proton);
  }
}

int main(int argc, char **argv)
{
  if(argc < 5){ Usage(); return -1; }

  bool isMC = (atoi(argv[1]) == 1);
  double Ebeam = atof(argv[2]);
  int nucleus_A = atoi(argv[4]);
  double targetMass = m_4He;
  int targetMassNumber = 4;
  int targetZ = 2;
  int targetN = 2;
  double targetSigmaCM = 0.10;
  const char *targetLabel = "He4";

  if(nucleus_A == 12) {
    targetMass = m_12C;
    targetMassNumber = 12;
    targetZ = 6;
    targetN = 6;
    targetSigmaCM = 0.15;
    targetLabel = "C12";
  } else {
    targetMass = m_4He;
    targetMassNumber = 4;
    targetZ = 2;
    targetN = 2;
    targetSigmaCM = 0.10;
    targetLabel = "He4";
  }

  TFile *outFile = new TFile(argv[3], "RECREATE");
  TTree *srcTree = new TTree("srcTree", "SRC Kinematics");

  // ---- output branches ----
  Float_t b_weight;
  Float_t b_weight_ep;
  Float_t b_weight_epp;

  Float_t b_xB, b_Q2, b_omega;

  // electron
  Float_t b_eP, b_eTheta, b_ePhi;

  // q vector (3-momentum transfer)
  Float_t b_qP, b_qTheta, b_qPhi;

  Int_t   b_nProtons;

  // lead proton (first candidate passing all SRC lead cuts)
  Float_t b_leadP, b_leadTheta, b_leadPhi;
  Float_t b_leadBeta, b_leadToF;
  Float_t b_leadVz;
  Float_t b_pMiss, b_pMissTheta, b_pMissPhi;
  Float_t b_mMiss;
  Float_t b_kMiss;
  Float_t b_EMiss, b_E0miss, b_E1miss;
  Float_t b_theta_PmQ;    // angle between pMiss and q
  Float_t b_theta_PmPlead; // angle between pMiss and lead proton momentum
  Float_t b_theta_PleadQ; // angle between lead proton momentum and q
  Float_t b_chi;          // dihedral angle about pMiss
  Bool_t  b_goodLead;     // lead passes all SRC lead cuts

  // pRel and pCM (filled only when a recoil partner exists)
  Float_t b_pRel, b_pRelTheta, b_pRelPhi;
  Float_t b_pRelx, b_pRely, b_pRelz;
  Float_t b_pLeadPlusRecOver2, b_pLeadMinusRec;
  Float_t b_pCM, b_pCMx, b_pCMy, b_pCMz;
  Float_t b_pcmx_lab, b_pcmy_lab, b_pcmz_lab;
  Float_t b_chi_frame;
  Float_t b_E2miss;
  Float_t b_alpha_1, b_alpha_2, b_alpha_CM, b_alpha_rel;
  Float_t b_p1_perp_x, b_p1_perp_y, b_p1_perp_mag;
  Float_t b_p2_perp_x, b_p2_perp_y, b_p2_perp_mag;
  Float_t b_pCM_perp_x, b_pCM_perp_y, b_pCM_perp_mag;
  Float_t b_prel_perp_x, b_prel_perp_y, b_prel_perp_mag;
  Float_t b_k, b_k2, b_k_z, b_m_bar;
  Bool_t  b_lc_quality;

  // event-level summary
  Int_t  b_nGoodLeads;
  Bool_t b_singleGoodLead;      // exactly one proton passes all SRC lead cuts

  // recoil kinematics (filled when a recoil is found)
  Float_t b_recP, b_recTheta, b_recPhi;
  Float_t b_recBeta, b_recToF;
  Float_t b_theta_PleadPrec; // angle between lead and recoil momentum
  Float_t b_theta_PmPrec; // angle between pMiss and recoil momentum
  Float_t b_theta_PrecQ;  // angle between recoil momentum and q
  Float_t b_phiTrento;

  // truth-level MC quantities, named to mirror reco branches
  Float_t b_xB_truth, b_Q2_truth, b_omega_truth;
  Float_t b_eP_truth, b_eTheta_truth, b_ePhi_truth;
  Float_t b_qP_truth, b_qTheta_truth, b_qPhi_truth;
  Float_t b_leadP_truth, b_leadTheta_truth, b_leadPhi_truth;
  Float_t b_pMiss_truth, b_pMissTheta_truth, b_pMissPhi_truth;
  Float_t b_mMiss_truth;
  Float_t b_kMiss_truth;
  Float_t b_EMiss_truth, b_E0miss_truth, b_E1miss_truth;
  Float_t b_theta_PmQ_truth;
  Float_t b_theta_PmPlead_truth;
  Float_t b_theta_PleadQ_truth;
  Float_t b_chi_truth;
  Float_t b_pRel_truth, b_pRelTheta_truth, b_pRelPhi_truth;
  Float_t b_pRelx_truth, b_pRely_truth, b_pRelz_truth;
  Float_t b_pLeadPlusRecOver2_truth, b_pLeadMinusRec_truth;
  Float_t b_pCM_truth, b_pCMx_truth, b_pCMy_truth, b_pCMz_truth;
  Float_t b_pcmx_lab_truth, b_pcmy_lab_truth, b_pcmz_lab_truth;
  Float_t b_chi_frame_truth;
  Float_t b_E2miss_truth;
  Float_t b_alpha_1_truth, b_alpha_2_truth, b_alpha_CM_truth, b_alpha_rel_truth;
  Float_t b_p1_perp_x_truth, b_p1_perp_y_truth, b_p1_perp_mag_truth;
  Float_t b_p2_perp_x_truth, b_p2_perp_y_truth, b_p2_perp_mag_truth;
  Float_t b_pCM_perp_x_truth, b_pCM_perp_y_truth, b_pCM_perp_mag_truth;
  Float_t b_prel_perp_x_truth, b_prel_perp_y_truth, b_prel_perp_mag_truth;
  Float_t b_k_truth, b_k2_truth, b_k_z_truth, b_m_bar_truth;
  Bool_t  b_lc_quality_truth;
  Float_t b_recP_truth, b_recTheta_truth, b_recPhi_truth;
  Float_t b_theta_PmPrec_truth;
  Float_t b_theta_PrecQ_truth;
  Float_t b_phiTrento_truth;

  srcTree->Branch("weight",      &b_weight,      "weight/F");
  srcTree->Branch("weight_ep",   &b_weight_ep,   "weight_ep/F");
  srcTree->Branch("weight_epp",  &b_weight_epp,  "weight_epp/F");

  srcTree->Branch("xB",          &b_xB,          "xB/F");
  srcTree->Branch("Q2",          &b_Q2,          "Q2/F");
  srcTree->Branch("omega",       &b_omega,       "omega/F");

  srcTree->Branch("eP",          &b_eP,          "eP/F");
  srcTree->Branch("eTheta",      &b_eTheta,      "eTheta/F");
  srcTree->Branch("ePhi",        &b_ePhi,        "ePhi/F");

  srcTree->Branch("qP",          &b_qP,          "qP/F");
  srcTree->Branch("qTheta",      &b_qTheta,      "qTheta/F");
  srcTree->Branch("qPhi",        &b_qPhi,        "qPhi/F");

  srcTree->Branch("nProtons",    &b_nProtons,    "nProtons/I");

  srcTree->Branch("leadP",       &b_leadP,       "leadP/F");
  srcTree->Branch("leadTheta",   &b_leadTheta,   "leadTheta/F");
  srcTree->Branch("leadPhi",     &b_leadPhi,     "leadPhi/F");
  srcTree->Branch("leadBeta",    &b_leadBeta,    "leadBeta/F");
  srcTree->Branch("leadToF",     &b_leadToF,     "leadToF/F");
  srcTree->Branch("leadVz",      &b_leadVz,      "leadVz/F");

  srcTree->Branch("pMiss",       &b_pMiss,       "pMiss/F");
  srcTree->Branch("pMissTheta",  &b_pMissTheta,  "pMissTheta/F");
  srcTree->Branch("pMissPhi",    &b_pMissPhi,    "pMissPhi/F");

  srcTree->Branch("mMiss",       &b_mMiss,       "mMiss/F");
  srcTree->Branch("kMiss",       &b_kMiss,       "kMiss/F");
  srcTree->Branch("EMiss",       &b_EMiss,       "EMiss/F");
  srcTree->Branch("E0miss",      &b_E0miss,      "E0miss/F");
  srcTree->Branch("E1miss",      &b_E1miss,      "E1miss/F");
  srcTree->Branch("theta_PmQ",   &b_theta_PmQ,   "theta_PmQ/F");
  srcTree->Branch("theta_PmPlead",&b_theta_PmPlead,"theta_PmPlead/F");
  srcTree->Branch("theta_PleadQ",&b_theta_PleadQ,"theta_PleadQ/F");
  srcTree->Branch("chi",         &b_chi,         "chi/F");
  srcTree->Branch("goodLead",    &b_goodLead,    "goodLead/O");

  srcTree->Branch("pRel",        &b_pRel,        "pRel/F");
  srcTree->Branch("pRelTheta",   &b_pRelTheta,   "pRelTheta/F");
  srcTree->Branch("pRelPhi",     &b_pRelPhi,     "pRelPhi/F");
  srcTree->Branch("pRelx",       &b_pRelx,       "pRelx/F");
  srcTree->Branch("pRely",       &b_pRely,       "pRely/F");
  srcTree->Branch("pRelz",       &b_pRelz,       "pRelz/F");
  srcTree->Branch("pLeadPlusRecOver2", &b_pLeadPlusRecOver2,
                  "pLeadPlusRecOver2/F");
  srcTree->Branch("pLeadMinusRec", &b_pLeadMinusRec,
                  "pLeadMinusRec/F");
  srcTree->Branch("pCM",         &b_pCM,         "pCM/F");
  srcTree->Branch("pCMx",        &b_pCMx,        "pCMx/F");
  srcTree->Branch("pCMy",        &b_pCMy,        "pCMy/F");
  srcTree->Branch("pCMz",        &b_pCMz,        "pCMz/F");
  srcTree->Branch("pcmx_lab",    &b_pcmx_lab,    "pcmx_lab/F");
  srcTree->Branch("pcmy_lab",    &b_pcmy_lab,    "pcmy_lab/F");
  srcTree->Branch("pcmz_lab",    &b_pcmz_lab,    "pcmz_lab/F");
  srcTree->Branch("chi_frame",   &b_chi_frame,   "chi_frame/F");
  srcTree->Branch("E2miss",      &b_E2miss,      "E2miss/F");
  srcTree->Branch("alpha_1",     &b_alpha_1,     "alpha_1/F");
  srcTree->Branch("alpha_2",     &b_alpha_2,     "alpha_2/F");
  srcTree->Branch("alpha_CM",    &b_alpha_CM,    "alpha_CM/F");
  srcTree->Branch("alpha_rel",   &b_alpha_rel,   "alpha_rel/F");
  srcTree->Branch("p1_perp_x",   &b_p1_perp_x,   "p1_perp_x/F");
  srcTree->Branch("p1_perp_y",   &b_p1_perp_y,   "p1_perp_y/F");
  srcTree->Branch("p1_perp_mag", &b_p1_perp_mag, "p1_perp_mag/F");
  srcTree->Branch("p2_perp_x",   &b_p2_perp_x,   "p2_perp_x/F");
  srcTree->Branch("p2_perp_y",   &b_p2_perp_y,   "p2_perp_y/F");
  srcTree->Branch("p2_perp_mag", &b_p2_perp_mag, "p2_perp_mag/F");
  srcTree->Branch("pCM_perp_x",  &b_pCM_perp_x,  "pCM_perp_x/F");
  srcTree->Branch("pCM_perp_y",  &b_pCM_perp_y,  "pCM_perp_y/F");
  srcTree->Branch("pCM_perp_mag",&b_pCM_perp_mag,"pCM_perp_mag/F");
  srcTree->Branch("prel_perp_x", &b_prel_perp_x, "prel_perp_x/F");
  srcTree->Branch("prel_perp_y", &b_prel_perp_y, "prel_perp_y/F");
  srcTree->Branch("prel_perp_mag", &b_prel_perp_mag, "prel_perp_mag/F");
  srcTree->Branch("k",           &b_k,           "k/F");
  srcTree->Branch("k2",          &b_k2,          "k2/F");
  srcTree->Branch("k_z",         &b_k_z,         "k_z/F");
  srcTree->Branch("m_bar",       &b_m_bar,       "m_bar/F");
  srcTree->Branch("lc_quality",  &b_lc_quality,  "lc_quality/O");

  srcTree->Branch("nGoodLeads",     &b_nGoodLeads,     "nGoodLeads/I");
  srcTree->Branch("singleGoodLead", &b_singleGoodLead, "singleGoodLead/O");

  srcTree->Branch("recP",        &b_recP,        "recP/F");
  srcTree->Branch("recTheta",    &b_recTheta,    "recTheta/F");
  srcTree->Branch("recPhi",      &b_recPhi,      "recPhi/F");
  srcTree->Branch("recBeta",     &b_recBeta,     "recBeta/F");
  srcTree->Branch("recToF",      &b_recToF,      "recToF/F");
  srcTree->Branch("theta_PleadPrec",&b_theta_PleadPrec,"theta_PleadPrec/F");
  srcTree->Branch("theta_PmPrec",&b_theta_PmPrec,"theta_PmPrec/F");
  srcTree->Branch("theta_PrecQ", &b_theta_PrecQ, "theta_PrecQ/F");
  srcTree->Branch("phiTrento",   &b_phiTrento,   "phiTrento/F");

  srcTree->Branch("xB_truth",          &b_xB_truth,          "xB_truth/F");
  srcTree->Branch("Q2_truth",          &b_Q2_truth,          "Q2_truth/F");
  srcTree->Branch("omega_truth",       &b_omega_truth,       "omega_truth/F");

  srcTree->Branch("eP_truth",          &b_eP_truth,          "eP_truth/F");
  srcTree->Branch("eTheta_truth",      &b_eTheta_truth,      "eTheta_truth/F");
  srcTree->Branch("ePhi_truth",        &b_ePhi_truth,        "ePhi_truth/F");

  srcTree->Branch("qP_truth",          &b_qP_truth,          "qP_truth/F");
  srcTree->Branch("qTheta_truth",      &b_qTheta_truth,      "qTheta_truth/F");
  srcTree->Branch("qPhi_truth",        &b_qPhi_truth,        "qPhi_truth/F");

  srcTree->Branch("leadP_truth",       &b_leadP_truth,       "leadP_truth/F");
  srcTree->Branch("leadTheta_truth",   &b_leadTheta_truth,   "leadTheta_truth/F");
  srcTree->Branch("leadPhi_truth",     &b_leadPhi_truth,     "leadPhi_truth/F");

  srcTree->Branch("pMiss_truth",       &b_pMiss_truth,       "pMiss_truth/F");
  srcTree->Branch("pMissTheta_truth",  &b_pMissTheta_truth,  "pMissTheta_truth/F");
  srcTree->Branch("pMissPhi_truth",    &b_pMissPhi_truth,    "pMissPhi_truth/F");

  srcTree->Branch("mMiss_truth",       &b_mMiss_truth,       "mMiss_truth/F");
  srcTree->Branch("kMiss_truth",       &b_kMiss_truth,       "kMiss_truth/F");
  srcTree->Branch("EMiss_truth",       &b_EMiss_truth,       "EMiss_truth/F");
  srcTree->Branch("E0miss_truth",      &b_E0miss_truth,      "E0miss_truth/F");
  srcTree->Branch("E1miss_truth",      &b_E1miss_truth,      "E1miss_truth/F");
  srcTree->Branch("theta_PmQ_truth",   &b_theta_PmQ_truth,   "theta_PmQ_truth/F");
  srcTree->Branch("theta_PmPlead_truth",&b_theta_PmPlead_truth,"theta_PmPlead_truth/F");
  srcTree->Branch("theta_PleadQ_truth",&b_theta_PleadQ_truth,"theta_PleadQ_truth/F");
  srcTree->Branch("chi_truth",         &b_chi_truth,         "chi_truth/F");

  srcTree->Branch("pRel_truth",        &b_pRel_truth,        "pRel_truth/F");
  srcTree->Branch("pRelTheta_truth",   &b_pRelTheta_truth,   "pRelTheta_truth/F");
  srcTree->Branch("pRelPhi_truth",     &b_pRelPhi_truth,     "pRelPhi_truth/F");
  srcTree->Branch("pRelx_truth",       &b_pRelx_truth,       "pRelx_truth/F");
  srcTree->Branch("pRely_truth",       &b_pRely_truth,       "pRely_truth/F");
  srcTree->Branch("pRelz_truth",       &b_pRelz_truth,       "pRelz_truth/F");
  srcTree->Branch("pLeadPlusRecOver2_truth", &b_pLeadPlusRecOver2_truth,
                  "pLeadPlusRecOver2_truth/F");
  srcTree->Branch("pLeadMinusRec_truth", &b_pLeadMinusRec_truth,
                  "pLeadMinusRec_truth/F");
  srcTree->Branch("pCM_truth",         &b_pCM_truth,         "pCM_truth/F");
  srcTree->Branch("pCMx_truth",        &b_pCMx_truth,        "pCMx_truth/F");
  srcTree->Branch("pCMy_truth",        &b_pCMy_truth,        "pCMy_truth/F");
  srcTree->Branch("pCMz_truth",        &b_pCMz_truth,        "pCMz_truth/F");
  srcTree->Branch("pcmx_lab_truth",    &b_pcmx_lab_truth,    "pcmx_lab_truth/F");
  srcTree->Branch("pcmy_lab_truth",    &b_pcmy_lab_truth,    "pcmy_lab_truth/F");
  srcTree->Branch("pcmz_lab_truth",    &b_pcmz_lab_truth,    "pcmz_lab_truth/F");
  srcTree->Branch("chi_frame_truth",   &b_chi_frame_truth,   "chi_frame_truth/F");
  srcTree->Branch("E2miss_truth",      &b_E2miss_truth,      "E2miss_truth/F");
  srcTree->Branch("alpha_1_truth",     &b_alpha_1_truth,     "alpha_1_truth/F");
  srcTree->Branch("alpha_2_truth",     &b_alpha_2_truth,     "alpha_2_truth/F");
  srcTree->Branch("alpha_CM_truth",    &b_alpha_CM_truth,    "alpha_CM_truth/F");
  srcTree->Branch("alpha_rel_truth",   &b_alpha_rel_truth,   "alpha_rel_truth/F");
  srcTree->Branch("p1_perp_x_truth",   &b_p1_perp_x_truth,   "p1_perp_x_truth/F");
  srcTree->Branch("p1_perp_y_truth",   &b_p1_perp_y_truth,   "p1_perp_y_truth/F");
  srcTree->Branch("p1_perp_mag_truth", &b_p1_perp_mag_truth, "p1_perp_mag_truth/F");
  srcTree->Branch("p2_perp_x_truth",   &b_p2_perp_x_truth,   "p2_perp_x_truth/F");
  srcTree->Branch("p2_perp_y_truth",   &b_p2_perp_y_truth,   "p2_perp_y_truth/F");
  srcTree->Branch("p2_perp_mag_truth", &b_p2_perp_mag_truth, "p2_perp_mag_truth/F");
  srcTree->Branch("pCM_perp_x_truth",  &b_pCM_perp_x_truth,  "pCM_perp_x_truth/F");
  srcTree->Branch("pCM_perp_y_truth",  &b_pCM_perp_y_truth,  "pCM_perp_y_truth/F");
  srcTree->Branch("pCM_perp_mag_truth",&b_pCM_perp_mag_truth,"pCM_perp_mag_truth/F");
  srcTree->Branch("prel_perp_x_truth", &b_prel_perp_x_truth, "prel_perp_x_truth/F");
  srcTree->Branch("prel_perp_y_truth", &b_prel_perp_y_truth, "prel_perp_y_truth/F");
  srcTree->Branch("prel_perp_mag_truth", &b_prel_perp_mag_truth, "prel_perp_mag_truth/F");
  srcTree->Branch("k_truth",           &b_k_truth,           "k_truth/F");
  srcTree->Branch("k2_truth",          &b_k2_truth,          "k2_truth/F");
  srcTree->Branch("k_z_truth",         &b_k_z_truth,         "k_z_truth/F");
  srcTree->Branch("m_bar_truth",       &b_m_bar_truth,       "m_bar_truth/F");
  srcTree->Branch("lc_quality_truth",  &b_lc_quality_truth,  "lc_quality_truth/O");

  srcTree->Branch("recP_truth",        &b_recP_truth,        "recP_truth/F");
  srcTree->Branch("recTheta_truth",    &b_recTheta_truth,    "recTheta_truth/F");
  srcTree->Branch("recPhi_truth",      &b_recPhi_truth,      "recPhi_truth/F");
  srcTree->Branch("theta_PmPrec_truth",&b_theta_PmPrec_truth,"theta_PmPrec_truth/F");
  srcTree->Branch("theta_PrecQ_truth", &b_theta_PrecQ_truth, "theta_PrecQ_truth/F");
  srcTree->Branch("phiTrento_truth",   &b_phiTrento_truth,   "phiTrento_truth/F");

  TH1D *h_alpha_1 = new TH1D("h_alpha_1", "alpha_{1};alpha_{1};Events", 80, 0., 2.);
  TH1D *h_alpha_2 = new TH1D("h_alpha_2", "alpha_{2};alpha_{2};Events", 80, 0., 2.);
  TH1D *h_alpha_CM = new TH1D("h_alpha_CM", "alpha_{CM};alpha_{CM};Events", 80, 0., 4.);
  TH1D *h_alpha_spectator = new TH1D("h_alpha_spectator",
                                     "A-alpha_{CM};A-alpha_{CM};Events", 80, 0., targetMassNumber);
  TH1D *h_k2 = new TH1D("h_k2", "k^{2};k^{2} [GeV^{2}];Events", 80, 0., 2.0);
  TH2D *h_k_vs_pRel = new TH2D("h_k_vs_pRel",
                               "k vs |p_{rel}^{lab}|;|p_{rel}^{lab}| [GeV/c];k [GeV/c]",
                               80, 0., 1.2, 80, 0., 1.2);

  // ---- chain setup ----
  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout << "Input file " << argv[k] << endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  clas12ana clasAna;

  // fixed 4-vectors
  TLorentzVector targP4(0., 0., 0., targetMass);
  TLorentzVector nucleusP4(0., 0., 0., targetMass);
  TLorentzVector beamP4(0., 0., Ebeam, Ebeam);

  TLorentzVector eP4(0., 0., 0., me);
  TLorentzVector leadP4(0., 0., 0., mP);
  TLorentzVector missP4;

  TVector3 p_beam(0., 0., Ebeam);

  runLightConeSelfTests();

  int counter = 0;

  char av18[] = "AV18";
  reweighter newWeight(Ebeam, targetZ, targetN, kelly, av18, targetSigmaCM);
  cout << "Target A = " << nucleus_A << " (" << targetLabel
       << ", sigmaCM = " << targetSigmaCM << ")" << endl;
  int ctr = 0;
  while(chain.Next())
  {
    if(counter % 100000 == 0)
      cout << "Processing event " << counter << "\t" << ctr << " saved." << endl;
    counter++;

    if(ctr > 1000000) break;  // for testing

    // ---- defaults ----
    b_weight = 1.f;
    b_weight_ep = 1.f;
    b_weight_epp = 1.f;
    b_xB = -9.f;  b_Q2 = -9.f;  b_omega = -9.f;
    b_eP = -9.f;  b_eTheta = -9.f;  b_ePhi = -9.f;
    b_qP = -9.f;  b_qTheta = -9.f;  b_qPhi = -9.f;

    b_nProtons    = 0;
    b_nGoodLeads  = 0;
    b_singleGoodLead = false;
    b_recP       = -9.f;  b_recTheta = -9.f;  b_recPhi = -9.f;
    b_recBeta    = -9.f;  b_recToF = -9.f;

    b_leadP = -9.f;  b_leadTheta = -9.f;  b_leadPhi = -9.f;
    b_leadBeta = -9.f;  b_leadToF = -9.f;  b_leadVz = -99.f;
    b_pMiss = -9.f;  b_pMissTheta = -9.f;  b_pMissPhi = -9.f;
    b_pRel = -9.f;   b_pRelTheta = -9.f;   b_pRelPhi = -9.f;
    b_pRelx = -9.f;  b_pRely = -9.f;       b_pRelz = -9.f;
    b_pLeadPlusRecOver2 = -9.f;
    b_pLeadMinusRec = -9.f;
    b_pCM = -9.f;    b_pCMx = -9.f;        b_pCMy = -9.f;    b_pCMz = -9.f;
    b_pcmx_lab = -9.f; b_pcmy_lab = -9.f;  b_pcmz_lab = -9.f;
    b_chi_frame = -9.f;
    b_E2miss = -9.f;
    b_alpha_1 = -9.f; b_alpha_2 = -9.f; b_alpha_CM = -9.f; b_alpha_rel = -9.f;
    b_p1_perp_x = -9.f; b_p1_perp_y = -9.f; b_p1_perp_mag = -9.f;
    b_p2_perp_x = -9.f; b_p2_perp_y = -9.f; b_p2_perp_mag = -9.f;
    b_pCM_perp_x = -9.f; b_pCM_perp_y = -9.f; b_pCM_perp_mag = -9.f;
    b_prel_perp_x = -9.f; b_prel_perp_y = -9.f; b_prel_perp_mag = -9.f;
    b_k = -9.f; b_k2 = -9.f; b_k_z = -9.f; b_m_bar = -9.f;
    b_lc_quality = false;
    b_mMiss = -9.f;  b_kMiss = -9.f;       b_EMiss = -9.f;
    b_E0miss = -9.f; b_E1miss = -9.f;
    b_theta_PmQ = -9.f;
    b_theta_PmPlead = -9.f;
    b_theta_PleadQ = -9.f;
    b_chi = -9.f;
    b_theta_PleadPrec = -9.f;
    b_theta_PmPrec = -9.f;
    b_theta_PrecQ  = -9.f;
    b_phiTrento = -9.f;
    b_goodLead = false;

    b_xB_truth = -9.f;  b_Q2_truth = -9.f;  b_omega_truth = -9.f;
    b_eP_truth = -9.f;  b_eTheta_truth = -9.f;  b_ePhi_truth = -9.f;
    b_qP_truth = -9.f;  b_qTheta_truth = -9.f;  b_qPhi_truth = -9.f;
    b_leadP_truth = -9.f;  b_leadTheta_truth = -9.f;  b_leadPhi_truth = -9.f;
    b_pMiss_truth = -9.f;  b_pMissTheta_truth = -9.f;  b_pMissPhi_truth = -9.f;
    b_mMiss_truth = -9.f;  b_kMiss_truth = -9.f;  b_EMiss_truth = -9.f;
    b_theta_PmQ_truth = -9.f;
    b_theta_PmPlead_truth = -9.f;
    b_theta_PleadQ_truth = -9.f;
    b_chi_truth = -9.f;
    b_pRel_truth = -9.f;   b_pRelTheta_truth = -9.f;   b_pRelPhi_truth = -9.f;
    b_pRelx_truth = -9.f;  b_pRely_truth = -9.f;       b_pRelz_truth = -9.f;
    b_pLeadPlusRecOver2_truth = -9.f;
    b_pLeadMinusRec_truth = -9.f;
    b_pCM_truth = -9.f;    b_pCMx_truth = -9.f;        b_pCMy_truth = -9.f;    b_pCMz_truth = -9.f;
    b_pcmx_lab_truth = -9.f; b_pcmy_lab_truth = -9.f;  b_pcmz_lab_truth = -9.f;
    b_chi_frame_truth = -9.f;
    b_E2miss_truth = -9.f;
    b_alpha_1_truth = -9.f; b_alpha_2_truth = -9.f; b_alpha_CM_truth = -9.f; b_alpha_rel_truth = -9.f;
    b_p1_perp_x_truth = -9.f; b_p1_perp_y_truth = -9.f; b_p1_perp_mag_truth = -9.f;
    b_p2_perp_x_truth = -9.f; b_p2_perp_y_truth = -9.f; b_p2_perp_mag_truth = -9.f;
    b_pCM_perp_x_truth = -9.f; b_pCM_perp_y_truth = -9.f; b_pCM_perp_mag_truth = -9.f;
    b_prel_perp_x_truth = -9.f; b_prel_perp_y_truth = -9.f; b_prel_perp_mag_truth = -9.f;
    b_k_truth = -9.f; b_k2_truth = -9.f; b_k_z_truth = -9.f; b_m_bar_truth = -9.f;
    b_lc_quality_truth = false;
    b_E0miss_truth = -9.f; b_E1miss_truth = -9.f;
    b_recP_truth = -9.f;   b_recTheta_truth = -9.f;    b_recPhi_truth = -9.f;
    b_theta_PmPrec_truth = -9.f;
    b_theta_PrecQ_truth  = -9.f;
    b_phiTrento_truth = -9.f;

    clasAna.Run(c12);

    auto electrons = clasAna.getByPid(11);
    auto protons   = clasAna.getByPid(2212);

    if(electrons.size() != 1) continue;
    if(protons.empty())       continue;
    if(!inFD(electrons[0]->getStatus())) continue;

    // ---- electron kinematics ----
    eP4.SetXYZM(0., 0., 0., me);
    setElectronP4(eP4, electrons[0], isMC);
    if(eP4.Vect().Mag() == 0.) continue;

    TVector3 eP3 = eP4.Vect();
    TVector3 qP3 = p_beam - eP3;

    double omega = Ebeam - eP4.E();
    double Q2    = qP3.Mag2() - omega * omega;
    double xB    = Q2 / (2. * mP * omega);

    if(xB < 1.2) continue;
    if(Q2 < 1.5) continue;

    if(isMC){
      b_weight     = c12->mcevent()->getWeight();
      b_weight_ep  = b_weight * newWeight.get_weight_ep(c12->mcparts());
      b_weight_epp = b_weight * newWeight.get_weight_epp(c12->mcparts());
    }

    b_eP     = eP3.Mag();
    b_eTheta = eP3.Theta();
    b_ePhi   = eP3.Phi();

    b_qP     = qP3.Mag();
    b_qTheta = qP3.Theta();
    b_qPhi   = qP3.Phi();

    b_Q2    = Q2;
    b_xB    = xB;
    b_omega = omega;

    TLorentzVector q(qP3, omega);

    // ---- pass 1: compute kinematics for every proton candidate ----
    int nFilled = 0;
    vector<TVector3> cand_p3;
    vector<float>    cand_beta, cand_tof;
    vector<float>    cand_vz;
    vector<TVector3> cand_pMissV;
    vector<float>    cand_mMiss, cand_kMiss, cand_EMiss, cand_E0miss,
             cand_E1miss, cand_theta_PmQ;
    vector<bool>     cand_goodLead;

    for(int pr = 0; pr < (int)protons.size(); pr++)
    {
      if(nFilled >= MAXP) break;

      setProtonP4(leadP4, protons[pr], isMC);
      TVector3 pLead3 = leadP4.Vect();

//      pLead3.SetMag(pLead3.Mag() - .03); // apply 1% momentum scale correction to lead proton
      leadP4.SetVectM(pLead3, mP);

      missP4 = targP4 + beamP4 - eP4 - leadP4;
      TLorentzVector lead_miss_p4 = leadP4 - (beamP4 - eP4);

      TVector3 pMissV = pLead3 - qP3;

      // kMiss (light-cone definition)
      TLorentzVector miss_LC = leadP4 - q;
      TVector3 u_ZQ    = q.Vect().Unit();
      double pmm_ZQ    = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
      double pmp_ZQ    = miss_LC.Vect().Perp(u_ZQ);
      double mass_p    = mP;
      double kmiss_ZQ  = sqrt(mass_p*mass_p*((pmp_ZQ*pmp_ZQ + mass_p*mass_p)/(pmm_ZQ*(2*mass_p - pmm_ZQ))) - mass_p*mass_p);
      double kMiss     = kmiss_ZQ;

      double EMiss = sqrt(pLead3.Mag2() + mP*mP) - omega;
      double E0miss = sqrt(pMissV.Mag2() + mN*mN) - mN;
      double TP = leadP4.E() - leadP4.M();
      double TB = missP4.E() - missP4.M();
      double E1miss = omega - TP - TB;

      // SRC lead cuts
      bool passCuts = true;
      if(pLead3.Mag() < 1.)                        passCuts = false;
      if(missP4.M() < 0.65 || missP4.M() > 1.1)   passCuts = false;
      if(kMiss < 0.3 || kMiss > 1.)                passCuts = false;
 //     if(pLead3.Angle(qP3) < 37.*M_PI/180.)             passCuts = false;

      cand_p3.push_back(pLead3);
      cand_beta.push_back(protons[pr]->par()->getBeta());
      cand_tof.push_back(protons[pr]->getTime() - c12->event()->getStartTime());
      cand_vz.push_back(protons[pr]->par()->getVz());
      cand_pMissV.push_back(pMissV);
      cand_mMiss.push_back(missP4.M());
      cand_kMiss.push_back(kMiss);
      cand_EMiss.push_back(EMiss);
      cand_E0miss.push_back(E0miss);
      cand_E1miss.push_back(E1miss);
      cand_theta_PmQ.push_back(pMissV.Angle(qP3));
      cand_goodLead.push_back(passCuts);

      if(passCuts) b_nGoodLeads++;

      nFilled++;
    }

    b_nProtons      = nFilled;
    b_singleGoodLead = (b_nGoodLeads == 1);

    // ---- identify the lead: first candidate passing all SRC lead cuts ----
    int leadIdx = -1;
    for(int i = 0; i < nFilled; i++){
      if(cand_goodLead[i]){ leadIdx = i; break; }
    }

    if(leadIdx >= 0)
    {
      b_leadP       = cand_p3[leadIdx].Mag();
      b_leadTheta   = cand_p3[leadIdx].Theta();
      b_leadPhi     = cand_p3[leadIdx].Phi();
      b_leadBeta    = cand_beta[leadIdx];
      b_leadToF     = cand_tof[leadIdx];
      b_leadVz      = cand_vz[leadIdx];

      b_pMiss       = cand_pMissV[leadIdx].Mag();
      b_pMissTheta  = cand_pMissV[leadIdx].Theta();
      b_pMissPhi    = cand_pMissV[leadIdx].Phi();
      b_mMiss       = cand_mMiss[leadIdx];
      b_kMiss       = cand_kMiss[leadIdx];
      b_EMiss       = cand_EMiss[leadIdx];
      b_E0miss      = cand_E0miss[leadIdx];
      b_E1miss      = cand_E1miss[leadIdx];
      b_theta_PmQ   = cand_theta_PmQ[leadIdx];
      b_theta_PmPlead = cand_pMissV[leadIdx].Angle(cand_p3[leadIdx]);
      b_goodLead    = cand_goodLead[leadIdx];
      b_theta_PleadQ = cand_p3[leadIdx].Angle(qP3);
      b_chi = getChi(cand_pMissV[leadIdx], qP3);
    }

    // ---- pass 2: find recoil and fill pRel / pCM (using the lead identified above) ----
    if(leadIdx >= 0 && nFilled >= 2)
    {
      // The lead proton's pMiss vector (miss_neg)
      TVector3 miss_neg = cand_pMissV[leadIdx];

      // Search for a recoil: any other proton with p > 0.3
      // Recoil angular cuts: open (no FD theta < 37 restriction); require CD (or simply p > 0.3).
      for(int j = 0; j < nFilled; j++)
      {
        if(j == leadIdx) continue;

        // cand_p3[j].SetMag(cand_p3[j].Mag() - .03); // apply nA potential correction 
        
        if(cand_p3[j].Mag() < 0.3) continue;   // recoil momentum cut

        // Optional recoil angular cut — looser than lead, accept both FD and CD.
        // Add detector-based cuts here if needed via protons[j]->getStatus().

        TVector3 recoil_p3 = cand_p3[j];

        // ---- pRel: relative momentum in the pair CM frame ----
        // pRel = (p_lead - p_recoil) / 2  in the pair rest frame; approximate as half the
        // difference of lab momenta (standard SRC convention).
        TVector3 lead_p3 = cand_p3[leadIdx];

        TVector3 pRelV = (miss_neg - recoil_p3) * 0.5;
        b_pRel      = pRelV.Mag();
        b_pRelTheta = pRelV.Theta();
        b_pRelPhi   = pRelV.Phi();
        b_pLeadPlusRecOver2 = 0.5f * static_cast<Float_t>((lead_p3 + recoil_p3).Mag());
        b_pLeadMinusRec = static_cast<Float_t>((lead_p3 - recoil_p3).Mag());

        // ---- pCM: pair CM momentum projected onto the (miss_neg, q) frame ----
        TVector3 v_rec = recoil_p3;
        TVector3 v_cm  = miss_neg + v_rec;

        TVector3 vz = miss_neg.Unit();
        TVector3 vy = miss_neg.Cross(qP3).Unit();
        TVector3 vx = vz.Cross(vy).Unit();

        b_pCM  = v_cm.Mag();
        b_pRelx = pRelV.Dot(vx);
        b_pRely = pRelV.Dot(vy);
        b_pRelz = pRelV.Dot(vz);
        b_pCMx = v_cm.Dot(vx);
        b_pCMy = v_cm.Dot(vy);
        b_pCMz = v_cm.Dot(vz);
        b_pcmx_lab = v_cm.X();
        b_pcmy_lab = v_cm.Y();
        b_pcmz_lab = v_cm.Z();
        b_chi_frame = getChiFrame(miss_neg, qP3, recoil_p3);

        // Same E2miss definition used by Main_Figs_Binned.cpp.
        TLorentzVector selectedLeadP4;
        selectedLeadP4.SetVectM(cand_p3[leadIdx], mP);
        TLorentzVector recoilP4;
        recoilP4.SetVectM(recoil_p3, mP);
        double TP = selectedLeadP4.E() - selectedLeadP4.M();
        double TP2 = recoilP4.E() - recoilP4.M();
        TLorentzVector miss_Am2 = q + nucleusP4 - selectedLeadP4 - recoilP4;
        double TB2 = miss_Am2.E() - miss_Am2.M();
        b_E2miss = q.E() - TP - TP2 - TB2;

        // recoil summary kinematics
        b_recP     = recoil_p3.Mag();
        b_recTheta = recoil_p3.Theta();
        b_recPhi   = recoil_p3.Phi();
        b_recBeta  = cand_beta[j];
        b_recToF   = cand_tof[j];

        // additional angles
        b_theta_PleadPrec = lead_p3.Angle(recoil_p3);
        b_theta_PmPrec = miss_neg.Angle(recoil_p3);
        b_theta_PrecQ  = recoil_p3.Angle(qP3);
        b_phiTrento = getPhiTrentoZeroCM(lead_p3, miss_neg, recoil_p3);

        const LightConeKinematics lc = computeLightConeKinematics(
          qP3, omega, lead_p3, recoil_p3, mP, mP, targetMass, targetMassNumber);
        if(lc.validBasis && lc.pairDefined){
          b_alpha_1 = lc.alpha1;
          b_alpha_2 = lc.alpha2;
          b_alpha_CM = lc.alphaCM;
          b_alpha_rel = lc.alphaRel;
          b_p1_perp_x = lc.p1PerpX;
          b_p1_perp_y = lc.p1PerpY;
          b_p1_perp_mag = lc.p1PerpMag;
          b_p2_perp_x = lc.p2PerpX;
          b_p2_perp_y = lc.p2PerpY;
          b_p2_perp_mag = lc.p2PerpMag;
          b_pCM_perp_x = lc.pCMPerpX;
          b_pCM_perp_y = lc.pCMPerpY;
          b_pCM_perp_mag = lc.pCMPerpMag;
          b_prel_perp_x = lc.pRelPerpX;
          b_prel_perp_y = lc.pRelPerpY;
          b_prel_perp_mag = lc.pRelPerpMag;
          b_k = lc.k;
          b_k2 = lc.k2;
          b_k_z = lc.kZ;
          b_m_bar = lc.mBar;
          b_lc_quality = lc.qualityBit;

          h_alpha_1->Fill(lc.alpha1);
          h_alpha_2->Fill(lc.alpha2);
          h_alpha_CM->Fill(lc.alphaCM);
          h_alpha_spectator->Fill(targetMassNumber - lc.alphaCM);
          h_k2->Fill(lc.k2);
          if(b_pRel > 0. && lc.k >= 0.) h_k_vs_pRel->Fill(b_pRel, lc.k);
        }

        break;  // one recoil per event
      }
    }

    if(isMC)
    {
      auto mcInfo = c12->mcparts();
      if(mcInfo && mcInfo->getRows() >= 2)
      {
        TVector3 e_truth(mcInfo->getPx(0), mcInfo->getPy(0), mcInfo->getPz(0));
        TVector3 lead_truth(mcInfo->getPx(1), mcInfo->getPy(1), mcInfo->getPz(1));
        TVector3 q_truth = p_beam - e_truth;

        double omega_truth = Ebeam - sqrt(e_truth.Mag2() + me*me);
        double Q2_truth = q_truth.Mag2() - omega_truth * omega_truth;
        double xB_truth = (omega_truth != 0.) ? Q2_truth / (2. * mP * omega_truth) : -9.;

        TLorentzVector eP4_truth;
        eP4_truth.SetVectM(e_truth, me);
        TLorentzVector leadP4_truth;
        leadP4_truth.SetVectM(lead_truth, mP);
        TLorentzVector q4_truth(q_truth, omega_truth);
        TLorentzVector missP4_truth = targP4 + beamP4 - eP4_truth - leadP4_truth;
        TVector3 pMiss_truth = lead_truth - q_truth;

        TLorentzVector miss_LC_truth = leadP4_truth - q4_truth;
        TVector3 u_ZQ_truth = q4_truth.Vect().Unit();
        double pmm_ZQ_truth = miss_LC_truth.E() - miss_LC_truth.Vect().Dot(u_ZQ_truth);
        double pmp_ZQ_truth = miss_LC_truth.Vect().Perp(u_ZQ_truth);
        double kmiss_denom_truth = pmm_ZQ_truth * (2*mP - pmm_ZQ_truth);
        double kMiss_truth = -9.;
        if(kmiss_denom_truth != 0.){
          double kmiss_arg_truth = mP*mP*((pmp_ZQ_truth*pmp_ZQ_truth + mP*mP)/kmiss_denom_truth) - mP*mP;
          kMiss_truth = (kmiss_arg_truth >= 0.) ? sqrt(kmiss_arg_truth) : -9.;
        }
        double EMiss_truth = sqrt(lead_truth.Mag2() + mP*mP) - omega_truth;
        double E0miss_truth = sqrt(pMiss_truth.Mag2() + mN*mN) - mN;
        double TP_truth = leadP4_truth.E() - leadP4_truth.M();
        double TB_truth = missP4_truth.E() - missP4_truth.M();
        double E1miss_truth = omega_truth - TP_truth - TB_truth;

        b_eP_truth     = e_truth.Mag();
        b_eTheta_truth = e_truth.Theta();
        b_ePhi_truth   = e_truth.Phi();

        b_qP_truth     = q_truth.Mag();
        b_qTheta_truth = q_truth.Theta();
        b_qPhi_truth   = q_truth.Phi();

        b_Q2_truth    = Q2_truth;
        b_xB_truth    = xB_truth;
        b_omega_truth = omega_truth;

        b_leadP_truth     = lead_truth.Mag();
        b_leadTheta_truth = lead_truth.Theta();
        b_leadPhi_truth   = lead_truth.Phi();

        b_pMiss_truth      = pMiss_truth.Mag();
        b_pMissTheta_truth = pMiss_truth.Theta();
        b_pMissPhi_truth   = pMiss_truth.Phi();
        b_mMiss_truth      = missP4_truth.M();
        b_kMiss_truth      = kMiss_truth;
        b_EMiss_truth      = EMiss_truth;
        b_E0miss_truth     = E0miss_truth;
        b_E1miss_truth     = E1miss_truth;
        b_theta_PmQ_truth  = pMiss_truth.Angle(q_truth);
        b_theta_PmPlead_truth = pMiss_truth.Angle(lead_truth);
        b_theta_PleadQ_truth = lead_truth.Angle(q_truth);
        b_chi_truth = getChi(pMiss_truth, q_truth);

        if(mcInfo->getRows() >= 3)
        {
          TVector3 rec_truth(mcInfo->getPx(2), mcInfo->getPy(2), mcInfo->getPz(2));
          TVector3 pRel_truth = (pMiss_truth - rec_truth) * 0.5;
          TVector3 pCM_truth = pMiss_truth + rec_truth;

          TVector3 vz_truth = pMiss_truth.Unit();
          TVector3 vy_truth = pMiss_truth.Cross(q_truth).Unit();
          TVector3 vx_truth = vz_truth.Cross(vy_truth).Unit();
          TVector3 v_cm_truth = pMiss_truth + rec_truth;

          b_pRel_truth      = pRel_truth.Mag();
          b_pRelTheta_truth = pRel_truth.Theta();
          b_pRelPhi_truth   = pRel_truth.Phi();
          b_pRelx_truth     = pRel_truth.Dot(vx_truth);
          b_pRely_truth     = pRel_truth.Dot(vy_truth);
          b_pRelz_truth     = pRel_truth.Dot(vz_truth);
          b_pLeadPlusRecOver2_truth = 0.5f * static_cast<Float_t>((lead_truth + rec_truth).Mag());
          b_pLeadMinusRec_truth = static_cast<Float_t>((lead_truth - rec_truth).Mag());

          b_pCM_truth  = pCM_truth.Mag();
          b_pCMx_truth = pCM_truth.Dot(vx_truth);
          b_pCMy_truth = pCM_truth.Dot(vy_truth);
          b_pCMz_truth = pCM_truth.Dot(vz_truth);
          b_pcmx_lab_truth = v_cm_truth.X();
          b_pcmy_lab_truth = v_cm_truth.Y();
          b_pcmz_lab_truth = v_cm_truth.Z();
          b_chi_frame_truth = getChiFrame(pMiss_truth, q_truth, rec_truth);

          TLorentzVector recP4_truth;
          recP4_truth.SetVectM(rec_truth, mP);
          double TP2_truth = recP4_truth.E() - recP4_truth.M();
          TLorentzVector miss_Am2_truth = q4_truth + nucleusP4 - leadP4_truth - recP4_truth;
          double TB2_truth = miss_Am2_truth.E() - miss_Am2_truth.M();
          b_E2miss_truth = q4_truth.E() - TP_truth - TP2_truth - TB2_truth;

          b_recP_truth     = rec_truth.Mag();
          b_recTheta_truth = rec_truth.Theta();
          b_recPhi_truth   = rec_truth.Phi();
          b_theta_PmPrec_truth = pMiss_truth.Angle(rec_truth);
          b_theta_PrecQ_truth  = rec_truth.Angle(q_truth);
          b_phiTrento_truth = getPhiTrentoZeroCM(lead_truth, pMiss_truth, rec_truth);

          const LightConeKinematics lcTruth = computeLightConeKinematics(
              q_truth, omega_truth, lead_truth, rec_truth, mP, mP, targetMass,
              targetMassNumber);
          if(lcTruth.validBasis && lcTruth.pairDefined){
            b_alpha_1_truth = lcTruth.alpha1;
            b_alpha_2_truth = lcTruth.alpha2;
            b_alpha_CM_truth = lcTruth.alphaCM;
            b_alpha_rel_truth = lcTruth.alphaRel;
            b_p1_perp_x_truth = lcTruth.p1PerpX;
            b_p1_perp_y_truth = lcTruth.p1PerpY;
            b_p1_perp_mag_truth = lcTruth.p1PerpMag;
            b_p2_perp_x_truth = lcTruth.p2PerpX;
            b_p2_perp_y_truth = lcTruth.p2PerpY;
            b_p2_perp_mag_truth = lcTruth.p2PerpMag;
            b_pCM_perp_x_truth = lcTruth.pCMPerpX;
            b_pCM_perp_y_truth = lcTruth.pCMPerpY;
            b_pCM_perp_mag_truth = lcTruth.pCMPerpMag;
            b_prel_perp_x_truth = lcTruth.pRelPerpX;
            b_prel_perp_y_truth = lcTruth.pRelPerpY;
            b_prel_perp_mag_truth = lcTruth.pRelPerpMag;
            b_k_truth = lcTruth.k;
            b_k2_truth = lcTruth.k2;
            b_k_z_truth = lcTruth.kZ;
            b_m_bar_truth = lcTruth.mBar;
            b_lc_quality_truth = lcTruth.qualityBit;
          }
        }
      }
    }

    srcTree->Fill();
    ctr++;
  }

  outFile->cd();
  srcTree->Write();
  h_alpha_1->Write();
  h_alpha_2->Write();
  h_alpha_CM->Write();
  h_alpha_spectator->Write();
  h_k2->Write();
  h_k_vs_pRel->Write();
  outFile->Close();

  cout << "Done. Processed " << counter << " events.\n";
  return 0;
}
