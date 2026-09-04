#include <algorithm>
#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TVector3.h>

namespace {

const double kMp = 0.9382720813;
const double kMn = 0.9395654133;
const double kMu = 0.9314941024;
const double kMe = 0.000510999;
const double kMd = 2.01410178 * kMu - kMe;
const double kLightConeEdgeEpsilon = 1e-6;
const double kLightConeRoundoffTolerance = 1e-9;

TRandom3 gSmearRand(0);

const double params_Smear_FD[6][6] = {
  {0.242269, 0.012285, 0.000722975, 0.229986, 0.00155551, 0.00088862},
  {-0.0331068, 0.0452348, -4.87551e-05, 0.557431, -0.0323144, 0.00163881},
  {-0.0781403, 0.0453914, 0.000252317, -0.465264, 0.0735487, -0.000816354},
  {0.0767921, 0.0289579, 0.000467256, 0.267525, -0.00692348, 0.00114717},
  {-0.0501808, 0.0404422, 0.000173173, 0.262579, -0.00198293, 0.000950558},
  {0.167938, 0.0297138, 0.000191341, -0.143651, 0.0362719, 0.000199702}};

const double params_Smear_CD[2] = {12.4637, 6.76708};

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
  double mBar = -9.;
  double alpha1 = -9.;
  double alpha2 = -9.;
  double alphaCM = -9.;
  double alphaRel = -9.;
  double qMinus = -9.;
  double p1Plus = -9.;
  double p2Plus = -9.;
  double k = -9.;
  double k2 = -9.;
  double kZ = -9.;
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

struct TargetConfig {
  int A = 12;
  double mass = 12.0 * kMu - 6.0 * kMe;
};

TargetConfig buildTargetConfig(int targetA) {
  TargetConfig cfg;
  cfg.A = targetA;
  switch (targetA) {
    case 2:
      cfg.mass = 2.01410178 * kMu - kMe;
      break;
    case 4:
      cfg.mass = 4.00260325415 * kMu - 2.0 * kMe;
      break;
    case 12:
      cfg.mass = 12.0 * kMu - 6.0 * kMe;
      break;
    case 40:
      cfg.mass = 39.962591 * kMu - 20.0 * kMe;
      break;
    case 48:
      cfg.mass = 47.952523 * kMu - 20.0 * kMe;
      break;
    default:
      cfg.A = 12;
      cfg.mass = 12.0 * kMu - 6.0 * kMe;
      std::cerr << "[convert_events2N_to_srcTree] Unsupported target A=" << targetA
                << ", defaulting to A=12.\n";
      break;
  }
  return cfg;
}

LightConeBasis makeLightConeBasis(const TVector3 &qVec) {
  LightConeBasis basis;
  if (qVec.Mag2() <= 0.) return basis;

  basis.zHat = qVec.Unit();

  const TVector3 beamAxis(0., 0., 1.);
  TVector3 xSeed = beamAxis - basis.zHat * beamAxis.Dot(basis.zHat);
  if (xSeed.Mag2() <= 1e-12) {
    const TVector3 xLab(1., 0., 0.);
    xSeed = xLab - basis.zHat * xLab.Dot(basis.zHat);
  }
  if (xSeed.Mag2() <= 1e-12) {
    const TVector3 yLab(0., 1., 0.);
    xSeed = yLab - basis.zHat * yLab.Dot(basis.zHat);
  }
  if (xSeed.Mag2() <= 1e-12) return basis;

  basis.xHat = xSeed.Unit();
  basis.yHat = basis.zHat.Cross(basis.xHat);
  if (basis.yHat.Mag2() <= 1e-12) return basis;

  basis.yHat = basis.yHat.Unit();
  basis.xHat = basis.yHat.Cross(basis.zHat).Unit();
  basis.valid = true;
  return basis;
}

LightConeKinematics computeLightConeKinematics(const TVector3 &qVec, double nu,
                                               const TVector3 &pLead,
                                               const TVector3 &pRec,
                                               double mLead,
                                               double mRec,
                                               double targetMass,
                                               int massNumber) {
  LightConeKinematics out;
  if (targetMass <= 0. || massNumber <= 0) return out;

  out.mBar = targetMass / static_cast<double>(massNumber);
  const LightConeBasis basis = makeLightConeBasis(qVec);
  if (!basis.valid) return out;

  out.validBasis = true;
  out.qualityBit = true;

  const double qMag = qVec.Mag();
  out.qMinus = nu - qMag;
  const double qPlus = nu + qMag;

  const double eLead = std::sqrt(pLead.Mag2() + mLead * mLead);
  const double eRec = std::sqrt(pRec.Mag2() + mRec * mRec);
  const double pLeadZ = pLead.Dot(basis.zHat);
  const double pRecZ = pRec.Dot(basis.zHat);

  out.p1Plus = eLead + pLeadZ - qPlus;
  out.p2Plus = eRec + pRecZ;

  out.alpha1 = (eLead - pLeadZ - out.qMinus) / out.mBar;
  out.alpha2 = (eRec - pRecZ) / out.mBar;

  const TVector3 p1Perp = pLead - basis.zHat * pLeadZ;
  const TVector3 p2Perp = pRec - basis.zHat * pRecZ;
  out.p1PerpX = p1Perp.Dot(basis.xHat);
  out.p1PerpY = p1Perp.Dot(basis.yHat);
  out.p1PerpMag = p1Perp.Mag();
  out.p2PerpX = p2Perp.Dot(basis.xHat);
  out.p2PerpY = p2Perp.Dot(basis.yHat);
  out.p2PerpMag = p2Perp.Mag();

  out.alphaCM = out.alpha1 + out.alpha2;
  if (std::abs(out.alphaCM) <= kLightConeRoundoffTolerance) {
    out.qualityBit = false;
    return out;
  }

  out.pairDefined = true;
  const TVector3 pCMPerp = p1Perp + p2Perp;
  out.pCMPerpX = pCMPerp.Dot(basis.xHat);
  out.pCMPerpY = pCMPerp.Dot(basis.yHat);
  out.pCMPerpMag = pCMPerp.Mag();

  out.alphaRel = 2. * out.alpha2 / out.alphaCM;
  const TVector3 pRelPerp = p2Perp - (out.alpha2 / out.alphaCM) * pCMPerp;
  out.pRelPerpX = pRelPerp.Dot(basis.xHat);
  out.pRelPerpY = pRelPerp.Dot(basis.yHat);
  out.pRelPerpMag = pRelPerp.Mag();

  const double alphaFactor = out.alphaRel * (2. - out.alphaRel);
  if (out.alphaRel <= kLightConeEdgeEpsilon ||
      out.alphaRel >= 2. - kLightConeEdgeEpsilon) {
    out.qualityBit = false;
  }
  if (alphaFactor <= 0.) {
    out.qualityBit = false;
    return out;
  }

  const double pairMass = 0.5 * (mLead + mRec);
  const double rawK2 =
      (pairMass * pairMass + pRelPerp.Mag2()) / alphaFactor - pairMass * pairMass;
  out.k2 = std::max(0., rawK2);
  out.k = std::sqrt(out.k2);
  out.kZ = (out.alphaRel - 1.) * std::sqrt(pairMass * pairMass + out.k2);
  return out;
}

double getMassFromType(int t) {
  const int at = std::abs(t);
  if (at == 2112) return kMn;
  return kMp;
}

void resetFloat(float &x) { x = -9.f; }

double functionAlgebraic5(double x, double a, double b, double c) {
  return a + b * x + c * x * x;
}

double subtractQuadratureFD(double x, double a, double b, double c,
                            double d, double e, double f) {
  const double p1 = functionAlgebraic5(x, a, b, c);
  const double p2 = functionAlgebraic5(x, d, e, f);
  const double val = p1 * p1 - p2 * p2;
  if (val > 0.01) return std::sqrt(val);
  return 0.1;
}

double subtractQuadratureCD(double a, double b) {
  const double val = a * a - b * b;
  if (val > 0.01) return std::sqrt(val);
  return 0.1;
}

double wrapPhiDeg(double phiDeg) {
  double out = phiDeg;
  while (out <= -180.) out += 360.;
  while (out > 180.) out -= 360.;
  return out;
}

int sectorFromPhiDeg(double phiDeg) {
  const double phi = wrapPhiDeg(phiDeg);
  if (phi > -30. && phi <= 30.) return 1;
  if (phi > 30. && phi <= 90.) return 2;
  if (phi > 90. && phi <= 150.) return 3;
  if (phi > 150. || phi <= -150.) return 4;
  if (phi > -150. && phi <= -90.) return 5;
  return 6;
}

void applyMomentumSimulationSmearApprox(TLorentzVector &p4) {
  TVector3 v3 = p4.Vect();
  const double mom = v3.Mag();
  if (mom <= 0.) return;

  const double thetaDeg = v3.Theta() * 180.0 / M_PI;
  const double phiDeg = v3.Phi() * 180.0 / M_PI;

  double smearFrac = 0.;
  if (thetaDeg < 37.0) {
    const int sector = sectorFromPhiDeg(phiDeg);
    const int idx = sector - 1;
    smearFrac = 0.01 * subtractQuadratureFD(
                            thetaDeg,
                            params_Smear_FD[idx][0], params_Smear_FD[idx][1],
                            params_Smear_FD[idx][2], params_Smear_FD[idx][3],
                            params_Smear_FD[idx][4], params_Smear_FD[idx][5]);
  } else if (thetaDeg > 45.0) {
    smearFrac = 0.01 * subtractQuadratureCD(params_Smear_CD[0], params_Smear_CD[1]);
  } else {
    return;
  }

  const double scale = gSmearRand.Gaus(1.0, smearFrac);
  v3.SetMag(mom * scale);
  p4.SetVectM(v3, p4.M());
}

}  // namespace

void convert_events2N_to_srcTree(
    const char *inputFileName = "~/data/RGM_DATA/events_2N.root",
    const char *outputFileName = "~/data/RGM_DATA/events_2N_srcTree.root",
    const char *inputTreeName = "events",
    bool useFSIAwareMomenta = true,
    int targetA = 12,
  double eBeamOverride = -1.0,
  bool requireEpp = true,
    bool applyMCSmearing = true,
  bool applyBasicSrcCuts = true,
  double minXB = 1.2,
  double minQ2 = 1.5,
  double minLeadP = 1.0,
  double minRecP = 0.3,
  double minKMiss = 0.3,
  double maxKMiss = 1.0,
  double minMMiss = 0.65,
  double maxMMiss = 1.1) {

  TFile *inFile = TFile::Open(inputFileName, "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "[convert_events2N_to_srcTree] Could not open input file: "
              << inputFileName << "\n";
    return;
  }

  TTree *inTree = dynamic_cast<TTree *>(inFile->Get(inputTreeName));
  if (!inTree) {
    std::cerr << "[convert_events2N_to_srcTree] Could not find tree \"" << inputTreeName
              << "\" in " << inputFileName << "\n";
    inFile->Close();
    return;
  }

  double weight = 1.0;
  int lead_type = 2212;
  int rec_type = 2212;
  int doFSI = 1;
  double electron[4] = {0., 0., 0., 0.};
  double lead_post[4] = {0., 0., 0., 0.};
  double recoil_post[4] = {0., 0., 0., 0.};
  double lead_pre[4] = {0., 0., 0., 0.};
  double recoil_pre[4] = {0., 0., 0., 0.};
  double qv[4] = {0., 0., 0., 0.};
  double Q2_in = -9.;
  double xB_in = -9.;
  double nu_in = -9.;

  inTree->SetBranchAddress("weight", &weight);
  inTree->SetBranchAddress("lead_type", &lead_type);
  inTree->SetBranchAddress("rec_type", &rec_type);
  inTree->SetBranchAddress("doFSI", &doFSI);
  inTree->SetBranchAddress("electron", electron);
  inTree->SetBranchAddress("lead_post", lead_post);
  inTree->SetBranchAddress("recoil_post", recoil_post);
  inTree->SetBranchAddress("lead_pre", lead_pre);
  inTree->SetBranchAddress("recoil_pre", recoil_pre);
  inTree->SetBranchAddress("q", qv);
  inTree->SetBranchAddress("Q2", &Q2_in);
  inTree->SetBranchAddress("xB", &xB_in);
  inTree->SetBranchAddress("nu", &nu_in);

  TFile *outFile = TFile::Open(outputFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "[convert_events2N_to_srcTree] Could not create output file: "
              << outputFileName << "\n";
    inFile->Close();
    return;
  }

  TTree *srcTree = new TTree("srcTree", "SRC kinematics converted from events_2N");

  std::cout << "[convert_events2N_to_srcTree] Settings"
            << " requireEpp=" << (requireEpp ? "true" : "false")
            << " applyMCSmearing=" << (applyMCSmearing ? "true" : "false")
            << " applyBasicSrcCuts=" << (applyBasicSrcCuts ? "true" : "false")
            << "\n";

  float weight_out, weight_ep, weight_epp;
  float xB, Q2, omega;
  float eP, eTheta, ePhi;
  float qP, qTheta, qPhi;
  int nProtons;
  float leadP, leadTheta, leadPhi;
  float recP, recTheta, recPhi;
  float pMiss, pMissTheta, pMissPhi;
  float kMiss, mMiss, EMiss, E0miss, E1miss, E2miss;
  float theta_PmQ, theta_PmPlead, theta_PleadQ;
  float theta_PleadPrec, theta_PmPrec, theta_PrecQ;
  float pRel, pRelTheta, pRelPhi, pRelx, pRely, pRelz;
  float pLeadPlusRecOver2, pLeadMinusRec;
  float pCM, pCMx, pCMy, pCMz;
  float pcmx_lab, pcmy_lab, pcmz_lab;
  float alpha_1, alpha_2, alpha_CM, alpha_rel;
  float alpha_q, alpha_pLead;
  float p1_plus, p2_plus;
  float p1_perp_x, p1_perp_y, p1_perp_mag;
  float p2_perp_x, p2_perp_y, p2_perp_mag;
  float pCM_perp_x, pCM_perp_y, pCM_perp_mag;
  float prel_perp_x, prel_perp_y, prel_perp_mag;
  float k, k2, k_z, m_bar;

  srcTree->Branch("weight", &weight_out, "weight/F");
  srcTree->Branch("weight_ep", &weight_ep, "weight_ep/F");
  srcTree->Branch("weight_epp", &weight_epp, "weight_epp/F");

  srcTree->Branch("xB", &xB, "xB/F");
  srcTree->Branch("Q2", &Q2, "Q2/F");
  srcTree->Branch("omega", &omega, "omega/F");

  srcTree->Branch("eP", &eP, "eP/F");
  srcTree->Branch("eTheta", &eTheta, "eTheta/F");
  srcTree->Branch("ePhi", &ePhi, "ePhi/F");

  srcTree->Branch("qP", &qP, "qP/F");
  srcTree->Branch("qTheta", &qTheta, "qTheta/F");
  srcTree->Branch("qPhi", &qPhi, "qPhi/F");

  srcTree->Branch("nProtons", &nProtons, "nProtons/I");

  srcTree->Branch("leadP", &leadP, "leadP/F");
  srcTree->Branch("leadTheta", &leadTheta, "leadTheta/F");
  srcTree->Branch("leadPhi", &leadPhi, "leadPhi/F");
  srcTree->Branch("recP", &recP, "recP/F");
  srcTree->Branch("recTheta", &recTheta, "recTheta/F");
  srcTree->Branch("recPhi", &recPhi, "recPhi/F");

  srcTree->Branch("pMiss", &pMiss, "pMiss/F");
  srcTree->Branch("pMissTheta", &pMissTheta, "pMissTheta/F");
  srcTree->Branch("pMissPhi", &pMissPhi, "pMissPhi/F");

  srcTree->Branch("kMiss", &kMiss, "kMiss/F");
  srcTree->Branch("mMiss", &mMiss, "mMiss/F");
  srcTree->Branch("EMiss", &EMiss, "EMiss/F");
  srcTree->Branch("E0miss", &E0miss, "E0miss/F");
  srcTree->Branch("E1miss", &E1miss, "E1miss/F");
  srcTree->Branch("E2miss", &E2miss, "E2miss/F");

  srcTree->Branch("theta_PmQ", &theta_PmQ, "theta_PmQ/F");
  srcTree->Branch("theta_PmPlead", &theta_PmPlead, "theta_PmPlead/F");
  srcTree->Branch("theta_PleadQ", &theta_PleadQ, "theta_PleadQ/F");
  srcTree->Branch("theta_PleadPrec", &theta_PleadPrec, "theta_PleadPrec/F");
  srcTree->Branch("theta_PmPrec", &theta_PmPrec, "theta_PmPrec/F");
  srcTree->Branch("theta_PrecQ", &theta_PrecQ, "theta_PrecQ/F");

  srcTree->Branch("pRel", &pRel, "pRel/F");
  srcTree->Branch("pRelTheta", &pRelTheta, "pRelTheta/F");
  srcTree->Branch("pRelPhi", &pRelPhi, "pRelPhi/F");
  srcTree->Branch("pRelx", &pRelx, "pRelx/F");
  srcTree->Branch("pRely", &pRely, "pRely/F");
  srcTree->Branch("pRelz", &pRelz, "pRelz/F");

  srcTree->Branch("pLeadPlusRecOver2", &pLeadPlusRecOver2, "pLeadPlusRecOver2/F");
  srcTree->Branch("pLeadMinusRec", &pLeadMinusRec, "pLeadMinusRec/F");

  srcTree->Branch("pCM", &pCM, "pCM/F");
  srcTree->Branch("pCMx", &pCMx, "pCMx/F");
  srcTree->Branch("pCMy", &pCMy, "pCMy/F");
  srcTree->Branch("pCMz", &pCMz, "pCMz/F");
  srcTree->Branch("pcmx_lab", &pcmx_lab, "pcmx_lab/F");
  srcTree->Branch("pcmy_lab", &pcmy_lab, "pcmy_lab/F");
  srcTree->Branch("pcmz_lab", &pcmz_lab, "pcmz_lab/F");

  srcTree->Branch("alpha_1", &alpha_1, "alpha_1/F");
  srcTree->Branch("alpha_2", &alpha_2, "alpha_2/F");
  srcTree->Branch("alpha_CM", &alpha_CM, "alpha_CM/F");
  srcTree->Branch("alpha_rel", &alpha_rel, "alpha_rel/F");
  srcTree->Branch("alpha_q", &alpha_q, "alpha_q/F");
  srcTree->Branch("alpha_pLead", &alpha_pLead, "alpha_pLead/F");

  srcTree->Branch("p1_plus", &p1_plus, "p1_plus/F");
  srcTree->Branch("p2_plus", &p2_plus, "p2_plus/F");

  srcTree->Branch("p1_perp_x", &p1_perp_x, "p1_perp_x/F");
  srcTree->Branch("p1_perp_y", &p1_perp_y, "p1_perp_y/F");
  srcTree->Branch("p1_perp_mag", &p1_perp_mag, "p1_perp_mag/F");
  srcTree->Branch("p2_perp_x", &p2_perp_x, "p2_perp_x/F");
  srcTree->Branch("p2_perp_y", &p2_perp_y, "p2_perp_y/F");
  srcTree->Branch("p2_perp_mag", &p2_perp_mag, "p2_perp_mag/F");

  srcTree->Branch("pCM_perp_x", &pCM_perp_x, "pCM_perp_x/F");
  srcTree->Branch("pCM_perp_y", &pCM_perp_y, "pCM_perp_y/F");
  srcTree->Branch("pCM_perp_mag", &pCM_perp_mag, "pCM_perp_mag/F");
  srcTree->Branch("prel_perp_x", &prel_perp_x, "prel_perp_x/F");
  srcTree->Branch("prel_perp_y", &prel_perp_y, "prel_perp_y/F");
  srcTree->Branch("prel_perp_mag", &prel_perp_mag, "prel_perp_mag/F");

  srcTree->Branch("k", &k, "k/F");
  srcTree->Branch("k2", &k2, "k2/F");
  srcTree->Branch("k_z", &k_z, "k_z/F");
  srcTree->Branch("m_bar", &m_bar, "m_bar/F");

  const TargetConfig targetCfg = buildTargetConfig(targetA);
  const TLorentzVector targP4(0., 0., 0., kMd);
  const TLorentzVector nucleusP4(0., 0., 0., targetCfg.mass);

  Long64_t nTotal = 0;
  Long64_t passEpp = 0;
  Long64_t passXB = 0;
  Long64_t passQ2 = 0;
  Long64_t passLeadP = 0;
  Long64_t passRecP = 0;
  Long64_t passMMiss = 0;
  Long64_t passKMiss = 0;
  Long64_t nWritten = 0;

  const Long64_t nEntries = inTree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    inTree->GetEntry(i);
    ++nTotal;

    if (requireEpp) {
      if (std::abs(lead_type) != 2212 || std::abs(rec_type) != 2212) continue;
      ++passEpp;
    }

    const bool usePost = (!useFSIAwareMomenta) || (doFSI != 0);
    const double *lead4 = usePost ? lead_post : lead_pre;
    const double *rec4 = usePost ? recoil_post : recoil_pre;

    const double leadMass = getMassFromType(lead_type);
    const double recMass = getMassFromType(rec_type);

    TLorentzVector eP4;
    eP4.SetPxPyPzE(electron[0], electron[1], electron[2], electron[3]);
    TLorentzVector leadP4;
    leadP4.SetPxPyPzE(lead4[0], lead4[1], lead4[2], lead4[3]);
    TLorentzVector recP4;
    recP4.SetPxPyPzE(rec4[0], rec4[1], rec4[2], rec4[3]);

    // Match simpleSRCSkim MC treatment by smearing reconstructed momenta.
    if (applyMCSmearing) {
      applyMomentumSimulationSmearApprox(eP4);
      applyMomentumSimulationSmearApprox(leadP4);
      applyMomentumSimulationSmearApprox(recP4);
    }

    const double eBeamEvent = nu_in + electron[3];
    const TVector3 beamVec(0., 0., eBeamEvent);
    const TVector3 eVec = eP4.Vect();
    const TVector3 qVec = beamVec - eVec;
    const double nu = eBeamEvent - eP4.E();
    TLorentzVector qP4;
    qP4.SetVectM(qVec, 0.0);
    qP4.SetE(nu);

    const TVector3 leadVec = leadP4.Vect();
    const TVector3 recVec = recP4.Vect();

    const double pLeadMag = leadVec.Mag();
    const double pRecMag = recVec.Mag();
    const TVector3 pMissVec = leadVec - qVec;
    const double pMissMag = pMissVec.Mag();

    const double pmm_ZQ = (leadP4.E() - qP4.E()) - (leadP4.Vect() - qP4.Vect()).Dot(qVec.Unit());
    const double pmp_ZQ = (leadP4.Vect() - qP4.Vect()).Perp(qVec.Unit());
    double kMissVal = -9.;
    const double kmissDen = pmm_ZQ * (2.0 * kMp - pmm_ZQ);
    if (std::abs(kmissDen) > 1e-12) {
      const double arg = kMp * kMp * ((pmp_ZQ * pmp_ZQ + kMp * kMp) / kmissDen) - kMp * kMp;
      if (arg >= 0.) kMissVal = std::sqrt(arg);
    }

    // Match simpleSRCSkim lead-cut convention: mMiss from deuteron target + q - pLead.
    const TLorentzVector missP4 = targP4 + qP4 - leadP4;
    const double mMissVal = missP4.M();
    const double EMissVal = std::sqrt(leadVec.Mag2() + leadMass * leadMass) - nu;
    const double E0missVal = std::sqrt(pMissVec.Mag2() + recMass * recMass) - recMass;

    const double TP = leadP4.E() - leadMass;
    const double TP2 = recP4.E() - recMass;
    const TLorentzVector missA2 = qP4 + nucleusP4 - leadP4 - recP4;
    const double TB2 = missA2.E() - missA2.M();
    const double E2missVal = nu - TP - TP2 - TB2;
    const double E1missVal = nu - TP - (missP4.E() - missP4.M());

    const TVector3 pRelVec = 0.5 * (pMissVec - recVec);
    const TVector3 pCMVec = pMissVec + recVec;

    TVector3 vz(0., 0., 1.);
    TVector3 vy(0., 1., 0.);
    TVector3 vx(1., 0., 0.);
    if (pMissVec.Mag2() > 0. && qVec.Mag2() > 0.) {
      vz = pMissVec.Unit();
      vy = pMissVec.Cross(qVec);
      if (vy.Mag2() > 1e-16) {
        vy = vy.Unit();
        vx = vz.Cross(vy).Unit();
      }
    }

    double omegaVal = nu;
    double q2Val = qVec.Mag2() - nu * nu;
    double xBVal = (std::abs(omegaVal) > 1e-12) ? q2Val / (2.0 * kMp * omegaVal) : -9.;

    // If requested, replace inferred beam-energy dependent quantities.
    if (eBeamOverride > 0.) {
      const double omegaFix = eBeamOverride - electron[3];
      const double q2Fix = qVec.Mag2() - omegaFix * omegaFix;
      const double xBFix =
          (std::abs(omegaFix) > 1e-12) ? q2Fix / (2.0 * kMp * omegaFix) : -9.;
      omegaVal = omegaFix;
      q2Val = q2Fix;
      xBVal = xBFix;
    }

    if (applyBasicSrcCuts) {
      if (xBVal < minXB) continue;
      ++passXB;
      if (q2Val < minQ2) continue;
      ++passQ2;
      if (pLeadMag < minLeadP) continue;
      ++passLeadP;
      if (pRecMag < minRecP) continue;
      ++passRecP;
      if (mMissVal < minMMiss || mMissVal > maxMMiss) continue;
      ++passMMiss;
      if (kMissVal < minKMiss || kMissVal > maxKMiss) continue;
      ++passKMiss;
    }

    weight_out = static_cast<float>(weight);
    weight_ep = static_cast<float>(weight);
    weight_epp = static_cast<float>(weight);

    xB = static_cast<float>(xBVal);
    Q2 = static_cast<float>(q2Val);
    omega = static_cast<float>(omegaVal);

    eP = static_cast<float>(eVec.Mag());
    eTheta = static_cast<float>(eVec.Theta());
    ePhi = static_cast<float>(eVec.Phi());

    qP = static_cast<float>(qVec.Mag());
    qTheta = static_cast<float>(qVec.Theta());
    qPhi = static_cast<float>(qVec.Phi());

    nProtons = 2;

    leadP = static_cast<float>(pLeadMag);
    leadTheta = static_cast<float>(leadVec.Theta());
    leadPhi = static_cast<float>(leadVec.Phi());
    recP = static_cast<float>(pRecMag);
    recTheta = static_cast<float>(recVec.Theta());
    recPhi = static_cast<float>(recVec.Phi());

    pMiss = static_cast<float>(pMissMag);
    pMissTheta = static_cast<float>(pMissVec.Theta());
    pMissPhi = static_cast<float>(pMissVec.Phi());

    kMiss = static_cast<float>(kMissVal);
    mMiss = static_cast<float>(mMissVal);
    EMiss = static_cast<float>(EMissVal);
    E0miss = static_cast<float>(E0missVal);
    E1miss = static_cast<float>(E1missVal);
    E2miss = static_cast<float>(E2missVal);

    theta_PmQ = static_cast<float>(pMissVec.Angle(qVec));
    theta_PmPlead = static_cast<float>(pMissVec.Angle(leadVec));
    theta_PleadQ = static_cast<float>(leadVec.Angle(qVec));
    theta_PleadPrec = static_cast<float>(leadVec.Angle(recVec));
    theta_PmPrec = static_cast<float>(pMissVec.Angle(recVec));
    theta_PrecQ = static_cast<float>(recVec.Angle(qVec));

    pRel = static_cast<float>(pRelVec.Mag());
    pRelTheta = static_cast<float>(pRelVec.Theta());
    pRelPhi = static_cast<float>(pRelVec.Phi());
    pRelx = static_cast<float>(pRelVec.Dot(vx));
    pRely = static_cast<float>(pRelVec.Dot(vy));
    pRelz = static_cast<float>(pRelVec.Dot(vz));

    pLeadPlusRecOver2 = static_cast<float>(0.5 * (leadVec + recVec).Mag());
    pLeadMinusRec = static_cast<float>((leadVec - recVec).Mag());

    pCM = static_cast<float>(pCMVec.Mag());
    pCMx = static_cast<float>(pCMVec.Dot(vx));
    pCMy = static_cast<float>(pCMVec.Dot(vy));
    pCMz = static_cast<float>(pCMVec.Dot(vz));

    pcmx_lab = static_cast<float>(pCMVec.X());
    pcmy_lab = static_cast<float>(pCMVec.Y());
    pcmz_lab = static_cast<float>(pCMVec.Z());

    resetFloat(alpha_1);
    resetFloat(alpha_2);
    resetFloat(alpha_CM);
    resetFloat(alpha_rel);
    resetFloat(alpha_q);
    resetFloat(alpha_pLead);
    resetFloat(p1_plus);
    resetFloat(p2_plus);
    resetFloat(p1_perp_x);
    resetFloat(p1_perp_y);
    resetFloat(p1_perp_mag);
    resetFloat(p2_perp_x);
    resetFloat(p2_perp_y);
    resetFloat(p2_perp_mag);
    resetFloat(pCM_perp_x);
    resetFloat(pCM_perp_y);
    resetFloat(pCM_perp_mag);
    resetFloat(prel_perp_x);
    resetFloat(prel_perp_y);
    resetFloat(prel_perp_mag);
    resetFloat(k);
    resetFloat(k2);
    resetFloat(k_z);
    resetFloat(m_bar);

    if (qVec.Mag2() > 0.) {
      const TVector3 qhat = qVec.Unit();
      alpha_q = static_cast<float>((qP4.E() - qVec.Dot(qhat)) / kMn);
      alpha_pLead = static_cast<float>((leadP4.E() - leadVec.Dot(qhat)) / kMn);
    }

    const LightConeKinematics lc = computeLightConeKinematics(
        qVec, nu, leadVec, recVec, leadMass, recMass, targetCfg.mass, targetCfg.A);
    if (lc.validBasis && lc.pairDefined) {
      alpha_1 = static_cast<float>(lc.alpha1);
      alpha_2 = static_cast<float>(lc.alpha2);
      alpha_CM = static_cast<float>(lc.alphaCM);
      alpha_rel = static_cast<float>(lc.alphaRel);
      p1_plus = static_cast<float>(lc.p1Plus);
      p2_plus = static_cast<float>(lc.p2Plus);
      p1_perp_x = static_cast<float>(lc.p1PerpX);
      p1_perp_y = static_cast<float>(lc.p1PerpY);
      p1_perp_mag = static_cast<float>(lc.p1PerpMag);
      p2_perp_x = static_cast<float>(lc.p2PerpX);
      p2_perp_y = static_cast<float>(lc.p2PerpY);
      p2_perp_mag = static_cast<float>(lc.p2PerpMag);
      pCM_perp_x = static_cast<float>(lc.pCMPerpX);
      pCM_perp_y = static_cast<float>(lc.pCMPerpY);
      pCM_perp_mag = static_cast<float>(lc.pCMPerpMag);
      prel_perp_x = static_cast<float>(lc.pRelPerpX);
      prel_perp_y = static_cast<float>(lc.pRelPerpY);
      prel_perp_mag = static_cast<float>(lc.pRelPerpMag);
      k = static_cast<float>(lc.k);
      k2 = static_cast<float>(lc.k2);
      k_z = static_cast<float>(lc.kZ);
      m_bar = static_cast<float>(lc.mBar);
    }

    srcTree->Fill();
    ++nWritten;
  }

  outFile->cd();
  srcTree->Write();
  outFile->Close();
  inFile->Close();

  std::cout << "[convert_events2N_to_srcTree] Wrote srcTree with " << nWritten
            << " entries to " << outputFileName << "\n";
  if (applyBasicSrcCuts) {
    std::cout << "[convert_events2N_to_srcTree] Cutflow\n"
              << "  total      : " << nTotal << "\n"
              << "  e'pp PID   : " << (requireEpp ? passEpp : nTotal) << "\n"
              << "  xB         : " << passXB << "\n"
              << "  Q2         : " << passQ2 << "\n"
              << "  leadP      : " << passLeadP << "\n"
              << "  recP       : " << passRecP << "\n"
              << "  mMiss      : " << passMMiss << "\n"
              << "  kMiss      : " << passKMiss << "\n"
              << "  written    : " << nWritten << "\n";
  } else if (requireEpp) {
    std::cout << "[convert_events2N_to_srcTree] Cutflow\n"
              << "  total      : " << nTotal << "\n"
              << "  e'pp PID   : " << passEpp << "\n"
              << "  written    : " << nWritten << "\n";
  }
}
