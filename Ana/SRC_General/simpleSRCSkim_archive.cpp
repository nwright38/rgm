#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>

#include "clas12reader.h"
#include "clas12ana.h"
#include "HipoChain.h"
#include "Corrections.h"
#include "reweighter.h"

using namespace std;
using namespace clas12;

const double mP = 0.938;
// const double mN = 0.939565;
const double mU = 0.9314941024;
const double me = 0.000511;
const double m_4He = 4.00260325415 * mU - 2 * me;

const int MAXP = 4;  // max number of proton candidates considered per event

void Usage()
{
  cerr << "Usage: ./code <MC=1,Data=0> <Ebeam(GeV)> <output.root> <input.hipo> [more hipo...]\n";
}

bool inFD(int status){ return (abs(status) >= 2000 && abs(status) < 4000); }
bool inCD(int status){ return (abs(status) >= 4000 && abs(status) < 8000); }

int main(int argc, char **argv)
{
  if(argc < 4){ Usage(); return -1; }

  bool isMC = (atoi(argv[1]) == 1);
  double Ebeam = atof(argv[2]);

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
  Float_t b_leadVz;
  Float_t b_pMiss, b_pMissTheta, b_pMissPhi;
  Float_t b_mMiss;
  Float_t b_kMiss;
  Float_t b_EMiss;
  Float_t b_theta_PmQ;    // angle between pMiss and q
  Float_t b_theta_PleadQ; // angle between lead proton momentum and q
  Bool_t  b_goodLead;     // lead passes all SRC lead cuts

  // pRel and pCM (filled only when a recoil partner exists)
  Float_t b_pRel, b_pRelTheta, b_pRelPhi;
  Float_t b_pCM, b_pCMx, b_pCMy, b_pCMz;
  Float_t b_E2miss;

  // event-level summary
  Int_t  b_nGoodLeads;
  Bool_t b_singleGoodLead;      // exactly one proton passes all SRC lead cuts

  // recoil kinematics (filled when a recoil is found)
  Float_t b_recP, b_recTheta, b_recPhi;
  Float_t b_theta_PmPrec; // angle between pMiss and recoil momentum
  Float_t b_theta_PrecQ;  // angle between recoil momentum and q

  // truth-level MC quantities, named to mirror reco branches
  Float_t b_xB_truth, b_Q2_truth, b_omega_truth;
  Float_t b_eP_truth, b_eTheta_truth, b_ePhi_truth;
  Float_t b_qP_truth, b_qTheta_truth, b_qPhi_truth;
  Float_t b_leadP_truth, b_leadTheta_truth, b_leadPhi_truth;
  Float_t b_pMiss_truth, b_pMissTheta_truth, b_pMissPhi_truth;
  Float_t b_mMiss_truth;
  Float_t b_kMiss_truth;
  Float_t b_EMiss_truth;
  Float_t b_theta_PmQ_truth;
  Float_t b_theta_PleadQ_truth;
  Float_t b_pRel_truth, b_pRelTheta_truth, b_pRelPhi_truth;
  Float_t b_pCM_truth, b_pCMx_truth, b_pCMy_truth, b_pCMz_truth;
  Float_t b_E2miss_truth;
  Float_t b_recP_truth, b_recTheta_truth, b_recPhi_truth;
  Float_t b_theta_PmPrec_truth;
  Float_t b_theta_PrecQ_truth;

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
  srcTree->Branch("leadVz",      &b_leadVz,      "leadVz/F");

  srcTree->Branch("pMiss",       &b_pMiss,       "pMiss/F");
  srcTree->Branch("pMissTheta",  &b_pMissTheta,  "pMissTheta/F");
  srcTree->Branch("pMissPhi",    &b_pMissPhi,    "pMissPhi/F");

  srcTree->Branch("mMiss",       &b_mMiss,       "mMiss/F");
  srcTree->Branch("kMiss",       &b_kMiss,       "kMiss/F");
  srcTree->Branch("EMiss",       &b_EMiss,       "EMiss/F");
  srcTree->Branch("theta_PmQ",   &b_theta_PmQ,   "theta_PmQ/F");
  srcTree->Branch("theta_PleadQ",&b_theta_PleadQ,"theta_PleadQ/F");
  srcTree->Branch("goodLead",    &b_goodLead,    "goodLead/O");

  srcTree->Branch("pRel",        &b_pRel,        "pRel/F");
  srcTree->Branch("pRelTheta",   &b_pRelTheta,   "pRelTheta/F");
  srcTree->Branch("pRelPhi",     &b_pRelPhi,     "pRelPhi/F");
  srcTree->Branch("pCM",         &b_pCM,         "pCM/F");
  srcTree->Branch("pCMx",        &b_pCMx,        "pCMx/F");
  srcTree->Branch("pCMy",        &b_pCMy,        "pCMy/F");
  srcTree->Branch("pCMz",        &b_pCMz,        "pCMz/F");
  srcTree->Branch("E2miss",      &b_E2miss,      "E2miss/F");

  srcTree->Branch("nGoodLeads",     &b_nGoodLeads,     "nGoodLeads/I");
  srcTree->Branch("singleGoodLead", &b_singleGoodLead, "singleGoodLead/O");

  srcTree->Branch("recP",        &b_recP,        "recP/F");
  srcTree->Branch("recTheta",    &b_recTheta,    "recTheta/F");
  srcTree->Branch("recPhi",      &b_recPhi,      "recPhi/F");
  srcTree->Branch("theta_PmPrec",&b_theta_PmPrec,"theta_PmPrec/F");
  srcTree->Branch("theta_PrecQ", &b_theta_PrecQ, "theta_PrecQ/F");

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
  srcTree->Branch("theta_PmQ_truth",   &b_theta_PmQ_truth,   "theta_PmQ_truth/F");
  srcTree->Branch("theta_PleadQ_truth",&b_theta_PleadQ_truth,"theta_PleadQ_truth/F");

  srcTree->Branch("pRel_truth",        &b_pRel_truth,        "pRel_truth/F");
  srcTree->Branch("pRelTheta_truth",   &b_pRelTheta_truth,   "pRelTheta_truth/F");
  srcTree->Branch("pRelPhi_truth",     &b_pRelPhi_truth,     "pRelPhi_truth/F");
  srcTree->Branch("pCM_truth",         &b_pCM_truth,         "pCM_truth/F");
  srcTree->Branch("pCMx_truth",        &b_pCMx_truth,        "pCMx_truth/F");
  srcTree->Branch("pCMy_truth",        &b_pCMy_truth,        "pCMy_truth/F");
  srcTree->Branch("pCMz_truth",        &b_pCMz_truth,        "pCMz_truth/F");
  srcTree->Branch("E2miss_truth",      &b_E2miss_truth,      "E2miss_truth/F");

  srcTree->Branch("recP_truth",        &b_recP_truth,        "recP_truth/F");
  srcTree->Branch("recTheta_truth",    &b_recTheta_truth,    "recTheta_truth/F");
  srcTree->Branch("recPhi_truth",      &b_recPhi_truth,      "recPhi_truth/F");
  srcTree->Branch("theta_PmPrec_truth",&b_theta_PmPrec_truth,"theta_PmPrec_truth/F");
  srcTree->Branch("theta_PrecQ_truth", &b_theta_PrecQ_truth, "theta_PrecQ_truth/F");

  // ---- chain setup ----
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout << "Input file " << argv[k] << endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  clas12ana clasAna;

  // fixed 4-vectors
  const double md = 2.01410178 * mU - me;
  TLorentzVector targP4(0., 0., 0., md);
  TLorentzVector nucleusP4(0., 0., 0., m_4He);
  TLorentzVector beamP4(0., 0., Ebeam, Ebeam);

  TLorentzVector eP4(0., 0., 0., me);
  TLorentzVector leadP4(0., 0., 0., mP);
  TLorentzVector missP4;

  TVector3 p_beam(0., 0., Ebeam);

  int counter = 0;

  reweighter newWeight(Ebeam,2,2,kelly,"AV18",.15);

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

    b_leadP = -9.f;  b_leadTheta = -9.f;  b_leadPhi = -9.f;  b_leadVz = -99.f;
    b_pMiss = -9.f;  b_pMissTheta = -9.f;  b_pMissPhi = -9.f;
    b_pRel = -9.f;   b_pRelTheta = -9.f;   b_pRelPhi = -9.f;
    b_pCM = -9.f;    b_pCMx = -9.f;        b_pCMy = -9.f;    b_pCMz = -9.f;
    b_E2miss = -9.f;
    b_mMiss = -9.f;  b_kMiss = -9.f;       b_EMiss = -9.f;
    b_theta_PmQ = -9.f;
    b_theta_PleadQ = -9.f;
    b_theta_PmPrec = -9.f;
    b_theta_PrecQ  = -9.f;
    b_goodLead = false;

    b_xB_truth = -9.f;  b_Q2_truth = -9.f;  b_omega_truth = -9.f;
    b_eP_truth = -9.f;  b_eTheta_truth = -9.f;  b_ePhi_truth = -9.f;
    b_qP_truth = -9.f;  b_qTheta_truth = -9.f;  b_qPhi_truth = -9.f;
    b_leadP_truth = -9.f;  b_leadTheta_truth = -9.f;  b_leadPhi_truth = -9.f;
    b_pMiss_truth = -9.f;  b_pMissTheta_truth = -9.f;  b_pMissPhi_truth = -9.f;
    b_mMiss_truth = -9.f;  b_kMiss_truth = -9.f;  b_EMiss_truth = -9.f;
    b_theta_PmQ_truth = -9.f;
    b_theta_PleadQ_truth = -9.f;
    b_pRel_truth = -9.f;   b_pRelTheta_truth = -9.f;   b_pRelPhi_truth = -9.f;
    b_pCM_truth = -9.f;    b_pCMx_truth = -9.f;        b_pCMy_truth = -9.f;    b_pCMz_truth = -9.f;
    b_E2miss_truth = -9.f;
    b_recP_truth = -9.f;   b_recTheta_truth = -9.f;    b_recPhi_truth = -9.f;
    b_theta_PmPrec_truth = -9.f;
    b_theta_PrecQ_truth  = -9.f;

    clasAna.Run(c12);

    auto electrons = clasAna.getByPid(11);
    auto protons   = clasAna.getByPid(2212);

    if(electrons.size() != 1) continue;
    if(protons.empty())       continue;
    if(!inFD(electrons[0]->getStatus())) continue;

    // ---- electron kinematics ----
    eP4.SetXYZM(0., 0., 0., me);
    GetLorentzVector_Corrected(eP4, electrons[0], 0);
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
    vector<float>    cand_vz;
    vector<TVector3> cand_pMissV;
    vector<float>    cand_mMiss, cand_kMiss, cand_EMiss, cand_theta_PmQ;
    vector<bool>     cand_goodLead;

    for(int pr = 0; pr < (int)protons.size(); pr++)
    {
      if(nFilled >= MAXP) break;

      GetLorentzVector_Corrected(leadP4, protons[pr], 0);
      TVector3 pLead3 = leadP4.Vect();

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

      // SRC lead cuts
      bool passCuts = true;
      if(pLead3.Mag() < 1.)                        passCuts = false;
      if(missP4.M() < 0.65 || missP4.M() > 1.1)   passCuts = false;
      if(kMiss < 0.3 || kMiss > 1.)                passCuts = false;
 //     if(pLead3.Angle(qP3) < 37.*M_PI/180.)             passCuts = false;

      cand_p3.push_back(pLead3);
      cand_vz.push_back(protons[pr]->par()->getVz());
      cand_pMissV.push_back(pMissV);
      cand_mMiss.push_back(missP4.M());
      cand_kMiss.push_back(kMiss);
      cand_EMiss.push_back(EMiss);
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
      b_leadVz      = cand_vz[leadIdx];

      b_pMiss       = cand_pMissV[leadIdx].Mag();
      b_pMissTheta  = cand_pMissV[leadIdx].Theta();
      b_pMissPhi    = cand_pMissV[leadIdx].Phi();
      b_mMiss       = cand_mMiss[leadIdx];
      b_kMiss       = cand_kMiss[leadIdx];
      b_EMiss       = cand_EMiss[leadIdx];
      b_theta_PmQ   = cand_theta_PmQ[leadIdx];
      b_goodLead    = cand_goodLead[leadIdx];
      b_theta_PleadQ = cand_p3[leadIdx].Angle(qP3);
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

        // ---- pCM: pair CM momentum projected onto the (miss_neg, q) frame ----
        TVector3 v_rec = recoil_p3;
        TVector3 v_cm  = miss_neg + v_rec;

        TVector3 vz = miss_neg.Unit();
        TVector3 vy = miss_neg.Cross(qP3).Unit();
        TVector3 vx = vz.Cross(vy).Unit();

        b_pCM  = v_cm.Mag();
        b_pCMx = v_cm.Dot(vx);
        b_pCMy = v_cm.Dot(vy);
        b_pCMz = v_cm.Dot(vz);

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

        // additional angles
        b_theta_PmPrec = miss_neg.Angle(recoil_p3);
        b_theta_PrecQ  = recoil_p3.Angle(qP3);

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
        b_theta_PmQ_truth  = pMiss_truth.Angle(q_truth);
        b_theta_PleadQ_truth = lead_truth.Angle(q_truth);

        if(mcInfo->getRows() >= 3)
        {
          TVector3 rec_truth(mcInfo->getPx(2), mcInfo->getPy(2), mcInfo->getPz(2));
          TVector3 pRel_truth = (pMiss_truth - rec_truth) * 0.5;
          TVector3 pCM_truth = pMiss_truth + rec_truth;

          TVector3 vz_truth = pMiss_truth.Unit();
          TVector3 vy_truth = pMiss_truth.Cross(q_truth).Unit();
          TVector3 vx_truth = vz_truth.Cross(vy_truth).Unit();

          b_pRel_truth      = pRel_truth.Mag();
          b_pRelTheta_truth = pRel_truth.Theta();
          b_pRelPhi_truth   = pRel_truth.Phi();

          b_pCM_truth  = pCM_truth.Mag();
          b_pCMx_truth = pCM_truth.Dot(vx_truth);
          b_pCMy_truth = pCM_truth.Dot(vy_truth);
          b_pCMz_truth = pCM_truth.Dot(vz_truth);

          TLorentzVector recP4_truth;
          recP4_truth.SetVectM(rec_truth, mP);
          double TP_truth = leadP4_truth.E() - leadP4_truth.M();
          double TP2_truth = recP4_truth.E() - recP4_truth.M();
          TLorentzVector miss_Am2_truth = q4_truth + nucleusP4 - leadP4_truth - recP4_truth;
          double TB2_truth = miss_Am2_truth.E() - miss_Am2_truth.M();
          b_E2miss_truth = q4_truth.E() - TP_truth - TP2_truth - TB2_truth;

          b_recP_truth     = rec_truth.Mag();
          b_recTheta_truth = rec_truth.Theta();
          b_recPhi_truth   = rec_truth.Phi();
          b_theta_PmPrec_truth = pMiss_truth.Angle(rec_truth);
          b_theta_PrecQ_truth  = rec_truth.Angle(q_truth);
        }
      }
    }

    srcTree->Fill();
    ctr++;
  }

  outFile->cd();
  srcTree->Write();
  outFile->Close();

  cout << "Done. Processed " << counter << " events.\n";
  return 0;
}
