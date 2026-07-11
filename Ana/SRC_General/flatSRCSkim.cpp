#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
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
const double mN = 0.939565;
const double mU = 0.9314941024;
const double me = 0.000511;
const double m_4He = 4.00260325415 * mU - 2 * me;

const int MAXP = 4;  // max number of proton candidates considered per event
const double rad2deg = 180. / M_PI;

string getCommandOutput(const char *cmd)
{
  string out;
  FILE *pipe = popen(cmd, "r");
  if(!pipe) return "unknown";
  char buffer[256];
  while(fgets(buffer, sizeof(buffer), pipe)) out += buffer;
  pclose(pipe);
  while(!out.empty() && (out.back() == '\n' || out.back() == ' ')) out.pop_back();
  return out.empty() ? "unknown" : out;
}

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

  Int_t   b_run;
  Long64_t b_event;

  Double_t b_xB, b_Q2, b_omega;

  // electron
  Double_t b_eP, b_eTheta, b_ePhi, b_vtzE;
  Int_t    b_sectorE;

  // q vector (3-momentum transfer)
  Double_t b_qP, b_qTheta, b_qPhi;

  Int_t   b_nProtons;

  // lead proton (first candidate passing all SRC lead cuts)
  Double_t b_leadP, b_leadTheta, b_leadPhi;
  Double_t b_leadVz;
  Int_t    b_leadRegion;
  Double_t b_pMiss, b_pMissTheta, b_pMissPhi;
  Double_t b_mMiss;
  Double_t b_kMiss;
  Double_t b_EMiss, b_E0miss, b_E1miss;
  Double_t b_theta_PmQ;    // angle between pMiss and q
  Double_t b_theta_PleadQ; // angle between lead proton momentum and q
  Bool_t  b_goodLead;     // lead passes all SRC lead cuts

  // pRel and pCM (filled only when a recoil partner exists)
  Double_t b_pRel, b_pRelTheta, b_pRelPhi;
  Double_t b_pCM, b_pCMx, b_pCMy, b_pCMz;
  Double_t b_E2miss;

  // event-level summary
  Int_t  b_nGoodLeads;
  Bool_t b_singleGoodLead;      // exactly one proton passes all SRC lead cuts

  // recoil kinematics (filled when a recoil is found)
  Bool_t   b_hasRecoil;
  Double_t b_recP, b_recTheta, b_recPhi, b_recVz;
  Double_t b_theta_PmPrec; // angle between pMiss and recoil momentum
  Double_t b_theta_PrecQ;  // angle between recoil momentum and q

  Bool_t   b_isMC;
  Double_t b_genWeight;
  Double_t b_pcmX_truth, b_pcmY_truth, b_pcmZ_truth;
  Double_t b_pcmX_truth_lab, b_pcmY_truth_lab, b_pcmZ_truth_lab;
  Double_t b_Q2_truth, b_pMiss_truth, b_pLead_truth, b_thetaLead_truth, b_pRec_truth;
  Double_t b_vtzDiffLeadE;

  srcTree->Branch("weight",      &b_weight,      "weight/F");
  srcTree->Branch("weight_ep",   &b_weight_ep,   "weight_ep/F");
  srcTree->Branch("weight_epp",  &b_weight_epp,  "weight_epp/F");

  srcTree->Branch("run",         &b_run,         "run/I");
  srcTree->Branch("event",       &b_event,       "event/L");

  srcTree->Branch("xB",          &b_xB,          "xB/D");
  srcTree->Branch("Q2",          &b_Q2,          "Q2/D");
  srcTree->Branch("omega",       &b_omega,       "omega/D");

  srcTree->Branch("eP",          &b_eP,          "eP/D");
  srcTree->Branch("eTheta",      &b_eTheta,      "eTheta/D");
  srcTree->Branch("ePhi",        &b_ePhi,        "ePhi/D");
  srcTree->Branch("pE",          &b_eP,          "pE/D");
  srcTree->Branch("thetaE",      &b_eTheta,      "thetaE/D");
  srcTree->Branch("phiE",        &b_ePhi,        "phiE/D");
  srcTree->Branch("vtzE",        &b_vtzE,        "vtzE/D");
  srcTree->Branch("sectorE",     &b_sectorE,     "sectorE/I");

  srcTree->Branch("qP",          &b_qP,          "qP/D");
  srcTree->Branch("qTheta",      &b_qTheta,      "qTheta/D");
  srcTree->Branch("qPhi",        &b_qPhi,        "qPhi/D");
  srcTree->Branch("q_mag",       &b_qP,          "q_mag/D");

  srcTree->Branch("nProtons",    &b_nProtons,    "nProtons/I");

  srcTree->Branch("leadP",       &b_leadP,       "leadP/D");
  srcTree->Branch("leadTheta",   &b_leadTheta,   "leadTheta/D");
  srcTree->Branch("leadPhi",     &b_leadPhi,     "leadPhi/D");
  srcTree->Branch("leadVz",      &b_leadVz,      "leadVz/D");
  srcTree->Branch("pLead",       &b_leadP,       "pLead/D");
  srcTree->Branch("thetaLead",   &b_leadTheta,   "thetaLead/D");
  srcTree->Branch("phiLead",     &b_leadPhi,     "phiLead/D");
  srcTree->Branch("vtzLead",     &b_leadVz,      "vtzLead/D");
  srcTree->Branch("leadRegion",  &b_leadRegion,  "leadRegion/I");

  srcTree->Branch("pMiss",       &b_pMiss,       "pMiss/D");
  srcTree->Branch("pMissTheta",  &b_pMissTheta,  "pMissTheta/D");
  srcTree->Branch("pMissPhi",    &b_pMissPhi,    "pMissPhi/D");

  srcTree->Branch("mMiss",       &b_mMiss,       "mMiss/D");
  srcTree->Branch("kMiss",       &b_kMiss,       "kMiss/D");
  srcTree->Branch("EMiss",       &b_EMiss,       "EMiss/D");
  srcTree->Branch("E0miss",      &b_E0miss,      "E0miss/D");
  srcTree->Branch("E1miss",      &b_E1miss,      "E1miss/D");
  srcTree->Branch("theta_PmQ",   &b_theta_PmQ,   "theta_PmQ/D");
  srcTree->Branch("theta_PleadQ",&b_theta_PleadQ,"theta_PleadQ/D");
  srcTree->Branch("thetaPmissQ", &b_theta_PmQ,   "thetaPmissQ/D");
  srcTree->Branch("thetaPLeadQ", &b_theta_PleadQ,"thetaPLeadQ/D");
  srcTree->Branch("goodLead",    &b_goodLead,    "goodLead/O");

  srcTree->Branch("pRel",        &b_pRel,        "pRel/D");
  srcTree->Branch("pRelTheta",   &b_pRelTheta,   "pRelTheta/D");
  srcTree->Branch("pRelPhi",     &b_pRelPhi,     "pRelPhi/D");
  srcTree->Branch("pCM",         &b_pCM,         "pCM/D");
  srcTree->Branch("pCMx",        &b_pCMx,        "pCMx/D");
  srcTree->Branch("pCMy",        &b_pCMy,        "pCMy/D");
  srcTree->Branch("pCMz",        &b_pCMz,        "pCMz/D");
  srcTree->Branch("pcmX",        &b_pCMx,        "pcmX/D");
  srcTree->Branch("pcmY",        &b_pCMy,        "pcmY/D");
  srcTree->Branch("pcmZ",        &b_pCMz,        "pcmZ/D");
  srcTree->Branch("E2miss",      &b_E2miss,      "E2miss/D");

  srcTree->Branch("nGoodLeads",     &b_nGoodLeads,     "nGoodLeads/I");
  srcTree->Branch("singleGoodLead", &b_singleGoodLead, "singleGoodLead/O");

  srcTree->Branch("hasRecoil",   &b_hasRecoil,   "hasRecoil/O");
  srcTree->Branch("recP",        &b_recP,        "recP/D");
  srcTree->Branch("recTheta",    &b_recTheta,    "recTheta/D");
  srcTree->Branch("recPhi",      &b_recPhi,      "recPhi/D");
  srcTree->Branch("pRec",        &b_recP,        "pRec/D");
  srcTree->Branch("thetaRec",    &b_recTheta,    "thetaRec/D");
  srcTree->Branch("phiRec",      &b_recPhi,      "phiRec/D");
  srcTree->Branch("vtzRec",      &b_recVz,       "vtzRec/D");
  srcTree->Branch("theta_PmPrec",&b_theta_PmPrec,"theta_PmPrec/D");
  srcTree->Branch("theta_PrecQ", &b_theta_PrecQ, "theta_PrecQ/D");
  srcTree->Branch("thetaPRecPmiss",&b_theta_PmPrec,"thetaPRecPmiss/D");
  srcTree->Branch("isMC",        &b_isMC,        "isMC/O");
  srcTree->Branch("genWeight",   &b_genWeight,   "genWeight/D");
  srcTree->Branch("pcmX_truth",  &b_pcmX_truth,  "pcmX_truth/D");
  srcTree->Branch("pcmY_truth",  &b_pcmY_truth,  "pcmY_truth/D");
  srcTree->Branch("pcmZ_truth",  &b_pcmZ_truth,  "pcmZ_truth/D");
  srcTree->Branch("pcmX_truth_lab",&b_pcmX_truth_lab,"pcmX_truth_lab/D");
  srcTree->Branch("pcmY_truth_lab",&b_pcmY_truth_lab,"pcmY_truth_lab/D");
  srcTree->Branch("pcmZ_truth_lab",&b_pcmZ_truth_lab,"pcmZ_truth_lab/D");
  srcTree->Branch("Q2_truth",    &b_Q2_truth,    "Q2_truth/D");
  srcTree->Branch("pMiss_truth", &b_pMiss_truth, "pMiss_truth/D");
  srcTree->Branch("pLead_truth", &b_pLead_truth, "pLead_truth/D");
  srcTree->Branch("thetaLead_truth",&b_thetaLead_truth,"thetaLead_truth/D");
  srcTree->Branch("pRec_truth",  &b_pRec_truth,  "pRec_truth/D");
  srcTree->Branch("vtzDiffLeadE",&b_vtzDiffLeadE,"vtzDiffLeadE/D");

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
  long long cut_read = 0;
  long long cut_one_electron = 0;
  long long cut_lead_found = 0;
  long long cut_kinematics = 0;
  long long cut_lead_region = 0;
  long long cut_written = 0;

  reweighter newWeight(Ebeam,2,2,kelly,"AV18",.15);

  int ctr = 0;
  while(chain.Next())
  {
    if(counter % 100000 == 0)
      cout << "Processing event " << counter << "\t" << ctr << " saved." << endl;
    counter++;
    cut_read++;

    // ---- defaults ----
    b_weight = 1.f;
    b_weight_ep = 1.f;
    b_weight_epp = 1.f;
    b_run = c12->runconfig()->getRun();
    b_event = c12->runconfig()->getEvent();
    b_xB = -9.;  b_Q2 = -9.;  b_omega = -9.;
    b_eP = -9.;  b_eTheta = -9.;  b_ePhi = -9.; b_vtzE = -99.;
    b_sectorE = -9;
    b_qP = -9.;  b_qTheta = -9.;  b_qPhi = -9.;

    b_nProtons    = 0;
    b_nGoodLeads  = 0;
    b_singleGoodLead = false;
    b_hasRecoil = false;
    b_recP       = -9.;  b_recTheta = -9.;  b_recPhi = -9.; b_recVz = -99.;

    b_leadP = -9.;  b_leadTheta = -9.;  b_leadPhi = -9.;  b_leadVz = -99.;
    b_leadRegion = -9;
    b_pMiss = -9.;  b_pMissTheta = -9.;  b_pMissPhi = -9.;
    b_pRel = -9.;   b_pRelTheta = -9.;   b_pRelPhi = -9.;
    b_pCM = -9.;    b_pCMx = -9.;        b_pCMy = -9.;    b_pCMz = -9.;
    b_E2miss = -9.;
    b_mMiss = -9.;  b_kMiss = -9.;       b_EMiss = -9.; b_E0miss = -9.; b_E1miss = -9.;
    b_theta_PmQ = -9.;
    b_theta_PleadQ = -9.;
    b_theta_PmPrec = -9.;
    b_theta_PrecQ  = -9.;
    b_goodLead = false;
    auto mcInfoForEvent = c12->mcparts();
    bool eventIsMC = isMC || (mcInfoForEvent && mcInfoForEvent->getRows() > 0);
    b_isMC = eventIsMC;
    b_genWeight = eventIsMC ? 1. : 0.;
    b_pcmX_truth = b_pcmY_truth = b_pcmZ_truth = -9.;
    b_pcmX_truth_lab = b_pcmY_truth_lab = b_pcmZ_truth_lab = -9.;
    b_Q2_truth = b_pMiss_truth = b_pLead_truth = b_thetaLead_truth = b_pRec_truth = -9.;
    b_vtzDiffLeadE = -99.;

    clasAna.Run(c12);

    auto electrons = clasAna.getByPid(11);
    auto protons   = clasAna.getByPid(2212);

    if(electrons.size() != 1) continue;
    cut_one_electron++;

    // ---- electron kinematics ----
    eP4.SetXYZM(0., 0., 0., me);
    GetLorentzVector_Corrected(eP4, electrons[0], eventIsMC);
    if(eP4.Vect().Mag() == 0.) continue;

    TVector3 eP3 = eP4.Vect();
    TVector3 qP3 = p_beam - eP3;

    double omega = Ebeam - eP4.E();
    double Q2    = qP3.Mag2() - omega * omega;
    double xB    = Q2 / (2. * mP * omega);

    if(eventIsMC){
      b_weight     = c12->mcevent()->getWeight();
      b_weight_ep  = b_weight * newWeight.get_weight_ep(c12->mcparts());
      b_weight_epp = b_weight * newWeight.get_weight_epp(c12->mcparts());
      b_genWeight  = b_weight;
    }

    b_eP     = eP3.Mag();
    b_eTheta = eP3.Theta() * rad2deg;
    b_ePhi   = eP3.Phi() * rad2deg;
    b_vtzE   = electrons[0]->par()->getVz();
    b_sectorE = electrons[0]->getSector();

    b_qP     = qP3.Mag();
    b_qTheta = qP3.Theta() * rad2deg;
    b_qPhi   = qP3.Phi() * rad2deg;

    b_Q2    = Q2;
    b_xB    = xB;
    b_omega = omega;

    TLorentzVector q(qP3, omega);

    clasAna.getLeadRecoilSRC(beamP4, targP4, eP4);
    auto lead   = clasAna.getLeadSRC();
    auto recoil = clasAna.getRecoilSRC();
    b_nProtons = protons.size() > MAXP ? MAXP : protons.size();
    b_nGoodLeads = lead.size();
    b_singleGoodLead = (lead.size() == 1);
    if(lead.empty()) continue;
    cut_lead_found++;

    GetLorentzVector_Corrected(leadP4, lead[0], eventIsMC);
    TVector3 pLead3 = leadP4.Vect();
    missP4 = targP4 + beamP4 - eP4 - leadP4;
    TVector3 miss_neg = pLead3 - qP3;

    TLorentzVector miss_LC = leadP4 - q;
    TVector3 u_ZQ    = q.Vect().Unit();
    double pmm_ZQ    = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
    double pmp_ZQ    = miss_LC.Vect().Perp(u_ZQ);
    double kMiss     = sqrt(mP*mP*((pmp_ZQ*pmp_ZQ + mP*mP)/(pmm_ZQ*(2*mP - pmm_ZQ))) - mP*mP);

    double EP = leadP4.E();
    double EB = omega + nucleusP4.M() - EP;
    double TB = EB - sqrt(EB*EB - miss_neg.Mag2());
    double TP = EP - leadP4.M();

    b_leadP       = pLead3.Mag();
    b_leadTheta   = pLead3.Theta() * rad2deg;
    b_leadPhi     = pLead3.Phi() * rad2deg;
    b_leadVz      = lead[0]->par()->getVz();
    b_leadRegion  = (int)lead[0]->getRegion();
    b_vtzDiffLeadE = b_leadVz - b_vtzE;

    b_pMiss       = miss_neg.Mag();
    b_pMissTheta  = miss_neg.Theta() * rad2deg;
    b_pMissPhi    = miss_neg.Phi() * rad2deg;
    b_mMiss       = missP4.M();
    b_kMiss       = kMiss;
    b_EMiss       = sqrt(pLead3.Mag2() + mP*mP) - omega;
    b_E0miss      = sqrt(b_pMiss*b_pMiss + mN*mN) - mN;
    b_E1miss      = omega - TP - TB;
    b_theta_PmQ   = miss_neg.Angle(qP3) * rad2deg;
    b_theta_PleadQ = pLead3.Angle(qP3) * rad2deg;
    b_goodLead    = true;

    if(xB <= 1.1) continue;
    if(Q2 <= 1.4) continue;
    if(b_mMiss < 0.50 || b_mMiss > 1.30) continue;
    if(b_kMiss <= 0.20) continue;
    if(b_leadP <= 0.85) continue;
    cut_kinematics++;

    bool leadFD = (lead[0]->getRegion() == FD) || inFD(lead[0]->getStatus());
    bool leadCD = (lead[0]->getRegion() == CD) || inCD(lead[0]->getStatus());
    if(!((leadFD && b_leadTheta < 45.) || (leadCD && b_leadTheta > 40.))) continue;
    cut_lead_region++;

    if(!recoil.empty())
    {
      TLorentzVector recoilP4;
      GetLorentzVector_Corrected(recoilP4, recoil[0], eventIsMC);
      TVector3 recoil_p3 = recoilP4.Vect();
      TVector3 pRelV = (miss_neg - recoil_p3) * 0.5;
      TVector3 v_cm  = miss_neg + recoil_p3;

      TVector3 vz = miss_neg.Unit();
      TVector3 vy = miss_neg.Cross(qP3).Unit();
      TVector3 vx = vz.Cross(vy).Unit();

      b_hasRecoil = true;
      b_pRel      = pRelV.Mag();
      b_pRelTheta = pRelV.Theta() * rad2deg;
      b_pRelPhi   = pRelV.Phi() * rad2deg;
      b_pCM       = v_cm.Mag();
      b_pCMx      = v_cm.Dot(vx);
      b_pCMy      = v_cm.Dot(vy);
      b_pCMz      = v_cm.Dot(vz);

      double TP2 = recoilP4.E() - recoilP4.M();
      TLorentzVector miss_Am2 = q + nucleusP4 - leadP4 - recoilP4;
      double TB2 = miss_Am2.E() - miss_Am2.M();
      b_E2miss = q.E() - TP - TP2 - TB2;

      b_recP     = recoil_p3.Mag();
      b_recTheta = recoil_p3.Theta() * rad2deg;
      b_recPhi   = recoil_p3.Phi() * rad2deg;
      b_recVz    = recoil[0]->par()->getVz();
      b_theta_PmPrec = miss_neg.Angle(recoil_p3) * rad2deg;
      b_theta_PrecQ  = recoil_p3.Angle(qP3) * rad2deg;
    }

    if(eventIsMC){
      auto mcInfo = c12->mcparts();
      if(mcInfo && mcInfo->getRows() >= 3){
        TVector3 ve_truth(mcInfo->getPx(0),mcInfo->getPy(0),mcInfo->getPz(0));
        TVector3 vlead_truth(mcInfo->getPx(1),mcInfo->getPy(1),mcInfo->getPz(1));
        TVector3 vrec_truth(mcInfo->getPx(2),mcInfo->getPy(2),mcInfo->getPz(2));
        TVector3 vq_truth = p_beam - ve_truth;
        TVector3 vmiss_truth = vlead_truth - vq_truth;
        TVector3 vcm_truth = vmiss_truth + vrec_truth;

        TVector3 tz = vmiss_truth.Unit();
        TVector3 ty = vmiss_truth.Cross(vq_truth).Unit();
        TVector3 tx = tz.Cross(ty).Unit();
        double omega_truth = Ebeam - sqrt(ve_truth.Mag2() + me*me);
        b_Q2_truth = vq_truth.Mag2() - omega_truth * omega_truth;
        b_pMiss_truth = vmiss_truth.Mag();
        b_pLead_truth = vlead_truth.Mag();
        b_thetaLead_truth = vlead_truth.Theta() * rad2deg;
        b_pRec_truth = vrec_truth.Mag();
        b_pcmX_truth = vcm_truth.Dot(tx);
        b_pcmY_truth = vcm_truth.Dot(ty);
        b_pcmZ_truth = vcm_truth.Dot(tz);
        b_pcmX_truth_lab = vcm_truth.X();
        b_pcmY_truth_lab = vcm_truth.Y();
        b_pcmZ_truth_lab = vcm_truth.Z();
      }
    }

    srcTree->Fill();
    ctr++;
    cut_written++;
  }

  outFile->cd();
  string inputList;
  for(int k = 4; k < argc; k++){
    if(k > 4) inputList += ",";
    inputList += argv[k];
  }
  time_t now = time(nullptr);
  string dateString = ctime(&now);
  while(!dateString.empty() && dateString.back() == '\n') dateString.pop_back();
  ostringstream meta;
  meta << "git_hash=" << getCommandOutput("git rev-parse HEAD") << "\n"
       << "date=" << dateString << "\n"
       << "inputs=" << inputList << "\n"
       << "cuts=xB>1.1,Q2>1.4,mMiss=[0.50,1.30],kMiss>0.20,pLead>0.85,leadFDtheta<45,leadCDtheta>40,recoil_not_required\n"
       << "clasAna_config=getLeadRecoilSRC default config; explicit flat-tree cuts recorded in cuts field\n"
       << "momentum_corrections=GetLorentzVector_Corrected\n"
       << "beam_energy=" << Ebeam << "\n"
       << "target=deuterium\n"
       << "run_group=RGM\n"
       << "QADB=OFF\n"
       << "isMC_arg=" << isMC << "\n"
       << "sigma_gen=0.2\n"
       << "NN_interaction_model=AV18\n"
       << "generator_version=unknown\n"
       << "cutflow_read=" << cut_read << "\n"
       << "cutflow_one_electron=" << cut_one_electron << "\n"
       << "cutflow_lead_found=" << cut_lead_found << "\n"
       << "cutflow_kinematics=" << cut_kinematics << "\n"
       << "cutflow_lead_region=" << cut_lead_region << "\n"
       << "cutflow_written=" << cut_written;
  TNamed metadata("skimmer_metadata", meta.str().c_str());
  metadata.Write();
  srcTree->Write();
  outFile->Close();

  cout << "Done. Processed " << counter << " events.\n";
  return 0;
}
