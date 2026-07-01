#include <cstdlib>
#include <iostream>

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

const int MAXP = 4;  // max number of proton candidates stored per event

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

  // per-proton arrays (all detected protons, indexed 0..nProtons-1)
  Int_t   b_nProtons;
  Float_t b_protonP[MAXP],      b_protonTheta[MAXP],     b_protonPhi[MAXP];
  Float_t b_protonVz[MAXP];
  Float_t b_pMiss[MAXP],        b_pMissTheta[MAXP],      b_pMissPhi[MAXP];
  Float_t b_mMiss[MAXP];
  Float_t b_kMiss[MAXP];
  Float_t b_EMiss[MAXP];
  Float_t b_theta_PmQ[MAXP];
  Bool_t  b_goodLead[MAXP];     // per-proton: passes all SRC lead cuts

  // pRel and pCM per proton (filled only when a recoil partner exists)
  Float_t b_pRel[MAXP],         b_pRelTheta[MAXP],       b_pRelPhi[MAXP];
  Float_t b_pCM[MAXP],          b_pCMx[MAXP],            b_pCMy[MAXP],     b_pCMz[MAXP];

  // event-level summary
  Int_t  b_nGoodLeads;
  Bool_t b_singleGoodLead;      // exactly one proton passes all SRC lead cuts

  // lead / recoil indices (-1 if not found)
  // leadIdx: index into proton arrays of the chosen lead (first goodLead)
  // recoilIdx: index of the recoil partner (p > 0.3, passes recoil angle cuts)
  Int_t  b_leadIdx;
  Int_t  b_recoilIdx;

  // recoil kinematics (filled when recoilIdx >= 0)
  Float_t b_recP, b_recTheta, b_recPhi;

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
  srcTree->Branch("protonP",     b_protonP,      Form("protonP[%d]/F", MAXP));
  srcTree->Branch("protonTheta", b_protonTheta,  Form("protonTheta[%d]/F", MAXP));
  srcTree->Branch("protonPhi",   b_protonPhi,    Form("protonPhi[%d]/F", MAXP));
  srcTree->Branch("protonVz",    b_protonVz,     Form("protonVz[%d]/F", MAXP));

  srcTree->Branch("pMiss",       b_pMiss,        Form("pMiss[%d]/F", MAXP));
  srcTree->Branch("pMissTheta",  b_pMissTheta,   Form("pMissTheta[%d]/F", MAXP));
  srcTree->Branch("pMissPhi",    b_pMissPhi,     Form("pMissPhi[%d]/F", MAXP));

  srcTree->Branch("mMiss",       b_mMiss,        Form("mMiss[%d]/F", MAXP));
  srcTree->Branch("kMiss",       b_kMiss,        Form("kMiss[%d]/F", MAXP));
  srcTree->Branch("EMiss",       b_EMiss,        Form("EMiss[%d]/F", MAXP));
  srcTree->Branch("theta_PmQ",   b_theta_PmQ,    Form("theta_PmQ[%d]/F", MAXP));
  srcTree->Branch("goodLead",    b_goodLead,     Form("goodLead[%d]/O", MAXP));

  srcTree->Branch("pRel",        b_pRel,         Form("pRel[%d]/F", MAXP));
  srcTree->Branch("pRelTheta",   b_pRelTheta,    Form("pRelTheta[%d]/F", MAXP));
  srcTree->Branch("pRelPhi",     b_pRelPhi,      Form("pRelPhi[%d]/F", MAXP));
  srcTree->Branch("pCM",         b_pCM,          Form("pCM[%d]/F", MAXP));
  srcTree->Branch("pCMx",        b_pCMx,         Form("pCMx[%d]/F", MAXP));
  srcTree->Branch("pCMy",        b_pCMy,         Form("pCMy[%d]/F", MAXP));
  srcTree->Branch("pCMz",        b_pCMz,         Form("pCMz[%d]/F", MAXP));

  srcTree->Branch("nGoodLeads",     &b_nGoodLeads,     "nGoodLeads/I");
  srcTree->Branch("singleGoodLead", &b_singleGoodLead, "singleGoodLead/O");

  srcTree->Branch("leadIdx",     &b_leadIdx,     "leadIdx/I");
  srcTree->Branch("recoilIdx",   &b_recoilIdx,   "recoilIdx/I");
  srcTree->Branch("recP",        &b_recP,        "recP/F");
  srcTree->Branch("recTheta",    &b_recTheta,    "recTheta/F");
  srcTree->Branch("recPhi",      &b_recPhi,      "recPhi/F");

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
    b_leadIdx    = -1;
    b_recoilIdx  = -1;
    b_recP       = -9.f;  b_recTheta = -9.f;  b_recPhi = -9.f;

    for(int i = 0; i < MAXP; i++){
      b_protonP[i] = -9.f;  b_protonTheta[i] = -9.f;  b_protonPhi[i] = -9.f;  b_protonVz[i] = -99.f;
      b_pMiss[i] = -9.f;    b_pMissTheta[i] = -9.f;   b_pMissPhi[i] = -9.f;
      b_pRel[i] = -9.f;     b_pRelTheta[i] = -9.f;    b_pRelPhi[i] = -9.f;
      b_pCM[i] = -9.f;      b_pCMx[i] = -9.f;         b_pCMy[i] = -9.f;    b_pCMz[i] = -9.f;
      b_mMiss[i] = -9.f;    b_kMiss[i] = -9.f;        b_EMiss[i] = -9.f;
      b_theta_PmQ[i] = -9.f;
      b_goodLead[i] = false;
    }

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

    if(xB < 1.2 || Q2 < 1.5) continue;

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

    // ---- pass 1: fill kinematics for every proton candidate ----
    int nFilled = 0;

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

      b_protonP[nFilled]     = pLead3.Mag();
      b_protonTheta[nFilled] = pLead3.Theta();
      b_protonPhi[nFilled]   = pLead3.Phi();
      b_protonVz[nFilled]    = protons[pr]->par()->getVz();

      b_pMiss[nFilled]       = pMissV.Mag();
      b_pMissTheta[nFilled]  = pMissV.Theta();
      b_pMissPhi[nFilled]    = pMissV.Phi();
      b_mMiss[nFilled]       = missP4.M();
      b_kMiss[nFilled]       = kMiss;
      b_EMiss[nFilled]       = EMiss;
      b_theta_PmQ[nFilled]   = pMissV.Angle(qP3);

      // SRC lead cuts
      bool passCuts = true;
      if(pLead3.Mag() < 1.)                        passCuts = false;
      if(missP4.M() < 0.65 || missP4.M() > 1.1)   passCuts = false;
      if(kMiss < 0.3 || kMiss > 1.)                passCuts = false;
 //     if(pLead3.Angle(qP3) < 37.*M_PI/180.)             passCuts = false;

      b_goodLead[nFilled] = passCuts;
      if(passCuts) b_nGoodLeads++;

      nFilled++;
    }

    b_nProtons      = nFilled;
    b_singleGoodLead = (b_nGoodLeads == 1);

    // ---- pass 2: identify lead + recoil and fill pRel / pCM ----
    // Choose the first (lowest-index) proton that passes the lead cuts.
    for(int i = 0; i < nFilled; i++){
      if(b_goodLead[i]){ b_leadIdx = i; break; }
    }

    if(b_leadIdx >= 0 && nFilled >= 2)
    {
      // The lead proton's pMiss is stored in b_pMiss[leadIdx]; use it as pn for rotation.
      TVector3 pn(b_pMiss[b_leadIdx] * sin(b_pMissTheta[b_leadIdx]) * cos(b_pMissPhi[b_leadIdx]),
                  b_pMiss[b_leadIdx] * sin(b_pMissTheta[b_leadIdx]) * sin(b_pMissPhi[b_leadIdx]),
                  b_pMiss[b_leadIdx] * cos(b_pMissTheta[b_leadIdx]));

      // Search for a recoil: any other proton with p > 0.3
      // Recoil angular cuts: open (no FD theta < 37 restriction); require CD (or simply p > 0.3).
      for(int j = 0; j < nFilled; j++)
      {
        if(j == b_leadIdx) continue;
        if(b_protonP[j] < 0.3) continue;   // recoil momentum cut

        // Optional recoil angular cut — looser than lead, accept both FD and CD.
        // Add detector-based cuts here if needed via protons[j]->getStatus().

        b_recoilIdx = j;  // take the first qualifying recoil

        TVector3 recoil_p3(b_protonP[j] * sin(b_protonTheta[j]) * cos(b_protonPhi[j]),
                           b_protonP[j] * sin(b_protonTheta[j]) * sin(b_protonPhi[j]),
                           b_protonP[j] * cos(b_protonTheta[j]));

        // ---- pRel: relative momentum in the pair CM frame ----
        // pRel = (p_lead - p_recoil) / 2  in the pair rest frame; approximate as half the
        // difference of lab momenta (standard SRC convention).
        TVector3 lead_p3(b_protonP[b_leadIdx] * sin(b_protonTheta[b_leadIdx]) * cos(b_protonPhi[b_leadIdx]),
                         b_protonP[b_leadIdx] * sin(b_protonTheta[b_leadIdx]) * sin(b_protonPhi[b_leadIdx]),
                         b_protonP[b_leadIdx] * cos(b_protonTheta[b_leadIdx]));

        TVector3 pRelV = (lead_p3 - recoil_p3) * 0.5;
        b_pRel[b_leadIdx]      = pRelV.Mag();
        b_pRelTheta[b_leadIdx] = pRelV.Theta();
        b_pRelPhi[b_leadIdx]   = pRelV.Phi();

        // ---- pCM: pair CM momentum via coordinate rotation ----
        // Rotate into the frame defined by pn (= pMiss of lead) and q.
        TVector3 vq      = qP3;               // q 3-vector
        TVector3 v1      = pn;                // miss momentum of lead (pn)
        TVector3 vq_copy = vq;

        double Pmiss_phi   = v1.Phi();
        double Pmiss_theta = v1.Theta();

        vq_copy.RotateZ(-Pmiss_phi);
        vq_copy.RotateY(-Pmiss_theta);
        double q_phi = vq_copy.Phi();
        vq_copy.RotateZ(-q_phi);

        TVector3 v2 = recoil_p3;
        v2.RotateZ(-Pmiss_phi);
        v2.RotateY(-Pmiss_theta);
        v2.RotateZ(-q_phi);

        v1.RotateZ(-Pmiss_phi);
        v1.RotateY(-Pmiss_theta);
        v1.RotateZ(-q_phi);

        TVector3 pcom = v1 + v2;

        b_pCM[b_leadIdx]  = pcom.Mag();
        b_pCMx[b_leadIdx] = pcom.X();
        b_pCMy[b_leadIdx] = pcom.Y();
        b_pCMz[b_leadIdx] = pcom.Z();

        // recoil summary kinematics
        b_recP     = recoil_p3.Mag();
        b_recTheta = recoil_p3.Theta();
        b_recPhi   = recoil_p3.Phi();

        break;  // one recoil per event
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