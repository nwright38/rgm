#include <cmath>
#include <cstdlib>
#include <iostream>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include "Corrections.h"
#include "HipoChain.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

const double mP = 0.9382720813;
const double mN = 0.9395654133;
const double mD = 1.8756;
const double me = 0.000511;

void Usage()
{
  cerr << "Usage: ./bandDISSkim <MC=1,Data=0> <Ebeam(GeV)> <output.root> <input.hipo> [more hipo...]\n";
}

bool inFD(int status){ return (abs(status) >= 2000 && abs(status) < 4000); }

void setParticleP4(TLorentzVector &p4, clas12::region_part_ptr part, int pid, bool isMC)
{
  p4.SetXYZM(0., 0., 0., (pid == 2112) ? mN : mP);
  GetLorentzVector_ReconVector(p4, part);

  if(pid == 2212){
    if(!isMC){ SetLorentzVector_ThetaCorrection(p4, part); }
    SetLorentzVector_EnergyLossCorrection(p4, part);
    if(!isMC){ SetLorentzVector_MomentumCorrection(p4, part); }
    if(isMC){ SetLorentzVector_MomentumSimulationSmear(p4, part); }
  }
}

int main(int argc, char **argv)
{
  if(argc < 5){ Usage(); return -1; }

  bool isMC = (atoi(argv[1]) == 1);
  double Ebeam = atof(argv[2]);

  TFile *outFile = new TFile(argv[3], "RECREATE");
  TTree *bandTree = new TTree("bandTree", "DIS spectator recoil kinematics");

  Int_t b_run, b_event;
  Int_t b_nProtons, b_nNeutrons;
  Float_t b_weight;

  Float_t b_Q2, b_W2, b_W, b_xB, b_omega, b_y;
  Float_t b_eP, b_eE, b_eTheta, b_ePhi, b_eTheta_deg, b_ePhi_deg;
  Float_t b_qP, b_qTheta, b_qPhi, b_qTheta_deg, b_qPhi_deg;

  Int_t b_recoilPid, b_recoilIndex, b_recoilStatus, b_recoilRegion;
  Float_t b_recP, b_recE, b_recMass, b_recTheta, b_recPhi, b_recTheta_deg, b_recPhi_deg;
  Float_t b_recPx, b_recPy, b_recPz, b_recPt;
  Float_t b_recVz, b_recBeta, b_recPath;
  Float_t b_theta_nq, b_theta_nq_deg, b_cos_theta_nq, b_aS;

  Float_t b_Xprime, b_Wprime2, b_Wprime, b_nuPrime;
  Float_t b_pStruckP, b_pStruckE, b_pStruckMass2;
  Float_t b_tSpectator, b_missingMass2;

  bandTree->Branch("run", &b_run, "run/I");
  bandTree->Branch("event", &b_event, "event/I");
  bandTree->Branch("nProtons", &b_nProtons, "nProtons/I");
  bandTree->Branch("nNeutrons", &b_nNeutrons, "nNeutrons/I");
  bandTree->Branch("weight", &b_weight, "weight/F");

  bandTree->Branch("Q2", &b_Q2, "Q2/F");
  bandTree->Branch("W2", &b_W2, "W2/F");
  bandTree->Branch("W", &b_W, "W/F");
  bandTree->Branch("xB", &b_xB, "xB/F");
  bandTree->Branch("omega", &b_omega, "omega/F");
  bandTree->Branch("y", &b_y, "y/F");

  bandTree->Branch("eP", &b_eP, "eP/F");
  bandTree->Branch("eE", &b_eE, "eE/F");
  bandTree->Branch("eTheta", &b_eTheta, "eTheta/F");
  bandTree->Branch("ePhi", &b_ePhi, "ePhi/F");
  bandTree->Branch("eTheta_deg", &b_eTheta_deg, "eTheta_deg/F");
  bandTree->Branch("ePhi_deg", &b_ePhi_deg, "ePhi_deg/F");

  bandTree->Branch("qP", &b_qP, "qP/F");
  bandTree->Branch("qTheta", &b_qTheta, "qTheta/F");
  bandTree->Branch("qPhi", &b_qPhi, "qPhi/F");
  bandTree->Branch("qTheta_deg", &b_qTheta_deg, "qTheta_deg/F");
  bandTree->Branch("qPhi_deg", &b_qPhi_deg, "qPhi_deg/F");

  bandTree->Branch("recoilPid", &b_recoilPid, "recoilPid/I");
  bandTree->Branch("recoilIndex", &b_recoilIndex, "recoilIndex/I");
  bandTree->Branch("recoilStatus", &b_recoilStatus, "recoilStatus/I");
  bandTree->Branch("recoilRegion", &b_recoilRegion, "recoilRegion/I");
  bandTree->Branch("recP", &b_recP, "recP/F");
  bandTree->Branch("recE", &b_recE, "recE/F");
  bandTree->Branch("recMass", &b_recMass, "recMass/F");
  bandTree->Branch("recTheta", &b_recTheta, "recTheta/F");
  bandTree->Branch("recPhi", &b_recPhi, "recPhi/F");
  bandTree->Branch("recTheta_deg", &b_recTheta_deg, "recTheta_deg/F");
  bandTree->Branch("recPhi_deg", &b_recPhi_deg, "recPhi_deg/F");
  bandTree->Branch("recPx", &b_recPx, "recPx/F");
  bandTree->Branch("recPy", &b_recPy, "recPy/F");
  bandTree->Branch("recPz", &b_recPz, "recPz/F");
  bandTree->Branch("recPt", &b_recPt, "recPt/F");
  bandTree->Branch("recVz", &b_recVz, "recVz/F");
  bandTree->Branch("recBeta", &b_recBeta, "recBeta/F");
  bandTree->Branch("recPath", &b_recPath, "recPath/F");

  bandTree->Branch("theta_nq", &b_theta_nq, "theta_nq/F");
  bandTree->Branch("theta_nq_deg", &b_theta_nq_deg, "theta_nq_deg/F");
  bandTree->Branch("cos_theta_nq", &b_cos_theta_nq, "cos_theta_nq/F");
  bandTree->Branch("aS", &b_aS, "aS/F");

  bandTree->Branch("Xprime", &b_Xprime, "Xprime/F");
  bandTree->Branch("Wprime2", &b_Wprime2, "Wprime2/F");
  bandTree->Branch("Wprime", &b_Wprime, "Wprime/F");
  bandTree->Branch("nuPrime", &b_nuPrime, "nuPrime/F");
  bandTree->Branch("pStruckP", &b_pStruckP, "pStruckP/F");
  bandTree->Branch("pStruckE", &b_pStruckE, "pStruckE/F");
  bandTree->Branch("pStruckMass2", &b_pStruckMass2, "pStruckMass2/F");
  bandTree->Branch("tSpectator", &b_tSpectator, "tSpectator/F");
  bandTree->Branch("missingMass2", &b_missingMass2, "missingMass2/F");

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

  TLorentzVector beam(0., 0., Ebeam, Ebeam);
  TLorentzVector target(0., 0., 0., mD);
  TLorentzVector freeNucleon(0., 0., 0., mN);
  TLorentzVector el(0., 0., 0., me);
  TLorentzVector recP4;

  long counter = 0;
  long saved = 0;

  while(chain.Next()){
    if(counter % 100000 == 0){
      cout << "Processing event " << counter << "\t" << saved << " saved." << endl;
    }
    counter++;

    clasAna.Run(c12);
    auto electrons = clasAna.getByPid(11);
    auto protons = clasAna.getByPid(2212);
    auto neutrons = c12->getByID(2112);

    if(electrons.size() != 1){ continue; }
    if(!inFD(electrons[0]->getStatus())){ continue; }

    el.SetXYZM(0., 0., 0., me);
    GetLorentzVector_Corrected(el, electrons[0], isMC);
    if(el.P() <= 0.){ continue; }

    TLorentzVector q = beam - el;
    double Q2 = -q.M2();
    double omega = q.E();
    if(omega <= 0.){ continue; }

    double xB = Q2 / (2. * mN * omega);
    double W2 = (freeNucleon + q).M2();
    if(Q2 <= 2.){ continue; }
    if(W2 <= 4.){ continue; }

    if(isMC){
      b_weight = c12->mcevent()->getWeight();
    }
    else{
      b_weight = 1.f;
    }

    b_run = c12->runconfig()->getRun();
    b_event = c12->runconfig()->getEvent();
    b_nProtons = (int)protons.size();
    b_nNeutrons = (int)neutrons.size();

    b_Q2 = Q2;
    b_W2 = W2;
    b_W = (W2 > 0.) ? sqrt(W2) : -9.f;
    b_xB = xB;
    b_omega = omega;
    b_y = omega / Ebeam;

    b_eP = el.P();
    b_eE = el.E();
    b_eTheta = el.Theta();
    b_ePhi = el.Phi();
    b_eTheta_deg = el.Theta() * 180. / M_PI;
    b_ePhi_deg = el.Phi() * 180. / M_PI;

    b_qP = q.P();
    b_qTheta = q.Theta();
    b_qPhi = q.Phi();
    b_qTheta_deg = q.Theta() * 180. / M_PI;
    b_qPhi_deg = q.Phi() * 180. / M_PI;

    auto fillRecoil = [&](clas12::region_part_ptr recoil, int pid, int idx){
      setParticleP4(recP4, recoil, pid, isMC);
      if(recP4.P() <= 0.2){ return; }

      double mass = (pid == 2112) ? mN : mP;
      TLorentzVector struck = target - recP4;
      TLorentzVector finalHadronic = struck + q;

      double denom = 2. * struck.Dot(q);
      double theta = recP4.Vect().Angle(q.Vect());
      double cosTheta = cos(theta);

      b_recoilPid = pid;
      b_recoilIndex = idx;
      b_recoilStatus = recoil->getStatus();
      b_recoilRegion = recoil->getRegion();
      b_recP = recP4.P();
      b_recE = recP4.E();
      b_recMass = mass;
      b_recTheta = recP4.Theta();
      b_recPhi = recP4.Phi();
      b_recTheta_deg = recP4.Theta() * 180. / M_PI;
      b_recPhi_deg = recP4.Phi() * 180. / M_PI;
      b_recPx = recP4.Px();
      b_recPy = recP4.Py();
      b_recPz = recP4.Pz();
      b_recPt = recP4.Pt();
      b_recVz = recoil->par()->getVz();
      b_recBeta = recoil->par()->getBeta();
      b_recPath = recoil->getPath();

      b_theta_nq = theta;
      b_theta_nq_deg = theta * 180. / M_PI;
      b_cos_theta_nq = cosTheta;
      b_aS = (recP4.E() - recP4.P() * cosTheta) / mass;

      b_Xprime = (denom != 0.) ? Q2 / denom : -9.f;
      b_Wprime2 = finalHadronic.M2();
      b_Wprime = (b_Wprime2 > 0.) ? sqrt(b_Wprime2) : -9.f;
      b_nuPrime = finalHadronic.E() - mass;
      b_pStruckP = struck.P();
      b_pStruckE = struck.E();
      b_pStruckMass2 = struck.M2();
      b_tSpectator = (target - recP4).M2();
      b_missingMass2 = (beam + target - el - recP4).M2();

      bandTree->Fill();
      saved++;
    };

    for(int i = 0; i < (int)protons.size(); i++){
      fillRecoil(protons[i], 2212, i);
    }
    for(int i = 0; i < (int)neutrons.size(); i++){
      fillRecoil(neutrons[i], 2112, i);
    }
  }

  outFile->cd();
  bandTree->Write();
  outFile->Close();

  cout << "Done. Processed " << counter << " events. Saved " << saved << " recoil candidates.\n";
  return 0;
}
