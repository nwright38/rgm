#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TTree.h>

#include "Corrections.h"
#include "HipoChain.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

const double mP = 0.9382720813;
const double mD = 1.8756;
const double mPi = 0.13957039;
const double me = 0.000511;
const double rad2deg = 180. / M_PI;

void Usage()
{
  cerr << "Usage: ./skim_deeppim <MC=1,Data=0> <output.root> <input.hipo> [more hipo...] [--beam Ebeam]\n"
       << "Default Ebeam is 5.98636 GeV, matching SRC_General/skim_ep.cpp.\n";
}

void setCorrectedP4(TLorentzVector &p4, clas12::region_part_ptr part, double mass, bool isMC)
{
  p4.SetXYZM(0., 0., 0., mass);
  GetLorentzVector_Corrected(p4, part, isMC);
}

void fillParticleBranches(TLorentzVector const &p4,
                          clas12::region_part_ptr part,
                          Double_t &p, Double_t &theta, Double_t &phi,
                          Double_t &px, Double_t &py, Double_t &pz,
                          Double_t &energy, Double_t &vz,
                          Int_t &status, Int_t &region)
{
  p = p4.P();
  theta = p4.Theta() * rad2deg;
  phi = p4.Phi() * rad2deg;
  px = p4.Px();
  py = p4.Py();
  pz = p4.Pz();
  energy = p4.E();
  vz = part->par()->getVz();
  status = part->getStatus();
  region = static_cast<int>(part->getRegion());
}

int main(int argc, char **argv)
{
  if(argc < 4){
    Usage();
    return -1;
  }

  bool isMC = (atoi(argv[1]) == 1);
  double Ebeam = 5.98636;
  vector<string> inputFiles;

  for(int k = 3; k < argc; k++){
    string arg = argv[k];
    if(arg == "--beam"){
      if(k + 1 >= argc){
        Usage();
        return -1;
      }
      Ebeam = atof(argv[++k]);
    }
    else{
      inputFiles.push_back(arg);
    }
  }

  if(inputFiles.empty()){
    Usage();
    return -1;
  }

  TFile *outFile = new TFile(argv[2], "RECREATE");
  TTree *tree = new TTree("deeppim", "d(e,e'pp pi-) missing-system skim");

  Int_t b_run = -9;
  Long64_t b_event = -9;
  Bool_t b_isMC = isMC;
  Float_t b_weight = 1.f;
  Int_t b_nElectrons = 0;
  Int_t b_nProtons = 0;
  Int_t b_nPiMinus = 0;
  Int_t b_p1Index = -1;
  Int_t b_p2Index = -1;
  Int_t b_pimIndex = -1;

  Double_t b_Q2 = -9.;
  Double_t b_xB = -9.;
  Double_t b_omega = -9.;

  Double_t b_mMiss = -9.;
  Double_t b_mMiss2 = -9.;
  Double_t b_pMiss = -9.;
  Double_t b_thetaMiss = -9.;
  Double_t b_phiMiss = -9.;
  Double_t b_pxMiss = -9.;
  Double_t b_pyMiss = -9.;
  Double_t b_pzMiss = -9.;

  Double_t b_eP = -9., b_eTheta = -9., b_ePhi = -9.;
  Double_t b_ePx = -9., b_ePy = -9., b_ePz = -9., b_eE = -9., b_eVz = -99.;
  Int_t b_eStatus = -9, b_eRegion = -9;

  Double_t b_pimP = -9., b_pimTheta = -9., b_pimPhi = -9.;
  Double_t b_pimPx = -9., b_pimPy = -9., b_pimPz = -9., b_pimE = -9., b_pimVz = -99.;
  Int_t b_pimStatus = -9, b_pimRegion = -9;

  Double_t b_p1P = -9., b_p1Theta = -9., b_p1Phi = -9.;
  Double_t b_p1Px = -9., b_p1Py = -9., b_p1Pz = -9., b_p1E = -9., b_p1Vz = -99.;
  Int_t b_p1Status = -9, b_p1Region = -9;

  Double_t b_p2P = -9., b_p2Theta = -9., b_p2Phi = -9.;
  Double_t b_p2Px = -9., b_p2Py = -9., b_p2Pz = -9., b_p2E = -9., b_p2Vz = -99.;
  Int_t b_p2Status = -9, b_p2Region = -9;

  tree->Branch("run", &b_run, "run/I");
  tree->Branch("event", &b_event, "event/L");
  tree->Branch("isMC", &b_isMC, "isMC/O");
  tree->Branch("weight", &b_weight, "weight/F");
  tree->Branch("nElectrons", &b_nElectrons, "nElectrons/I");
  tree->Branch("nProtons", &b_nProtons, "nProtons/I");
  tree->Branch("nPiMinus", &b_nPiMinus, "nPiMinus/I");
  tree->Branch("p1Index", &b_p1Index, "p1Index/I");
  tree->Branch("p2Index", &b_p2Index, "p2Index/I");
  tree->Branch("pimIndex", &b_pimIndex, "pimIndex/I");

  tree->Branch("Q2", &b_Q2, "Q2/D");
  tree->Branch("xB", &b_xB, "xB/D");
  tree->Branch("omega", &b_omega, "omega/D");

  tree->Branch("mMiss", &b_mMiss, "mMiss/D");
  tree->Branch("mMiss2", &b_mMiss2, "mMiss2/D");
  tree->Branch("pMiss", &b_pMiss, "pMiss/D");
  tree->Branch("thetaMiss", &b_thetaMiss, "thetaMiss/D");
  tree->Branch("phiMiss", &b_phiMiss, "phiMiss/D");
  tree->Branch("pxMiss", &b_pxMiss, "pxMiss/D");
  tree->Branch("pyMiss", &b_pyMiss, "pyMiss/D");
  tree->Branch("pzMiss", &b_pzMiss, "pzMiss/D");

  tree->Branch("eP", &b_eP, "eP/D");
  tree->Branch("eTheta", &b_eTheta, "eTheta/D");
  tree->Branch("ePhi", &b_ePhi, "ePhi/D");
  tree->Branch("ePx", &b_ePx, "ePx/D");
  tree->Branch("ePy", &b_ePy, "ePy/D");
  tree->Branch("ePz", &b_ePz, "ePz/D");
  tree->Branch("eE", &b_eE, "eE/D");
  tree->Branch("eVz", &b_eVz, "eVz/D");
  tree->Branch("eStatus", &b_eStatus, "eStatus/I");
  tree->Branch("eRegion", &b_eRegion, "eRegion/I");

  tree->Branch("pimP", &b_pimP, "pimP/D");
  tree->Branch("pimTheta", &b_pimTheta, "pimTheta/D");
  tree->Branch("pimPhi", &b_pimPhi, "pimPhi/D");
  tree->Branch("pimPx", &b_pimPx, "pimPx/D");
  tree->Branch("pimPy", &b_pimPy, "pimPy/D");
  tree->Branch("pimPz", &b_pimPz, "pimPz/D");
  tree->Branch("pimE", &b_pimE, "pimE/D");
  tree->Branch("pimVz", &b_pimVz, "pimVz/D");
  tree->Branch("pimStatus", &b_pimStatus, "pimStatus/I");
  tree->Branch("pimRegion", &b_pimRegion, "pimRegion/I");

  tree->Branch("p1P", &b_p1P, "p1P/D");
  tree->Branch("p1Theta", &b_p1Theta, "p1Theta/D");
  tree->Branch("p1Phi", &b_p1Phi, "p1Phi/D");
  tree->Branch("p1Px", &b_p1Px, "p1Px/D");
  tree->Branch("p1Py", &b_p1Py, "p1Py/D");
  tree->Branch("p1Pz", &b_p1Pz, "p1Pz/D");
  tree->Branch("p1E", &b_p1E, "p1E/D");
  tree->Branch("p1Vz", &b_p1Vz, "p1Vz/D");
  tree->Branch("p1Status", &b_p1Status, "p1Status/I");
  tree->Branch("p1Region", &b_p1Region, "p1Region/I");

  tree->Branch("p2P", &b_p2P, "p2P/D");
  tree->Branch("p2Theta", &b_p2Theta, "p2Theta/D");
  tree->Branch("p2Phi", &b_p2Phi, "p2Phi/D");
  tree->Branch("p2Px", &b_p2Px, "p2Px/D");
  tree->Branch("p2Py", &b_p2Py, "p2Py/D");
  tree->Branch("p2Pz", &b_p2Pz, "p2Pz/D");
  tree->Branch("p2E", &b_p2E, "p2E/D");
  tree->Branch("p2Vz", &b_p2Vz, "p2Vz/D");
  tree->Branch("p2Status", &b_p2Status, "p2Status/I");
  tree->Branch("p2Region", &b_p2Region, "p2Region/I");

  clas12root::HipoChain chain;
  for(const auto &fname : inputFiles){
    cout << "Input file " << fname << endl;
    chain.Add(fname.c_str());
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  clas12ana clasAna;

  TLorentzVector beam(0., 0., Ebeam, Ebeam);
  TLorentzVector target(0., 0., 0., mP);
  TLorentzVector el;
  TLorentzVector pim;
  TLorentzVector p1;
  TLorentzVector p2;

  long long nRead = 0;
  long long nOneElectron = 0;
  long long nTwoProtons = 0;
  long long nPiMinus = 0;
  long long nTopology = 0;
  long long nGoodElectronKinematics = 0;
  long long nWritten = 0;

  while(chain.Next()){
    if(nRead % 100000 == 0){
      cout << "Processing event " << nRead << "\t" << nWritten << " rows saved." << endl;
    }
    nRead++;

   if(nWritten > 100000) break;

    clasAna.Run(c12);
    auto electrons = clasAna.getByPid(11);
    auto protons = clasAna.getByPid(2212);
    auto piminus = clasAna.getByPid(-211);

    b_nElectrons = static_cast<int>(electrons.size());
    b_nProtons = static_cast<int>(protons.size());
    b_nPiMinus = static_cast<int>(piminus.size());

    if(electrons.size() != 1){ continue; }
    nOneElectron++;
    if(protons.size() < 2){ continue; }
    nTwoProtons++;
    if(piminus.empty()){ continue; }
    nPiMinus++;
    nTopology++;

    setCorrectedP4(el, electrons[0], me, isMC);
    if(el.P() <= 0.){ continue; }

    TLorentzVector q = beam - el;
    b_Q2 = -q.M2();
    b_omega = q.E();
    if(b_omega <= 0.){ continue; }
    b_xB = b_Q2 / (2. * mP * b_omega);
    nGoodElectronKinematics++;

    b_run = c12->runconfig()->getRun();
    b_event = c12->runconfig()->getEvent();
    b_weight = isMC ? c12->mcevent()->getWeight() : 1.f;
    fillParticleBranches(el, electrons[0],
                         b_eP, b_eTheta, b_ePhi, b_ePx, b_ePy, b_ePz,
                         b_eE, b_eVz, b_eStatus, b_eRegion);

    // for(int ipim = 0; ipim < static_cast<int>(piminus.size()); ipim++){
    //   setCorrectedP4(pim, piminus[ipim], mPi, isMC);
    //   fillParticleBranches(pim, piminus[ipim],
    //                        b_pimP, b_pimTheta, b_pimPhi, b_pimPx, b_pimPy, b_pimPz,
    //                        b_pimE, b_pimVz, b_pimStatus, b_pimRegion);
    //   b_pimIndex = ipim;

      for(int ip1 = 0; ip1 < static_cast<int>(protons.size()); ip1++){
        setCorrectedP4(p1, protons[ip1], mP, isMC);
        fillParticleBranches(p1, protons[ip1],
                             b_p1P, b_p1Theta, b_p1Phi, b_p1Px, b_p1Py, b_p1Pz,
                             b_p1E, b_p1Vz, b_p1Status, b_p1Region);
        b_p1Index = ip1;

        // for(int ip2 = ip1 + 1; ip2 < static_cast<int>(protons.size()); ip2++){
        //   setCorrectedP4(p2, protons[ip2], mP, isMC);
        //   fillParticleBranches(p2, protons[ip2],
        //                        b_p2P, b_p2Theta, b_p2Phi, b_p2Px, b_p2Py, b_p2Pz,
        //                        b_p2E, b_p2Vz, b_p2Status, b_p2Region);
        //   b_p2Index = ip2;

          TLorentzVector miss = beam + target - el - p1; // - p2 - pim;
          b_mMiss = miss.M();
          b_mMiss2 = miss.M2();
          b_pMiss = miss.P();
          b_thetaMiss = miss.Theta() * rad2deg;
          b_phiMiss = miss.Phi() * rad2deg;
          b_pxMiss = miss.Px();
          b_pyMiss = miss.Py();
          b_pzMiss = miss.Pz();

          tree->Fill();
          nWritten++;
        //}
      }
    //}
  }

  outFile->cd();
  ostringstream meta;
  time_t now = time(nullptr);
  string dateString = ctime(&now);
  while(!dateString.empty() && dateString.back() == '\n'){ dateString.pop_back(); }
  meta << "date=" << dateString << "\n"
       << "topology=d(e,e'pp pi-); one row per pi-/proton-pair combination\n"
       << "cuts=exactly_one_electron,nProtons>=2,nPiMinus>=1\n"
       << "corrections=GetLorentzVector_Corrected from Corrections.h, matching the SRC_General correction prescription\n"
       << "beam_energy=" << Ebeam << "\n"
       << "target=deuterium\n"
       << "QADB=OFF\n"
       << "isMC_arg=" << isMC << "\n"
       << "events_read=" << nRead << "\n"
       << "cutflow_one_electron=" << nOneElectron << "\n"
       << "cutflow_two_protons=" << nTwoProtons << "\n"
       << "cutflow_pi_minus=" << nPiMinus << "\n"
       << "events_with_topology=" << nTopology << "\n"
       << "cutflow_good_electron_kinematics=" << nGoodElectronKinematics << "\n"
       << "rows_written=" << nWritten;
  TNamed metadata("skimmer_metadata", meta.str().c_str());
  metadata.Write();
  tree->Write();
  outFile->Close();

  cout << "Done. Processed " << nRead << " events. Saved " << nWritten << " rows.\n"
       << "Cutflow: one electron " << nOneElectron
       << ", >=2 protons " << nTwoProtons
       << ", >=1 pi- " << nPiMinus
       << ", electron kinematics " << nGoodElectronKinematics << ".\n";
  return 0;
}
