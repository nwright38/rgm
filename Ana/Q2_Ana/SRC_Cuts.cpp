#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "Functions.h"

using namespace std;
using namespace clas12;

const double mN = 0.938272;

void Usage()
{
  std::cerr << "Usage: ./code isMC A outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";
}

int main(int argc, char ** argv)
{
  if(argc < 5)
    {
      Usage();
      return -1;
    }

  int isMC      = atoi(argv[1]);
  int A         = atoi(argv[2]);
  TString outFile = argv[3];
  char * pdfFile  = argv[4];

  cout << "Output file "     << outFile  << endl;
  cout << "Output PDF file " << pdfFile  << endl;

  clas12ana clasAna;
  clasAna.printParams();

  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout << "Input file " << argv[k] << endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  ////////////////////////////////////////////////////////////////////
  // Particle masses and beam
  ////////////////////////////////////////////////////////////////////
  auto db      = TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD     = 1.8756;
  double mU     = 0.9315;
  double me     = 0.000510998;
  double mA     = 0;
  switch(A){
    case 4:  mA = 4.002603   * mU - 2*me; break;
    case 12: mA = 12.0       * mU - 6*me; break;
    case 40: mA = 39.962591  * mU - 20*me; break;
    case 48: mA = 47.952523  * mU - 20*me; break;
  }

  double beam_E = 5.98636;
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());

  ////////////////////////////////////////////////////////////////////
  // Output ROOT file and TTree
  ////////////////////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  TTree *tree = new TTree("T","myTree");
  double t_e_mom, t_e_theta, t_e_phi, t_mom, t_theta, t_phi, t_xB, t_Q2;
  double t_thetapq, t_pq, t_mmiss, t_pmiss_mom, t_pmiss_theta, t_pmiss_phi, t_theta_pmq;
  tree->Branch("e_mom",      &t_e_mom,      "e_mom/D");
  tree->Branch("e_theta",    &t_e_theta,    "e_theta/D");
  tree->Branch("e_phi",      &t_e_phi,      "e_phi/D");
  tree->Branch("xB",         &t_xB,         "xB/D");
  tree->Branch("Q2",         &t_Q2,         "Q2/D");
  tree->Branch("mom",        &t_mom,        "mom/D");
  tree->Branch("theta",      &t_theta,      "theta/D");
  tree->Branch("phi",        &t_phi,        "phi/D");
  tree->Branch("thetapq",    &t_thetapq,    "thetapq/D");
  tree->Branch("pq",         &t_pq,         "pq/D");
  tree->Branch("mmiss",      &t_mmiss,      "mmiss/D");
  tree->Branch("pmiss_mom",  &t_pmiss_mom,  "pmiss_mom/D");
  tree->Branch("pmiss_theta",&t_pmiss_theta,"pmiss_theta/D");
  tree->Branch("pmiss_phi",  &t_pmiss_phi,  "pmiss_phi/D");
  tree->Branch("theta_pmq",  &t_theta_pmq,  "theta_pmq/D");

  ////////////////////////////////////////////////////////////////////
  // Histograms — identical to the original file
  ////////////////////////////////////////////////////////////////////
  TH1D * h_Q2_bc        = new TH1D("Q2_bc","Q^{2}",1000,0,5);
  TH1D * h_xB_bc        = new TH1D("xB_bc","x_{B}",1000,0,2);
  TH2D * h_phi_theta_bc = new TH2D("phi_theta_bc","#phi_{e} vs. #theta_{e};#phi_{e};#theta_{e}",100,-180,180,100,5,40);
  TH2D * h_qmag_qtheta  = new TH2D("qmag_qtheta","Q Magnitude vs Theta;q Magnitude (GeV);#theta_{q} (deg)",100,0,4,100,0,100);

  // FD
  TH2D * h_thetae_q2_fd    = new TH2D("thetae_q2_fd","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_mmiss_xb_fd     = new TH2D("mmiss_xb_fd","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_fd     = new TH2D("mmiss_q2_fd","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_fd = new TH1D("mmiss_nocuts_fd","Missing Mass (before SRC cuts);M_{miss} (GeV/c^{2})",100,0,2);
  TH1D * h_xb_fd           = new TH1D("xb_fd","x_{B};x_{B};Counts",100,0.5,2.5);
  TH1D * h_kmiss_fd        = new TH1D("kmiss_fd","k_{miss};k_{miss} (GeV/c);Counts",25,0,1.2);
  TH1D * h_q2_fd           = new TH1D("q2_fd","Q^{2};Q^{2} (GeV^{2});Counts",25,0,4);
  TH2D * h_mmiss_thetapq_fd= new TH2D("mmiss_thetapq_fd","M_{miss} vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_thetalead_fd = new TH2D("mmiss_thetalead_fd","M_{miss} vs #theta_{Lead};#theta_{Lead} (deg);M_{miss} (GeV/c^{2})",60,0,50,60,0,2);
  TH2D * h_thetapq_thetalead_fd = new TH2D("thetapq_thetalead_fd","#theta_{pq} vs #theta_{Lead};#theta_{Lead} (deg);#theta_{pq} (deg)",50,0,50,50,0,60);
  TH1D * h_mmiss_fd        = new TH1D("mmiss_fd","M_{miss};M_{miss} (GeV/c^{2});Counts",50,0,2);
  TH2D * h_p_theta_fd      = new TH2D("p_theta_fd","Phase Space (Lead Protons);#theta_{p} (deg);p_{L} (GeV/c)",45,0,180,50,0,2.5);

  // CD
  TH2D * h_thetae_q2_cd    = new TH2D("thetae_q2_cd","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_mmiss_xb_cd     = new TH2D("mmiss_xb_cd","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_cd     = new TH2D("mmiss_q2_cd","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_cd = new TH1D("mmiss_nocuts_cd","Missing Mass (before SRC cuts);M_{miss} (GeV/c^{2})",100,0,2);
  TH1D * h_xb_cd           = new TH1D("xb_cd","x_{B};x_{B};Counts",100,0.5,2.5);
  TH1D * h_kmiss_cd        = new TH1D("kmiss_cd","k_{miss};k_{miss} (GeV/c);Counts",25,0,1.2);
  TH1D * h_q2_cd           = new TH1D("q2_cd","Q^{2};Q^{2} (GeV^{2});Counts",25,0,4);
  TH2D * h_mmiss_thetapq_cd= new TH2D("mmiss_thetapq_cd","M_{miss} vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_thetalead_cd = new TH2D("mmiss_thetalead_cd","M_{miss} vs #theta_{Lead};#theta_{Lead} (deg);M_{miss} (GeV/c^{2})",60,0,50,60,0,2);
  TH2D * h_thetapq_thetalead_cd = new TH2D("thetapq_thetalead_cd","#theta_{pq} vs #theta_{Lead};#theta_{Lead} (deg);#theta_{pq} (deg)",50,0,50,50,0,60);
  TH1D * h_mmiss_cd        = new TH1D("mmiss_cd","M_{miss};M_{miss} (GeV/c^{2});Counts",50,0,2);
  TH2D * h_p_theta_cd      = new TH2D("p_theta_cd","Phase Space (Lead Protons);#theta_{p} (deg);p_{L} (GeV/c)",45,0,180,50,0,2.5);

  // All
  TH2D * h_thetae_q2_all   = new TH2D("thetae_q2_all","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_thetap_q2_all   = new TH2D("thetap_q2_all","Proton Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{p} (deg)",100,0,4,60,0,60);
  TH1D * h_xb_all          = new TH1D("xb_all","x_{B};x_{B};Counts",100,0.5,2.5);
  TH2D * h_mmiss_xb_all    = new TH2D("mmiss_xb_all","M_{miss} vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_all    = new TH2D("mmiss_q2_all","M_{miss} vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_all= new TH1D("mmiss_nocuts_all","Missing Mass (before SRC cuts);M_{miss} (GeV/c^{2})",100,0,2);
  TH1D * h_kmiss_all       = new TH1D("kmiss_all","k_{miss};k_{miss} (GeV/c);Counts",25,0,1.2);
  TH1D * h_q2_all          = new TH1D("q2_all","Q^{2};Q^{2} (GeV^{2});Counts",25,0,4);
  TH2D * h_mmiss_thetapq_all= new TH2D("mmiss_thetapq_all","M_{miss} vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_thetalead_all = new TH2D("mmiss_thetalead_all","M_{miss} vs #theta_{Lead};#theta_{Lead} (deg);M_{miss} (GeV/c^{2})",60,0,50,60,0,2);
  TH2D * h_thetapq_thetalead_all = new TH2D("thetapq_thetalead_all","#theta_{pq} vs #theta_{Lead};#theta_{Lead} (deg);#theta_{pq} (deg)",50,0,50,50,0,60);
  TH1D * h_mmiss_all       = new TH1D("mmiss_all","M_{miss};M_{miss} (GeV/c^{2});Counts",50,0,2);
  TH2D * h_p_theta_all     = new TH2D("p_theta_all","Phase Space (Lead Protons);#theta_{p} (deg);p_{L} (GeV/c)",45,0,180,50,0,2.5);
  TH1D * h_kmiss_src_all   = new TH1D("kmiss_src_all","k_{miss} (SRC lead);k_{miss} (GeV/c);Counts",50,0.3,1);
  TH2D * h_xb_q2_src_all   = new TH2D("xb_q2_src_all","x_{B} vs Q^{2} (SRC lead);Q^{2} (GeV^{2});x_{B}",50,1.5,4,50,1.2,2.5);
  TH2D * h_thetapmq_thetalead_all = new TH2D("thetapmq_thetalead_all","#theta_{pmiss,q} vs #theta_{Lead};#theta_{Lead} (deg);#theta_{pmiss,q}",50,0,50,90,0,180);
  TH2D * h_thetapmq_thetalead_fd  = new TH2D("thetapmq_thetalead_fd","#theta_{pmiss,q} vs #theta_{Lead} (FD);#theta_{Lead} (deg);#theta_{pmiss,q}",50,0,50,90,0,180);
  TH2D * h_thetapmq_thetalead_cd  = new TH2D("thetapmq_thetalead_cd","#theta_{pmiss,q} vs #theta_{Lead} (CD);#theta_{Lead} (deg);#theta_{pmiss,q}",50,0,50,90,0,180);

  // Final e'p
  TH1D * h_xb_final        = new TH1D("xb_final","x_{B};x_{B};Counts",100,1.2,2.5);
  TH1D * h_kmiss_final     = new TH1D("kmiss_final","k_{miss};k_{miss} (GeV/c);Counts",120,0.3,1.2);
  TH1D * h_q2_final        = new TH1D("q2_final","Q^{2};Q^{2} (GeV^{2});Counts",100,1.4,4);
  TH2D * h_thetapq_thetalead_final = new TH2D("thetapq_thetalead_final","#theta_{pq} vs #theta_{Lead};#theta_{Lead} (deg);#theta_{pq} (deg)",100,0,50,100,0,60);
  TH1D * h_mmiss_final     = new TH1D("mmiss_final","M_{miss};M_{miss} (GeV/c^{2});Counts",100,0,2);

  // e'p kinematics
  TH1D * h_emiss           = new TH1D("emiss","E_{miss};E_{miss} (GeV);Counts",50,-0.2,0.5);
  TH2D * h_emiss_omega     = new TH2D("emiss_omega","E_{miss} vs #omega;#omega (GeV);E_{miss} (GeV)",50,0,2.5,50,-0.2,0.6);
  TH2D * h_q2_omega        = new TH2D("q2_omega","Q^{2} vs #omega;#omega (GeV);Q^{2} (GeV/c^{2})",50,0,2.5,50,1,6);
  TH2D * h_mmiss_emiss     = new TH2D("mmiss_emiss","M_{miss} vs E_{miss};E_{miss} (GeV);M_{miss}",50,-0.3,0.6,50,0.4,1.3);

  // Recoil selection
  TH2D * h_lead_rec        = new TH2D("lead_rec","Lead #theta vs Recoil #theta;#theta_{rec} (deg);#theta_{L} (deg)",90,0,180,90,0,180);
  TH1D * h_prec            = new TH1D("prec","Recoil Momentum;p_{rec} (GeV/c);Counts",100,0,1.5);
  TH1D * h_cos0            = new TH1D("cos0","Angle between Lead and Recoil",100,-1,1);

  // e'pp kinematics
  TH1D * h_emiss_pp        = new TH1D("emiss_pp","E_{miss};E_{miss} (GeV)",100,-0.1,0.5);
  TH2D * h_emiss_omega_pp  = new TH2D("emiss_omega_pp","E_{miss} vs #omega;#omega (GeV);E_{miss} (GeV)",50,0,2.5,50,-0.1,0.5);
  TH1D * h_xb_pp           = new TH1D("xb_pp","x_{B};x_{B};Counts",100,1.2,2.5);
  TH1D * h_kmiss_pp        = new TH1D("kmiss_pp","k_{miss};k_{miss} (GeV/c);Counts",120,0.3,1.2);
  TH1D * h_q2_pp           = new TH1D("q2_pp","Q^{2};Q^{2} (GeV^{2});Counts",100,1,4);
  TH2D * h_thetapq_thetalead_pp = new TH2D("thetapq_thetalead_pp","#theta_{pq} vs #theta_{Lead};#theta_{Lead} (deg);#theta_{pq} (deg)",100,0,50,100,0,60);
  TH1D * h_mmiss_pp        = new TH1D("mmiss_pp","M_{miss};M_{miss} (GeV/c^{2});Counts",100,0,1.2);
  TH2D * h_plead_prec      = new TH2D("plead_prec","p_{L} vs p_{rec};p_{rec} (GeV/c);p_{L} (GeV/c)",50,0,1,50,0,3);
  TH2D * h_tlead_trec      = new TH2D("tlead_trec","#theta_{L} vs #theta_{rec};#theta_{rec} (deg);#theta_{L} (deg)",90,0,180,90,0,180);
  TH2D * h_kmiss_prec      = new TH2D("kmiss_prec","k_{miss} vs p_{rec};p_{rec} (GeV/c);k_{miss} (GeV/c)",50,0.3,1,50,0.3,1.2);
  TH1D * h_thetapmq_pp     = new TH1D("thetapmq_pp","#theta_{pmiss,q};#theta_{pmiss,q};Counts",50,90,180);

  // Additional
  TH1D * h_plead_all       = new TH1D("plead_all","p_{L};p_{L} (GeV/c)",100,0,3);
  TH1D * h_plead_fd        = new TH1D("plead_fd","p_{L} (FD);p_{L} (GeV/c)",100,0,3);
  TH1D * h_plead_cd        = new TH1D("plead_cd","p_{L} (CD);p_{L} (GeV/c)",100,0,3);
  TH1D * h_doublelead      = new TH1D("doublelead","N protons passing SRC cuts;N",5,1,6);

  ////////////////////////////////////////////////////////////////////
  // Event loop
  ////////////////////////////////////////////////////////////////////
  int counter = 0;

  while(chain.Next())
    {
      counter++;
      if((counter%1000000) == 0){ cerr << "\n" << counter/1000000 << " million completed"; }
      if((counter%100000)  == 0){ cerr << "."; }

      double weight = 1.0;
      if(isMC){ weight = c12->mcevent()->getWeight(); }

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons   = clasAna.getByPid(2212);

      if(electrons.size() != 1){ continue; }

      // --- Electron with corrections ---
      GetLorentzVector_ReconVector(el,electrons[0]);
      if(!isMC){
        SetLorentzVector_ThetaCorrection(el,electrons[0]);
        SetLorentzVector_MomentumCorrection(el,electrons[0]);
      }
      if(isMC){ SetLorentzVector_MomentumSimulationSmear(el,electrons[0]); }

      TLorentzVector q = beam - el;
      double Q2    = -q.M2();
      double xB    = Q2 / (2 * mass_p * (beam.E() - el.E()));
      double omega = beam.E() - el.E();

      h_Q2_bc->Fill(Q2,weight);
      h_xB_bc->Fill(xB,weight);
      h_phi_theta_bc->Fill(el.Phi()*180/M_PI, el.Theta()*180/M_PI, weight);
      h_qmag_qtheta->Fill(q.Vect().Mag(), q.Vect().Theta()*180/M_PI, weight);

      // --- Find SRC lead and recoil candidates via clasAna ---
      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead   = clasAna.getLeadSRC();
      auto recoil = clasAna.getRecoilSRC();

      // --------------------------------------------------------
      // "Before SRC cuts" histograms: loop over all clasAna
      // protons (PID/fiducial/vertex already applied by clasAna)
      // and fill the staged pre-cut histograms, mirroring what
      // the original code did in its proton loop before applying
      // the SRC cut sequence.
      // --------------------------------------------------------
      int num_lead = 0;
      for(auto p = protons.begin(); p != protons.end(); ++p){
        if((*p)->par()->getCharge() < 1){ continue; }

        GetLorentzVector_ReconVector(lead_ptr,(*p));
        if(!isMC){ SetLorentzVector_ThetaCorrection(lead_ptr,(*p)); }
        SetLorentzVector_EnergyLossCorrection(lead_ptr,(*p));
        if(!isMC){ SetLorentzVector_MomentumCorrection(lead_ptr,(*p)); }
        if(isMC){ SetLorentzVector_MomentumSimulationSmear(lead_ptr,(*p)); }

        double mom_p   = lead_ptr.P();
        double theta_p = lead_ptr.Theta() * 180 / M_PI;
        double phi_p   = lead_ptr.Phi()   * 180 / M_PI;
        double momT_p  = lead_ptr.Vect().Perp();
        double beta_p  = (*p)->par()->getBeta();

        if(beta_p < 0.2){ continue; }

        // raw kinematics for staged histos (using corrected 4-vector)
        TLorentzVector miss_p  = q + deut_ptr - lead_ptr;
        double mmiss_p         = miss_p.M();
        TVector3 pmiss_v       = miss_p.Vect();
        double thetapq_p       = lead_ptr.Vect().Angle(q.Vect()) * 180 / M_PI;
        double theta_pmq       = pmiss_v.Angle(q.Vect()) * 180 / M_PI;

        // kmiss_ZQ (light-cone missing momentum along q)
        TLorentzVector miss_LC_p = lead_ptr - q;
        TVector3 u_ZQ_p          = q.Vect().Unit();
        double pmm_ZQ_p          = miss_LC_p.E() - miss_LC_p.Vect().Dot(u_ZQ_p);
        double pmp_ZQ_p          = miss_LC_p.Vect().Perp(u_ZQ_p);
        double kmiss_ZQ_p        = sqrt(mass_p*mass_p*((pmp_ZQ_p*pmp_ZQ_p + mass_p*mass_p)/(pmm_ZQ_p*(2*mass_p - pmm_ZQ_p))) - mass_p*mass_p);

        // Fill TTree (all protons passing beta cut, same as original)
        t_e_mom       = el.P();
        t_e_theta     = el.Theta()*180/M_PI;
        t_e_phi       = el.Phi()*180/M_PI;
        t_xB          = xB;
        t_Q2          = Q2;
        t_mom         = mom_p;
        t_theta       = theta_p;
        t_phi         = phi_p;
        t_thetapq     = thetapq_p;
        t_pq          = kmiss_ZQ_p;   // store kmiss_ZQ in pq branch for continuity
        t_mmiss       = mmiss_p;
        t_pmiss_mom   = pmiss_v.Mag();
        t_pmiss_theta = pmiss_v.Theta()*180/M_PI;
        t_pmiss_phi   = pmiss_v.Phi()*180/M_PI;
        t_theta_pmq   = theta_pmq;
        tree->Fill();

        bool is_FD = ((*p)->getRegion() == FD);
        bool is_CD = ((*p)->getRegion() == CD);

        if(is_FD){
          h_thetae_q2_fd->Fill(Q2, el.Theta()*180/M_PI, weight);
          h_thetae_q2_all->Fill(Q2, el.Theta()*180/M_PI, weight);
          h_thetap_q2_all->Fill(Q2, theta_p, weight);
          h_plead_all->Fill(mom_p, weight);
          h_plead_fd->Fill(mom_p, weight);
          if(mom_p < 1.0){ continue; }
          h_mmiss_nocuts_fd->Fill(mmiss_p, weight);
          h_mmiss_xb_fd->Fill(xB, mmiss_p, weight);
          h_mmiss_nocuts_all->Fill(mmiss_p, weight);
          h_mmiss_xb_all->Fill(xB, mmiss_p, weight);
          h_xb_fd->Fill(xB, weight);
          h_xb_all->Fill(xB, weight);
          if(xB < 1.2){ continue; }
          h_kmiss_fd->Fill(kmiss_ZQ_p, weight);
          h_kmiss_all->Fill(kmiss_ZQ_p, weight);
          if(kmiss_ZQ_p < 0.3){ continue; }
          h_mmiss_q2_fd->Fill(Q2, mmiss_p, weight);
          h_mmiss_q2_all->Fill(Q2, mmiss_p, weight);
          h_q2_fd->Fill(Q2, weight);
          h_q2_all->Fill(Q2, weight);
          if(Q2 < 1.5){ continue; }
          h_mmiss_thetapq_fd->Fill(thetapq_p, mmiss_p, weight);
          h_mmiss_thetalead_fd->Fill(theta_p, mmiss_p, weight);
          h_thetapq_thetalead_fd->Fill(theta_p, thetapq_p, weight);
          h_thetapq_thetalead_all->Fill(theta_p, thetapq_p, weight);
          h_mmiss_thetapq_all->Fill(thetapq_p, mmiss_p, weight);
          h_mmiss_thetalead_all->Fill(theta_p, mmiss_p, weight);
          h_thetapmq_thetalead_fd->Fill(theta_p, theta_pmq, weight);
          h_thetapmq_thetalead_all->Fill(theta_p, theta_pmq, weight);
          if(theta_p > 37){ continue; }
          h_mmiss_fd->Fill(mmiss_p, weight);
          h_mmiss_all->Fill(mmiss_p, weight);
          if(mmiss_p < 0.65 || mmiss_p > 1.1){ continue; }
          h_p_theta_fd->Fill(theta_p, mom_p, weight);
          h_p_theta_all->Fill(theta_p, mom_p, weight);
          h_kmiss_src_all->Fill(kmiss_ZQ_p, weight);
          h_xb_q2_src_all->Fill(Q2, xB, weight);
          num_lead++;
        }
        else if(is_CD){
          h_thetae_q2_cd->Fill(Q2, el.Theta()*180/M_PI, weight);
          h_thetae_q2_all->Fill(Q2, el.Theta()*180/M_PI, weight);
          h_plead_all->Fill(mom_p, weight);
          h_plead_cd->Fill(mom_p, weight);
          if(mom_p < 1.0){ continue; }
          h_mmiss_nocuts_cd->Fill(mmiss_p, weight);
          h_mmiss_xb_cd->Fill(xB, mmiss_p, weight);
          h_mmiss_nocuts_all->Fill(mmiss_p, weight);
          h_mmiss_xb_all->Fill(xB, mmiss_p, weight);
          h_xb_cd->Fill(xB, weight);
          h_xb_all->Fill(xB, weight);
          if(xB < 1.2){ continue; }
          h_kmiss_cd->Fill(kmiss_ZQ_p, weight);
          h_kmiss_all->Fill(kmiss_ZQ_p, weight);
          if(kmiss_ZQ_p < 0.3){ continue; }
          h_mmiss_q2_cd->Fill(Q2, mmiss_p, weight);
          h_mmiss_q2_all->Fill(Q2, mmiss_p, weight);
          h_q2_cd->Fill(Q2, weight);
          h_q2_all->Fill(Q2, weight);
          if(Q2 < 1.5){ continue; }
          h_mmiss_thetapq_cd->Fill(thetapq_p, mmiss_p, weight);
          h_mmiss_thetalead_cd->Fill(theta_p, mmiss_p, weight);
          h_thetapq_thetalead_cd->Fill(theta_p, thetapq_p, weight);
          h_thetapq_thetalead_all->Fill(theta_p, thetapq_p, weight);
          h_mmiss_thetapq_all->Fill(thetapq_p, mmiss_p, weight);
          h_mmiss_thetalead_all->Fill(theta_p, mmiss_p, weight);
          h_thetapmq_thetalead_cd->Fill(theta_p, theta_pmq, weight);
          h_thetapmq_thetalead_all->Fill(theta_p, theta_pmq, weight);
          if(theta_p > 37){ continue; }
          h_mmiss_cd->Fill(mmiss_p, weight);
          h_mmiss_all->Fill(mmiss_p, weight);
          if(mmiss_p < 0.65 || mmiss_p > 1.1){ continue; }
          h_p_theta_cd->Fill(theta_p, mom_p, weight);
          h_p_theta_all->Fill(theta_p, mom_p, weight);
          h_kmiss_src_all->Fill(kmiss_ZQ_p, weight);
          h_xb_q2_src_all->Fill(Q2, xB, weight);
          num_lead++;
        }
      } // end proton loop

      if(num_lead > 0){ h_doublelead->Fill(num_lead, weight); }

      // --------------------------------------------------------
      // Now use clasAna's SRC-selected lead for the final e'p
      // distributions and the recoil search.
      // --------------------------------------------------------
      if(lead.size() != 1){ continue; }

      GetLorentzVector_ReconVector(lead_ptr,lead[0]);
      if(!isMC){ SetLorentzVector_ThetaCorrection(lead_ptr,lead[0]); }
      SetLorentzVector_EnergyLossCorrection(lead_ptr,lead[0]);
      if(!isMC){ SetLorentzVector_MomentumCorrection(lead_ptr,lead[0]); }
      if(isMC){ SetLorentzVector_MomentumSimulationSmear(lead_ptr,protons[0]); }

      TLorentzVector miss_lead = q + deut_ptr - lead_ptr;
      TVector3 lead_pmiss  = miss_lead.Vect();
      double lead_mom      = lead_ptr.P();
      double lead_theta    = lead_ptr.Theta() * 180 / M_PI;
      double lead_thetapq  = lead_ptr.Vect().Angle(q.Vect()) * 180 / M_PI;
      double lead_mmiss    = miss_lead.M();

      TLorentzVector miss_LC   = lead_ptr - q;
      TVector3 u_ZQ            = q.Vect().Unit();
      double pmm_ZQ            = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
      double pmp_ZQ            = miss_LC.Vect().Perp(u_ZQ);
      double lead_kmiss        = sqrt(mass_p*mass_p*((pmp_ZQ*pmp_ZQ + mass_p*mass_p)/(pmm_ZQ*(2*mass_p - pmm_ZQ))) - mass_p*mass_p);

      double Tp           = lead_ptr.E() - mN;
      double Tb           = omega + mA - lead_ptr.E()
                            - sqrt(pow(omega + mA - lead_ptr.E(), 2) - lead_pmiss.Mag2());
      double lead_emiss   = omega - Tp - Tb;

      h_xb_final->Fill(xB, weight);
      h_kmiss_final->Fill(lead_kmiss, weight);
      h_q2_final->Fill(Q2, weight);
      h_thetapq_thetalead_final->Fill(lead_theta, lead_thetapq, weight);
      h_mmiss_final->Fill(lead_mmiss, weight);
      h_emiss->Fill(lead_emiss, weight);
      h_emiss_omega->Fill(omega, lead_emiss, weight);
      h_q2_omega->Fill(omega, Q2, weight);
      h_mmiss_emiss->Fill(lead_emiss, lead_mmiss, weight);

      // --------------------------------------------------------
      // Recoil (e'pp): look at all clasAna protons in the CD
      // that are not the lead candidate, mirroring the original.
      // --------------------------------------------------------
      for(auto p = protons.begin(); p != protons.end(); ++p){
        if((*p)->par()->getCharge() < 1){ continue; }
        if((*p)->getRegion() != CD){ continue; }

        // Exclude the lead: compare by CVT pindex if lead is in CD,
        // or just confirm it's a different object if lead is in FD.
        if(lead[0]->getRegion() == CD){
          if((*p)->trk(CVT)->getPindex() == lead[0]->trk(CVT)->getPindex()){ continue; }
        }

        GetLorentzVector_ReconVector(recoil_ptr,(*p));
        SetLorentzVector_ThetaCorrection(recoil_ptr,(*p));
        if(!isMC){ SetLorentzVector_EnergyLossCorrection(recoil_ptr,(*p)); }
        SetLorentzVector_MomentumCorrection(recoil_ptr,(*p));

        double mom_rec   = recoil_ptr.P();
        double theta_rec = recoil_ptr.Theta() * 180 / M_PI;

        h_lead_rec->Fill(theta_rec, lead_theta, weight);

        if(mom_rec < 0.35){ continue; }

        h_prec->Fill(mom_rec, weight);
        h_cos0->Fill(cos(lead_pmiss.Angle(recoil_ptr.Vect())), weight);

        h_emiss_pp->Fill(lead_emiss, weight);
        h_emiss_omega_pp->Fill(omega, lead_emiss, weight);
        h_xb_pp->Fill(xB, weight);
        h_kmiss_pp->Fill(lead_kmiss, weight);
        h_q2_pp->Fill(Q2, weight);
        h_thetapq_thetalead_pp->Fill(lead_theta, lead_thetapq, weight);
        h_mmiss_pp->Fill(lead_mmiss, weight);
        h_plead_prec->Fill(mom_rec, lead_mom, weight);
        h_tlead_trec->Fill(theta_rec, lead_theta, weight);
        h_kmiss_prec->Fill(mom_rec, lead_kmiss, weight);
        h_thetapmq_pp->Fill(lead_pmiss.Angle(q.Vect())*180/M_PI, weight);
      } // end recoil loop

    } // end event loop

  ////////////////////////////////////////////////////////////////////
  // Write and plot
  ////////////////////////////////////////////////////////////////////
  f->cd();
  tree->Write();
  h_Q2_bc->Write();
  h_xB_bc->Write();

  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText   = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);

  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);

  auto page1 = [&](TH1* h){ myCanvas->Divide(1,1); myCanvas->cd(1); h->Draw(); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); };
  auto page1logy = [&](TH1* h){ myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogy(); h->Draw(); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); };
  auto page2 = [&](TH2* h){ myCanvas->Divide(1,1); myCanvas->cd(1); h->Draw("colz"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); };
  auto page2logz = [&](TH2* h){ myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogz(); h->Draw("colz"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); };

  auto vline = [](TH1* h, double x) -> TLine* {
    TLine * l = new TLine(x,0,x,h->GetMaximum());
    l->SetLineColor(kRed); l->SetLineWidth(3); return l;
  };
  auto vline2d = [](double x, double ylo, double yhi) -> TLine* {
    TLine * l = new TLine(x,ylo,x,yhi);
    l->SetLineColor(kRed); l->SetLineWidth(3); return l;
  };

  page1(h_Q2_bc);
  page1logy(h_xB_bc);
  page2(h_thetap_q2_all);
  page2(h_qmag_qtheta);

  // All
  page1(h_mmiss_nocuts_all);
  page2(h_thetae_q2_all);
  { myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogz(); h_mmiss_xb_all->Draw("colz"); vline2d(1.2,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogy(); h_xb_all->Draw(); vline(h_xb_all,1.2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_kmiss_all->Draw(); vline(h_kmiss_all,0.3)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page2(h_mmiss_q2_all);
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_q2_all->Draw(); vline(h_q2_all,1.5)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_thetapq_all->Draw("colz"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_thetalead_all->Draw("colz"); vline2d(37,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_thetapq_thetalead_all->Draw("colz"); vline2d(37,0,60)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_all->Draw(); vline(h_mmiss_all,0.65)->Draw("same"); vline(h_mmiss_all,1.1)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page2(h_p_theta_all);
  page1(h_kmiss_src_all);
  page2(h_xb_q2_src_all);

  // FD
  page1(h_mmiss_nocuts_fd);
  page2(h_thetae_q2_fd);
  { myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogz(); h_mmiss_xb_fd->Draw("colz"); vline2d(1.2,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogy(); h_xb_fd->Draw(); vline(h_xb_fd,1.2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_kmiss_fd->Draw(); vline(h_kmiss_fd,0.3)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page2(h_mmiss_q2_fd);
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_q2_fd->Draw(); vline(h_q2_fd,1.5)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_thetapq_fd->Draw("colz"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_thetalead_fd->Draw("colz"); vline2d(37,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_thetapq_thetalead_fd->Draw("colz"); vline2d(37,0,60)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_fd->Draw(); vline(h_mmiss_fd,0.65)->Draw("same"); vline(h_mmiss_fd,1.1)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page2(h_p_theta_fd);

  // CD
  page1(h_mmiss_nocuts_cd);
  page2(h_thetae_q2_cd);
  { myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogz(); h_mmiss_xb_cd->Draw("colz"); vline2d(1.2,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); myCanvas->cd(1)->SetLogy(); h_xb_cd->Draw(); vline(h_xb_cd,1.2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_kmiss_cd->Draw(); vline(h_kmiss_cd,0.3)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_q2_cd->Draw(); vline(h_q2_cd,1.5)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_thetapq_cd->Draw("colz"); vline2d(25,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_thetalead_cd->Draw("colz"); vline2d(37,0,2)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_thetapq_thetalead_cd->Draw("colz"); vline2d(37,0,60)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_mmiss_cd->Draw(); vline(h_mmiss_cd,0.65)->Draw("same"); vline(h_mmiss_cd,1.1)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page2(h_p_theta_cd);

  // Final e'p
  page1(h_xb_final);
  page1(h_kmiss_final);
  page1(h_q2_final);
  page2(h_thetapq_thetalead_final);
  page1(h_mmiss_final);

  // e'p kinematics
  page1(h_emiss);
  page2(h_emiss_omega);
  page2(h_q2_omega);
  page2(h_mmiss_emiss);

  // Lead/recoil
  page2(h_lead_rec);
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_prec->Draw(); vline(h_prec,0.35)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page1(h_cos0);

  // e'pp kinematics
  page1(h_emiss_pp);
  page2(h_emiss_omega_pp);
  page1(h_xb_pp);
  page1(h_kmiss_pp);
  page1(h_q2_pp);
  page2(h_thetapq_thetalead_pp);
  page1(h_mmiss_pp);
  page2(h_plead_prec);
  page2(h_tlead_trec);
  page2(h_kmiss_prec);
  page1(h_thetapmq_pp);

  // Additional
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_plead_all->Draw(); vline(h_plead_all,1.0)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_plead_fd->Draw();  vline(h_plead_fd,1.0)->Draw("same");  myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_plead_cd->Draw();  vline(h_plead_cd,1.0)->Draw("same");  myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_thetapmq_thetalead_all->Draw("colz"); vline2d(37,0,180)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_thetapmq_thetalead_fd->Draw("colz");  vline2d(37,0,180)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  { myCanvas->Divide(1,1); myCanvas->cd(1); h_thetapmq_thetalead_cd->Draw("colz");  vline2d(37,0,180)->Draw("same"); myCanvas->Print(fileName,"pdf"); myCanvas->Clear(); }
  page1(h_doublelead);

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();
  return 0;
}