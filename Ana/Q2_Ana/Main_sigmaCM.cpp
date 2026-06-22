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
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "many_plots.h"
#include "reweighter.h"
#include "TGraphErrors.h"
#include "Corrections.h"


// For Fitting
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/FitResult.h"
#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"

// For defining the functions
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include <vector>

using namespace std;
using namespace clas12;

const int linbin = 100;
const double min_sigma = 0.050;
const double max_sigma = 0.250;


const double c = 29.9792458;

auto db=TDatabasePDG::Instance();
const double beam_E = 5.98636;
const double mass_n = db->GetParticle(2112)->Mass();
const double mass_p = db->GetParticle(2212)->Mass();
const double mass_pi = db->GetParticle(-211)->Mass();
const double mD = 1.8756;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

vector<double> bE_Q2 = {1.5,1.80,2.10,2.40,2.70,3.00,3.50,5.0}; 
const int bQ2 = bE_Q2.size()-1;

int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

void getG(TFile *f, TCanvas * myCanvas, char fileName[100], string objectName, TH1D * h_myhist[7], double min, double max);
void runEvent(const std::unique_ptr<clas12::clas12reader>& c12, clas12ana & clasAna, bool isMC ,double & q, double & x, double & y, double & z, bool & passepp, int & ctr);
void getChi2(TH1D * h_d, TH1D * h_s, double min, double max, double & final_scale, double & min_chi2);
void Usage()
{
  std::cerr << "Usage: ./code A outputfile.root outputfile.pdf inputdatafile.hipo inputsimfiles.hipo \n\n\n";
}

int main(int argc, char ** argv)
{
  
  if(argc < 5)
    {
      Usage();
      return -1;
    }
  
  int nucleus_A = atoi(argv[1]);
  TString outFile = argv[2];
  cout<<"Ouput file "<< outFile <<endl;
    
  ////////////////////////////
  clas12root::HipoChain chain;
  cout<<"Input file "<<argv[3]<<endl;
  chain.Add(argv[3]);
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader(); 
  auto &c12=chain.C12ref();
  ////////////////////////////
  clas12root::HipoChain chain_Sim;
  cout<<"Input file "<<argv[4]<<endl;
  chain_Sim.Add(argv[4]);
  chain_Sim.SetReaderTags({0});
  chain_Sim.db()->turnOffQADB();
  auto config_c12_Sim=chain_Sim.GetC12Reader();
  auto &c12_Sim=chain_Sim.C12ref();

  clas12ana clasAna;
  clasAna.printParams();

  int Z=nucleus_A/2;
  int N=nucleus_A/2;
  if(nucleus_A==48){
    Z=20;
    N=20;
  }
  if(nucleus_A==120){
    Z=50;
    N=70;
  }
  
  char temp_name[100];
  char temp_title[100];

  TH1D * h_Q2[bQ2];
  for(int i=0; i<(bQ2); i++){
    sprintf(temp_title,"%f - %f",bE_Q2[i],bE_Q2[i+1]);

    sprintf(temp_name,"h_Q2_%d",i);
    h_Q2[i] = new TH1D(temp_name,temp_title,200,1.5,5.0);
  }
  
  TH1D * h_pcmx_epp = new TH1D("pcmx_epp","pcmx_epp",25,-0.55,0.55);
  TH1D * h_pcmy_epp = new TH1D("pcmy_epp","pcmy_epp",25,-0.55,0.55);
  TH1D * h_pcmz_epp = new TH1D("pcmz_epp","pcmz_epp",25,-0.5,1.0);
  TH1D * h_pcmT_epp = new TH1D("pcmT_epp","pcmT_epp",25,0.0,0.75);
  TH1D * h_pcmyP_epp = new TH1D("pcmyP_epp","pcmyP_epp",25,-0.55,0.55);

  TH1D * h_pcmx_epp_SRC_Q2[bQ2];
  TH1D * h_pcmy_epp_SRC_Q2[bQ2];
  TH1D * h_pcmz_epp_SRC_Q2[bQ2];
  TH1D * h_pcmT_epp_SRC_Q2[bQ2];
  TH1D * h_pcmyP_epp_SRC_Q2[bQ2];
  for(int i=0; i<(bQ2); i++){
    sprintf(temp_title,"%f - %f",bE_Q2[i],bE_Q2[i+1]);

    sprintf(temp_name,"h_pcmx_epp_SRC_Q2_%d",i);
    h_pcmx_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.55,0.55);

    sprintf(temp_name,"h_pcmy_epp_SRC_Q2_%d",i);
    h_pcmy_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.55,0.55);
    
    sprintf(temp_name,"h_pcmz_epp_SRC_Q2_%d",i);
    h_pcmz_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.5,1.0);
    
    sprintf(temp_name,"h_pcmT_epp_SRC_Q2_%d",i);
    h_pcmT_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,0.0,0.75);

    sprintf(temp_name,"h_pcmyP_epp_SRC_Q2_%d",i);
    h_pcmyP_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.55,0.55);
  }



  std::vector<reweighter> myWeights;
  for(int j=0; j<linbin; j++){
    double sCM= ((double)j/(double)linbin)*(max_sigma-min_sigma) + min_sigma;
    reweighter newWeight(beam_E,Z,N,kelly,"AV18",sCM);
    myWeights.push_back(newWeight);
  }  
  //////////////////////////////////

  TH1D * h_pcmx_epp_simSCM[linbin];
  TH1D * h_pcmy_epp_simSCM[linbin];
  TH1D * h_pcmz_epp_simSCM[linbin];
  TH1D * h_pcmT_epp_simSCM[linbin];
  TH1D * h_pcmyP_epp_simSCM[linbin];
  TH1D * h_pcmx_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmy_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmz_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmT_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmyP_epp_SRC_simSCM_Q2[linbin][bQ2];
  for(int j=0; j<linbin; j++){

    sprintf(temp_name,"h_pcmx_epp_simSCM_%d",j);
    h_pcmx_epp_simSCM[j] = new TH1D(temp_name,"pcmx_epp",25,-0.55,0.55);

    sprintf(temp_name,"h_pcmy_epp_simSCM_%d",j);
    h_pcmy_epp_simSCM[j] = new TH1D(temp_name,"pcmy_epp",25,-0.55,0.55);

    sprintf(temp_name,"h_pcmz_epp_simSCM_%d",j);
    h_pcmz_epp_simSCM[j] = new TH1D(temp_name,"pcmz_epp",25,-0.5,1.0);

    sprintf(temp_name,"h_pcmT_epp_simSCM_%d",j);
    h_pcmT_epp_simSCM[j] = new TH1D(temp_name,"pcmT_epp",25,0.0,0.75);

    sprintf(temp_name,"h_pcmyP_epp_simSCM_%d",j);
    h_pcmyP_epp_simSCM[j] = new TH1D(temp_name,"pcmyP_epp",25,-0.55,0.55);

    for(int i=0; i<(bQ2); i++){
      int sCM= ((double)j/(double)linbin)*(max_sigma-min_sigma) + min_sigma;
      sprintf(temp_title,"sCM=%f %f - %f",sCM,bE_Q2[i],bE_Q2[i+1]);
      
      sprintf(temp_name,"h_pcmx_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmx_epp_SRC_simSCM_Q2[j][i] = new TH1D(temp_name,temp_title,15,-0.55,0.55);
      
      sprintf(temp_name,"h_pcmy_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmy_epp_SRC_simSCM_Q2[j][i] = new TH1D(temp_name,temp_title,15,-0.55,0.55);
      
      sprintf(temp_name,"h_pcmz_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmz_epp_SRC_simSCM_Q2[j][i] = new TH1D(temp_name,temp_title,15,-0.5,1.0);

      sprintf(temp_name,"h_pcmT_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmT_epp_SRC_simSCM_Q2[j][i] = new TH1D(temp_name,temp_title,15,0.0,0.75);

      sprintf(temp_name,"h_pcmyP_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmyP_epp_SRC_simSCM_Q2[j][i] = new TH1D(temp_name,temp_title,15,-0.55,0.55);
    }
  }

  
  ////////////////////////////////////////
  int ctr = 0;
  while(chain.Next() && ctr <1000000)
    {
      double q,x,y,z,yP;
      bool pass = false;
      runEvent(c12,clasAna,false,q,x,y,z,pass,ctr); 
      double T = sqrt(x*x+y*y);
      if(pass){
	h_pcmx_epp->Fill(x,1.0);
	h_pcmy_epp->Fill(y,1.0);
	h_pcmz_epp->Fill(z,1.0);
	h_pcmT_epp->Fill(T,1.0);
	h_pcmyP_epp->Fill(yP,1.0);
      }      
      if((binX(bE_Q2,q)!=-1)&&pass){
	h_Q2[binX(bE_Q2,q)]->Fill(q,1.0);
	h_pcmx_epp_SRC_Q2[binX(bE_Q2,q)]->Fill(x,1.0);
	h_pcmy_epp_SRC_Q2[binX(bE_Q2,q)]->Fill(y,1.0);
	h_pcmz_epp_SRC_Q2[binX(bE_Q2,q)]->Fill(z,1.0);
	h_pcmT_epp_SRC_Q2[binX(bE_Q2,q)]->Fill(T,1.0); 
	h_pcmyP_epp_SRC_Q2[binX(bE_Q2,q)]->Fill(yP,1.0);
     }
    }
  ////////////////////////////////////////
  ctr = 0;
  while(chain_Sim.Next() && ctr <1000000)
    {
      double q,x,y,z,yP;
      bool pass = false;
      double original_weight = c12_Sim->mcevent()->getWeight();
      runEvent(c12_Sim,clasAna,true,q,x,y,z,pass,ctr); 
      double T = sqrt(x*x+y*y);
      for(int i=0; i<linbin; i++){
	double wepp = original_weight * myWeights[i].get_weight_epp(c12_Sim->mcparts());
	if(pass){
	  h_pcmx_epp_simSCM[i]->Fill(x,wepp);
	  h_pcmy_epp_simSCM[i]->Fill(y,wepp);
	  h_pcmz_epp_simSCM[i]->Fill(z,wepp);
	  h_pcmT_epp_simSCM[i]->Fill(T,wepp);
	  h_pcmyP_epp_simSCM[i]->Fill(yP,wepp);
	}      
	if((binX(bE_Q2,q)!=-1)&&pass){
	  h_pcmx_epp_SRC_simSCM_Q2[i][binX(bE_Q2,q)]->Fill(x,wepp);
	  h_pcmy_epp_SRC_simSCM_Q2[i][binX(bE_Q2,q)]->Fill(y,wepp);
	  h_pcmz_epp_SRC_simSCM_Q2[i][binX(bE_Q2,q)]->Fill(z,wepp);
	  h_pcmT_epp_SRC_simSCM_Q2[i][binX(bE_Q2,q)]->Fill(T,wepp);
	  h_pcmyP_epp_SRC_simSCM_Q2[i][binX(bE_Q2,q)]->Fill(yP,wepp);
	}	
      }
    }

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<(bQ2); i++){
    h_Q2[i]->Write();
  }
  h_pcmx_epp->Write();
  h_pcmy_epp->Write();
  h_pcmz_epp->Write();
  h_pcmT_epp->Write();
  h_pcmyP_epp->Write();
  for(int i=0; i<(bQ2); i++){
    h_pcmx_epp_SRC_Q2[i]->Write();
    h_pcmy_epp_SRC_Q2[i]->Write();
    h_pcmz_epp_SRC_Q2[i]->Write();
    h_pcmT_epp_SRC_Q2[i]->Write();
    h_pcmyP_epp_SRC_Q2[i]->Write();
  }
  for(int j=0; j<linbin; j++){
    h_pcmx_epp_simSCM[j]->Write();
    h_pcmy_epp_simSCM[j]->Write();
    h_pcmz_epp_simSCM[j]->Write();
    h_pcmT_epp_simSCM[j]->Write();
    h_pcmyP_epp_simSCM[j]->Write();
    for(int i=0; i<(bQ2); i++){
      h_pcmx_epp_SRC_simSCM_Q2[j][i]->Write();
      h_pcmy_epp_SRC_simSCM_Q2[j][i]->Write();
      h_pcmz_epp_SRC_simSCM_Q2[j][i]->Write();
      h_pcmT_epp_SRC_simSCM_Q2[j][i]->Write();
      h_pcmyP_epp_SRC_simSCM_Q2[j][i]->Write();
    }
  }
  f->Close();

  return 0;
}

void runEvent(const std::unique_ptr<clas12::clas12reader>& c12, clas12ana & clasAna, bool isMC ,double & qSq, double & x, double & y, double & z, bool & passepp, int & ctr){

  
  passepp=false;
  
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector elcpy(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptrcpy(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptrcpy(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  
  //Display completed  
  ctr++;
  if((ctr%100000) == 0){
    cerr << "\n" <<ctr/100000 <<" hundred thousand completed";
  }    
  if((ctr%10000) == 0){
    cerr << ".";
  }    
  
  clasAna.Run(c12);
  auto electrons = clasAna.getByPid(11);
  auto protons = clasAna.getByPid(2212);
  auto pims = clasAna.getByPid(-211);
  auto pips = clasAna.getByPid(211);
  if(electrons.size() != 1){return;}

  GetLorentzVector_Corrected(el,electrons[0],isMC);
  TLorentzVector q = beam - el;
  double Q2        = -q.M2();
  double omega = q.E();
  double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
  double WSq = (mN*mN) - Q2 + (2*omega*mN);
  double W = sqrt(WSq);
  double vtz_e = electrons[0]->par()->getVz();
	  
  int sector_e = electrons[0]->getSector()-1;
  double mom_e = el.P();
  double theta_e = el.Theta() * 180 / M_PI;
  double phi_e = el.Vect().Phi() * 180 / M_PI;
  double shift_e = 7.5;
  clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
  auto lead    = clasAna.getLeadSRC();
  auto recoil  = clasAna.getRecoilSRC();      

  if(lead.size()!=1){return;}
  GetLorentzVector_Corrected(lead_ptr,lead[0],isMC);
  
  TLorentzVector miss = q + deut_ptr - lead_ptr;
  double mmiss2 = miss.M2();
  double mmiss= sqrt(mmiss2);
  double alphamiss = (miss.E() - miss.Vect().Dot(q.Vect().Unit()))/mN;
  TVector3 miss_neg = -miss.Vect();
  double mom_miss = miss.P();
      
  double E0miss = sqrt(mom_miss*mom_miss + mN*mN)-mN;
  double beta_lead = lead[0]->par()->getBeta();
  double mom_lead = lead_ptr.P();
  double momT_lead = lead_ptr.Vect().Perp();
  double theta_lead = lead_ptr.Theta() * 180 / M_PI;
  double phi_lead = lead_ptr.Phi() * 180 / M_PI;
  double vtz_lead = lead[0]->par()->getVz();
  double EP = lead_ptr.E();
  double EB = omega + nucleus_ptr.M() - EP;
  double TB = EB - sqrt(EB*EB - mom_miss*mom_miss);
  double TP = EP - sqrt(EP*EP - mom_lead*mom_lead);
  double E1miss = omega - TP - TB;
  double thetamissq = miss_neg.Angle(q.Vect())*180/M_PI;
  double thetapq = lead_ptr.Vect().Angle(q.Vect())*180/M_PI;
  double poq = mom_lead/q.P();

  TLorentzVector miss_LC = lead_ptr - q;

  TVector3 u = q.Vect().Unit();
  double pmm = miss_LC.E() - miss_LC.Vect().Dot(u);
  double pmp = miss_LC.Vect().Perp(u);
  double pmiss = miss.P();
  double kmiss = sqrt(mN*mN*((pmp*pmp+mN*mN)/(pmm*(2*mN-pmm))) - mN*mN);

  double phidiff = (q.Vect().Phi()*180/M_PI)-phi_lead;
  if(phidiff<-180){phidiff+=360;}
  if(phidiff>180){phidiff-=360;}
            
  bool rec = false;
  if(recoil.size()==1){rec = true;}

  if(xB<1.2){return;}
  if(Q2<1.5){return;}
  if(mmiss<0.65){return;}
  if(mmiss>1.10){return;}
  if(lead[0]->getRegion()!=FD){return;}
  if(kmiss<0.3){return;}            
  if(kmiss>1.0){return;}            
  if(mom_lead<1.0){return;}
  if(theta_lead>37){return;}
  
  if(!rec){return;}
  GetLorentzVector_Corrected(recoil_ptr,recoil[0],isMC);

  double vtz_recoil = recoil[0]->par()->getVz();
  double phi_rec = recoil_ptr.Phi() * 180 / M_PI;
  double TP2 = recoil_ptr.E() - recoil_ptr.M();
  TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
  double TB2 = miss_Am2.E() - miss_Am2.M();
  double E2miss = q.E() - TP - TP2 - TB2;
      
  TVector3 v_rec  = recoil_ptr.Vect();
  TVector3 v_rel  = (miss_neg - v_rec) * 0.5;
  TVector3 v_cm   = miss_neg + v_rec;
      
  TVector3 vz = miss_neg.Unit();
  TVector3 vy = miss_neg.Cross(q.Vect()).Unit();
  TVector3 vx = vz.Cross(vy).Unit();

  double phidiffpp = phi_rec-phi_lead;
  if(phidiffpp<-180){phidiffpp+=360;}
  if(phidiffpp>180){phidiffpp-=360;}

  if(v_rec.Mag()<0.3){return;}
  if(v_rec.Mag()>1.0){return;}

  //std::cout<<mmiss<<endl;
  
  qSq=Q2;
  x=xB;
  x=v_cm.Dot(vx);
  y=v_cm.Dot(vy);
  z=v_cm.Dot(vz);
  passepp=true;
}


void runEvent_oldSave(const std::unique_ptr<clas12::clas12reader>& c12, clas12ana & clasAna, bool isMC ,double & qSq, double & x, double & y, double & z, bool & passepp, int & ctr){

  passepp=false;
  
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector elcpy(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptrcpy(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptrcpy(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  
  //Display completed  
  ctr++;
  if((ctr%100000) == 0){
    cerr << "\n" <<ctr/100000 <<" hundred thousand completed";
  }    
  if((ctr%10000) == 0){
    cerr << ".";
  }    

  clasAna.Run(c12);
  auto electrons = clasAna.getByPid(11);
  auto protons = clasAna.getByPid(2212);
  auto pims = clasAna.getByPid(-211);
  auto pips = clasAna.getByPid(211);
  if(electrons.size() != 1){return;}

  GetLorentzVector_Corrected(el,electrons[0],isMC);
  TLorentzVector q = beam - el;
  double Q2        = -q.M2();
  double omega = q.E();
  double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
  double WSq = (mN*mN) - Q2 + (2*omega*mN);
  double W = sqrt(WSq);
  double vtz_e = electrons[0]->par()->getVz();
	  
  int sector_e = electrons[0]->getSector()-1;
  double mom_e = el.P();
  double theta_e = el.Theta() * 180 / M_PI;
  double phi_e = el.Vect().Phi() * 180 / M_PI;
  double shift_e = 7.5;

  clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
  auto lead    = clasAna.getLeadSRC();
  auto recoil  = clasAna.getRecoilSRC();
      
  if(lead.size()!=1){return;}
  GetLorentzVector_Corrected(lead_ptr,lead[0],isMC);

  TLorentzVector miss = q + deut_ptr - lead_ptr;
  double mmiss2 = miss.M2();
  double mmiss= sqrt(mmiss2);
  double alphamiss = (miss.E() - miss.Vect().Dot(q.Vect().Unit()))/mN;
  TVector3 miss_neg = -miss.Vect();
  double mom_miss = miss.P();
      
  double E0miss = sqrt(mom_miss*mom_miss + mN*mN)-mN;
  double beta_lead = lead[0]->par()->getBeta();
  double mom_lead = lead_ptr.P();
  double momT_lead = lead_ptr.Vect().Perp();
  double theta_lead = lead_ptr.Theta() * 180 / M_PI;
  double phi_lead = lead_ptr.Phi() * 180 / M_PI;
  double vtz_lead = lead[0]->par()->getVz();
  double EP = lead_ptr.E();
  double EB = omega + nucleus_ptr.M() - EP;
  double TB = EB - sqrt(EB*EB - mom_miss*mom_miss);
  double TP = EP - sqrt(EP*EP - mom_lead*mom_lead);
  double E1miss = omega - TP - TB;
  double thetamissq = miss_neg.Angle(q.Vect())*180/M_PI;
  double thetapq = lead_ptr.Vect().Angle(q.Vect())*180/M_PI;
  double poq = mom_lead/q.P();

  TLorentzVector miss_LC = lead_ptr - q;

  TVector3 u = q.Vect().Unit();
  double pmm = miss_LC.E() - miss_LC.Vect().Dot(u);
  double pmp = miss_LC.Vect().Perp(u);
  double pmiss = miss.P();
  double kmiss = sqrt(mN*mN*((pmp*pmp+mN*mN)/(pmm*(2*mN-pmm))) - mN*mN);

  double phidiff = (q.Vect().Phi()*180/M_PI)-phi_lead;
  if(phidiff<-180){phidiff+=360;}
  if(phidiff>180){phidiff-=360;}
            
  bool rec = false;
  if(recoil.size()==1){rec = true;}

  if(xB<1.2){return;}
  if(Q2<1.5){return;}
  if(Q2>5){return;}
  if(mmiss<0.65){return;}
  if(mmiss>1.10){return;}
  if(lead[0]->getRegion()!=FD){return;}

  if(kmiss<0.3){return;}            

  if(pmiss<0.3){return;}
  if(pmiss>0.6){return;}
  //if(kmiss>0.7){return;}            

  if(mom_lead<1.0){return;}
  if(theta_lead>37){return;}
      
  if(!rec){return;}
  GetLorentzVector_Corrected(recoil_ptr,recoil[0],isMC);

  double vtz_recoil = recoil[0]->par()->getVz();
  double phi_rec = recoil_ptr.Phi() * 180 / M_PI;
  double TP2 = recoil_ptr.E() - recoil_ptr.M();
  TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
  double TB2 = miss_Am2.E() - miss_Am2.M();
  double E2miss = q.E() - TP - TP2 - TB2;
      
  TVector3 v_rec  = recoil_ptr.Vect();
  TVector3 v_rel  = (miss_neg - v_rec) * 0.5;
  TVector3 v_cm   = miss_neg + v_rec;
      
  TVector3 vz = miss_neg.Unit();
  TVector3 vy = miss_neg.Cross(q.Vect()).Unit();
  TVector3 vx = vz.Cross(vy).Unit();

  double phidiffpp = phi_rec-phi_lead;
  if(phidiffpp<-180){phidiffpp+=360;}
  if(phidiffpp>180){phidiffpp-=360;}
      
  qSq=Q2;
  x=v_cm.Dot(vx);
  y=v_cm.Dot(vy);
  z=v_cm.Dot(vz);
  passepp=true;

}

void getChi2(TH1D * h_d, TH1D * h_s, double min, double max, double & final_scale, double & min_chi2){
  /////////////////////////////////////  
  //Look first then
  //Calculate the Chi2
  /////////////////////////////////////  
  double start_scale=0.5*h_d->GetMaximum()/h_s->GetMaximum();
  double   end_scale=1.5*h_d->GetMaximum()/h_s->GetMaximum();
  final_scale=1;
  min_chi2 = 100000000000000;
  for(int i = 0; i < 100; i++){
    double scale = ((double)i/100.0)*(end_scale-start_scale) + start_scale;
    double chi2 = 0;
    for(int j = h_d->FindBin(min); j < h_d->FindBin(max)+1; j++){
      TH1D * h_s_clone = (TH1D*) h_s->Clone();
      h_s_clone->Scale(scale);
      double ob = h_d->GetBinContent(j);
      double ob_err = h_d->GetBinError(j);
      double ex = h_s_clone->GetBinContent(j);
      double ex_err = h_s_clone->GetBinError(j);
      double chi2_point = 0.5 * sq(ob - ex) / sq(ob_err + ex_err);   
      chi2+=chi2_point;
    }
    if(chi2<min_chi2){
      final_scale=scale;
      min_chi2=chi2;}
  }
  return;
}


void getG(TFile *f, TCanvas * myCanvas, char fileName[100], string objectName, TH1D * h_myhist[7], double min, double max){

  TGraphErrors * g_mu = new TGraphErrors();
  g_mu->SetName(("g_mu_"+objectName).c_str());
  TGraphErrors * g_sigma = new TGraphErrors();
  g_sigma->SetName(("g_sigma_"+objectName).c_str());

  int ctr = 0;
  //Now project the histogram    
  for(int j = 0; j < bQ2; j++){
    //Define x and y(1D histogram)
    double x = (bE_Q2[j]+bE_Q2[j+1])/2.0;
    ctr++;
    //Now preform a guassian fit
    if(h_myhist[j]->GetEntries()<15){continue;}
    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },min,max,3);
    gFit->SetParameter(0,h_myhist[j]->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,(max+min)/2);
    gFit->SetParLimits(1,min,max);
    gFit->SetParameter(2,(max-min)/4);
    gFit->SetParLimits(2,0.00,max-min);
    
    TFitResultPtr gPoint = h_myhist[j]->Fit(gFit,"SrBeqn","",min,max);
    if(gPoint == 0){
      g_mu->SetPoint(g_mu->GetN(),x,gPoint->Parameter(1));
      g_mu->SetPointError(g_mu->GetN()-1,0,gPoint->ParError(1));

      g_sigma->SetPoint(g_sigma->GetN(),x,gPoint->Parameter(2));
      g_sigma->SetPointError(g_sigma->GetN()-1,0,gPoint->ParError(2));
    }
    
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_myhist[j]->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_mu->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  g_mu->Write();
  g_sigma->Write();
  
}


void getGraph(TFile *f, TCanvas * myCanvas, char fileName[100], string objectName, TH2D * h_myhist, double min, double max){

  TGraphErrors * g_mu = new TGraphErrors();
  g_mu->SetName(("g_mu_"+objectName).c_str());
  TGraphErrors * g_sigma = new TGraphErrors();
  g_sigma->SetName(("g_sigma_"+objectName).c_str());

  int ctr = 0;
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    TH1D * proj = h_myhist->ProjectionY(("h_proj_"+objectName+"_"+to_string(ctr)).c_str(),j+1,j+1);

    //Now preform a guassian fit
    if(proj->GetEntries()<15){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },min,max,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,(max+min)/2);
    gFit->SetParLimits(1,min,max);
    gFit->SetParameter(2,(max-min)/4);
    gFit->SetParLimits(2,0.00,max-min);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",min,max);
    if(gPoint == 0){
      g_mu->SetPoint(g_mu->GetN(),x,gPoint->Parameter(1));
      g_mu->SetPointError(g_mu->GetN()-1,0,gPoint->ParError(1));

      g_sigma->SetPoint(g_sigma->GetN(),x,gPoint->Parameter(2));
      g_sigma->SetPointError(g_sigma->GetN()-1,0,gPoint->ParError(2));
    }
    /*
    proj->Write();

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    proj->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
    */
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_mu->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  g_mu->Write();
  g_sigma->Write();
  
}



  /*
    std::vector<double> x ={1.0, 2.0 , 2.1, 3.3, 3.1, 3.};

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range(0,5); 
  ROOT::Fit::UnBinData data(opt, range, x.size());
  for( auto it:x){data.Add(it);}
  TF1 * myFunction = new TF1("myFunction","ROOT::Math::normal_pdf(x,[1],[0])", 0, 5);
  ROOT::Math::WrappedMultiTF1 fitFunction( *myFunction, myFunction->GetNdim() );
  ROOT::Fit::Fitter fitter;
  fitter.SetFunction( fitFunction, false);
  double initialParams[] = { 2,1};
  fitter.Config().SetParamsSettings(2,initialParams);
  fitter.Config().SetUpdateAfterFit();
  fitter.LikelihoodFit(data);
  ROOT::Fit::FitResult r=fitter.Result();
  r.Print(std::cout);
  return -1;
  */
