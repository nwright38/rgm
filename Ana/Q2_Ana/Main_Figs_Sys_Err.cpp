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

void getMeanStddev(vector<double> list, double & mean, double & stddev){
  mean = 0;
  stddev = 0;

  double sum_mean = 0;
  double sum_stddev = 0;

  for(int i = 0; i < list.size(); i++){
    sum_mean+=list[i];
  }
  mean = sum_mean/list.size();
  
  for(int i = 0; i < list.size(); i++){
    sum_stddev+=((list[i]-mean)*(list[i]-mean));
  }
  stddev = sqrt(sum_stddev/list.size());
}

vector<double> bE_Q2 = {1.5,1.80,2.10,2.40,2.70,3.00,3.50,5.0}; 
static const int bQ2 = bE_Q2.size()-1;
vector<double> bE_pmiss = {0.4,0.55,0.7,0.85,1.0};
vector<double> bE_kmiss = {0.3,0.45,0.6,0.75,0.9};
vector<double> bE_pmiss_long = {0.4,0.45,0.5,0.55,0.6,0.65,0.75,0.85,0.95,1.05,1.2};
vector<double> bE_kmiss_long = {0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1};

struct HistGroup{

  vector<TH1*> hist_list;

  TH1D * h_Q2_ep;
  TH1D * h_Q2_epp;

  TH1D * h_Q2_ep_SRC_Q2[7];
  TH1D * h_Q2_epp_SRC_Q2[7];
  
  TH1D * h_pcmx_epp;
  TH1D * h_pcmy_epp;
  TH1D * h_pcmz_epp;
  
  TH1D * h_pcmx_epp_SRC_Q2[7];
  TH1D * h_pcmy_epp_SRC_Q2[7];
  TH1D * h_pcmz_epp_SRC_Q2[7];

  TH1D * h_pMiss_ep;
  TH1D * h_pMiss_epp;

  TH1D * h_kMiss_ep;
  TH1D * h_kMiss_epp;
 
  TH1D * h_Q2_ep_SRC_pmiss[4];
  TH1D * h_Q2_epp_SRC_pmiss[4];
  TH1D * h_Q2_ep_SRC_kmiss[4];
  TH1D * h_Q2_epp_SRC_kmiss[4];

  TH1D * h_E0miss_ep_SRC_pmiss[4];
  TH1D * h_E1miss_ep_SRC_pmiss[4];
  TH1D * h_E1miss_epp_SRC_pmiss[4];
  TH1D * h_E2miss_epp_SRC_pmiss[4];
    
  TH1D * h_E0miss_ep_SRC_kmiss[4];
  TH1D * h_E1miss_ep_SRC_kmiss[4];
  TH1D * h_E1miss_epp_SRC_kmiss[4];
  TH1D * h_E2miss_epp_SRC_kmiss[4];

  TH1D * h_E1miss_ep_SRC_pmiss_Q2[4][7];
  TH1D * h_E1miss_epp_SRC_pmiss_Q2[4][7];
  TH1D * h_E2miss_epp_SRC_pmiss_Q2[4][7];
    
  TH1D * h_E1miss_ep_SRC_kmiss_Q2[4][7];
  TH1D * h_E1miss_epp_SRC_kmiss_Q2[4][7];
  TH1D * h_E2miss_epp_SRC_kmiss_Q2[4][7];

};

int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

void getG(TFile *f, TCanvas * myCanvas, char fileName[100], string objectName, TH1D * h_myhist[7], double min, double max);
double getSigma(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * h_myhist, double min, double max);
double getStdDev(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * h_myhist, double min, double max);
double getMean(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * h_myhist, double min, double max);
void runEvent(const std::unique_ptr<clas12::clas12reader>& c12, clas12ana & clasAna, bool isMC, double & qSq, double & x, double & pL, double & tL, double & pR, double & px, double & py, double & pz, double & mM, double & pM, double & kM, double & E0, double & E1, double & E2, bool & passep, bool & passepp, int & ctr);
void CutRandom(double x, double qSq, double mM, double kM, double pL, double tL, double pR, bool & passep, bool & passepp, bool randomize);
void fillUpHistGroup(HistGroup & myGroup,double qSq,double px,double py,double pz,double pM,double kM,double E0,double E1,double E2,bool passep,bool passepp,double wep,double wepp);
TGraphErrors * getGraphWithError(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * nominal, vector<TH1*> hist_Errors);

void setUpHistGroup(HistGroup & myGroup);

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
  
  int isMC = atoi(argv[1]);
  char*uType="AV18";
  if(isMC>0){
    if(isMC==1){
      uType="AV18";
    }
    if(isMC==2){
      uType="AV4";
    }
    if(isMC==3){
      uType="N2LO10";
    }
    if(isMC==4){
      uType="N2LO12";
    }
    if(isMC==5){
      uType="NV";
    }
    isMC=1;
  }
  
  int nucleus_A = atoi(argv[2]);
  TString outFile = argv[3];
  char * pdfFile = argv[4];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  //clas12root::HipoChain chain;
  /*
  */
  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  
  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();
  
  int Z=2;
  int N=2;
  if(isMC){
    Z=nucleus_A/2;
    N=nucleus_A/2;
  }
  reweighter newWeight(beam_E,Z,N,kelly,uType,.15);
  vector<reweighter> Weight_List;
  for(int i = 0; i < 100; i++){
    reweighter Random_Weight(beam_E,Z,N,kelly,uType,.15);
    Random_Weight.randomize_Config();
    Weight_List.push_back(Random_Weight);
  }
  //reweighter newWeight(beam_E,Z,N,kelly,"AV4");
  //reweighter newWeight(beam_E,Z,N,kelly,"N2LO10");
  //reweighter newWeight(beam_E,Z,N,kelly,"N2LO12");
  //reweighter newWeight(beam_E,Z,N,kelly,"NV");
  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];


  HistGroup dataGroup;
  setUpHistGroup(dataGroup);

  HistGroup dataGroup_SysErr[100];
  for(int i = 0; i <100; i++){
    setUpHistGroup(dataGroup_SysErr[i]);
  }
  
  int ctr = 0;
  while(chain.Next() && ctr < 10000000)
    {

      if(ctr%1000 == 0){
        cout << "Event " << ctr << endl;
      }
      double wep = 1;
      double wepp = 1;
      double original_weight; 
      if(isMC){
        original_weight = c12->mcevent()->getWeight(); 
        wep = original_weight * newWeight.get_weight_ep(c12->mcparts());
        wepp = original_weight * newWeight.get_weight_epp(c12->mcparts());
      }

      double qSq,x,pL,tL,pR,px,py,pz,mM,pM,kM,E0,E1,E2;      
      bool passep = false;
      bool passepp = false;
      runEvent(c12,clasAna,false,qSq,x,pL,tL,pR,px,py,pz,mM,pM,kM,E0,E1,E2,passep,passepp,ctr);

      bool passep_cut = passep;
      bool passepp_cut = passepp;
      CutRandom(x,qSq,mM,kM,pL,tL,pR,passep_cut,passepp_cut,false);
      fillUpHistGroup(dataGroup,qSq,px,py,pz,pM,kM,E0,E1,E2,passep_cut,passepp_cut,wep,wepp);      

      for(int i = 0; i < 100; i++){
        double wep_Sys,wepp_Sys;
            if(isMC){
              wep_Sys = original_weight * Weight_List[i].get_weight_ep(c12->mcparts());
              wepp_Sys = original_weight * Weight_List[i].get_weight_epp(c12->mcparts());	  
        }
        bool passep_cut_SysErr = passep;
        bool passepp_cut_SysErr = passepp;
        CutRandom(x,qSq,mM,kM,pL,tL,pR,passep_cut_SysErr,passepp_cut_SysErr,true);
        fillUpHistGroup(dataGroup_SysErr[i],qSq,px,py,pz,pM,kM,E0,E1,E2,passep_cut_SysErr,passepp_cut_SysErr,wep_Sys,wepp_Sys);	
      }
    }

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  double Q2_mean[bQ2];
  for(int j=0; j<(bQ2); j++){
    Q2_mean[j]=dataGroup.h_Q2_epp_SRC_Q2[j]->GetMean();
  }

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  dataGroup.h_Q2_ep->Write();
  dataGroup.h_Q2_epp->Write();
  for(int j = 0; j < 7; j++){
    dataGroup.h_Q2_ep_SRC_Q2[j]->Write();
    dataGroup.h_Q2_epp_SRC_Q2[j]->Write();
  }
  
  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);

  vector<TH1*> hist_pcmx_epp;
  for(int i = 0; i < 100; i++){hist_pcmx_epp.push_back(dataGroup_SysErr[i].h_pcmx_epp);}
  TGraphErrors * g_pcmx_epp = getGraphWithError(f,myCanvas,fileName,dataGroup.h_pcmx_epp,hist_pcmx_epp);

  vector<TH1*> hist_pcmy_epp;
  for(int i = 0; i < 100; i++){hist_pcmy_epp.push_back(dataGroup_SysErr[i].h_pcmy_epp);}
  TGraphErrors * g_pcmy_epp = getGraphWithError(f,myCanvas,fileName,dataGroup.h_pcmy_epp,hist_pcmy_epp);

  vector<TH1*> hist_pcmz_epp;
  for(int i = 0; i < 100; i++){hist_pcmz_epp.push_back(dataGroup_SysErr[i].h_pcmz_epp);}
  TGraphErrors * g_pcmz_epp = getGraphWithError(f,myCanvas,fileName,dataGroup.h_pcmz_epp,hist_pcmz_epp);

  for(int j=0; j<(bQ2); j++){
    vector<TH1*> hist_pcmx_epp_SRC_Q2;
    for(int i = 0; i < 100; i++){hist_pcmx_epp_SRC_Q2.push_back(dataGroup_SysErr[i].h_pcmx_epp_SRC_Q2[j]);}
    TGraphErrors * g_pcmx_epp_SRC_Q2 = getGraphWithError(f,myCanvas,fileName,dataGroup.h_pcmx_epp_SRC_Q2[j],hist_pcmx_epp_SRC_Q2);
  }
  
  for(int j=0; j<(bQ2); j++){
    vector<TH1*> hist_pcmy_epp_SRC_Q2;
    for(int i = 0; i < 100; i++){hist_pcmy_epp_SRC_Q2.push_back(dataGroup_SysErr[i].h_pcmy_epp_SRC_Q2[j]);}
    TGraphErrors * g_pcmy_epp_SRC_Q2 = getGraphWithError(f,myCanvas,fileName,dataGroup.h_pcmy_epp_SRC_Q2[j],hist_pcmy_epp_SRC_Q2);
  }

  ////////////
  
  TGraphErrors * g_sigma_pcmx = new TGraphErrors();
  g_sigma_pcmx->SetName("g_sigma_pcmx");
  for(int j=0; j<(bQ2); j++){
    vector<double> sigma_pcmx_epp_SRC_Q2;
    for(int i = 0; i < 100; i++){sigma_pcmx_epp_SRC_Q2.push_back(getSigma(f,myCanvas,fileName,dataGroup_SysErr[i].h_pcmx_epp_SRC_Q2[j],-0.2,0.2));}
    double mean, sigma;
    getMeanStddev(sigma_pcmx_epp_SRC_Q2,mean,sigma);
    g_sigma_pcmx->SetPoint(g_sigma_pcmx->GetN(),Q2_mean[j],mean);
    g_sigma_pcmx->SetPointError(g_sigma_pcmx->GetN()-1,0,sigma);
  }
  g_sigma_pcmx->Write();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma_pcmx->Draw("");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  TGraphErrors * g_sigma_pcmy = new TGraphErrors();
  g_sigma_pcmy->SetName("g_sigma_pcmy");
  for(int j=0; j<(bQ2); j++){
    vector<double> sigma_pcmy_epp_SRC_Q2;
    for(int i = 0; i < 100; i++){sigma_pcmy_epp_SRC_Q2.push_back(getSigma(f,myCanvas,fileName,dataGroup_SysErr[i].h_pcmy_epp_SRC_Q2[j],-0.2,0.2));}
    double mean, sigma;
    getMeanStddev(sigma_pcmy_epp_SRC_Q2,mean,sigma);
    g_sigma_pcmy->SetPoint(g_sigma_pcmy->GetN(),Q2_mean[j],mean);
    g_sigma_pcmy->SetPointError(g_sigma_pcmy->GetN()-1,0,sigma);
  }
  g_sigma_pcmy->Write();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma_pcmy->Draw("");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  ////////////

  vector<TH1*> hist_rat_pMiss;
  dataGroup.h_pMiss_epp->Divide(dataGroup.h_pMiss_ep);
  for(int i = 0; i < 100; i++){
    dataGroup_SysErr[i].h_pMiss_epp->Divide(dataGroup_SysErr[i].h_pMiss_ep);
    hist_rat_pMiss.push_back(dataGroup_SysErr[i].h_pMiss_epp);
  }
  TGraphErrors * g_rat_pMiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_pMiss_epp,hist_rat_pMiss);

  vector<TH1*> hist_rat_kMiss;
  dataGroup.h_kMiss_epp->Divide(dataGroup.h_kMiss_ep);
  for(int i = 0; i < 100; i++){
    dataGroup_SysErr[i].h_kMiss_epp->Divide(dataGroup_SysErr[i].h_kMiss_ep);
    hist_rat_kMiss.push_back(dataGroup_SysErr[i].h_kMiss_epp);
  }
  TGraphErrors * g_rat_kMiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_kMiss_epp,hist_rat_kMiss);

  for(int j=0; j<4; j++){
    vector<TH1*> hist_rat_Q2_pMiss;
    dataGroup.h_Q2_epp_SRC_pmiss[j]->Divide(dataGroup.h_Q2_ep_SRC_pmiss[j]);
    for(int i = 0; i < 100; i++){
      dataGroup_SysErr[i].h_Q2_epp_SRC_pmiss[j]->Divide(dataGroup_SysErr[i].h_Q2_ep_SRC_pmiss[j]);
      hist_rat_Q2_pMiss.push_back(dataGroup_SysErr[i].h_Q2_epp_SRC_pmiss[j]);
    }
    TGraphErrors * g_rat_Q2_pMiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_Q2_epp_SRC_pmiss[j],hist_rat_Q2_pMiss);    
  }

  for(int j=0; j<4; j++){
    vector<TH1*> hist_rat_Q2_kMiss;
    dataGroup.h_Q2_epp_SRC_kmiss[j]->Divide(dataGroup.h_Q2_ep_SRC_kmiss[j]);
    for(int i = 0; i < 100; i++){
      dataGroup_SysErr[i].h_Q2_epp_SRC_kmiss[j]->Divide(dataGroup_SysErr[i].h_Q2_ep_SRC_kmiss[j]);
      hist_rat_Q2_kMiss.push_back(dataGroup_SysErr[i].h_Q2_epp_SRC_kmiss[j]);
    }
    TGraphErrors * g_rat_Q2_kMiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_Q2_epp_SRC_kmiss[j],hist_rat_Q2_kMiss);    
  }

  ////////////

  for(int j=0; j<4; j++){
    vector<TH1*> hist_E1miss_ep_pmiss;
    for(int i = 0; i < 100; i++){hist_E1miss_ep_pmiss.push_back(dataGroup_SysErr[i].h_E1miss_ep_SRC_pmiss[j]);}
    TGraphErrors * g_E1miss_ep_pmiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_E1miss_ep_SRC_pmiss[j],hist_E1miss_ep_pmiss);
  }

  for(int j=0; j<4; j++){
    vector<TH1*> hist_E1miss_epp_pmiss;
    for(int i = 0; i < 100; i++){hist_E1miss_epp_pmiss.push_back(dataGroup_SysErr[i].h_E1miss_epp_SRC_pmiss[j]);}
    TGraphErrors * g_E1miss_epp_pmiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_E1miss_epp_SRC_pmiss[j],hist_E1miss_epp_pmiss);
  }

  for(int j=0; j<4; j++){
    vector<TH1*> hist_E2miss_epp_pmiss;
    for(int i = 0; i < 100; i++){hist_E2miss_epp_pmiss.push_back(dataGroup_SysErr[i].h_E2miss_epp_SRC_pmiss[j]);}
    TGraphErrors * g_E2miss_epp_pmiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_E2miss_epp_SRC_pmiss[j],hist_E2miss_epp_pmiss);
  }

  for(int j=0; j<4; j++){
    vector<TH1*> hist_E1miss_ep_kmiss;
    for(int i = 0; i < 100; i++){hist_E1miss_ep_kmiss.push_back(dataGroup_SysErr[i].h_E1miss_ep_SRC_kmiss[j]);}
    TGraphErrors * g_E1miss_ep_kmiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_E1miss_ep_SRC_kmiss[j],hist_E1miss_ep_kmiss);
  }

  for(int j=0; j<4; j++){
    vector<TH1*> hist_E1miss_epp_kmiss;
    for(int i = 0; i < 100; i++){hist_E1miss_epp_kmiss.push_back(dataGroup_SysErr[i].h_E1miss_epp_SRC_kmiss[j]);}
    TGraphErrors * g_E1miss_epp_kmiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_E1miss_epp_SRC_kmiss[j],hist_E1miss_epp_kmiss);
  }

  for(int j=0; j<4; j++){
    vector<TH1*> hist_E2miss_epp_kmiss;
    for(int i = 0; i < 100; i++){hist_E2miss_epp_kmiss.push_back(dataGroup_SysErr[i].h_E2miss_epp_SRC_kmiss[j]);}
    TGraphErrors * g_E2miss_epp_kmiss = getGraphWithError(f,myCanvas,fileName,dataGroup.h_E2miss_epp_SRC_kmiss[j],hist_E2miss_epp_kmiss);
  }

  ////////////  
  
  TGraphErrors * g_sigma_E1miss_ep_pmiss[4];
  TGraphErrors * g_sigma_E1miss_epp_pmiss[4];
  TGraphErrors * g_sigma_E2miss_epp_pmiss[4];
  TGraphErrors * g_sigma_E1miss_ep_kmiss[4];
  TGraphErrors * g_sigma_E1miss_epp_kmiss[4];
  TGraphErrors * g_sigma_E2miss_epp_kmiss[4];
  
  TGraphErrors * g_mean_E1miss_ep_pmiss[4];
  TGraphErrors * g_mean_E1miss_epp_pmiss[4];
  TGraphErrors * g_mean_E2miss_epp_pmiss[4];
  TGraphErrors * g_mean_E1miss_ep_kmiss[4];
  TGraphErrors * g_mean_E1miss_epp_kmiss[4];
  TGraphErrors * g_mean_E2miss_epp_kmiss[4];
  for(int k = 0; k < 4; k++){
    g_sigma_E1miss_ep_pmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_sigma_E1miss_ep_pmiss_%d",k+1);
    g_sigma_E1miss_ep_pmiss[k]->SetName(temp_name);

    g_sigma_E1miss_epp_pmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_sigma_E1miss_epp_pmiss_%d",k+1);
    g_sigma_E1miss_epp_pmiss[k]->SetName(temp_name);

    g_sigma_E2miss_epp_pmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_sigma_E2miss_epp_pmiss_%d",k+1);
    g_sigma_E2miss_epp_pmiss[k]->SetName(temp_name);

    g_sigma_E1miss_ep_kmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_sigma_E1miss_ep_kmiss_%d",k+1);
    g_sigma_E1miss_ep_kmiss[k]->SetName(temp_name);

    g_sigma_E1miss_epp_kmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_sigma_E1miss_epp_kmiss_%d",k+1);
    g_sigma_E1miss_epp_kmiss[k]->SetName(temp_name);

    g_sigma_E2miss_epp_kmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_sigma_E2miss_epp_kmiss_%d",k+1);
    g_sigma_E2miss_epp_kmiss[k]->SetName(temp_name);

    g_mean_E1miss_ep_pmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_mean_E1miss_ep_pmiss_%d",k+1);
    g_mean_E1miss_ep_pmiss[k]->SetName(temp_name);

    g_mean_E1miss_epp_pmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_mean_E1miss_epp_pmiss_%d",k+1);
    g_mean_E1miss_epp_pmiss[k]->SetName(temp_name);

    g_mean_E2miss_epp_pmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_mean_E2miss_epp_pmiss_%d",k+1);
    g_mean_E2miss_epp_pmiss[k]->SetName(temp_name);

    g_mean_E1miss_ep_kmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_mean_E1miss_ep_kmiss_%d",k+1);
    g_mean_E1miss_ep_kmiss[k]->SetName(temp_name);

    g_mean_E1miss_epp_kmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_mean_E1miss_epp_kmiss_%d",k+1);
    g_mean_E1miss_epp_kmiss[k]->SetName(temp_name);

    g_mean_E2miss_epp_kmiss[k] = new TGraphErrors();
    sprintf(temp_name,"g_mean_E2miss_epp_kmiss_%d",k+1);
    g_mean_E2miss_epp_kmiss[k]->SetName(temp_name);
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> sigma_list;
      for(int i = 0; i < 100; i++){sigma_list.push_back(getStdDev(f,myCanvas,fileName,dataGroup_SysErr[i].h_E1miss_ep_SRC_pmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(sigma_list,mean,sigma);
      g_sigma_E1miss_ep_pmiss[k]->SetPoint(g_sigma_E1miss_ep_pmiss[k]->GetN(),Q2_mean[j],mean);
      g_sigma_E1miss_ep_pmiss[k]->SetPointError(g_sigma_E1miss_ep_pmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> sigma_list;
      for(int i = 0; i < 100; i++){sigma_list.push_back(getStdDev(f,myCanvas,fileName,dataGroup_SysErr[i].h_E1miss_epp_SRC_pmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(sigma_list,mean,sigma);
      g_sigma_E1miss_epp_pmiss[k]->SetPoint(g_sigma_E1miss_epp_pmiss[k]->GetN(),Q2_mean[j],mean);
      g_sigma_E1miss_epp_pmiss[k]->SetPointError(g_sigma_E1miss_epp_pmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> sigma_list;
      for(int i = 0; i < 100; i++){sigma_list.push_back(getStdDev(f,myCanvas,fileName,dataGroup_SysErr[i].h_E2miss_epp_SRC_pmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(sigma_list,mean,sigma);
      g_sigma_E2miss_epp_pmiss[k]->SetPoint(g_sigma_E2miss_epp_pmiss[k]->GetN(),Q2_mean[j],mean);
      g_sigma_E2miss_epp_pmiss[k]->SetPointError(g_sigma_E2miss_epp_pmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> sigma_list;
      for(int i = 0; i < 100; i++){sigma_list.push_back(getStdDev(f,myCanvas,fileName,dataGroup_SysErr[i].h_E1miss_ep_SRC_kmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(sigma_list,mean,sigma);
      g_sigma_E1miss_ep_kmiss[k]->SetPoint(g_sigma_E1miss_ep_kmiss[k]->GetN(),Q2_mean[j],mean);
      g_sigma_E1miss_ep_kmiss[k]->SetPointError(g_sigma_E1miss_ep_kmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> sigma_list;
      for(int i = 0; i < 100; i++){sigma_list.push_back(getStdDev(f,myCanvas,fileName,dataGroup_SysErr[i].h_E1miss_epp_SRC_kmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(sigma_list,mean,sigma);
      g_sigma_E1miss_epp_kmiss[k]->SetPoint(g_sigma_E1miss_epp_kmiss[k]->GetN(),Q2_mean[j],mean);
      g_sigma_E1miss_epp_kmiss[k]->SetPointError(g_sigma_E1miss_epp_kmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> sigma_list;
      for(int i = 0; i < 100; i++){sigma_list.push_back(getStdDev(f,myCanvas,fileName,dataGroup_SysErr[i].h_E2miss_epp_SRC_kmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(sigma_list,mean,sigma);
      g_sigma_E2miss_epp_kmiss[k]->SetPoint(g_sigma_E2miss_epp_kmiss[k]->GetN(),Q2_mean[j],mean);
      g_sigma_E2miss_epp_kmiss[k]->SetPointError(g_sigma_E2miss_epp_kmiss[k]->GetN()-1,0,sigma);
    }
  }

    // --- Mean of E_miss vs Q^2 (parallels the sigma loops above) ---
  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> mean_list;
      for(int i = 0; i < 100; i++){mean_list.push_back(getMean(f,myCanvas,fileName,dataGroup_SysErr[i].h_E1miss_ep_SRC_pmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(mean_list,mean,sigma);
      g_mean_E1miss_ep_pmiss[k]->SetPoint(g_mean_E1miss_ep_pmiss[k]->GetN(),Q2_mean[j],mean);
      g_mean_E1miss_ep_pmiss[k]->SetPointError(g_mean_E1miss_ep_pmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> mean_list;
      for(int i = 0; i < 100; i++){mean_list.push_back(getMean(f,myCanvas,fileName,dataGroup_SysErr[i].h_E1miss_epp_SRC_pmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(mean_list,mean,sigma);
      g_mean_E1miss_epp_pmiss[k]->SetPoint(g_mean_E1miss_epp_pmiss[k]->GetN(),Q2_mean[j],mean);
      g_mean_E1miss_epp_pmiss[k]->SetPointError(g_mean_E1miss_epp_pmiss[k]->GetN()-1,0,sigma);
    }
  }

  for(int k = 0; k < 4; k++){
    for(int j=0; j<(bQ2); j++){
      vector<double> mean_list;
      for(int i = 0; i < 100; i++){mean_list.push_back(getMean(f,myCanvas,fileName,dataGroup_SysErr[i].h_E2miss_epp_SRC_pmiss_Q2[k][j],-0.2,0.4));}
      double mean, sigma;
      getMeanStddev(mean_list,mean,sigma);
      g_mean_E2miss_epp_pmiss[k]->SetPoint(g_mean_E2miss_epp_pmiss[k]->GetN(),Q2_mean[j],mean);
      g_mean_E2miss_epp_pmiss[k]->SetPointError(g_mean_E2miss_epp_pmiss[k]->GetN()-1,0,sigma);
    }
  }

  
  for(int k = 0; k < 4; k++){
    g_sigma_E1miss_ep_pmiss[k]->Write();
    g_sigma_E1miss_epp_pmiss[k]->Write();
    g_sigma_E2miss_epp_pmiss[k]->Write();
    g_sigma_E1miss_ep_kmiss[k]->Write();
    g_sigma_E1miss_epp_kmiss[k]->Write();
    g_sigma_E2miss_epp_kmiss[k]->Write();
    
    myCanvas->Divide(2,3);
    myCanvas->cd(1);    
    g_sigma_E1miss_ep_pmiss[k]->Draw("");
    myCanvas->cd(2);    
    g_sigma_E1miss_epp_pmiss[k]->Draw("");
    myCanvas->cd(3);    
    g_sigma_E2miss_epp_pmiss[k]->Draw("");
    myCanvas->cd(4);    
    g_sigma_E1miss_ep_kmiss[k]->Draw("");
    myCanvas->cd(5);    
    g_sigma_E1miss_epp_kmiss[k]->Draw("");
    myCanvas->cd(6);    
    g_sigma_E2miss_epp_kmiss[k]->Draw("");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  for(int k = 0; k < 4; k++){
    g_mean_E1miss_ep_pmiss[k]->Write();
    g_mean_E1miss_epp_pmiss[k]->Write();
    g_mean_E2miss_epp_pmiss[k]->Write();
    g_mean_E1miss_ep_kmiss[k]->Write();
    g_mean_E1miss_epp_kmiss[k]->Write();
    g_mean_E2miss_epp_kmiss[k]->Write();
    
    myCanvas->Divide(2,3);
    myCanvas->cd(1);    
    g_mean_E1miss_ep_pmiss[k]->Draw("");
    myCanvas->cd(2);    
    g_mean_E1miss_epp_pmiss[k]->Draw("");
    myCanvas->cd(3);    
    g_mean_E2miss_epp_pmiss[k]->Draw("");
    myCanvas->cd(4);    
    g_mean_E1miss_ep_kmiss[k]->Draw("");
    myCanvas->cd(5);    
    g_mean_E1miss_epp_kmiss[k]->Draw("");
    myCanvas->cd(6);    
    g_mean_E2miss_epp_kmiss[k]->Draw("");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }
  
  ////////////  
  
  for(int i=0; i<dataGroup.hist_list.size(); i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    dataGroup.hist_list[i]->Draw("colz");    
    for(int j=0; j<100; j++){
      dataGroup_SysErr[j].hist_list[i]->SetLineColor(2);
      dataGroup_SysErr[j].hist_list[i]->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }
  
  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();

  return 0;
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

double getSigma(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * h_myhist, double min, double max){

  if(h_myhist->GetEntries()==0){
    return 0;
  }

  if(h_myhist->GetEntries()<20){
    return h_myhist->GetStdDev();
  }
  
  TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },min,max,3);
    gFit->SetParameter(0,h_myhist->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,(max+min)/2);
    gFit->SetParLimits(1,min,max);
    gFit->SetParameter(2,(max-min)/4);
    gFit->SetParLimits(2,0.00,max-min);    
    TFitResultPtr gPoint = h_myhist->Fit(gFit,"SrBeqn","",min,max);
    /*
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_myhist->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
    */
    return gPoint->Parameter(2);    
}

double getStdDev(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * h_myhist, double min, double max){

  if(h_myhist->GetEntries()==0){
    return 0;
  }  
  return h_myhist->GetStdDev();
}

double getMean(TFile *f, TCanvas * myCanvas, char fileName[100],TH1D * h_myhist, double min, double max){

  if(h_myhist->GetEntries()==0){
    return 0;
  }  
  return h_myhist->GetMean();
}


void runEvent(const std::unique_ptr<clas12::clas12reader>& c12, clas12ana & clasAna, bool isMC, double & qSq, double & x, double & pL, double & tL, double & pR, double & px, double & py, double & pz, double & mM, double & pM, double & kM, double & E0, double & E1, double & E2, bool & passep, bool & passepp, int & ctr){

  passep=false;
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
  if(lead[0]->getRegion()!=FD){return;}

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
  
  passep=true;

  qSq=Q2;
  x = xB;
  pL = mom_lead;
  tL = theta_lead;
  mM = mmiss;
  pM = pmiss;
  kM = kmiss;
  E0 = E0miss;
  E1 = E1miss;
  
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

  passepp=true;

  pR = v_rec.Mag();  
  px = v_cm.Dot(vx);
  py = v_cm.Dot(vy);
  pz = v_cm.Dot(vz);
  E2 = E2miss;

}

void CutRandom(double x, double qSq, double mM, double kM, double pL, double tL, double pR, bool & passep, bool & passepp, bool randomize){
  double xB_lower = 1.2;
  double Q2_lower = 1.5;
  double mMiss_lower = 0.65;
  double mMiss_upper = 1.10;
  double kMiss_lower = 0.3;
  double pLead_lower = 1.0;
  double tLead_upper = 37;
  double pRecoil_lower = 0.3;

  if(randomize){
    TRandom3 * thisRand = new TRandom3(0);;
    xB_lower += thisRand->Gaus(0.0,0.01);
    Q2_lower += thisRand->Gaus(0.0,0.01);

    mMiss_lower += thisRand->Gaus(0.0,0.03);
    mMiss_upper += thisRand->Gaus(0.0,0.03);
    kMiss_lower += thisRand->Gaus(0.0,0.03);
    pLead_lower += thisRand->Gaus(0.0,0.03);
    tLead_upper += thisRand->Gaus(0.0,1.0);

    pRecoil_lower += thisRand->Gaus(0.0,0.045);
  }
  
  if(passep==false){return;}  
  
  if(x<xB_lower){passep=false;}
  if(qSq<Q2_lower){passep=false;}
  if(qSq>5.0){passep=false;}
  if(mM<mMiss_lower){passep=false;}
  if(mM>mMiss_upper){passep=false;}
  if(kM<kMiss_lower){passep=false;}            
  if(pL<pLead_lower){passep=false;}
  if(tL>tLead_upper){passep=false;}
      
  if(passep==false){return;}

  if(pR<pRecoil_lower){passepp=false;}

  
}


//void setUpHistGroup(HistGroup & myGroup,std::string temp_name){
void setUpHistGroup(HistGroup & myGroup){

  char temp_name[100];
  char temp_title[100];

  myGroup.h_Q2_ep = new TH1D("Q2","Q2 ep;Q2;Counts",50,1.5,5);
  myGroup.hist_list.push_back(myGroup.h_Q2_ep);
  myGroup.h_Q2_epp = new TH1D("Q2","Q2 epp;Q2;Counts",50,1.5,5);
  myGroup.hist_list.push_back(myGroup.h_Q2_epp);

  for(int i=0; i<(bQ2); i++){
    sprintf(temp_title,"%f - %f",bE_Q2[i],bE_Q2[i+1]);

    sprintf(temp_name,"h_Q2_ep_SRC_Q2_%d",i);
    myGroup.h_Q2_ep_SRC_Q2[i] = new TH1D(temp_name,temp_title,100,1.5,5.0);
    myGroup.hist_list.push_back(myGroup.h_Q2_ep_SRC_Q2[i]);

    sprintf(temp_name,"h_Q2_epp_SRC_Q2_%d",i);
    myGroup.h_Q2_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,100,1.5,5.0);
    myGroup.hist_list.push_back(myGroup.h_Q2_epp_SRC_Q2[i]);
  }
  
  myGroup.h_pcmx_epp = new TH1D("h_pcmx_epp","pcmx",15,-0.75,0.75);
  myGroup.hist_list.push_back(myGroup.h_pcmx_epp);
  myGroup.h_pcmy_epp = new TH1D("h_pcmy_epp","pcmy",15,-0.75,0.75);
  myGroup.hist_list.push_back(myGroup.h_pcmy_epp);
  myGroup.h_pcmz_epp = new TH1D("h_pcmz_epp","pcmz",15,-0.75,0.75);
  myGroup.hist_list.push_back(myGroup.h_pcmz_epp);

  
  for(int i=0; i<(bQ2); i++){
    sprintf(temp_title,"%f - %f",bE_Q2[i],bE_Q2[i+1]);

    sprintf(temp_name,"h_pcmx_epp_SRC_Q2_%d",i);
    myGroup.h_pcmx_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.75,0.75);
    myGroup.hist_list.push_back(myGroup.h_pcmx_epp_SRC_Q2[i]);

    sprintf(temp_name,"h_pcmy_epp_SRC_Q2_%d",i);
    myGroup.h_pcmy_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.75,0.75);
    myGroup.hist_list.push_back(myGroup.h_pcmy_epp_SRC_Q2[i]);
    
    sprintf(temp_name,"h_pcmz_epp_SRC_Q2_%d",i);
    myGroup.h_pcmz_epp_SRC_Q2[i] = new TH1D(temp_name,temp_title,15,-0.5,1.0);
    myGroup.hist_list.push_back(myGroup.h_pcmz_epp_SRC_Q2[i]);
  }
  

  
  myGroup.h_pMiss_ep = new TH1D("pMiss_ep","p_{Miss} ep;p_{Miss};Counts",bE_pmiss_long.size()-1,&bE_pmiss_long[0]);
  myGroup.hist_list.push_back(myGroup.h_pMiss_ep);
  myGroup.h_pMiss_epp = new TH1D("pMiss_epp","p_{Miss} epp;p_{Miss};Counts",bE_pmiss_long.size()-1,&bE_pmiss_long[0]);
  myGroup.hist_list.push_back(myGroup.h_pMiss_epp);

  myGroup.h_kMiss_ep = new TH1D("kMiss_ep","k_{Miss} ep;k_{Miss};Counts",bE_kmiss_long.size()-1,&bE_kmiss_long[0]);
  myGroup.hist_list.push_back(myGroup.h_kMiss_ep);
  myGroup.h_kMiss_epp = new TH1D("kMiss_epp","k_{Miss} epp;k_{Miss};Counts",bE_kmiss_long.size()-1,&bE_kmiss_long[0]);
  myGroup.hist_list.push_back(myGroup.h_kMiss_epp);
 
  for(int i=0; i<4; i++){
    sprintf(temp_title,"%f - %f",bE_pmiss[i],bE_pmiss[i+1]);
    sprintf(temp_name,"h_Q2_ep_SRC_pmiss_%d",i+1);
    myGroup.h_Q2_ep_SRC_pmiss[i] = new TH1D(temp_name,temp_title,bQ2,&bE_Q2[0]);
    myGroup.hist_list.push_back(myGroup.h_Q2_ep_SRC_pmiss[i]);
    sprintf(temp_name,"h_Q2_epp_SRC_pmiss_%d",i+1);
    myGroup.h_Q2_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,bQ2,&bE_Q2[0]);
    myGroup.hist_list.push_back(myGroup.h_Q2_epp_SRC_pmiss[i]);

    sprintf(temp_title,"%f - %f",bE_kmiss[i],bE_kmiss[i+1]);
    sprintf(temp_name,"h_Q2_ep_SRC_kmiss_%d",i+1);
    myGroup.h_Q2_ep_SRC_kmiss[i] = new TH1D(temp_name,temp_title,bQ2,&bE_Q2[0]);
    myGroup.hist_list.push_back(myGroup.h_Q2_ep_SRC_kmiss[i]);
    sprintf(temp_name,"h_Q2_epp_SRC_kmiss_%d",i+1);
    myGroup.h_Q2_epp_SRC_kmiss[i] = new TH1D(temp_name,temp_title,bQ2,&bE_Q2[0]);
    myGroup.hist_list.push_back(myGroup.h_Q2_epp_SRC_kmiss[i]);
  }

  for(int i=0; i<4; i++){
    sprintf(temp_title,"%f - %f",bE_pmiss[i],bE_pmiss[i+1]);

    sprintf(temp_name,"h_E0miss_ep_SRC_pmiss_%d",i+1);
    myGroup.h_E0miss_ep_SRC_pmiss[i] = new TH1D(temp_name,temp_title,30,-0.1,0.85);
    myGroup.hist_list.push_back(myGroup.h_E0miss_ep_SRC_pmiss[i]);
    sprintf(temp_name,"h_E1miss_ep_SRC_pmiss_%d",i+1);
    myGroup.h_E1miss_ep_SRC_pmiss[i] = new TH1D(temp_name,temp_title,15,-0.25,0.55);
    myGroup.hist_list.push_back(myGroup.h_E1miss_ep_SRC_pmiss[i]);
    sprintf(temp_name,"h_E1miss_epp_SRC_pmiss_%d",i+1);
    myGroup.h_E1miss_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,15,-0.2,0.5);
    myGroup.hist_list.push_back(myGroup.h_E1miss_epp_SRC_pmiss[i]);
    sprintf(temp_name,"h_E2miss_epp_SRC_pmiss_%d",i+1);
    myGroup.h_E2miss_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,15,-0.25,0.5);
    myGroup.hist_list.push_back(myGroup.h_E2miss_epp_SRC_pmiss[i]);


    sprintf(temp_title,"%f - %f",bE_kmiss[i],bE_kmiss[i+1]);

    sprintf(temp_name,"h_E0miss_ep_SRC_kmiss_%d",i+1);
    myGroup.h_E0miss_ep_SRC_kmiss[i] = new TH1D(temp_name,temp_title,30,-0.1,0.85);
    myGroup.hist_list.push_back(myGroup.h_E0miss_ep_SRC_kmiss[i]);
    sprintf(temp_name,"h_E1miss_ep_SRC_kmiss_%d",i+1);
    myGroup.h_E1miss_ep_SRC_kmiss[i] = new TH1D(temp_name,temp_title,15,-0.25,0.55);
    myGroup.hist_list.push_back(myGroup.h_E1miss_ep_SRC_kmiss[i]);
    sprintf(temp_name,"h_E1miss_epp_SRC_kmiss_%d",i+1);
    myGroup.h_E1miss_epp_SRC_kmiss[i] = new TH1D(temp_name,temp_title,15,-0.2,0.5);
    myGroup.hist_list.push_back(myGroup.h_E1miss_epp_SRC_kmiss[i]);
    sprintf(temp_name,"h_E2miss_epp_SRC_kmiss_%d",i+1);
    myGroup.h_E2miss_epp_SRC_kmiss[i] = new TH1D(temp_name,temp_title,15,-0.25,0.5);
    myGroup.hist_list.push_back(myGroup.h_E2miss_epp_SRC_kmiss[i]);
  }


  for(int i=0; i<4; i++){
    for(int j=0; j<7; j++){
      sprintf(temp_title,"%f - %f",bE_pmiss[i],bE_pmiss[i+1]);
      
      sprintf(temp_name,"h_E1miss_ep_SRC_pmiss_%d_Q2_%d",i+1,j+1);
      myGroup.h_E1miss_ep_SRC_pmiss_Q2[i][j] = new TH1D(temp_name,temp_title,15,-0.25,0.55);
      myGroup.hist_list.push_back(myGroup.h_E1miss_ep_SRC_pmiss_Q2[i][j]);
      sprintf(temp_name,"h_E1miss_epp_SRC_pmiss_%d_Q2_%d",i+1,j+1);
      myGroup.h_E1miss_epp_SRC_pmiss_Q2[i][j] = new TH1D(temp_name,temp_title,15,-0.2,0.5);
      myGroup.hist_list.push_back(myGroup.h_E1miss_epp_SRC_pmiss_Q2[i][j]);
      sprintf(temp_name,"h_E2miss_epp_SRC_pmiss_%d_Q2_%d",i+1,j+1);
      myGroup.h_E2miss_epp_SRC_pmiss_Q2[i][j] = new TH1D(temp_name,temp_title,15,-0.25,0.5);
      myGroup.hist_list.push_back(myGroup.h_E2miss_epp_SRC_pmiss_Q2[i][j]);
      

      sprintf(temp_title,"%f - %f",bE_kmiss[i],bE_kmiss[i+1]);
      
      sprintf(temp_name,"h_E1miss_ep_SRC_kmiss_%d_Q2_%d",i+1,j+1);
      myGroup.h_E1miss_ep_SRC_kmiss_Q2[i][j] = new TH1D(temp_name,temp_title,15,-0.25,0.55);
      myGroup.hist_list.push_back(myGroup.h_E1miss_ep_SRC_kmiss_Q2[i][j]);
      sprintf(temp_name,"h_E1miss_epp_SRC_kmiss_%d_Q2_%d",i+1,j+1);
      myGroup.h_E1miss_epp_SRC_kmiss_Q2[i][j] = new TH1D(temp_name,temp_title,15,-0.2,0.5);
      myGroup.hist_list.push_back(myGroup.h_E1miss_epp_SRC_kmiss_Q2[i][j]);
      sprintf(temp_name,"h_E2miss_epp_SRC_kmiss_%d_Q2_%d",i+1,j+1);
      myGroup.h_E2miss_epp_SRC_kmiss_Q2[i][j] = new TH1D(temp_name,temp_title,15,-0.25,0.5);
      myGroup.hist_list.push_back(myGroup.h_E2miss_epp_SRC_kmiss_Q2[i][j]);
    }
  }
  
  for(int i=0; i<myGroup.hist_list.size(); i++){
    myGroup.hist_list[i]->Sumw2();
    myGroup.hist_list[i]->GetXaxis()->CenterTitle();
    myGroup.hist_list[i]->GetYaxis()->CenterTitle();
  }

  
}


void fillUpHistGroup(HistGroup & myGroup,double qSq,double px,double py,double pz,double pM,double kM,double E0,double E1,double E2,bool passep,bool passepp,double wep,double wepp){
        
  if(!passep){return;}
  
  myGroup.h_Q2_ep->Fill(qSq,wep);
  myGroup.h_pMiss_ep->Fill(pM,wep);
  myGroup.h_kMiss_ep->Fill(kM,wep);
  if(binX(bE_pmiss,pM)!=-1){
    myGroup.h_Q2_ep_SRC_pmiss[binX(bE_pmiss,pM)]->Fill(qSq,wep);
    myGroup.h_E0miss_ep_SRC_pmiss[binX(bE_pmiss,pM)]->Fill(E0,wep);
    myGroup.h_E1miss_ep_SRC_pmiss[binX(bE_pmiss,pM)]->Fill(E1,wep);
  }
  if(binX(bE_kmiss,kM)!=-1){
    myGroup.h_Q2_ep_SRC_kmiss[binX(bE_kmiss,kM)]->Fill(qSq,wep);
    myGroup.h_E0miss_ep_SRC_kmiss[binX(bE_kmiss,kM)]->Fill(E0,wep);
    myGroup.h_E1miss_ep_SRC_kmiss[binX(bE_kmiss,kM)]->Fill(E1,wep);
  }
  if(binX(bE_Q2,qSq)!=-1){
    myGroup.h_Q2_ep_SRC_Q2[binX(bE_Q2,qSq)]->Fill(qSq,wepp);
  }
  if(binX(bE_Q2,qSq)!=-1){
    if(binX(bE_pmiss,pM)!=-1){
      myGroup.h_E1miss_ep_SRC_pmiss_Q2[binX(bE_pmiss,pM)][binX(bE_Q2,qSq)]->Fill(E1,wepp);
    }
    if(binX(bE_kmiss,kM)!=-1){
      myGroup.h_E1miss_ep_SRC_kmiss_Q2[binX(bE_kmiss,kM)][binX(bE_Q2,qSq)]->Fill(E1,wepp);
    }
  }


  if(!passepp){return;}

  myGroup.h_Q2_epp->Fill(qSq,wepp);
  myGroup.h_pMiss_epp->Fill(pM,wepp);
  myGroup.h_kMiss_epp->Fill(kM,wepp);      
  if(binX(bE_pmiss,pM)!=-1){
    myGroup.h_Q2_epp_SRC_pmiss[binX(bE_pmiss,pM)]->Fill(qSq,wepp);
    myGroup.h_E1miss_epp_SRC_pmiss[binX(bE_pmiss,pM)]->Fill(E1,wepp);
    myGroup.h_E2miss_epp_SRC_pmiss[binX(bE_pmiss,pM)]->Fill(E2,wepp);
  }
  if(binX(bE_kmiss,kM)!=-1){
    myGroup.h_Q2_epp_SRC_kmiss[binX(bE_kmiss,kM)]->Fill(qSq,wepp);
    myGroup.h_E1miss_epp_SRC_kmiss[binX(bE_kmiss,kM)]->Fill(E1,wepp);
    myGroup.h_E2miss_epp_SRC_kmiss[binX(bE_kmiss,kM)]->Fill(E2,wepp);
  }
  if(binX(bE_Q2,qSq)!=-1){
    myGroup.h_Q2_epp_SRC_Q2[binX(bE_Q2,qSq)]->Fill(qSq,wepp);
  }

  myGroup.h_pcmx_epp->Fill(px,wepp);
  myGroup.h_pcmy_epp->Fill(py,wepp);
  myGroup.h_pcmz_epp->Fill(pz,wepp);
  if(binX(bE_Q2,qSq)!=-1){
    myGroup.h_pcmx_epp_SRC_Q2[binX(bE_Q2,qSq)]->Fill(px,wepp);
    myGroup.h_pcmy_epp_SRC_Q2[binX(bE_Q2,qSq)]->Fill(py,wepp);
    myGroup.h_pcmz_epp_SRC_Q2[binX(bE_Q2,qSq)]->Fill(pz,wepp);
  }

  if(binX(bE_Q2,qSq)!=-1){
    if(binX(bE_pmiss,pM)!=-1){
      myGroup.h_E1miss_epp_SRC_pmiss_Q2[binX(bE_pmiss,pM)][binX(bE_Q2,qSq)]->Fill(E1,wepp);
      myGroup.h_E2miss_epp_SRC_pmiss_Q2[binX(bE_pmiss,pM)][binX(bE_Q2,qSq)]->Fill(E2,wepp);
    }
    if(binX(bE_kmiss,kM)!=-1){
      myGroup.h_E1miss_epp_SRC_kmiss_Q2[binX(bE_kmiss,kM)][binX(bE_Q2,qSq)]->Fill(E1,wepp);
      myGroup.h_E2miss_epp_SRC_kmiss_Q2[binX(bE_kmiss,kM)][binX(bE_Q2,qSq)]->Fill(E2,wepp);
    }
  }
  
}

TGraphErrors * getGraphWithError(TFile *f, TCanvas * myCanvas, char fileName[100], TH1D * nominal, vector<TH1*> hist_Errors){
  
  TGraphErrors * g_mu = new TGraphErrors();
  g_mu->SetName(nominal->GetName());
  for(int j = 0; j < hist_Errors[0]->GetXaxis()->GetNbins(); j++){
    double mean, stddev;
    double x = hist_Errors[0]->GetXaxis()->GetBinCenter(j+1);
    vector<double> list;
    for(int i = 0; i < 100; i++){
      double y = hist_Errors[i]->GetBinContent(j+1);
      list.push_back(y);
    }
    getMeanStddev(list,mean,stddev);
    double stat_err = nominal->GetBinError(j);
    g_mu->SetPoint(g_mu->GetN(),x,mean);
    g_mu->SetPointError(g_mu->GetN()-1,0,sqrt(stddev*stddev + stat_err*stat_err));
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  //nominal->Draw();
  g_mu->Draw("");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  g_mu->Write();
  
  return g_mu;
}
