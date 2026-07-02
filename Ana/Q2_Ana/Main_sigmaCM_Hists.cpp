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
#include "TGraphAsymmErrors.h"
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
const double min_sigma = 0.08;
const double max_sigma = 0.25;
// .05, .2 were orig
//////////////////////////////////

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
void getValue(TGraph * thisGraph,double min_sigma, double max_sigma, double & chi2NDF, double & center, double & lower, double & upper);
int getMinBin(double sigma){
  // double bin_double = (double)linbin*(sigma-min_sigma)/(max_sigma-min_sigma);
  // int bin_truncated = bin_double;
  // if(bin_double-bin_truncated<0.5){return bin_truncated;}
  // return bin_truncated+1;

  double bin_double = (double)linbin*(sigma-min_sigma)/(max_sigma-min_sigma);
  int bin_truncated = bin_double;
  int result;
  if(bin_double - bin_truncated < 0.5){ result = bin_truncated; }
  else { result = bin_truncated + 1; }
  // Clamp to valid range
  if(result < 0) result = 0;
  if(result >= linbin) result = linbin - 1;
  return result;
}

void Usage()
{
  std::cerr << "Usage: ./code outputfile.root outputfile.pdf [CutRangeUpper] inputfile.root  \n\n\n";
}

int main(int argc, char ** argv)
{
  
  if(argc < 4)
    {
      Usage();
      return -1;
    }
  
  TString outFile = argv[1];
  char * pdfFile = argv[2];
  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;
  int intcutRange = atoi(argv[3]);
  double cutRange = (double)intcutRange/100.0; // 
  TFile * inFile = new TFile(argv[4]);

  
  char temp_name[100];
  char temp_title[100];

  TH1D * h_Q2[bQ2];
  for(int i=0; i<(bQ2); i++){
    sprintf(temp_title,"%f - %f",bE_Q2[i],bE_Q2[i+1]);

    sprintf(temp_name,"h_Q2_%d",i);
    h_Q2[i] = (TH1D*)inFile->Get(temp_name);
  }

  TH1D * h_pcmx_epp = (TH1D*)inFile->Get("pcmx_epp");
  TH1D * h_pcmy_epp = (TH1D*)inFile->Get("pcmy_epp");
  TH1D * h_pcmz_epp = (TH1D*)inFile->Get("pcmz_epp");
  TH1D * h_pcmT_epp = (TH1D*)inFile->Get("pcmT_epp");
  
  TH1D * h_pcmx_epp_SRC_Q2[bQ2];
  TH1D * h_pcmy_epp_SRC_Q2[bQ2];
  TH1D * h_pcmz_epp_SRC_Q2[bQ2];
  TH1D * h_pcmT_epp_SRC_Q2[bQ2];
  for(int i=0; i<(bQ2); i++){
    sprintf(temp_title,"%f - %f",bE_Q2[i],bE_Q2[i+1]);

    sprintf(temp_name,"h_pcmx_epp_SRC_Q2_%d",i);
    h_pcmx_epp_SRC_Q2[i] = (TH1D*)inFile->Get(temp_name);

    sprintf(temp_name,"h_pcmy_epp_SRC_Q2_%d",i);
    h_pcmy_epp_SRC_Q2[i] = (TH1D*)inFile->Get(temp_name);
    
    sprintf(temp_name,"h_pcmz_epp_SRC_Q2_%d",i);
    h_pcmz_epp_SRC_Q2[i] = (TH1D*)inFile->Get(temp_name);

    sprintf(temp_name,"h_pcmT_epp_SRC_Q2_%d",i);
    h_pcmT_epp_SRC_Q2[i] = (TH1D*)inFile->Get(temp_name);
  }



  double h_pcmx_epp_simSCM_scale[linbin];
  double h_pcmy_epp_simSCM_scale[linbin];
  double h_pcmz_epp_simSCM_scale[linbin];
  double h_pcmT_epp_simSCM_scale[linbin];
  TH1D * h_pcmx_epp_simSCM[linbin];
  TH1D * h_pcmy_epp_simSCM[linbin];
  TH1D * h_pcmz_epp_simSCM[linbin];
  TH1D * h_pcmT_epp_simSCM[linbin];
  
  double h_pcmx_epp_SRC_simSCM_Q2_scale[linbin][bQ2];
  double h_pcmy_epp_SRC_simSCM_Q2_scale[linbin][bQ2];
  double h_pcmz_epp_SRC_simSCM_Q2_scale[linbin][bQ2];
  double h_pcmT_epp_SRC_simSCM_Q2_scale[linbin][bQ2];
  TH1D * h_pcmx_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmy_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmz_epp_SRC_simSCM_Q2[linbin][bQ2];
  TH1D * h_pcmT_epp_SRC_simSCM_Q2[linbin][bQ2];
  for(int j=0; j<linbin; j++){
    sprintf(temp_name,"h_pcmx_epp_simSCM_%d",j);
    h_pcmx_epp_simSCM[j] = (TH1D*)inFile->Get(temp_name);

    sprintf(temp_name,"h_pcmy_epp_simSCM_%d",j);
    h_pcmy_epp_simSCM[j] = (TH1D*)inFile->Get(temp_name);

    sprintf(temp_name,"h_pcmz_epp_simSCM_%d",j);
    h_pcmz_epp_simSCM[j] = (TH1D*)inFile->Get(temp_name);

    sprintf(temp_name,"h_pcmT_epp_simSCM_%d",j);
    h_pcmT_epp_simSCM[j] = (TH1D*)inFile->Get(temp_name);

    for(int i=0; i<(bQ2); i++){
      int sCM= ((double)j/(double)linbin)*(max_sigma-min_sigma) + min_sigma;
      sprintf(temp_title,"sCM=%f %f - %f",sCM,bE_Q2[i],bE_Q2[i+1]);
      
      sprintf(temp_name,"h_pcmx_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmx_epp_SRC_simSCM_Q2[j][i] = (TH1D*)inFile->Get(temp_name);
      
      sprintf(temp_name,"h_pcmy_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmy_epp_SRC_simSCM_Q2[j][i] = (TH1D*)inFile->Get(temp_name);
      
      sprintf(temp_name,"h_pcmz_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmz_epp_SRC_simSCM_Q2[j][i] = (TH1D*)inFile->Get(temp_name);

      sprintf(temp_name,"h_pcmT_epp_SRC_simSCM_Q2_%d_%d",j,i);
      h_pcmT_epp_SRC_simSCM_Q2[j][i] = (TH1D*)inFile->Get(temp_name);
    }
  }


  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  TGraph * g_chi2_pcmx_epp = new TGraph;
  g_chi2_pcmx_epp->SetName("g_chi2_pcmx_epp");
  TGraph * g_chi2_pcmy_epp = new TGraph;
  g_chi2_pcmy_epp->SetName("g_chi2_pcmy_epp");
  TGraph * g_chi2_pcmz_epp = new TGraph;
  g_chi2_pcmz_epp->SetName("g_chi2_pcmz_epp");
  TGraph * g_chi2_pcmT_epp = new TGraph;
  g_chi2_pcmT_epp->SetName("g_chi2_pcmT_epp");
  TGraph * g_scale_pcmx_epp = new TGraph;
  g_scale_pcmx_epp->SetName("g_scale_pcmx_epp");
  TGraph * g_scale_pcmy_epp = new TGraph;
  g_scale_pcmy_epp->SetName("g_scale_pcmy_epp");
  TGraph * g_scale_pcmz_epp = new TGraph;
  g_scale_pcmz_epp->SetName("g_scale_pcmz_epp");
  TGraph * g_scale_pcmT_epp = new TGraph;
  g_scale_pcmT_epp->SetName("g_scale_pcmT_epp");
  
  TGraph * g_chi2_pcmx_epp_SRC_Q2[bQ2];
  TGraph * g_chi2_pcmy_epp_SRC_Q2[bQ2];
  TGraph * g_chi2_pcmz_epp_SRC_Q2[bQ2];
  TGraph * g_chi2_pcmT_epp_SRC_Q2[bQ2];

  TGraph * g_scale_pcmx_epp_SRC_Q2[bQ2];
  TGraph * g_scale_pcmy_epp_SRC_Q2[bQ2];
  TGraph * g_scale_pcmz_epp_SRC_Q2[bQ2];
  TGraph * g_scale_pcmT_epp_SRC_Q2[bQ2];
  for(int i=0; i<(bQ2); i++){
    g_chi2_pcmx_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_chi2_pcmx_epp_SRC_Q2_%d",i);
    g_chi2_pcmx_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM X Q2bin=%d;#sigma_{cmx};#chi^2",i);
    g_chi2_pcmx_epp_SRC_Q2[i]->SetTitle(temp_title);

    g_chi2_pcmy_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_chi2_pcmy_epp_SRC_Q2_%d",i);
    g_chi2_pcmy_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM Y Q2bin=%d;#sigma_{cmy};#chi^2",i);
    g_chi2_pcmy_epp_SRC_Q2[i]->SetTitle(temp_title);
    
    g_chi2_pcmz_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_chi2_pcmz_epp_SRC_Q2_%d",i);
    g_chi2_pcmz_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM Z Q2bin=%d;#sigma_{cmz};#chi^2",i);
    g_chi2_pcmz_epp_SRC_Q2[i]->SetTitle(temp_title);

    g_chi2_pcmT_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_chi2_pcmT_epp_SRC_Q2_%d",i);
    g_chi2_pcmT_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM T Q2bin=%d;#sigma_{cmT};#chi^2",i);
    g_chi2_pcmT_epp_SRC_Q2[i]->SetTitle(temp_title);

    /////////
    
    g_scale_pcmx_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_scale_pcmx_epp_SRC_Q2_%d",i);
    g_scale_pcmx_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM X Q2bin=%d;#sigma_{cmx};Scale",i);
    g_scale_pcmx_epp_SRC_Q2[i]->SetTitle(temp_title);

    g_scale_pcmy_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_scale_pcmy_epp_SRC_Q2_%d",i);
    g_scale_pcmy_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM Y Q2bin=%d;#sigma_{cmy};Scale",i);
    g_scale_pcmy_epp_SRC_Q2[i]->SetTitle(temp_title);

    g_scale_pcmz_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_scale_pcmz_epp_SRC_Q2_%d",i);
    g_scale_pcmz_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM Z Q2bin=%d;#sigma_{cmz};Scale",i);
    g_scale_pcmz_epp_SRC_Q2[i]->SetTitle(temp_title);

    g_scale_pcmT_epp_SRC_Q2[i] = new TGraph;
    sprintf(temp_name,"g_scale_pcmT_epp_SRC_Q2_%d",i);
    g_scale_pcmT_epp_SRC_Q2[i]->SetName(temp_name);
    sprintf(temp_title,"pCM T Q2bin=%d;#sigma_{cmT};Scale",i);
    g_scale_pcmT_epp_SRC_Q2[i]->SetTitle(temp_title);

  }

  for(int j=0; j<linbin; j++){
    cout<<"Working on Bin="<<j<<endl;
    double sCM= ((double)j/(double)linbin)*(max_sigma-min_sigma) + min_sigma;

    double chi2_x_int, chi2_y_int, chi2_z_int, chi2_T_int, scale_x_int, scale_y_int, scale_z_int, scale_T_int;
    
    getChi2(h_pcmx_epp,h_pcmx_epp_simSCM[j],-cutRange/sqrt(2.0),cutRange/sqrt(2.0),scale_x_int,chi2_x_int);
    h_pcmx_epp_simSCM_scale[j]=scale_x_int;
    g_chi2_pcmx_epp->SetPoint(g_chi2_pcmx_epp->GetN(),sCM,chi2_x_int);
    g_scale_pcmx_epp->SetPoint(g_scale_pcmx_epp->GetN(),sCM,scale_x_int);

    getChi2(h_pcmy_epp,h_pcmy_epp_simSCM[j],-cutRange/sqrt(2.0),cutRange/sqrt(2.0),scale_y_int,chi2_y_int);
    h_pcmy_epp_simSCM_scale[j]=scale_y_int;
    g_chi2_pcmy_epp->SetPoint(g_chi2_pcmy_epp->GetN(),sCM,chi2_y_int);
    g_scale_pcmy_epp->SetPoint(g_scale_pcmy_epp->GetN(),sCM,scale_y_int);

    getChi2(h_pcmz_epp,h_pcmz_epp_simSCM[j],-cutRange/sqrt(2.0),cutRange/sqrt(2.0),scale_z_int,chi2_z_int);
    h_pcmz_epp_simSCM_scale[j]=scale_z_int;
    g_chi2_pcmz_epp->SetPoint(g_chi2_pcmz_epp->GetN(),sCM,chi2_z_int);
    g_scale_pcmz_epp->SetPoint(g_scale_pcmz_epp->GetN(),sCM,scale_z_int);
    
    getChi2(h_pcmT_epp,h_pcmT_epp_simSCM[j],0.0,cutRange,scale_T_int,chi2_T_int);
    h_pcmT_epp_simSCM_scale[j]=scale_T_int;     
    g_chi2_pcmT_epp->SetPoint(g_chi2_pcmT_epp->GetN(),sCM,chi2_T_int);    
    g_scale_pcmT_epp->SetPoint(g_scale_pcmT_epp->GetN(),sCM,scale_T_int);
    
    for(int i=0; i<(bQ2); i++){
      double chi2_x, chi2_y, chi2_z, chi2_T, scale_x, scale_y, scale_z, scale_T;
      
      getChi2(h_pcmx_epp_SRC_Q2[i],h_pcmx_epp_SRC_simSCM_Q2[j][i],-0.2,0.2,scale_x,chi2_x);
      h_pcmx_epp_SRC_simSCM_Q2_scale[j][i]=scale_x;
      g_chi2_pcmx_epp_SRC_Q2[i]->SetPoint(g_chi2_pcmx_epp_SRC_Q2[i]->GetN(),sCM,chi2_x);
      g_scale_pcmx_epp_SRC_Q2[i]->SetPoint(g_scale_pcmx_epp_SRC_Q2[i]->GetN(),sCM,scale_x);

      getChi2(h_pcmy_epp_SRC_Q2[i],h_pcmy_epp_SRC_simSCM_Q2[j][i],-0.2,0.2,scale_y,chi2_y);
      h_pcmy_epp_SRC_simSCM_Q2_scale[j][i]=scale_y;
      g_chi2_pcmy_epp_SRC_Q2[i]->SetPoint(g_chi2_pcmy_epp_SRC_Q2[i]->GetN(),sCM,chi2_y);
      g_scale_pcmy_epp_SRC_Q2[i]->SetPoint(g_scale_pcmy_epp_SRC_Q2[i]->GetN(),sCM,scale_y);

      getChi2(h_pcmz_epp_SRC_Q2[i],h_pcmz_epp_SRC_simSCM_Q2[j][i],-0.2,0.2,scale_z,chi2_z);
      h_pcmz_epp_SRC_simSCM_Q2_scale[j][i]=scale_z;
      g_chi2_pcmz_epp_SRC_Q2[i]->SetPoint(g_chi2_pcmz_epp_SRC_Q2[i]->GetN(),sCM,chi2_z);
      g_scale_pcmz_epp_SRC_Q2[i]->SetPoint(g_scale_pcmz_epp_SRC_Q2[i]->GetN(),sCM,scale_z);
      
      getChi2(h_pcmT_epp_SRC_Q2[i],h_pcmT_epp_SRC_simSCM_Q2[j][i],0.0,cutRange,scale_T,chi2_T);
      h_pcmT_epp_SRC_simSCM_Q2_scale[j][i]=scale_T;      
      g_chi2_pcmT_epp_SRC_Q2[i]->SetPoint(g_chi2_pcmT_epp_SRC_Q2[i]->GetN(),sCM,chi2_T);
      g_scale_pcmT_epp_SRC_Q2[i]->SetPoint(g_scale_pcmT_epp_SRC_Q2[i]->GetN(),sCM,scale_T);
      
    }
  }

  double pcmx_chiNDF_int, pcmx_cent_int, pcmx_lower_err_int, pcmx_upper_err_int;
  getValue(g_chi2_pcmx_epp,min_sigma,max_sigma,pcmx_chiNDF_int,pcmx_cent_int,pcmx_lower_err_int,pcmx_upper_err_int);  
  TGraphAsymmErrors * g_sigmacmx_int = new TGraphAsymmErrors;
  g_sigmacmx_int->SetName("sigmacmx_int");
  g_sigmacmx_int->SetPoint(g_sigmacmx_int->GetN(),1.0,pcmx_cent_int);
  g_sigmacmx_int->SetPointError(g_sigmacmx_int->GetN()-1,0.0,2.0,pcmx_cent_int-pcmx_lower_err_int,pcmx_upper_err_int-pcmx_cent_int);

  double pcmy_chiNDF_int, pcmy_cent_int, pcmy_lower_err_int, pcmy_upper_err_int;
  getValue(g_chi2_pcmy_epp,min_sigma,max_sigma,pcmy_chiNDF_int,pcmy_cent_int,pcmy_lower_err_int,pcmy_upper_err_int);  
  TGraphAsymmErrors * g_sigmacmy_int = new TGraphAsymmErrors;
  g_sigmacmy_int->SetName("sigmacmy_int");
  g_sigmacmy_int->SetPoint(g_sigmacmy_int->GetN(),1.0,pcmy_cent_int);
  g_sigmacmy_int->SetPointError(g_sigmacmy_int->GetN()-1,0.0,2.0,pcmy_cent_int-pcmy_lower_err_int,pcmy_upper_err_int-pcmy_cent_int);

  double pcmz_chiNDF_int, pcmz_cent_int, pcmz_lower_err_int, pcmz_upper_err_int;
  getValue(g_chi2_pcmz_epp,min_sigma,max_sigma,pcmz_chiNDF_int,pcmz_cent_int,pcmz_lower_err_int,pcmz_upper_err_int);  
  TGraphAsymmErrors * g_sigmacmz_int = new TGraphAsymmErrors;
  g_sigmacmz_int->SetName("sigmacmz_int");
  g_sigmacmz_int->SetPoint(g_sigmacmz_int->GetN(),1.0,pcmz_cent_int);
  g_sigmacmz_int->SetPointError(g_sigmacmz_int->GetN()-1,0.0,2.0,pcmz_cent_int-pcmz_lower_err_int,pcmz_upper_err_int-pcmz_cent_int);

  double pcmT_chiNDF_int, pcmT_cent_int, pcmT_lower_err_int, pcmT_upper_err_int;
  getValue(g_chi2_pcmT_epp,min_sigma,max_sigma,pcmT_chiNDF_int,pcmT_cent_int,pcmT_lower_err_int,pcmT_upper_err_int);  
  TGraphAsymmErrors * g_sigmacmT_int = new TGraphAsymmErrors;
  g_sigmacmT_int->SetName("sigmacmT_int");
  g_sigmacmT_int->SetPoint(g_sigmacmT_int->GetN(),1.0,pcmT_cent_int);
  g_sigmacmT_int->SetPointError(g_sigmacmT_int->GetN()-1,0.0,2.0,pcmT_cent_int-pcmT_lower_err_int,pcmT_upper_err_int-pcmT_cent_int);
  
  ////For Q2 bins///
  TGraphAsymmErrors * g_sigmacmT_Q2 = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_sigmacmx_Q2 = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_sigmacmy_Q2 = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_sigmacmz_Q2 = new TGraphAsymmErrors;


  g_sigmacmT_Q2->SetName("sigmacmT_Q2");
  g_sigmacmx_Q2->SetName("sigmacmx_Q2");
  g_sigmacmy_Q2->SetName("sigmacmy_Q2");
  g_sigmacmz_Q2->SetName("sigmacmz_Q2");
  double pcmT_cent_Q2bin[bQ2];
  double pcmx_cent_Q2bin[bQ2];
  double pcmy_cent_Q2bin[bQ2];
  double pcmz_cent_Q2bin[bQ2];
  for(int i=0; i<(bQ2); i++){    
    double thisQ2_center = h_Q2[i]->GetMean();
    double thisQ2_lower = thisQ2_center - bE_Q2[i];
    double thisQ2_upper = bE_Q2[i+1]-thisQ2_center;
    
    double pcmT_chiNDF, pcmT_cent, pcmT_lower_err, pcmT_upper_err;
    getValue(g_chi2_pcmT_epp_SRC_Q2[i],min_sigma,max_sigma,pcmT_chiNDF,pcmT_cent,pcmT_lower_err,pcmT_upper_err);      
    pcmT_cent_Q2bin[i]=pcmT_cent;
    g_sigmacmT_Q2->SetPoint(g_sigmacmT_Q2->GetN(),thisQ2_center,pcmT_cent);
    g_sigmacmT_Q2->SetPointError(g_sigmacmT_Q2->GetN()-1,thisQ2_lower,thisQ2_upper,pcmT_cent-pcmT_lower_err,pcmT_upper_err-pcmT_cent);
 
    double pcmx_chiNDF, pcmx_cent, pcmx_lower_err, pcmx_upper_err;
    getValue(g_chi2_pcmx_epp_SRC_Q2[i],min_sigma,max_sigma,pcmx_chiNDF,pcmx_cent,pcmx_lower_err,pcmx_upper_err);      
    pcmx_cent_Q2bin[i]=pcmx_cent;
    g_sigmacmx_Q2->SetPoint(g_sigmacmx_Q2->GetN(),thisQ2_center,pcmx_cent);
    g_sigmacmx_Q2->SetPointError(g_sigmacmx_Q2->GetN()-1,thisQ2_lower,thisQ2_upper,pcmx_cent-pcmx_lower_err,pcmx_upper_err-pcmx_cent);

    double pcmy_chiNDF, pcmy_cent, pcmy_lower_err, pcmy_upper_err;
    getValue(g_chi2_pcmy_epp_SRC_Q2[i],min_sigma,max_sigma,pcmy_chiNDF,pcmy_cent,pcmy_lower_err,pcmy_upper_err);      
    pcmy_cent_Q2bin[i]=pcmy_cent;
    g_sigmacmy_Q2->SetPoint(g_sigmacmy_Q2->GetN(),thisQ2_center,pcmy_cent);
    g_sigmacmy_Q2->SetPointError(g_sigmacmy_Q2->GetN()-1,thisQ2_lower,thisQ2_upper,pcmy_cent-pcmy_lower_err,pcmy_upper_err-pcmy_cent);    

    double pcmz_chiNDF, pcmz_cent, pcmz_lower_err, pcmz_upper_err;
    getValue(g_chi2_pcmz_epp_SRC_Q2[i],min_sigma,max_sigma,pcmz_chiNDF,pcmz_cent,pcmz_lower_err,pcmz_upper_err);
    pcmz_cent_Q2bin[i]=pcmz_cent;
    g_sigmacmz_Q2->SetPoint(g_sigmacmz_Q2->GetN(),thisQ2_center,pcmz_cent);
    g_sigmacmz_Q2->SetPointError(g_sigmacmz_Q2->GetN()-1,thisQ2_lower,thisQ2_upper,pcmz_cent-pcmz_lower_err,pcmz_upper_err-pcmz_cent);
 
  }

  
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();

  h_pcmx_epp->Write();
  TH1D * h_x_clone = (TH1D*) h_pcmx_epp_simSCM[getMinBin(pcmx_cent_int)]->Clone();
  h_x_clone->Scale(h_pcmx_epp_simSCM_scale[getMinBin(pcmx_cent_int)]);
  h_x_clone->SetName("pcmx_epp_fit");
  h_x_clone->Write();
  g_sigmacmx_int->Write();

  h_pcmy_epp->Write();
  TH1D * h_y_clone = (TH1D*) h_pcmy_epp_simSCM[getMinBin(pcmy_cent_int)]->Clone();
  h_y_clone->Scale(h_pcmy_epp_simSCM_scale[getMinBin(pcmy_cent_int)]);
  h_y_clone->SetName("pcmy_epp_fit");
  h_y_clone->Write();
  g_sigmacmy_int->Write();
  
  h_pcmz_epp->Write();
  TH1D * h_z_clone = (TH1D*) h_pcmz_epp_simSCM[getMinBin(pcmz_cent_int)]->Clone();
  h_z_clone->Scale(h_pcmz_epp_simSCM_scale[getMinBin(pcmz_cent_int)]);
  h_z_clone->SetName("pcmz_epp_fit");
  h_z_clone->Write();
  g_sigmacmz_int->Write();

  h_pcmT_epp->Write();
  TH1D * h_T_clone = (TH1D*) h_pcmT_epp_simSCM[getMinBin(pcmT_cent_int)]->Clone();
  h_T_clone->Scale(h_pcmT_epp_simSCM_scale[getMinBin(pcmT_cent_int)]);
  h_T_clone->SetName("pcmT_epp_fit");
  h_T_clone->Write();
  g_sigmacmT_int->Write();

  g_sigmacmT_Q2->Write();
  g_sigmacmx_Q2->Write();
  g_sigmacmy_Q2->Write();
  g_sigmacmz_Q2->Write();
  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  // char fileName[100];
  // sprintf(fileName,"%s[",pdfFile);
  // myText->SaveAs(fileName);
  // sprintf(fileName,"%s",pdfFile);

  // helper lambda: write a named canvas snapshot to the ROOT file
  // (data hist already drawn on myCanvas; clone it into a dedicated TCanvas)
  auto writeCanvas = [&](const char* cname, const char* ctitle){
    f->cd();
    TCanvas * cOut = new TCanvas(cname, ctitle, pixelx, pixely);
    myCanvas->DrawClonePad();
    cOut->Write();
    delete cOut;
    f->cd();
  };

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_chi2_pcmx_epp->Draw();
  //myCanvas->Print(fileName,"pdf");
  writeCanvas("c_chi2_pcmx_epp","chi2 vs sigma - pcmX integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_scale_pcmx_epp->Draw();
 // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_scale_pcmx_epp","scale vs sigma - pcmX integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmx_epp->Draw();
  {
    int thisbin=getMinBin(pcmx_cent_int);
    double scale = h_pcmx_epp_simSCM_scale[thisbin];
    h_pcmx_epp_simSCM[thisbin]->SetLineColor(2);
    h_pcmx_epp_simSCM[thisbin]->Scale(scale);
    h_pcmx_epp_simSCM[thisbin]->Draw("SAME");
  }
//  myCanvas->Print(fileName,"pdf");
  writeCanvas("c_overlay_pcmx_epp","data+sim overlay - pcmX integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_chi2_pcmy_epp->Draw();
// myCanvas->Print(fileName,"pdf");
  writeCanvas("c_chi2_pcmy_epp","chi2 vs sigma - pcmY integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_scale_pcmy_epp->Draw();
 // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_scale_pcmy_epp","scale vs sigma - pcmY integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmy_epp->Draw();
  {
    int thisbin=getMinBin(pcmy_cent_int);
    double scale = h_pcmy_epp_simSCM_scale[thisbin];
    h_pcmy_epp_simSCM[thisbin]->SetLineColor(2);
    h_pcmy_epp_simSCM[thisbin]->Scale(scale);
    h_pcmy_epp_simSCM[thisbin]->Draw("SAME");
  }
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_overlay_pcmy_epp","data+sim overlay - pcmY integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_chi2_pcmz_epp->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_chi2_pcmz_epp","chi2 vs sigma - pcmZ integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_scale_pcmz_epp->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_scale_pcmz_epp","scale vs sigma - pcmZ integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmz_epp->Draw();
  {
    int thisbin=getMinBin(pcmz_cent_int);
    double scale = h_pcmz_epp_simSCM_scale[thisbin];
    h_pcmz_epp_simSCM[thisbin]->SetLineColor(2);
    h_pcmz_epp_simSCM[thisbin]->Scale(scale);
    h_pcmz_epp_simSCM[thisbin]->Draw("SAME");
  }
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_overlay_pcmz_epp","data+sim overlay - pcmZ integrated");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_chi2_pcmT_epp->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_chi2_pcmT_epp","chi2 vs sigma - pcmT integrated");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_scale_pcmT_epp->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_scale_pcmT_epp","scale vs sigma - pcmT integrated");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmT_epp->Draw();
  {
    int thisbin=getMinBin(pcmT_cent_int);
    double scale = h_pcmT_epp_simSCM_scale[thisbin];
    h_pcmT_epp_simSCM[thisbin]->SetLineColor(2);
    h_pcmT_epp_simSCM[thisbin]->Scale(scale);
    h_pcmT_epp_simSCM[thisbin]->Draw("SAME");
  }
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_overlay_pcmT_epp","data+sim overlay - pcmT integrated");
  myCanvas->Clear();  

  //////////////////////////////////////
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigmacmx_Q2->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_sigmacmx_Q2","sigma_cmx vs Q2");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigmacmy_Q2->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_sigmacmy_Q2","sigma_cmy vs Q2");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigmacmz_Q2->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_sigmacmz_Q2","sigma_cmz vs Q2");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigmacmT_Q2->Draw();
  // myCanvas->Print(fileName,"pdf");
  writeCanvas("c_sigmacmT_Q2","sigma_cmT vs Q2");
  myCanvas->Clear();  

  for(int i=0; i<(bQ2); i++){    

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_chi2_pcmx_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_chi2_pcmx_Q2_%d",i),Form("chi2 pcmX Q2bin=%d",i));
    myCanvas->Clear();  
  
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_scale_pcmx_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_scale_pcmx_Q2_%d",i),Form("scale pcmX Q2bin=%d",i));
    myCanvas->Clear();  
    
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_pcmx_epp_SRC_Q2[i]->Draw();
    {
      int thisbin=getMinBin(pcmx_cent_Q2bin[i]);
      double scale = h_pcmx_epp_SRC_simSCM_Q2_scale[thisbin][i];
      h_pcmx_epp_SRC_simSCM_Q2[thisbin][i]->SetLineColor(2);
      h_pcmx_epp_SRC_simSCM_Q2[thisbin][i]->Scale(scale);
      h_pcmx_epp_SRC_simSCM_Q2[thisbin][i]->Draw("SAME");
    }
  //  myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_overlay_pcmx_Q2_%d",i),Form("data+sim pcmX Q2bin=%d",i));
    myCanvas->Clear();  

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_chi2_pcmy_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_chi2_pcmy_Q2_%d",i),Form("chi2 pcmY Q2bin=%d",i));
    myCanvas->Clear();  
  
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_scale_pcmy_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_scale_pcmy_Q2_%d",i),Form("scale pcmY Q2bin=%d",i));
    myCanvas->Clear();  
    
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_pcmy_epp_SRC_Q2[i]->Draw();
    {
      int thisbin=getMinBin(pcmy_cent_Q2bin[i]);
      double scale = h_pcmy_epp_SRC_simSCM_Q2_scale[thisbin][i];
      h_pcmy_epp_SRC_simSCM_Q2[thisbin][i]->SetLineColor(2);
      h_pcmy_epp_SRC_simSCM_Q2[thisbin][i]->Scale(scale);
      h_pcmy_epp_SRC_simSCM_Q2[thisbin][i]->Draw("SAME");
    }
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_overlay_pcmy_Q2_%d",i),Form("data+sim pcmY Q2bin=%d",i));
    myCanvas->Clear();  

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_chi2_pcmz_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_chi2_pcmz_Q2_%d",i),Form("chi2 pcmZ Q2bin=%d",i));
    myCanvas->Clear();  
  
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_scale_pcmz_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_scale_pcmz_Q2_%d",i),Form("scale pcmZ Q2bin=%d",i));
    myCanvas->Clear();  
    
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_pcmz_epp_SRC_Q2[i]->Draw();
    {
      int thisbin=getMinBin(pcmz_cent_Q2bin[i]);
      double scale = h_pcmz_epp_SRC_simSCM_Q2_scale[thisbin][i];
      h_pcmz_epp_SRC_simSCM_Q2[thisbin][i]->SetLineColor(2);
      h_pcmz_epp_SRC_simSCM_Q2[thisbin][i]->Scale(scale);
      h_pcmz_epp_SRC_simSCM_Q2[thisbin][i]->Draw("SAME");
    }
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_overlay_pcmz_Q2_%d",i),Form("data+sim pcmZ Q2bin=%d",i));
    myCanvas->Clear();  

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_chi2_pcmT_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_chi2_pcmT_Q2_%d",i),Form("chi2 pcmT Q2bin=%d",i));
    myCanvas->Clear();  
  
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    g_scale_pcmT_epp_SRC_Q2[i]->Draw();
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_scale_pcmT_Q2_%d",i),Form("scale pcmT Q2bin=%d",i));
    myCanvas->Clear();  
    
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_pcmT_epp_SRC_Q2[i]->Draw();
    {
      int thisbin=getMinBin(pcmT_cent_Q2bin[i]);
      double scale = h_pcmT_epp_SRC_simSCM_Q2_scale[thisbin][i];
      h_pcmT_epp_SRC_simSCM_Q2[thisbin][i]->SetLineColor(2);
      h_pcmT_epp_SRC_simSCM_Q2[thisbin][i]->Scale(scale);
      h_pcmT_epp_SRC_simSCM_Q2[thisbin][i]->Draw("SAME");
    }
    // myCanvas->Print(fileName,"pdf");
    writeCanvas(Form("c_overlay_pcmT_Q2_%d",i),Form("data+sim pcmT Q2bin=%d",i));
    myCanvas->Clear();  
  }
  
  //sprintf(fileName,"%s]",pdfFile);
  // myCanvas->Print(fileName,"pdf");

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
  if(Q2>5){return;}
  if(mmiss<0.65){return;}
  if(mmiss>1.10){return;}
  if(lead[0]->getRegion()!=FD){return;}

  if(kmiss<0.3){return;}            
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
  double int_scale = h_d->Integral()/h_s->Integral();
  double start_scale=0.1*int_scale;
  double   end_scale=3.0*int_scale;
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
   //   double chi2_point = 0.5 * sq(ob - ex) / sq(ob_err + ex_err);   
      double chi2_point = sq(ob - ex) / (sq(ob_err));// + sq(ex_err));

      chi2+=chi2_point;
    }
    if(chi2<min_chi2){
      final_scale=scale;
      min_chi2=chi2;
    }
  }
  //double NDF = h_d->FindBin(max)+1-h_d->FindBin(min);
  //cout<<NDF<<endl;
  //min_chi2=min_chi2/NDF;
  //return chi2/NDF instead of chi2
  
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
  //  myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_mu->Draw();
  //myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma->Draw();
 // myCanvas->Print(fileName,"pdf");
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
  //(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma->Draw();
 // myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  g_mu->Write();
  g_sigma->Write();
  
}



void getValue(TGraph * thisGraph,double min_sigma, double max_sigma, double & chi2NDF, double & center, double & lower, double & upper){
  double center_int_y = 10000000000;
  double center_int_x = 0;
  double lower_int_y = 10000000000;
  double lower_int_x = 0;
  double upper_int_y = 10000000000;
  double upper_int_x = 0;

  for(int k = 0; k < 200; k++){
    double x = ((double)k/(double)200)*(max_sigma-min_sigma) + min_sigma;
    if((thisGraph->Eval(x) < (center_int_y - .1))){
      center_int_x = x;
      center_int_y = thisGraph->Eval(x);
    }
  }
  for(int k = 0; k < 200; k++){
    double x = ((double)k/(double)200)*(center_int_x-min_sigma) + min_sigma;
    double new_y = fabs(thisGraph->Eval(x)-(center_int_y+1));
      if(new_y<lower_int_y){
	lower_int_x = x;
	lower_int_y = new_y;
      }
  }
  for(int k = 0; k < 200; k++){
    double x = ((double)k/(double)200)*(max_sigma-center_int_x) + center_int_x;
    double new_y = fabs(thisGraph->Eval(x)-(center_int_y+1));
    if(new_y<upper_int_y){
      upper_int_x = x;
      upper_int_y = new_y;
    }
  }

  center = center_int_x;
  lower = lower_int_x;
  upper = upper_int_x;
  chi2NDF = center_int_y;
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