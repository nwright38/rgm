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

#include "Functions.h"

using namespace std;
using namespace clas12;


void Usage()
{
  std::cerr << "Usage: ./code isMC(data=0,mc=1) A(4He=4) outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";
}

int main(int argc, char ** argv)
{
  
  if(argc < 5)
    {
      Usage();
      return -1;
    }

  //Start by defining the inputs
  int isMC = atoi(argv[1]);
  int nucleus_A = atoi(argv[2]);
  TString outFile = argv[3];
  char * pdfFile = argv[4];
  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  //This is the clas12ana class that helps us
  //cut on detector level and SRC variables.
  clas12ana clasAna;
  clasAna.printParams();
    
  //make hipochain for input
  clas12root::HipoChain chain;
  //Now add the input files to the chain
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  //Some necessary tags for the chain
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  //And this is the object we use to get event
  //details.
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();

  ////////////////////////////////////////////////
  
  //Define some useful variables
  auto db=TDatabasePDG::Instance();
  const double beam_E = 5.98636;
  const double mass_n = db->GetParticle(2112)->Mass();
  const double mass_p = db->GetParticle(2212)->Mass();
  const double mD = 1.8756;
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector el;
  TLorentzVector lead_ptr;

  //This is the reweighter which you can use to define
  //new GCF variables on the weight such as sigma_CM
  //and the potential.
  int Z=2;
  int N=2;
  if(isMC){
    Z=nucleus_A/2;
    N=nucleus_A/2;
  }
  reweighter newWeight(beam_E,Z,N,kelly,"AV18");
  //reweighter newWeight(beam_E,Z,N,kelly,"AV4");
  //reweighter newWeight(beam_E,Z,N,kelly,"N2LO10");
  //reweighter newWeight(beam_E,Z,N,kelly,"N2LO12");
  //reweighter newWeight(beam_E,Z,N,kelly,"NV");
  //////////////////////////////////

  //Define some histograms and a vector to put them
  //in so that we can manipulate them quickly.
  vector<TH1*> hist_list;
  TH1D * h_xB = new TH1D("xB","xB",100,1.0,2.0);
  hist_list.push_back(h_xB);
  TH1D * h_Q2 = new TH1D("Q2","Q2",200,1.0,5.0);
  hist_list.push_back(h_Q2);
  TH1D * h_plead = new TH1D("plead","plead",100,0.9,4.0);
  hist_list.push_back(h_plead);
  TH1D * h_thetalead = new TH1D("thetalead","thetalead",100,0,45);
  hist_list.push_back(h_thetalead);
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }
  
  //Now loop over all events in the input files.
  int counter = 0;
  while(chain.Next() && counter <100000000000000)
    {

      //Define the weight which is 1 for data.
      double weight = 1;
      if(isMC){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	weight = original_weight * newWeight.get_weight_ep(c12->mcparts());
      }

      //Display number of completed events
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    
      if((counter%10000) == 0){
	cerr << ".";
      }    

      //Here is where you run the clas12ana class.
      //It does the pid, fiducial, and vertex cuts
      //on all particles and returns particles
      //by PID number.
      //Right underneath we do the SRC cuts.
      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      clasAna.getLeadRecoilSRCwithCorrections(beam,isMC);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      auto lead_vectors = clasAna.getLeadSRC4Vectors();
      auto recoil_vectors = clasAna.getRecoilSRC4Vectors();

      if(electrons.size() != 1){continue;}
      el = clasAna.getElectronSRC4Vector();

      TLorentzVector q = beam - el;
      double Q2        = -q.M2();
      double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
      
      if(lead.size()!=1){continue;}
      lead_ptr = lead_vectors[0];
      double mom_lead = lead_ptr.P();
      double theta_lead = lead_ptr.Theta() * 180 / M_PI;

      h_xB->Fill(xB,weight);
      h_Q2->Fill(Q2,weight);
      h_plead->Fill(mom_lead,weight);
      h_thetalead->Fill(theta_lead,weight);
       
    }

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }

  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);
  
  myCanvas->Divide(2,2);
  myCanvas->cd(1);    
  h_xB->Draw();
  myCanvas->cd(2);    
  h_Q2->Draw();
  myCanvas->cd(3);    
  h_plead->Draw();
  myCanvas->cd(4);    
  h_thetalead->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
    
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();

  return 0;
}
