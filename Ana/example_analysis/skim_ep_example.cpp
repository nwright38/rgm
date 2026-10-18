#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vector>

#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "HipoChainWriter.h"
#include "clas12ana.h"


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

using namespace std;
using namespace clas12;

auto db = TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD     = 1.8756;

void Usage()
{
  std::cerr << "Usage: ./skim_ep isMC Ebeam outputfile.hipo inputfiles.hipo \n\n";
}

int main(int argc, char ** argv)
{
  if(argc < 4)
    {
      Usage();
      return -1;
    }

  int isMC = atoi(argv[1]);

  // This is the clas12ana class that helps us cut on detector level
  // and SRC variables (PID, fiducial, vertex cuts).
  clas12ana clasAna;
  clasAna.printParams();

  double beam_E = atof(argv[2]); // GeV


  // Output hipo file
  char * outName = argv[3];
  cout << "Output file " << outName << endl;


  // Writer chain
  clas12root::HipoChainWriter chain(outName);
  int n_inputs_added = 0;
  for(int k = 4; k < argc; k++){
    std::string fname = argv[k];
    cout << "Input file " << fname << endl;
    chain.Add(fname.c_str());
    n_inputs_added++;
  }

  if(n_inputs_added == 0){
    cerr << "ERROR: no valid input files were added." << endl;
    return -1;
  }


  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());

  long counter      = 0;
  long n_ep_written  = 0;
  long n_epp_written = 0;

  int acc_counter = 0;
  while(chain.Next())
    {
      counter++;
      if((counter % 100000) == 0){
	      cout << "\n" << counter << " completed";
      } 

      // Run PID/fiducial/vertex cuts, get particles by PID
      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons   = clasAna.getByPid(2212);

      if(electrons.size() != 1){continue;}

      // --- Electron 4-vector + corrections (data only; MC uses smearing) ---
      GetLorentzVector_ReconVector(el,electrons[0]);
      if(!isMC){
	      SetLorentzVector_ThetaCorrection(el,electrons[0]);
	      SetLorentzVector_MomentumCorrection(el,electrons[0]);
      }
      if(isMC){SetLorentzVector_MomentumSimulationSmear(el,electrons[0]);}


      TLorentzVector q = beam - el;
      double Q2  = -q.M2();
      double xB  = Q2 / (2 * mass_p * (beam.E() - el.E()));

      if(xB < 1.){continue;}
      if(Q2 < 1.){continue;}

      // --- Build lead/recoil directly from reconstructed protons ---
      int lead_idx = -1;
      double max_lead_mom = -1.;
      double min_lead_mom = 0.7; //Skimming for protons with momentum > .7

      std::vector<TLorentzVector> proton_p4;
      proton_p4.reserve(protons.size());

      for(int ip = 0; ip < (int)protons.size(); ip++){

        TLorentzVector p4(0,0,0,mass_p);
        GetLorentzVector_ReconVector(p4,protons[ip]);
        if(!isMC){SetLorentzVector_ThetaCorrection(p4,protons[ip]);}
        SetLorentzVector_EnergyLossCorrection(p4,protons[ip]);
        if(!isMC){SetLorentzVector_MomentumCorrection(p4,protons[ip]);}
        if(isMC){SetLorentzVector_MomentumSimulationSmear(p4,protons[ip]);}
        proton_p4.push_back(p4);

        const double pMag = p4.P();
        if(pMag > min_lead_mom && pMag > max_lead_mom){
          max_lead_mom = pMag;
          lead_idx = ip;
        }
      }

      if(lead_idx < 0){continue;}
      lead_ptr = proton_p4[lead_idx];

      TLorentzVector miss = q + deut_ptr - lead_ptr;
      double mmiss = sqrt(miss.M2());
      double mom_lead = lead_ptr.P();
      double theta_lead = lead_ptr.Theta() * 180 / M_PI;

      TLorentzVector miss_LC = lead_ptr - q;
      TVector3 u_ZQ = q.Vect().Unit();
      double pmm_ZQ = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
      double pmp_ZQ = miss_LC.Vect().Perp(u_ZQ);
      double kmiss_ZQ = sqrt(mass_p*mass_p*((pmp_ZQ*pmp_ZQ + mass_p*mass_p)/(pmm_ZQ*(2*mass_p - pmm_ZQ))) - mass_p*mass_p);

      if(kmiss_ZQ < 0.3){continue;}
      if(mmiss < 0.65){continue;}
      if(mmiss > 1.1){continue;}

      bool is_epp = (protons.size() >= 2);
      if(is_epp){ n_epp_written++; }
      else      { n_ep_written++;  }

      chain.WriteEvent();
      acc_counter++;
    }

  cerr << "\n\nDone.\n";
  cerr << "Total events processed : " << counter << "\n";
  cerr << "e'p events written      : " << n_ep_written << "\n";
  cerr << "e'pp events written     : " << n_epp_written << "\n";
  cerr << "Total written           : " << (n_ep_written + n_epp_written) << "\n";

  return 0;
}