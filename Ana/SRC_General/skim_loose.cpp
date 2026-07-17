#include <cstdlib>
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

// ---------------------------------------------------------------------
// Constants (kept identical to the full analysis code so the skim cuts
// stay consistent with downstream analysis cuts)
// ---------------------------------------------------------------------
auto db = TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD     = 1.8756;
double beam_E = 5.98636;

void Usage()
{
  std::cerr << "Usage: ./skim_ep isMC outputfile.hipo inputfiles.hipo \n\n";
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

  // Output hipo file
  char * outName = argv[2];
  cout << "Output file " << outName << endl;

  // Writer chain
  clas12root::HipoChainWriter chain(outName);
  for(int k = 3; k < argc; k++){
    std::string fname = argv[k];

    bool ok = true;
    try {
      hipo::reader testReader;
      testReader.open(fname.c_str());   // void return; bad files throw or leave reader empty
      if(testReader.getEntries() <= 0){
        ok = false;
      }
    } catch (const std::exception &e) {
      cerr << "WARNING: exception opening " << fname << ": " << e.what() << endl;
      ok = false;
    } catch (...) {
      cerr << "WARNING: unknown error opening " << fname << endl;
      ok = false;
    }

    if(!ok){
      cerr << "SKIPPING bad/corrupt file: " << fname << endl;
      continue;
    }

    cout << "Input file " << fname << endl;
    chain.Add(fname.c_str());
  }


  // for(int k = 3; k < argc; k++){
  //   cout << "Input file " << argv[k] << endl;
  //   chain.Add(argv[k]);
  // }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  ////////////////////////////////////////////////////////////////////
  // Lorentz vectors needed for the e'p selection
  ////////////////////////////////////////////////////////////////////
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
	cerr << "\n" << counter << " completed";
      }
    //  if(isMC && acc_counter > 10000000) break;
  //     if((counter % 100000) == 0){
	// cerr << ".";
  //     }

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

      if(xB < 1.1){continue;}
      if(Q2 < 1.){continue;}

      chain.WriteEvent();
    }

  cerr << "\n\nDone.\n";
  cerr << "Total events processed : " << counter << "\n";
  cerr << "e'p events written      : " << n_ep_written << "\n";
  cerr << "e'pp events written     : " << n_epp_written << "\n";
  cerr << "Total written           : " << (n_ep_written + n_epp_written) << "\n";

  return 0;
}