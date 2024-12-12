#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <algorithm> 

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "clas12reader.h"
#include "clas12ana.h"
#include "HipoChain.h"
//#include "functions.h"
//#include "/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/NeutronVeto/include/veto_functions.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;

#include "neutronSkim.hh"

void printProgress(double percentage);

void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root>  <path/to/input.hipo> \n";
}

void init_arrays(TTree *nTree);
void reset_arrays();

int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  /////////////////////////////////////
  
  bool isMC = false;
  if(atoi(argv[1]) == 1){isMC=true;}

  double Ebeam = atof(argv[2]);
  
  TFile * outFile = new TFile(argv[3],"RECREATE");
  TTree *nTree = new TTree("nTree","Isotropic Neutrons");

  init_arrays(nTree);

 // eventcut myCut(Ebeam,argv[4]);
 //myCut.print_cuts();
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
 
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  auto& c12 = chain.C12ref(); // region particle  


  // CND::hits
  auto cnd_hits = config_c12->addBank("CND::hits");
  auto cnd_hit_id = config_c12->getBankOrder(cnd_hits,"id");
  auto cnd_hit_status = config_c12->getBankOrder(cnd_hits,"status");
  auto cnd_hit_trkID = config_c12->getBankOrder(cnd_hits,"trkID");
  auto cnd_hit_sector = config_c12->getBankOrder(cnd_hits,"sector");
  auto cnd_hit_layer = config_c12->getBankOrder(cnd_hits,"layer");
  auto cnd_hit_component = config_c12->getBankOrder(cnd_hits,"component");
  auto cnd_hit_energy = config_c12->getBankOrder(cnd_hits,"energy");
  auto cnd_hit_time = config_c12->getBankOrder(cnd_hits,"time");
  auto cnd_hit_energy_unc = config_c12->getBankOrder(cnd_hits,"energy_unc");
  auto cnd_hit_time_unc = config_c12->getBankOrder(cnd_hits,"time_unc");
  auto cnd_hit_x = config_c12->getBankOrder(cnd_hits,"x");
  auto cnd_hit_y = config_c12->getBankOrder(cnd_hits,"y");
  auto cnd_hit_z = config_c12->getBankOrder(cnd_hits,"z");
  auto cnd_hit_x_unc = config_c12->getBankOrder(cnd_hits,"x_unc");
  auto cnd_hit_y_unc = config_c12->getBankOrder(cnd_hits,"y_unc");
  auto cnd_hit_z_unc = config_c12->getBankOrder(cnd_hits,"z_unc");
  auto cnd_hit_tx = config_c12->getBankOrder(cnd_hits,"tx");
  auto cnd_hit_ty = config_c12->getBankOrder(cnd_hits,"ty");
  auto cnd_hit_tz = config_c12->getBankOrder(cnd_hits,"tz");
  auto cnd_hit_tlength = config_c12->getBankOrder(cnd_hits,"tlength");
  auto cnd_hit_pathlength = config_c12->getBankOrder(cnd_hits,"pathlength");
  auto cnd_hit_indexLadc = config_c12->getBankOrder(cnd_hits,"indexLadc");
  auto cnd_hit_indexRadc = config_c12->getBankOrder(cnd_hits,"indexRadc");
  auto cnd_hit_indexLtdc = config_c12->getBankOrder(cnd_hits,"indexLtdc");
  auto cnd_hit_indexRtdc = config_c12->getBankOrder(cnd_hits,"indexRtdc");

  // CND::clusters
  auto cnd_clusters = config_c12->addBank("CND::clusters");
  auto cnd_clusters_id = config_c12->getBankOrder(cnd_clusters,"id");
  auto cnd_clusters_sector = config_c12->getBankOrder(cnd_clusters,"sector");
  auto cnd_clusters_layer = config_c12->getBankOrder(cnd_clusters,"layer");
  auto cnd_clusters_component = config_c12->getBankOrder(cnd_clusters,"component");
  auto cnd_clusters_nhits = config_c12->getBankOrder(cnd_clusters,"nhits");
  auto cnd_clusters_energy = config_c12->getBankOrder(cnd_clusters,"energy");
  auto cnd_clusters_x = config_c12->getBankOrder(cnd_clusters,"x");
  auto cnd_clusters_y = config_c12->getBankOrder(cnd_clusters,"y");
  auto cnd_clusters_z = config_c12->getBankOrder(cnd_clusters,"z");
  auto cnd_clusters_time = config_c12->getBankOrder(cnd_clusters,"time");
  auto cnd_clusters_status = config_c12->getBankOrder(cnd_clusters,"status");

  // MC banks
  auto mc_p = config_c12->addBank("MC::Particle");
  auto mc_px = config_c12->getBankOrder(mc_p,"px");
  auto mc_py = config_c12->getBankOrder(mc_p,"py");
  auto mc_pz = config_c12->getBankOrder(mc_p,"pz");

  auto mc_rec_match = config_c12->addBank("MC::RecMatch");
  auto mc_rec_idx = config_c12->getBankOrder(mc_rec_match, "pindex");
  auto mc_gen_idx = config_c12->getBankOrder(mc_rec_match, "mcindex");


  TVector3 p_mc(0.,0.,0.); //  MC particle momentum
  TVector3 p_beam(0., 0., Ebeam);

  clas12ana clasAna;

  int counter = 0;
  double weight;
  bool any_cnd;
  //Define cut class
  while(chain.Next()==true){


    reset_arrays();
    clasAna.Run(c12); // set parameters in include/clas12ana.h - using defaults for electrons

    weight = 1;
    if(isMC) weight=c12->mcevent()->getWeight();

    // get particles by type
    auto electrons = clasAna.getByPid(11);
    auto allParticles = c12->getDetParticles();

    if(electrons.size() != 1) continue;
    

    // reading in CND hit information without particle information
    if (c12->getBank(cnd_hits)->getRows() == 0) continue;
    for (auto iRow=0; iRow < c12->getBank(cnd_hits)->getRows(); iRow++){

      cnd_hit_id_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_id, iRow);
      cnd_hit_status_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_status, iRow);
      cnd_hit_trkID_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_trkID, iRow);
      cnd_hit_sector_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_sector, iRow);
      cnd_hit_layer_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_layer, iRow);
      cnd_hit_component_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_component, iRow);

      cnd_hit_energy_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy, iRow);
      cnd_hit_time_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_time, iRow);
      cnd_hit_energy_unc_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy_unc, iRow);
      cnd_hit_time_unc_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_time_unc, iRow);
      cnd_hit_x_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_x, iRow);
      cnd_hit_y_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_y, iRow);
      cnd_hit_z_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_z, iRow);
      cnd_hit_x_unc_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_x_unc, iRow);
      cnd_hit_y_unc_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_y_unc, iRow);
      cnd_hit_z_unc_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_z_unc, iRow);
      cnd_hit_tx_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_tx, iRow);
      cnd_hit_ty_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_ty, iRow);
      cnd_hit_tz_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_tz, iRow);
      cnd_hit_tlength_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_tlength, iRow);
      cnd_hit_pathlength_arr[iRow] = c12->getBank(cnd_hits)->getFloat(cnd_hit_pathlength, iRow);

      cnd_hit_indexLadc_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_indexLadc, iRow);
      cnd_hit_indexRadc_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_indexRadc, iRow);
      cnd_hit_indexLtdc_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_indexLtdc, iRow);
      cnd_hit_indexRtdc_arr[iRow] = c12->getBank(cnd_hits)->getInt(cnd_hit_indexRtdc, iRow);


    }

    any_cnd = false;
    for(int particle = 0; particle < allParticles.size(); particle++){
      
      int mcparticle = -1;

      inCND1[particle] = (allParticles[particle]->sci(clas12::CND1)->getDetector() == 3);
      inCND2[particle] = (allParticles[particle]->sci(clas12::CND2)->getDetector() == 3);
      inCND3[particle] = (allParticles[particle]->sci(clas12::CND3)->getDetector() == 3);

      inCND[particle] = (inCND1 || inCND2 || inCND3);
      if(inCND[particle]) any_cnd = true;


      pid[particle] = allParticles[particle]->getPid();
      momentum[particle] = allParticles[particle]->getP();
      theta[particle] = allParticles[particle]->getTheta();
      phi[particle] = allParticles[particle]->getPhi();
      sector[particle] = allParticles[particle]->getSector();
      path_length[particle] = allParticles[particle]->getPath();

      if(inCND[particle]){

        de_dx[particle] = allParticles[particle]->sci(clas12::CND1)->getDedx() + allParticles[particle]->sci(clas12::CND2)->getDedx() + allParticles[particle]->sci(clas12::CND3)->getDedx();
        layermult[particle] = allParticles[particle]->sci(clas12::CND1)->getLayermulti() + allParticles[particle]->sci(clas12::CND2)->getLayermulti() + allParticles[particle]->sci(clas12::CND3)->getLayermulti();
      }

      // read Monte Carlo nucleon PID and momentum
      /*if(particle == c12->getBank(mc_rec_match)->getInt(mc_rec_idx)) // if there is an MC particle associated with this one
      {
        mcparticle = c12->getBank(mc_rec_match)->getInt(mc_gen_idx);
        p_mc.SetX(c12->getBank(mc_p)->getFloat(mc_px,mcparticle));
        p_mc.SetY(c12->getBank(mc_p)->getFloat(mc_py,mcparticle));
        p_mc.SetZ(c12->getBank(mc_p)->getFloat(mc_pz,mcparticle));
        momentum_mc[mcparticle] = p_mc.Mag();
        theta_mc[mcparticle] = p_mc.Theta();
        phi_mc[mcparticle] = p_mc.Phi();

        idx_recon[mcparticle] = particle;
      }*/


      


    }

    if(!any_cnd) continue;


/*
  
  /////////////////////////////////////
  //Neutral Particle Information
  /////////////////////////////////////
      for(int j = 0; j < allParticles.size(); j ++){
	if(allParticles[j]->par()->getCharge()!=0){continue;}
	//cout<<"PID="<<allParticles[j]->getPid()<<endl;
	bool PCAL = (allParticles[j]->cal(clas12::PCAL)->getDetector() == 7);
	bool ECin = (allParticles[j]->cal(clas12::ECIN)->getDetector() == 7);
  bool ECout = (allParticles[j]->cal(clas12::ECOUT)->getDetector() == 7);
      

	double path_n = neutrons[j]->getPath();
	double beta_frommom_n = p_n.Mag()/E_n;
	double time_frommom_n = path_n / (c*beta_frommom_n);
	double time_frombeta_n = path_n / (c*beta_n);


  // CND::clusters banks info
for (auto iRow=0; iRow < c12->getBank(cnd_clusters)->getRows(); iRow++){
  int nhits = c12->getBank(cnd_clusters)->getInt(cnd_clusters_nhits,iRow);
}

  // CND cluster size 
  double csize = neutrons[j]->sci(clas12::CND1)->getSize() + neutrons[j]->sci(clas12::CND2)->getSize() + neutrons[j]->sci(clas12::CND3)->getSize();
	

  */

    nTree->Fill();
    
  }



  outFile->cd();
  nTree->Write();

  outFile->Close();
}

void init_arrays(TTree *nTree){


  

  nTree->Branch("cnd_hit_id", &cnd_hit_id_arr, Form("cnd_hit_id[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_status", &cnd_hit_status_arr, Form("cnd_hit_status[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_trkID", &cnd_hit_trkID_arr, Form("cnd_hit_trkID[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_sector", &cnd_hit_sector_arr, Form("cnd_hit_sector[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_layer", &cnd_hit_layer_arr, Form("cnd_hit_layer[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_component", &cnd_hit_component_arr, Form("cnd_hit_component[%d]/I", n_allowed_hits));

  nTree->Branch("cnd_hit_energy", &cnd_hit_energy_arr, Form("cnd_hit_energy[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_time", &cnd_hit_time_arr, Form("cnd_hit_time[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_energy_unc", &cnd_hit_energy_unc_arr, Form("cnd_hit_energy_unc[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_time_unc", &cnd_hit_time_unc_arr, Form("cnd_hit_time_unc[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_x", &cnd_hit_x_arr, Form("cnd_hit_x[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_y", &cnd_hit_y_arr, Form("cnd_hit_y[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_z", &cnd_hit_z_arr, Form("cnd_hit_z[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_x_unc", &cnd_hit_x_unc_arr, Form("cnd_hit_x_unc[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_y_unc", &cnd_hit_y_unc_arr, Form("cnd_hit_y_unc[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_z_unc", &cnd_hit_z_unc_arr, Form("cnd_hit_z_unc[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_tx", &cnd_hit_tx_arr, Form("cnd_hit_tx[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_ty", &cnd_hit_ty_arr, Form("cnd_hit_ty[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_tz", &cnd_hit_tz_arr, Form("cnd_hit_tz[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_tlength", &cnd_hit_tlength_arr, Form("cnd_hit_tlength[%d]/F", n_allowed_hits));
  nTree->Branch("cnd_hit_pathlength", &cnd_hit_pathlength_arr, Form("cnd_hit_pathlength[%d]/F", n_allowed_hits));

  nTree->Branch("cnd_hit_indexLadc", &cnd_hit_indexLadc_arr, Form("cnd_hit_indexLadc[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_indexRadc", &cnd_hit_indexRadc_arr, Form("cnd_hit_indexRadc[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_indexLtdc", &cnd_hit_indexLtdc_arr, Form("cnd_hit_indexLtdc[%d]/I", n_allowed_hits));
  nTree->Branch("cnd_hit_indexRtdc", &cnd_hit_indexRtdc_arr, Form("cnd_hit_indexRtdc[%d]/I", n_allowed_hits));

  nTree->Branch("idx_recon", &idx_recon, Form("idx_recon[%d]/I", n_allowed_particles));

  nTree->Branch("pid", &pid, Form("pid[%d]/I", n_allowed_particles));
  nTree->Branch("momentum", &momentum, Form("momentum[%d]/F", n_allowed_particles));
  nTree->Branch("theta", &theta, Form("theta[%d]/F", n_allowed_particles));
  nTree->Branch("phi", &phi, Form("phi[%d]/F", n_allowed_particles));
  nTree->Branch("sector", &sector, Form("sector[%d]/I", n_allowed_particles));
  nTree->Branch("charge", &charge, Form("charge[%d]/I", n_allowed_particles));
  nTree->Branch("beta", &beta_arr, Form("beta[%d]/F", n_allowed_particles));
  nTree->Branch("path_length", &path_length, Form("path_length[%d]/F", n_allowed_particles));

  nTree->Branch("pid_mc", &pid_mc, Form("pid_mc[%d]/I", n_allowed_particles));
  nTree->Branch("momentum_mc", &momentum_mc, Form("momentum_mc[%d]/F", n_allowed_particles));
  nTree->Branch("theta_mc", &theta_mc, Form("theta_mc[%d]/F", n_allowed_particles));
  nTree->Branch("phi_mc", &phi_mc, Form("phi_mc[%d]/F", n_allowed_particles));

  nTree->Branch("de_dx", &de_dx, Form("de_dx[%d]/F", n_allowed_particles));
  nTree->Branch("layermult", &layermult, Form("layermult[%d]/I", n_allowed_particles));

  nTree->Branch("inCND", &inCND, Form("inCND[%d]/O", n_allowed_particles));
  nTree->Branch("inCND1", &inCND1, Form("inCND1[%d]/O", n_allowed_particles));
  nTree->Branch("inCND2", &inCND2, Form("inCND2[%d]/O", n_allowed_particles));
  nTree->Branch("inCND3", &inCND3, Form("inCND3[%d]/O", n_allowed_particles));



}
void reset_arrays(){


  std::fill_n(idx_recon,n_allowed_particles, 0);

    // PARTICLE VARIABLES
  std::fill_n(pid,n_allowed_particles, 0);
  std::fill_n(momentum,n_allowed_particles, 0.0f);
  std::fill_n(theta,n_allowed_particles, 0.0f);
  std::fill_n(phi,n_allowed_particles, 0.0f);
  std::fill_n(sector,n_allowed_particles, 0);
  std::fill_n(charge,n_allowed_particles, 0);
  std::fill_n(beta_arr,n_allowed_particles, 0.0f);
  std::fill_n(path_length,n_allowed_particles, 0.0f);

  std::fill_n(pid_mc,n_allowed_particles, 0);
  std::fill_n(momentum_mc,n_allowed_particles, 0.0f);
  std::fill_n(theta_mc,n_allowed_particles, 0.0f);
  std::fill_n(phi_mc,n_allowed_particles, 0.0f);

  std::fill_n(inCND,n_allowed_particles, false);
  std::fill_n(inCND1,n_allowed_particles, false);
  std::fill_n(inCND2,n_allowed_particles, false);

  std::fill_n(de_dx,n_allowed_particles, 0.0f);
  std::fill_n(layermult,n_allowed_particles, 0);

  // RAW HIT VARIABLES
  std::fill_n(cnd_hit_id_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_status_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_trkID_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_sector_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_layer_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_component_arr,n_allowed_hits, 0);

  std::fill_n(cnd_hit_energy_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_time_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_energy_unc_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_time_unc_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_x_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_y_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_z_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_x_unc_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_y_unc_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_z_unc_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_tx_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_ty_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_tz_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_tlength_arr,n_allowed_hits, 0.0f);
  std::fill_n(cnd_hit_pathlength_arr,n_allowed_hits, 0.0f);

  std::fill_n(cnd_hit_indexLadc_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_indexRadc_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_indexLtdc_arr,n_allowed_hits, 0);
  std::fill_n(cnd_hit_indexRtdc_arr,n_allowed_hits, 0);

  return;


}
