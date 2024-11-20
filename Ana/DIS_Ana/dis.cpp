#include <cstdlib>
#include <iostream>
#include <chrono>

#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TString.h>

#include "clas12reader.h"
#include "clas12writer.h"
#include "HipoChain.h"

#include "clas12ana.h"

#include "parameters.h"
#include "e_pid/e_pid.h"
#include "e_pid/e_pid.cpp"
//#include "p_pid/p_pid.h"
//#include "p_pid/p_pid.cpp"
#include "kinematics/kinematics.h"
#include "kinematics/kinematics.cpp"


using namespace std;
using namespace clas12;

//#include "/group/clas12/packages/qadb/1.2.0/srcC/include/QADB.h"
//using namespace QA;

//#define HE 
//#define H 
//#define D 
//#define CA40_B //to be used: 80nA
//#define CA40 //150nA
#define CA48
//#define C12 //6GeV 
//#define SIM 
//#define SN


void Usage(){
  std::cerr << "Usage: ./testAna outputfile runnumber... \n\n\n";
}

int sectorCD(clas12::region_part_ptr pHit);


int main(int argc, char ** argv){

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	//gROOT->ProcessLine("#include <vector>");
	//gROOT->ProcessLine("#include <array>");


	/////////////////////////////////////
	clas12root::HipoChain chain;
	//chain.db()->turnOffQADB();
	
	if(argc < 2){
		Usage();
		return -1;
	}

	TString outFile = argv[1];
	cout<<"Ouput file "<< outFile <<endl;

	TString fileName = "";
	if(argc >= 2){
		for(int i = 2; i != argc; ++i){
			TString inFile(argv[i]);
			
			//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.00040-00044.hipo",inFile.Data(),inFile.Data()));
			//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.00045-00049.hipo",inFile.Data(),inFile.Data()));
			//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.00050-00054.hipo",inFile.Data(),inFile.Data()));
			//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",inFile.Data(),inFile.Data()));
			cout<<"Input file "<< inFile << "\n";
			fileName = inFile;
		}
	}
	
	//tree maker
	TFile fout(Form("%s",outFile.Data()), "RECREATE");
	TTree tout("clas", "clas");
	TVector3 chargeVec(0,0,0);
	#include "ca_tree_dis.hh"

	
	int tgtA;
	Int_t cal_period = 2; // ~ 40Ca, after run 15542
#if defined(CA40_B)
	tgtA = 402;
	tgt_no = 40;
	ntarget = "40Ca";
	cal_period = 1;
		chain.Add(Form("/lustre/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));
	//	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));
	//cout << Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/train/meanField/meanField_0%s.hipo",fileName.Data()) << endl;
		cout << "adding " << Form("/lustre/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()) << endl;
		
#elif defined(CA40)
	tgtA = 40;
	tgt_no = 40;
	ntarget = "40Ca";
	cal_period = 1;
	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/40Ca/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));

#elif defined(CA48)
        tgtA = 48;
	tgt_no = 48;
	ntarget = "48Ca";
	cal_period = 2;
	//Reco
	chain.Add(Form("/lustre/expphy/cache/clas12/rg-m/production/pass1/6gev/48Ca/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));
	//Train
	//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/48Ca/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));

	//cout << Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/48Ca/dst/train/meanField/meanField_0%s.hipo",fileName.Data()) << endl;
#elif defined(HE)
        tgtA = 4;
	tgt_no = 4;
	ntarget = "He";
	cal_period = 1;
	//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/He/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));
	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/He/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));
	//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/He/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));

#elif defined(H)
        tgtA = 1;
	tgt_no = 1;
	ntarget = "1H";
	cal_period = 1;
	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/H/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));

#elif defined(D)
        tgtA = 2;
	tgt_no = 2;
	ntarget = "2D";
	cal_period = 1;
	//	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/D/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));
	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));


#elif defined(C12)
	//!!!!!
	// 4 GeV setting !!
	// change pbeam in parameters !!
	//
        tgtA = 12;
	tgt_no = 12;
	ntarget = "12C";
	cal_period = 1;
	//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/4gev/C/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));
	//chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/Cx4/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));

	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/Cx4/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));
	
#elif defined(SN)
        //!!!!!                                                                                                                                                           
        // 4 GeV setting !!                                                                                                    
        //                                                                                                                                                                
        tgtA = 120;
        tgt_no = 120;
        ntarget = "120Sn";
        cal_period = 1;
        //chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/4gev/C/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));                           /lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/Snx4/dst/train/meanField/meanField
	//	        chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/Snx4/dst/train/meanField/meanField_0%s.hipo",fileName.Data()));
	chain.Add(Form("/lustre19/expphy/cache/clas12/rg-m/production/pass1/6gev/Snx4/dst/recon/0%s/rec_clas_0%s.evio.*.hipo",fileName.Data(),fileName.Data()));

#elif defined(SIM)
	tgtA = 40;
	tgt_no = 40;
	ntarget = "40Ca";
	cal_period = 2;
	/*	for(Int_t i=610; i<710; i++){
		chain.Add(Form("/lustre19/expphy/volatile/clas12/rg-m/mc/reconhipo/recon_qe_ca_598636gev_%d_torus-1.0.hipo",i));
	}*/	
		/*        tgtA = 4;
	tgt_no = 4;
	ntarget = "He";
	cal_period = 1;*/
	//for(Int_t i=890; i<940; i++){
	//	//chain.Add(Form("/lustre19/expphy/volatile/clas12/rg-m/mc/reconhipo/recon_qe_he_6gev_sigmacm_100_%d_torus-1.0.hipo",i));	
	//}	
	/*	for(Int_t i=875; i<1000; i++){ //950
		//chain.Add(Form("/lustre19/expphy/volatile/clas12/rg-m/mc/reconhipo_old/recon_qe_he_6gev_sigmacm_140_%d_torus-1.0.hipo",i));	
	chain.Add(Form("/lustre19/expphy/volatile/clas12/rg-m/mc/reconhipo/bkgmerged_recon_qe_he_6gev_sigmacm_140_prel_200_%d_torus-1.0.hipo",i));
	}
	*/
	chain.Add("/volatile/clas12/rg-m/nwright/mc/mchipo/bkgmerged_mc_dis_40ca_6gev__torus-1.0.hipo");
		
#endif
	

	
	/////////////////////////////////////
	//Set masses
	/*setPara(tgtA);
	
	if(ntarget == "40Ca") tgt_no = 40;
	else if(ntarget == "48Ca") tgt_no = 48;
	else if(ntarget == "12C") tgt_no = 12;
	else if(ntarget == "120Sn") tgt_no = 120;
	else if(ntarget == "He") tgt_no = 4;
	else if(ntarget == "D") tgt_no = 2;
	else tgt_no = -1;
	*/
	
	//Ana class
	
 	//Only needed for simulation anymore
	/*       clas12ana clasAna;
	clasAna.readInputParam("ana.par");
		TString cal_file = "40Ca";
	//if(cal_period == 1) cal_file = "LD2";
	//else cal_file = "40Ca";
	clasAna.readEcalSFPar(Form("paramsSF_%s_x2.dat",cal_file));
	clasAna.readEcalPPar(Form("paramsPI_%s_x2.dat",cal_file));
	*/
	

	clas12ana clasAna;	

//#if defined(SIM)
        TString cal_file = "40Ca";
        TString anapar_file = "ca40";
        if(tgtA == 4){
                cal_file = "LD2";
                anapar_file = "h";
        }else if(tgtA == 1){
                cal_file = "LD2";
                anapar_file = "ca_48";
        }else if(tgtA == 48){
	  //cal_file = "40Ca";
	   cal_file = "LD2";
                anapar_file = "ca_48";
	}else if(tgtA == 40){
                cal_file = "LD2";
                anapar_file = "ca_40";
	}else if(tgtA == 12){
                cal_file = "12C";
                anapar_file = "c_12";
	}else if(tgtA == 120){
	  cal_file = "120Sn";
	  anapar_file = "sn_120";
        }

	//clasAna.readInputParam(Form("/u/home/jkahlbow/software/clas12root/macros/justin_software/rgm/Ana/cutFiles/ana_%s.par",anapar_file));
	//    clasAna.readInputParam(Form("/u/home/jkahlbow/software/clas12root/macros/justin_software/rgm/Ana/cutFiles/ana_he4.par"));
	//clasAna.readEcalSFPar(Form("/u/home/jkahlbow/software/clas12root/macros/justin_software/rgm/Ana/cutFiles/paramsSF_%s_x2.dat",cal_file));
	//	clasAna.readEcalPPar(Form("/u/home/jkahlbow/software/clas12root/macros/justin_software/rgm/Ana/cutFiles/paramsPI_%s_x2.dat",cal_file));
//#endif
        clasAna.printParams();

	

	//e relevant
	bool useCuts = true;
	clasAna.setPidCuts(useCuts); //!!!!
	clasAna.setEcalPCuts(useCuts);
	clasAna.setEcalSFCuts(useCuts);
	clasAna.setEcalDiagCuts(useCuts);
	clasAna.setEcalEdgeCuts(useCuts);
	clasAna.setDCEdgeCuts(useCuts);
	clasAna.setVertexCuts(useCuts);

	//clasAna.printParams();
	
	
	//p relevant
	/*	clasAna.setPidCuts();
	//clasAna.setProtonPidCuts();
	clasAna.setProtonPidCuts(true);
	//clasAna.setCDEdgeCuts();
	//clasAna.setCDEdgeCuts();
	//clasAna.setGhostTrackCuts();
	*/
	//clasAna.setPidCuts();
	//clasAna.setProtonPidCuts(true);



	//Ana electrons
	e_pid ePID;
        ePID.setEbeam(pbeam);
        ePID.setTarget(ntarget);
	e_pid ePIDtrue;
        ePIDtrue.setEbeam(pbeam);
        ePIDtrue.setTarget(ntarget);
       		
	//std::vector<p_pid> pPID;
	//std::vector<p_pid> pPIDtrue;

		        
	//clas12reader *c12;
	//for(int f = 0; f < chain.GetNFiles(); f++){

	auto config_c12 = chain.GetC12Reader();
	  //	  c12 = new clas12reader(chain.GetFileName(f).Data());
	  auto& c12 = chain.C12ref(); // region particle
	  auto db=TDatabasePDG::Instance();


	/////////////////////////////////////
	//Charge
	Double_t totBeamCharge = 0;
	Double_t totCharge = 0;
	Double_t charge = 0;
	/*auto config_c12=chain.GetC12Reader();
	chain.SetReaderTags({1});
	*/
	clas12reader* c12ptr = nullptr;
	
	
	

	/////////////////////////////////////

	int count = 0.;
	bool loopEnd = false; 
	double prevrunCharge =0;
	double prevCharge = 0;
	double addedCharge = 0;
	
	cout << "before chain " << endl;
	auto scal=c12->scalerReader();

	//	QADB * qa = new QADB();
	int tot = 1;
	while(chain.Next()){
       

	  //  break;
	  //	  count++;
	  // qa->AccumulateCharge();
	  
	  if(count % 100000 == 0) cout << count << " of " << tot << endl;
	  // if (count == 1000000) break;
	  //	  if(loopEnd) break;
		//Charge scaler
	  	   if(c12ptr!=c12.get()){//check if file changed
			c12ptr=c12.get();
			//create the scaler reader
			scal=c12->scalerReader();
			cout<<"Charge "<< c12->getRunBeamCharge()<<endl;
			totBeamCharge+=c12->getRunBeamCharge();
			if(loopEnd) break;
			}
		   prevCharge = charge;
		charge = c12->event()->getBeamCharge();
		if(charge!=0 && charge>totCharge) totCharge = charge;
		if(charge != prevCharge) addedCharge+=prevCharge;

		auto scal=c12->scalerReader(); 
		
		//	cout << "beam charge: " << c12->getRunBeamCharge() << endl;
		//continue; 
		prevrunCharge = runCharge;
		runCharge = -999.;
		runCharge = c12->getRunBeamCharge();
		if(runCharge != prevrunCharge) totBeamCharge+=runCharge;


		//		continue;
		Pe.SetXYZM(0,0,0,0);
		xB = -999.;
		isGoodElectron = false;
		e_dE_pcal = -999.;
		e_dE_ecin = -999.;
		e_dE_ecout = -999.;
		SF = -999.;
		theta_e = -999.;
                phi_e = -999;
                sector_e = -1;
                detector_e = -1;
		e_sec_dc = -1;
		vtx_e.SetXYZ(-999.,-999.,-999.);
		for(int n = 0; n < 3; n++){
		  e_dc_x[n] = -999.;
		  e_dc_y[n] = -999.;
		  e_dc_z[n] = -999.;
		}
		e_edge_dc.SetXYZ(-999.,-999.,-999.);
		q.SetXYZM(0,0,0,0);
                theta_q = -999.;
		Q2 = -999.;
		xB = -999.;
		W2 = -999.;
		e_V = -999.;
		e_W = -999.;
		Yscale = -999.;
		/* //proton
		tcorr.clear();
		pbeta.clear();
		ptof.clear();
		ptofp.clear();
		ppath.clear();
		p_dE_cnd1.clear();
		p_dE_cnd2.clear();
		p_dE_cnd3.clear();
		p_dE_ctof.clear();
		p_tof_cnd1.clear();
		p_tof_cnd2.clear();
		p_tof_cnd3.clear();
		p_tof_ctof.clear();
		p_path_cnd1.clear();
		p_path_cnd2.clear();
		p_path_cnd3.clear();
		p_path_ctof.clear();
		theta_p.clear();
		phi_p.clear();
		vtxp.clear();
		sector_p_cd.clear();
		p_cd.clear();
		p_fd.clear();
		mult_p_cd = 0;
		mult_p_fd = 0;
                Pp.clear();				
                pp.clear();				
		theta_pq.clear();
		theta_pmq.clear();
		pq.clear();
		p_no = 0;
		p_no_cd = 0;
		p_no_fd = 0;
             	Pmiss.clear();
             	pmiss.clear();
		Mmiss.clear();
		Mmiss2.clear();
		Emiss0.clear();
		Emiss.clear();
		plead_no = 0;
		plead_no_cd = 0;
		plead_no_fd = 0;
		plead_index = -1;
		plead_index_reco = -1;
		double lead_reco_true = 0;	
		double lead_reco_true_min = 9999.;	
		
		pim_dE_cnd1.clear();
		pim_dE_cnd2.clear();
		pim_dE_cnd3.clear();
		pim_dE_ctof.clear();
		pim_p.clear();
		mult_pim = 0;
		*/
		
		//Sim
		weight = 1.;
				
		theta_e_true = -999.;
		phi_e_true = -999.;
		q_true.SetXYZM(0,0,0,0);
		theta_q_true = -999.;
		Q2_true = -999.;
		xB_true = -999.;
                /*
		Pp_true.clear();				
		theta_p_true.clear();
		phi_p_true.clear();
		theta_pq_true.clear();
		pq_true.clear();
             	Pmiss_true.clear();
		Mmiss_true.clear();
		Emiss_true.clear();
             	pp_true.clear();
             	pmiss_true.clear();
		plead_index_true = -1;

		double q_prev = 0;
		*/

		/////////////////////////////////////
		//	clasAna.setEcalSFCuts();
		//	clasAna.setEcalEdgeCuts();
		//	clasAna.setPidCuts();
		//clasAna.setVertexCuts();
		
		clasAna.Run(c12);
	       
		double start_time = c12->event()->getStartTime();
		//c12->event()->getRFTime();

#if defined (SIM)
		weight = c12->mcevent()->getWeight();

		//mult_p_true = 0;
		//pPIDtrue.clear();

		//cout << "-------------" << endl;
		//cout << c12->mcevent()->getNpart() << endl;

		for(int i=0; i<c12->mcevent()->getNpart(); i++){
			
			c12->mcparts()->setEntry(i);

			if(c12->mcparts()->getPid(i)==11){
			
				Double_t electron_true[4] = {c12->mcparts()->getPx(),c12->mcparts()->getPy(),c12->mcparts()->getPz(),c12->mcparts()->getMass()};
				SetLorentzVectorTrue(Pe_true,electron_true);
				ePIDtrue.setParticle(Pe_true);

				theta_e_true = c12->mcparts()->getTheta();
				phi_e_true = c12->mcparts()->getPhi();

				// electron kinematics
				ePIDtrue.calcInclusive();
				q_true = ePIDtrue.getQ();
				theta_q_true = ePIDtrue.getThetaq();
				Q2_true = ePIDtrue.getQ2();
				xB_true = ePIDtrue.getXB();
			}


			/*if(c12->mcparts()->getPid(i)==2212){

				mult_p_true++;
				
				TLorentzVector Pp_true_cand;
				Double_t pp_e_true = sqrt(pow(c12->mcparts()->getP(),2)+pow(c12->mcparts()->getMass(),2));
				Pp_true_cand.SetPxPyPzE(c12->mcparts()->getPx(), c12->mcparts()->getPy(), c12->mcparts()->getPz(),pp_e_true);

				std::vector<double> iPp_true;
				iPp_true.push_back(Pp_true_cand.Px()); iPp_true.push_back(Pp_true_cand.Py()); 
				iPp_true.push_back(Pp_true_cand.Pz()); iPp_true.push_back(Pp_true_cand.E()); 
				Pp_true.push_back(iPp_true);				
				pp_true.push_back(Pp_true_cand.P());				
				
				
				theta_p_true.push_back(Pp_true_cand.Theta());
				phi_p_true.push_back(Pp_true_cand.Phi());
				p_pid pPIDcand_true;
                                pPIDtrue.push_back(pPIDcand_true);
				pPIDtrue[mult_p_true-1].setTarget(ntarget);
				pPIDtrue[mult_p_true-1].setParticle(Pp_true_cand);
				pPIDtrue[mult_p_true-1].calcPmiss(ePIDtrue.getQ());

				theta_pq_true.push_back(pPIDtrue[mult_p_true-1].getThetapq());
                                pq_true.push_back(pPIDtrue[mult_p_true-1].getPQ());

				TVector3 ipmiss_true(0,0,0);
                                ipmiss_true = pPIDtrue[mult_p_true-1].getPmiss();
                                std::vector<Double_t> iPmiss_true;
                                iPmiss_true.push_back(ipmiss_true.X()); iPmiss_true.push_back(ipmiss_true.Y()); iPmiss_true.push_back(ipmiss_true.Z());
                                Mmiss_true.push_back(pPIDtrue[mult_p_true-1].getMmiss());
                                Emiss_true.push_back(pPIDtrue[mult_p_true-1].getEmiss());

                                Pmiss_true.push_back(iPmiss_true);
                                pmiss_true.push_back(sqrt(ipmiss_true.X()*ipmiss_true.X()+ipmiss_true.Y()*ipmiss_true.Y()+ipmiss_true.Z()*ipmiss_true.Z()));
			
				//sorting in simulation: 0=electron; 1=lead p; 2=recoil nucleon
				if(i==1) plead_index_true = pmiss_true.size()-1;
				//cout << plead_index_true << endl;

				//if(theta_pq_true<25 && pq_true<0.96 && pq_true>0.62){
				//	plead_true = true;
				//	plead_no_true++;
				//	Plead_true = Pp_true;
				//}
				//
			}*/
		}		

#endif


		//Good ID selection done
		auto electrons = clasAna.getByPid(11);
		//auto protons = clasAna.getByPid(2212);
		//auto pims = c12->getByID(-211);


		//cout << "nelectorn: " << electrons.size() << endl;
		/*if(electrons.size()!=1){
		  cout << "bad electrons" << endl;
		  continue;
		  }*/
		if(electrons.size()==1){
                	SetLorentzVector(Pe,electrons[0]);
			ePID.setParticle(Pe);
			ePID.calcInclusive();


			xB = ePID.getXB();
			
			//cout << "xb: " << xB << endl;

			//if(xB < minXB) continue;

			isGoodElectron = true;

			ft_x = electrons[0]->ft(FTCAL)->getX();
			ft_y = electrons[0]->ft(FTCAL)->getY();
			ft_r = electrons[0]->ft(FTCAL)->getRadius();
			nphe = electrons[0]->che(HTCC)->getNphe();
			e_htcc_x = electrons[0]->che(HTCC)->getX();
			e_htcc_y = electrons[0]->che(HTCC)->getY();
			e_htcc_z = electrons[0]->che(HTCC)->getZ();
			e_dE_pcal = electrons[0]->cal(PCAL)->getEnergy();
			e_dE_ecin = electrons[0]->cal(ECIN)->getEnergy();
			e_dE_ecout = electrons[0]->cal(ECOUT)->getEnergy();
			e_V = electrons[0]->cal(PCAL)->getLv();
			e_W = electrons[0]->cal(PCAL)->getLw();
			e_Chi2DoF = electrons[0]->trk(DC)->getChi2()/electrons[0]->trk(DC)->getNDF();
			SF = (e_dE_ecin+e_dE_pcal+e_dE_ecout)/electrons[0]->getP();

			e_x_pcal = electrons[0]->cal(PCAL)->getX();
			e_y_pcal = electrons[0]->cal(PCAL)->getY();

			e_x_ecin = electrons[0]->cal(ECIN)->getX();
			e_y_ecin = electrons[0]->cal(ECIN)->getY();

			e_x_ecout = electrons[0]->cal(ECOUT)->getX();
			e_y_ecout = electrons[0]->cal(ECOUT)->getY();

			theta_e = electrons[0]->getTheta();
			phi_e = electrons[0]->getPhi();
			sector_e = electrons[0]->getSector();
			detector_e = electrons[0]->getRegion();

			vtx_e.SetXYZ(electrons[0]->par()->getVx(),electrons[0]->par()->getVy(),electrons[0]->par()->getVz());

			auto traj_index_1 = electrons[0]->traj(DC,6)->getIndex();  //layer 1 
			auto traj_index_2 = electrons[0]->traj(DC,18)->getIndex(); //layer 2
			auto traj_index_3 = electrons[0]->traj(DC,36)->getIndex(); //layer 3

			auto traj_edge_1  = electrons[0]->traj(DC,6)->getFloat("edge",traj_index_1);
			auto traj_edge_2  = electrons[0]->traj(DC,18)->getFloat("edge",traj_index_2);
			auto traj_edge_3  = electrons[0]->traj(DC,36)->getFloat("edge",traj_index_3);

			auto traj_x_1  = electrons[0]->traj(DC,6)->getFloat("x",traj_index_1);
                        auto traj_x_2  = electrons[0]->traj(DC,18)->getFloat("x",traj_index_2);
			auto traj_x_3  = electrons[0]->traj(DC,36)->getFloat("x",traj_index_3);

			auto traj_y_1  = electrons[0]->traj(DC,6)->getFloat("y",traj_index_1);
                        auto traj_y_2  = electrons[0]->traj(DC,18)->getFloat("y",traj_index_2);
                        auto traj_y_3  = electrons[0]->traj(DC,36)->getFloat("y",traj_index_3);
			
			auto traj_z_1  = electrons[0]->traj(DC,6)->getFloat("z",traj_index_1);
                        auto traj_z_2  = electrons[0]->traj(DC,18)->getFloat("z",traj_index_2);
                        auto traj_z_3  = electrons[0]->traj(DC,36)->getFloat("z",traj_index_3);
			
			
			int traj_idx[3] = {traj_index_1,traj_index_2,traj_index_3};
			int layerNum[3] = {6,18,36};
			for(int n = 0; n < 3; n++){

			  e_dc_x[n] = electrons[0]->traj(DC,layerNum[n])->getFloat("x",traj_idx[n]);
			  e_dc_y[n] = electrons[0]->traj(DC,layerNum[n])->getFloat("y",traj_idx[n]);
                          e_dc_z[n] = electrons[0]->traj(DC,layerNum[n])->getFloat("z",traj_idx[n]);

			}
			
		
		

			e_edge_dc.SetXYZ(traj_edge_1,traj_edge_2,traj_edge_3);
			//		e_pos_dc.SetXYZ(electrons[0]->traj(DC,6 )->getFloat("edge",electrons[0]->traj(DC,6 )->getIndex()),electrons[0]->traj(DC,18)->getFloat("edge",electrons[0]->traj(DC,18)->getIndex()),electrons[0]->traj(DC,36)->getFloat("edge",electrons[0]->traj(DC,36)->getIndex()));

			q = ePID.getQ();
			theta_q = ePID.getThetaq();
			Q2 = ePID.getQ2();
			W2 = ePID.getW2();
			Yscale = ePID.getY();

			if(Q2 < 1.) continue;
			if(xB < .2 || xB > .7) continue;
			if(sqrt(W2) < 1.2) continue;
			
		}
		

		/*mult_pim = pims.size();
		if(pims.size()>0 && isGoodElectron){
			for(Int_t pi=0; pi<pims.size(); pi++){

				pim_dE_cnd1.push_back(pims[pi]->sci(CND1)->getEnergy());
				pim_dE_cnd2.push_back(pims[pi]->sci(CND2)->getEnergy());
				pim_dE_cnd3.push_back(pims[pi]->sci(CND3)->getEnergy());
				pim_dE_ctof.push_back(pims[pi]->sci(CTOF)->getEnergy());
				
				TLorentzVector Pim(0,0,0,db->GetParticle(-211)->Mass());
				SetLorentzVector(Pim,pims[pi]);
				pim_p.push_back(Pim.P());
				Pim.Delete();
			}
		}
		*/
		

		/*if(protons.size()>=1 && isGoodElectron){
			pPID.clear();
			TVector3 ipmiss;
                      	ipmiss.SetXYZ(0,0,0);

			for(Int_t p=0; p<protons.size(); p++){
				Double_t itcorr = -999.;
				Double_t ipbeta = -999.;
				Double_t iptof = -999.;
				Double_t iptofp = -999.;
				Double_t ippath = -999.;
				std::vector<Double_t> ivtx_p;
				Int_t isector_cd = -1;
				Bool_t ip_cd = false;
                                Bool_t ip_fd = false;			
				Double_t itheta_pq = -999.;
				Double_t itheta_pmq = -999.;
				Double_t ipq = -999.;

				p_pid pPIDcand;
                                pPID.push_back(pPIDcand);
                                pPID[p].setTarget(ntarget);

				pPID[p].calcTCorr(protons[p],electrons[0]);
				itcorr = pPID[p].getTCorr();
				
				ipbeta = protons[p]->par()->getBeta();
				//iptof = protons[p]->getPath()/(clight*ipbeta);
				iptof = protons[p]->getTime() - start_time;
				ippath = protons[p]->getPath();
				pbeta.push_back(ipbeta);
				ptof.push_back(iptof);
				ppath.push_back(ippath);
				
				//Tof from momentum
				double imom = protons[p]->par()->getP();
				iptofp = (ippath/clight) / (imom/sqrt(pow(imom,2) + pow(mp,2)));
				ptofp.push_back(iptofp);

				p_dE_cnd1.push_back(protons[p]->sci(CND1)->getEnergy());
				p_dE_cnd2.push_back(protons[p]->sci(CND2)->getEnergy());
				p_dE_cnd3.push_back(protons[p]->sci(CND3)->getEnergy());
				p_dE_ctof.push_back(protons[p]->sci(CTOF)->getEnergy());
				
				p_tof_cnd1.push_back(protons[p]->sci(CND1)->getTime());
				p_tof_cnd2.push_back(protons[p]->sci(CND2)->getTime());
				p_tof_cnd3.push_back(protons[p]->sci(CND3)->getTime());
				p_tof_ctof.push_back(protons[p]->sci(CTOF)->getTime());
				p_path_cnd1.push_back(protons[p]->sci(CND1)->getPath());
				p_path_cnd2.push_back(protons[p]->sci(CND2)->getPath());
				p_path_cnd3.push_back(protons[p]->sci(CND3)->getPath());
				p_path_ctof.push_back(protons[p]->sci(CTOF)->getPath());
	
                                ivtx_p.push_back(protons[p]->par()->getVx());
                                ivtx_p.push_back(protons[p]->par()->getVy());
                                ivtx_p.push_back(protons[p]->par()->getVz());

				theta_p.push_back(protons[p]->getTheta());
                                phi_p.push_back(protons[p]->getPhi());

				p_no++;
				if(protons[p]->getRegion() == CD){
					mult_p_cd++; ip_cd = true; p_no_cd++;
					//isector_cd = sectorCD(protons[p]);
					isector_cd = clasAna.getCDRegion(protons[p]);
				}else if(protons[p]->getRegion() == FD){ mult_p_fd++; ip_fd = true; p_no_fd++;}
	
				
				//p momenta
				TLorentzVector Pp_cand;
				//Pp_cand.SetPxPyPzE(0,0,0,db->GetParticle(2212)->Mass());
				Pp_cand.SetPxPyPzE(0,0,0,mp);
				SetLorentzVector(Pp_cand,protons[p]);
				std::vector<Double_t> iPp;
                                iPp.push_back(Pp_cand.Px()); iPp.push_back(Pp_cand.Py()); iPp.push_back(Pp_cand.Pz()); 
				iPp.push_back(Pp_cand.E());
                                Pp.push_back(iPp);				
                                pp.push_back(Pp_cand.P());				
				
				//pmiss
				pPID[p].setParticle(Pp_cand);
				pPID[p].calcPmiss(ePID.getQ());
				//std::cout << protons.size() << "   " << (ePID.getQ()).P() << "  " << Pp_cand.P() << std::endl;
				///double q_new = (ePID.getQ()).P();
				///if(q_prev == q_new) std::cout << protons.size() << "   " << (ePID.getQ()).P() << "  " << Pp_cand.P() << std::endl;
				///q_prev = q_new;

				itheta_pq = pPID[p].getThetapq();
				itheta_pmq = pPID[p].getThetaPmissQ();
				ipq = pPID[p].getPQ();

				ipmiss = pPID[p].getPmiss();
				std::vector<Double_t> iPmiss;
                                iPmiss.push_back(ipmiss.X()); iPmiss.push_back(ipmiss.Y()); iPmiss.push_back(ipmiss.Z());
                                Mmiss.push_back(pPID[p].getMmiss());
				Mmiss2.push_back(pPID[p].getMmiss2());
				Emiss.push_back(pPID[p].getEmiss());
				Emiss0.push_back(pPID[p].getEmiss0());
				TLorentzVector Pp_miss(ipmiss.X(),ipmiss.Y(),ipmiss.Z(),pPID[p].getEmiss());

				//plead
				//compare to MC
				if(plead_index_true>-1){ 
					lead_reco_true = abs(ipmiss.Mag() - pmiss_true.at(plead_index_true));	
					if(lead_reco_true < lead_reco_true_min){
						plead_index_reco = p;
						lead_reco_true_min = lead_reco_true;
					}
				}
				


				tcorr.push_back(itcorr);
				vtxp.push_back(ivtx_p);
				sector_p_cd.push_back(isector_cd);
				p_cd.push_back(ip_cd);
                                p_fd.push_back(ip_fd);
				theta_pq.push_back(itheta_pq);
				theta_pmq.push_back(itheta_pmq);
                                pq.push_back(ipq);
                                Pmiss.push_back(iPmiss);
                                pmiss.push_back(ipmiss.Mag());
				//cout << "L  " << plead_index_reco << endl;
			}
			*/	

			//data
			/*clasAna.getLeadRecoilSRC(xB,q);
			auto lead = clasAna.getLeadSRC();

			if(lead.size()>0){
				//iplead = true;
				for(int pl=0; pl<lead.size(); pl++){

					plead_index = clasAna.getLeadSRCInd().at(pl);

					plead_no++;
					if(lead[pl]->getRegion() == CD) plead_no_cd++;
					else if(lead[pl]->getRegion() == FD) plead_no_fd++;

				}
			}
		}*/


		
		//	cout << "xB " << xB << endl;
		//	cout << "Q2 " << Q2 << endl;
		if(isGoodElectron && xB>0.15 && Q2>1.2){
			tout.Fill();
			//cout << "good event" << endl;
			count++;
		}

		//	if(count > 10000) break;
	}

	///	if(count > 10000) break;
	




	/////////////////////////////////////
	
	fout.cd();
	chargeVec.SetX(totCharge);
        chargeVec.SetY(totBeamCharge);
	chargeVec.Write("chargeVec");
	fout.Write();
	//	fout.Close();


	cout << "Charge sum  " << addedCharge << endl;


	cout << "Beam Charge  " << totBeamCharge << endl;
	cout << "Charge       " << totCharge << endl;

	
	cout << "charge ?? " << chain.TotalBeamCharge() << endl;

	cout << "\ntotal accumulated charge analyzed: " << endl;
	//	cout << "run=" << "  charge=" <<
	// qa->GetAccumulatedCharge() << " nC" << endl;


	//gBenchmark->Stop("timer");
	//gBenchmark->Print("timer");

	//auto finish = std::chrono::high_resolution_clock::now();

}


int sectorCD(clas12::region_part_ptr pHit){
	Int_t isector_p_cd = -1;

	if(pHit->getPhi()*TMath::RadToDeg()>30 && pHit->getPhi()*TMath::RadToDeg()<=150 ) isector_p_cd = 1;
	else if(pHit->getPhi()*TMath::RadToDeg()>150 || pHit->getPhi()*TMath::RadToDeg()<=-90 ) isector_p_cd = 2;
	else if(pHit->getPhi()*TMath::RadToDeg()>-90 && pHit->getPhi()*TMath::RadToDeg()<=30 ) isector_p_cd = 3;

	return isector_p_cd;
}
