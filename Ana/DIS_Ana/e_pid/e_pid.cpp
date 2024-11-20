#include "e_pid.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace std;
using namespace clas12;


bool e_pid::isElectronFull(clas12::region_part_ptr eHit){

	//if any of the energy values is 0, cut away the electron
	double eop = (eHit->cal(PCAL)->getEnergy()+eHit->cal(ECIN)->getEnergy()+eHit->cal(ECOUT)->getEnergy()) / eHit->getP();
	if (eHit->cal(PCAL)->getEnergy() == 0 || eop == 0 || eHit->cal(ECIN)->getEnergy() == 0) {
		return false;
	}

	bool passSFEpcal = SF_Epcal_Cut(eHit->getSector()-1,eHit->cal(PCAL)->getEnergy(),eop);
	bool passSFMom = SF_Mom_Cut(eHit->getSector()-1,eHit->getP(),eop);
	bool passSFpi = SFpcal_SFecin_Cut(eHit->getP(),eHit->cal(PCAL)->getEnergy(),eHit->cal(ECIN)->getEnergy());
	bool passDC = DC_Track_Cut(eHit->getSector()-1, eHit);
	//For Debugging
	//std::cout << "ePID: passSFEpcal = " << passSFEpcal << " , passSFMom = " << passSFMom << " , passSFpi = " << passSFpi << std::endl;

	if(eHit->getPid() != 11){ return false; }
	if(eHit->par()->getCharge() != -1){ return false; }
	if(eHit->getRegion() != FD){ return false; }
	if(eHit->che(HTCC)->getNphe() <= 2){ return false; }	
	if(eHit->cal(PCAL)->getLv() < 14){ return false; }
	if(eHit->cal(PCAL)->getLw() < 14){ return false; }

	if(!passSFEpcal){ return false; }
	if(!passSFMom){ return false; }
	//!!! include once calibraiton is final !!!
	/////if(!passSFpi){ return false; }
	if(!passDC){ return false; } 

	return true;

}




void e_pid::setParamsRGB(double Ebeam){

	//Determine correct file name
	if(std::fabs(Ebeam-6.0) < 0.015){
		//new CLAS calibrations
		//40Ca, 12C, LD2
		//sprintf(paramFileNameEpcal,"cal_par/SFcut_test_new40Ca_bettersigma_optimizedSF1_r15392_19.dat");
		//sprintf(paramFileNameMom,"cal_par/PIcut_test_40Ca_new.dat");
		//sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/e_pid/cal_par/SFcut_new40Ca_final.dat");
                //sprintf(paramFileNameMom,  "/u/home/jkahlbow/software/clas12root/macros/e_pid/cal_par/PIcut_new40Ca_final.dat");
		//FINAL March 23 -- a+b/sqrt(x)+c/x
		//sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/SFcut_new40Ca_final.dat");
		//sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/PIcut_new40Ca_final.dat");
		//48Ca (and simulation ?)
		//sprintf(paramFileNameEpcal,"cal_par/SFcut_test_new40Ca_bettersigma_optimizedSF1_shiftto48Ca.dat");
		//sprintf(paramFileNameMom,"cal_par/PIcut_test_48Ca_new.dat");
		//FINAL March 23 -- a+b/sqrt(x)+c/x
		//sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/SFcut_new48Ca_final.dat");
		//sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/PIcut_new48Ca_final.dat");
		
		//FINAL USE March 23 -- a+b/x+c/x2
		//40Ca
		sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/SFcut_new40Ca_final_x2.dat");
		sprintf(paramFileNameMom,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/PIcut_new40Ca_final_x2.dat");
		//LD2
		//sprintf(paramFileNameEpcal,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/SFcut_newLD2_final_x2.dat");
		//sprintf(paramFileNameMom,"/u/home/jkahlbow/software/clas12root/macros/unpacker/e_pid/cal_par/PIcut_newLD2_final_x2.dat");

	}
	else{
		std::cout<<"Attempting to set a beam energy "<< Ebeam <<" GeV\n"
			<<"without a defined RGB fiducal cut\n"
			<<"\t Using E=10.6 GeV as Default. \n\n\n";
	}

	sprintf(paramFileNameDC1, "e_pid/cal_par/DC_linear_cuts_region1.dat");
	sprintf(paramFileNameDC2, "e_pid/cal_par/DC_linear_cuts_region2.dat");
	sprintf(paramFileNameDC3, "e_pid/cal_par/DC_linear_cuts_region3.dat");

	//Load file data
	fillParams();
}

void e_pid::fillParams(){

	//Load file data
	ifstream paramFileEpcal(paramFileNameEpcal);
	ifstream paramFileMom(paramFileNameMom);

	ifstream paramFileDC1(paramFileNameDC1);
	ifstream paramFileDC2(paramFileNameDC2);
	ifstream paramFileDC3(paramFileNameDC3);

	if(!paramFileEpcal.is_open() || !paramFileMom.is_open() || !paramFileDC1.is_open() || !paramFileDC2.is_open() || !paramFileDC3.is_open()){

		cout << paramFileNameEpcal << endl;
		cout << paramFileNameMom << endl;
		cout << paramFileNameDC1 << endl; 
		cout << paramFileNameDC2 << endl; 
		cout << paramFileNameDC3 << endl; 
		cout<<"Fiducial cut parameter files failed to load.\n"
			<<"Aborting...\n\n\n";
		exit(-2);
	}

	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			paramsEpcal[i][j] = 0;
			paramFileEpcal >> paramsEpcal[i][j];
			paramsMom[i][j] = 0;
			paramFileMom >> paramsMom[i][j];
			cout << paramsEpcal[i][j] << "  " << paramsMom[i][j] << endl;
		}
	}

	for(int i = 0; i < 12; i++){
		for(int j = 0; j < 2; j++){
			paramsDC[0][i][j] = 0;
			paramFileDC1 >> paramsDC[0][i][j];
			paramsDC[1][i][j] = 0;
			paramFileDC2 >> paramsDC[1][i][j];
			paramsDC[2][i][j] = 0;
			paramFileDC3 >> paramsDC[2][i][j];
		}
	}
	cout << "region1:" << endl;
	for(int i = 0; i < 12; i++){
		for(int j = 0; j < 2; j++){
			cout << paramsDC[0][i][j];
		}
		cout << endl;
	}
	cout << "region2:" << endl;
	for(int i = 0; i < 12; i++){
		for(int j = 0; j < 2; j++){
			cout << paramsDC[1][i][j];
		}
		cout << endl;
	}
	cout << "region3:" << endl;
	for(int i = 0; i < 12; i++){
		for(int j = 0; j < 2; j++){
			cout << paramsDC[2][i][j];
		}
		cout << endl;
	}


	paramFileEpcal.close();
	paramFileMom.close();
	paramFileDC1.close();
	paramFileDC2.close();
	paramFileDC3.close();

}

void e_pid::setIntervalEpcal(double newInterval){
	intervalEpcal = newInterval;
}

void e_pid::setIntervalMom(double newInterval){
	intervalMom = newInterval;
}

void e_pid::setColor(int newColor){
	color = newColor;
}

void e_pid::drawEpcal(int sector, TCanvas * myCanvas){

	sector--;

	TF1 * meanFunction = new TF1("Mean",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]); },0.06,1.6,3);
	meanFunction->SetLineColor(color);
	meanFunction->SetParameters(paramsEpcal[sector][0],paramsEpcal[sector][1],paramsEpcal[sector][2]);

	TF1 * maxFunction = new TF1("Max",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) + intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.06,1.6,6);
	maxFunction->SetLineColor(color);
	maxFunction->SetParameters(paramsEpcal[sector][0],paramsEpcal[sector][1],paramsEpcal[sector][2],paramsEpcal[sector][3],paramsEpcal[sector][4],paramsEpcal[sector][5]);

	TF1 * minFunction = new TF1("Min",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) - intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.06,1.6,6);
	minFunction->SetLineColor(color);
	minFunction->SetParameters(paramsEpcal[sector][0],paramsEpcal[sector][1],paramsEpcal[sector][2],paramsEpcal[sector][3],paramsEpcal[sector][4],paramsEpcal[sector][5]);


	//myCanvas->cd();
	meanFunction->Draw("SAME");
	maxFunction->Draw("SAME");
	minFunction->Draw("SAME");

}

void e_pid::drawMom(int sector, TCanvas * myCanvas){

	sector--;

	TF1 * meanFunction = new TF1("Mean",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]); },0.5,6.5,3);
	meanFunction->SetLineColor(color);
	meanFunction->SetParameters(paramsMom[sector][0],paramsMom[sector][1],paramsMom[sector][2]);

	TF1 * maxFunction = new TF1("Max",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) + intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.5,6.5,6);
	maxFunction->SetLineColor(color);
	maxFunction->SetParameters(paramsMom[sector][0],paramsMom[sector][1],paramsMom[sector][2],paramsMom[sector][3],paramsMom[sector][4],paramsMom[sector][5]);

	TF1 * minFunction = new TF1("Min",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) - intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.5,6.5,6);
	minFunction->SetLineColor(color);
	minFunction->SetParameters(paramsMom[sector][0],paramsMom[sector][1],paramsMom[sector][2],paramsMom[sector][3],paramsMom[sector][4],paramsMom[sector][5]);


	//myCanvas->cd();
	meanFunction->Draw("SAME");
	maxFunction->Draw("SAME");
	minFunction->Draw("SAME");

}

bool e_pid::SF_Epcal_Cut(int sector, double Epcal, double SF){
	return SF_Cut(sector,Epcal,SF,paramsEpcal,intervalEpcal);
}

bool e_pid::SF_Mom_Cut(int sector, double p, double SF){
	return SF_Cut(sector,p,SF,paramsMom,intervalMom);
}

bool e_pid::SFpcal_SFecin_Cut(double p, double Epcal, double Eecin){

	//This cut is only for high momentum particles
	if(p < 4.5){ return true; }

	double SFpi = (Epcal + Eecin) / p;
	if(SFpi < 0.2){ return false; }
	return true;

}

bool e_pid::SF_Cut(int sector, double x, double SF, double params[6][6], double interval){
	//cout << "in SF_cut: SF =  " << SF << " , minSF = " <<minSF(sector,x,params,interval) << " , maxSF " << maxSF(sector,x,params,interval) << endl;
	if(SF < minSF(sector,x,params,interval)){
		return false;
	}
	if(SF > maxSF(sector,x,params,interval)){
		return false;
	}
	return true;
}

bool e_pid::DC_Track_Cut(int sector, clas12::region_part_ptr eHit){
	bool DC_cuts[3];
	TVector3 v[3];
	double rot_angles[6] = {0, -TMath::Pi()/3, -2*TMath::Pi()/3,-3*TMath::Pi()/3, -4*TMath::Pi()/3, -5*TMath::Pi()/3};

	if((sector<0) || (sector>5)){return false;}

	auto track = eHit->trk(DC);

	v[0].SetXYZ(eHit->traj(DC,DC1)->getX(), eHit->traj(DC,DC1)->getY(), eHit->traj(DC,DC1)->getZ());
	v[1].SetXYZ(eHit->traj(DC,DC3)->getX(), eHit->traj(DC,DC3)->getY(), eHit->traj(DC,DC3)->getZ());
	v[2].SetXYZ(eHit->traj(DC,DC6)->getX(), eHit->traj(DC,DC6)->getY(), eHit->traj(DC,DC6)->getZ());

	for(int region = 0; region < 3; region++){
		v[region].RotateZ(rot_angles[sector]);
		if(!DC_Cut(sector, region, v[region], paramsDC)){
			return false;
		}
	}

	return true;

}

bool e_pid::DC_Cut(int sector, int region, TVector3 v, double params[3][12][2]){
	if(v.Y() < minDC(sector, region, v, params)){
		return false;
	}
	if(v.Y() > maxDC(sector, region, v, params)){
		return false;
	}
	return true;
}


double e_pid::meanSF(int sector, double x, double params[6][6]){
	return FF(x,params[sector][0],params[sector][1],params[sector][2]);
}

double e_pid::sigmaSF(int sector, double x, double params[6][6]){
	return FF(x,params[sector][3],params[sector][4],params[sector][5]);
}

double e_pid::maxSF(int sector, double x, double params[6][6], double interval){
	return meanSF(sector,x,params) + interval * sigmaSF(sector,x,params);
}

double e_pid::minSF(int sector, double x, double params[6][6], double interval){
	return meanSF(sector,x,params) - interval * sigmaSF(sector,x,params);
}

double e_pid::FF(double x, double sf1, double sf2, double sf3){
	//return sf1 + (sf2/sqrt(x)) + (sf3/x);
	return sf1 + (sf2/x) + (sf3/pow(x,2));
}

double e_pid::minDC(int sector, int region, TVector3 v, double params[3][12][2]){
	return params[region][2*sector+1][0] + params[region][2*sector+1][1]*v.X();
}

double e_pid::maxDC(int sector, int region, TVector3 v, double params[3][12][2]){
	return params[region][2*sector][0] + params[region][2*sector][1]*v.X();
}

