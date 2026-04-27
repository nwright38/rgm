#ifndef CORRECTIONS_HH
#define CORRECTIONS_HH
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

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
#include <TGraph.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "TRandom3.h"

inline TRandom3 * thisRand = new TRandom3(0);;

inline const double params_Theta_FD[6][4][4]={{{-0.191151,-0.669629,10,1.00023},
			    {0.113741,-0.0599251,19.0221,31.6623},
			    {350.715,460.172,17.1609,1},
			    {3.32022,2.7311,29.6054,25.1635}},
			   {{-0.0891744,-0.14685,7.50659,25.2638},
			    {0.0450352,-0.109521,6.0011,19.0888},
			    {350.357,140.077,6,15.8616},
			    {-3.25284,-1.24789,10.6854,24.956}},
			   {{-0.0849098,-0.0651043,11.7406,19.6928},
			    {0.0209555,-0.0273879,6.00308,20.1268},
			    {98.1241,20.8179,6,29.528},
			    {-5.06788,-4.9607,20.1145,1.00211}},
			   {{-0.0739777,0.0963433,6.00048,17.6548},
			    {0.014477,-0.156077,9.29576,12.1747},
			    {341.693,154.122,4.66104,22.2281},
			    {-6.38361,-1.95148,6,29.7998}},
			   {{-0.0528101,-0.100172,6.00982,12.4738},
			    {0.010489,-0.0339853,6.09838,16.0642},
			    {348.057,161.215,6,21.4269},
			    {-4.07125,-10,6,2.57348}},
			   {{-0.0623622,-0.179411,7.77701,9.66338},
			    {0.040617,-0.185013,6.00574,9.96181},
			    {199.345,74.073,6,28.6864},
			    {-3.83103,-1.41067,6,11.1274}}};
  
inline const double params_Theta_CD[3][2]={{-0.19764,8.14947},
			      {-0.0207092,4.52151},
			      {-1.45427,0}};

inline const double params_EnergyLoss_FD[3][3]={{-0.000695124,-0.000355869,0},
				   {0.00182181,7.77933e-05,0},
				   {0.000266541,0.424055,49.068}};

inline const double params_EnergyLoss_CD[3][2]={{-0.00555986,-6.06201e-05},
				   {0.00695634,8.24535e-05},
				   {0.00155103,1.74283e-05}};

inline const double params_Momentum_negparts_FD[6][4][4]={{{0.257158,-0.136574,4,25},
						    {1.79326,0.729487,4,25},
						    {121.024,40.091,5.00025,15},
						    {0.15236,-1.30429,1.00004,19.1085}},
						   {{0.781339,0.526124,4.35603,21.3139},
						    {1.22716,0.448927,4,24.2824},
						    {169.925,53.0067,4,17.9088},
						    {-0.139585,-1.28625,1.08963,19.6235}},
						   {{1.68681,0.96684,5.33823,23.4936},
						    {2.07862,1.13561,9.31451,24.8375},
						    {91.215,11.2227,4.00005,15},
						    {-3.31088,-0.431251,5.92834,15.0001}},
						   {{1.63967,0.683461,4,23.4088},
						    {0.833452,0.668007,4,20.6145},
						    {129.332,42.812,4.85292,18.6698},
						    {-3.10738,-0.605016,1,20.4228}},
						   {{1.34151,0.990672,10,22.6629},
						    {0.144001,-0.0533304,9.12876,15.0012},
						    {56.678,-7.6276,4,21.8985},
						    {-1.03723,-0.831154,3.09479,15}},
						   {{1.29883,0.888824,6.47301,21.2259},
						    {0.194675,-0.127275,4.00001,15},
						    {44.0007,-8.33546,4,15.1237},
						    {-0.439147,-0.599292,1,16.0641}}};

inline const double params_Momentum_posparts_FD[6][4][4]={{{-0.705061,1.01168,5,26.0005},
						    {3.48797,0.754098,1.92478,30.9666},
						    {114.46,8.69429,3,26},
						    {-1.7963,1.01494,5,26.5376}},
						   {{-0.0388399,0.831155,2.26375,27.1916},
						    {3.49646,1.03748,1.80573,27.729},
						    {161.625,21.6904,3,27.6331},
						    {-1.78998,0.611645,5,28.0376}},
						   {{-0.698553,1.3504,1.24212,26.4306},
						    {-0.442019,1.66664,4.99999,26.0007},
						    {200,-16.5106,4.97897,34},
						    {2.21317,-1.01685,1.93538,28.1037}},
						   {{0.0864055,1.33218,5,28.3018},
						    {1.80606,1.79593,4.23224,27.6494},
						    {167.639,19.6429,3,31.0384},
						    {-1.15901,0.465383,1.00009,27.7387}},
						   {{-1.00385,0.387042,2.4097,26.0041},
						    {-1.7593,-0.514224,1.00013,26.0001},
						    {146.135,20.5291,3,28.1081},
						    {1.28355,-1.12501,1.42483,26.0018}},
						   {{0.902805,1.11944,3,34},
						    {-2.32694,-0.909035,1,31.8645},
						    {108.115,22.0863,4.88316,26.0001},
						    {2.16539,0.371047,1.00015,26}}};

inline const double params_Momentum_CD[5][4]={{-0.250757,1.68519,14.261,74.5813},
				 {-0.95819,0.234407,6.15029,45.6566},
				 {-1.06401,2.89001,18.8068,45},
				 {0.669503,-1.43197,6.56704,45},
				 {0.285612,0.561206,5.00003,58.4323}};


inline const double params_Smear_FD[6][6]={{0.242269,0.012285,0.000722975,0.229986,0.00155551,0.00088862},
{-0.0331068,0.0452348,-4.87551e-05,0.557431,-0.0323144,0.00163881},
{-0.0781403,0.0453914,0.000252317,-0.465264,0.0735487,-0.000816354},
{0.0767921,0.0289579,0.000467256,0.267525,-0.00692348,0.00114717},
{-0.0501808,0.0404422,0.000173173,0.262579,-0.00198293,0.000950558},
  {0.167938,0.0297138,0.000191341,-0.143651,0.0362719,0.000199702}};

inline const double params_Smear_CD[2]={12.4637,6.76708};

/*
//Vertex Smearing

Electrons=
0.309387&1.85845
0.27904&2.80253
0.336598&1.9909
0.310611&2.25277
0.245308&2.59497
0.232689&2.9493

Protons=
0.763939&16.1967
0.629776&18.9645
0.588898&19.4205
0.500838&21.3636
0.729782&15.4736
0.643492&16.808

Protons CD 
0.82702&-0.016491&9.69232e-05
*/
inline double Function_Erf(double x, double A, double B, double C, double D){
  return A - B*(1+erf(((-x+D)/C))); 
}

inline double Function_Trig(double x, double A, double B, double C, double D){
  return A + B*sin((x*2*M_PI/C)+D); 
}

inline double Function_Cosine(double x, double A, double B, double C, double D){
  return A + (B*cos((M_PI/2)*(1/C)*(x-D)));
}

inline double Function_TrigFixedPeriodFixedOffset(double x, double A, double B){
  return A*sin((x*2*M_PI/180)+B); 
}

inline double Function_TrigFixedPeriod(double x, double A, double B, double C){
  return A + B*sin((x*2*M_PI/180)+C); 
}

inline double Function_TrigTripleFixedPeriod(double x, double A, double B, double C, double D, double E, double F, double G){
  return A + B*sin((x*1*M_PI/180)+C) + D*sin((x*2*M_PI/180)+E) + F*sin((x*3*M_PI/180)+G); 
}

inline double Function_TrigQuadFixedPeriod(double x, double A, double B, double C, double D, double E){
  return A +
    C*sin((x*1*M_PI/180)+ B) +
    D*sin((x*2*M_PI/180)+ B) +
    E*sin((x*3*M_PI/180)+ B);
}

inline double Function_Algebraic1(double x, double A, double B){
  return A + B*x; 
}

inline double Function_Algebraic2(double x, double A, double B){
  return A + B/(x); 
}

inline double Function_Algebraic3(double x, double A, double B, double C){
  return A + B*x + (C/x); 
}

inline double Function_Algebraic4(double x, double A, double B, double C){
  return A - B/(x-C); 
}

inline double Function_Algebraic5(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}

inline double SubtractQuadriture_FD(double x, double A, double B, double C, double D, double E, double F){
  double X = Function_Algebraic5(x,A,B,C)*Function_Algebraic5(x,A,B,C) - Function_Algebraic5(x,D,E,F)*Function_Algebraic5(x,D,E,F);
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

inline double SubtractQuadriture_CD(double x, double A, double B){
  double X = A*A - B*B;
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

inline void GetLorentzVector_ReconVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

inline void SetLorentzVector_ThetaCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==clas12::FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      params[i] = Function_Erf(theta,params_Theta_FD[sector-1][i][0],params_Theta_FD[sector-1][i][1],params_Theta_FD[sector-1][i][2],params_Theta_FD[sector-1][i][3]);
    }
    theta+=Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==clas12::CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic2(theta,params_Theta_CD[i][0],params_Theta_CD[i][1]);
    }
    theta+=Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);        
  }
  else{
    std::cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetTheta(theta*M_PI/180);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

inline void SetLorentzVector_EnergyLossCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->par()->getCharge()<=0){return;}
  if(p->getRegion()==clas12::FD){
    double params0 = Function_Algebraic1(theta,params_EnergyLoss_FD[0][0],params_EnergyLoss_FD[0][1]);
    double params1 = Function_Algebraic1(theta,params_EnergyLoss_FD[1][0],params_EnergyLoss_FD[1][1]);
    double params2 = Function_Algebraic4(theta,params_EnergyLoss_FD[2][0],params_EnergyLoss_FD[2][1],params_EnergyLoss_FD[2][2]);
    mom+=Function_Algebraic3(mom,params0,params1,params2);
  }
  else if(p->getRegion()==clas12::CD){
    //CD Momentum Corrections were so small that we have turned them off
    /*
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic1(theta,params_EnergyLoss_CD[i][0],params_EnergyLoss_CD[i][1]);
    }
    mom+=Function_Algebraic3(mom,params[0],params[1],params[2]);    
    */
    }
  else{
    std::cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
} 

inline void SetLorentzVector_MomentumCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  double correction = 0;
  if(p->getRegion()==clas12::FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      if(p->par()->getCharge()<0){
	params[i] = Function_Erf(theta,params_Momentum_negparts_FD[sector-1][i][0],params_Momentum_negparts_FD[sector-1][i][1],params_Momentum_negparts_FD[sector-1][i][2],params_Momentum_negparts_FD[sector-1][i][3]);
      }
      else{
	params[i] = Function_Erf(theta,params_Momentum_posparts_FD[sector-1][i][0],params_Momentum_posparts_FD[sector-1][i][1],params_Momentum_posparts_FD[sector-1][i][2],params_Momentum_posparts_FD[sector-1][i][3]);
      }
    }
    correction=Function_Trig(phi,params[0],params[1],params[2],params[3]);
    mom/=(1-(correction/100));
  }
  else if(p->getRegion()==clas12::CD){
    double params[5];
    for(int i = 0; i < 5; i++){
      if(i==0){
	params[i] = Function_Cosine(theta,params_Momentum_CD[i][0],params_Momentum_CD[i][1],params_Momentum_CD[i][2],params_Momentum_CD[i][3]);
      }
      else{
	params[i] = Function_Erf(theta,params_Momentum_CD[i][0],params_Momentum_CD[i][1],params_Momentum_CD[i][2],params_Momentum_CD[i][3]);
      }
    }
    correction=Function_TrigQuadFixedPeriod(phi,params[0],params[1],params[2],params[3],params[4]);
    mom/=(1-(correction/100));   
  }
  else{
    std::cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

inline void SetLorentzVector_MomentumSimulationSmear(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==clas12::FD){
    int sector = p->getSector();
    double smear = 0.01*SubtractQuadriture_FD(theta,params_Smear_FD[sector-1][0],params_Smear_FD[sector-1][1],params_Smear_FD[sector-1][2],params_Smear_FD[sector-1][3],params_Smear_FD[sector-1][4],params_Smear_FD[sector-1][5]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else if(p->getRegion()==clas12::CD){
    double smear = 0.01*SubtractQuadriture_CD(phi,params_Smear_CD[0],params_Smear_CD[1]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else{
    std::cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}


inline void GetLorentzVector_Corrected(TLorentzVector &p4, clas12::region_particle* p, bool isMC){
  GetLorentzVector_ReconVector(p4,p);
  if(!isMC){SetLorentzVector_ThetaCorrection(p4,p);}
  SetLorentzVector_EnergyLossCorrection(p4,p);  
  if(!isMC){SetLorentzVector_MomentumCorrection(p4,p);}
  if(isMC){SetLorentzVector_MomentumSimulationSmear(p4,p);}
}


#endif

