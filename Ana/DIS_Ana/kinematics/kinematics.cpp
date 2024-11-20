#include "kinematics.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <unistd.h>


using namespace std;
//using namespace clas12;


void kinematics::setEbeam(double Ebeam){

	ebeam = Ebeam;

	Pbeam.SetPxPyPzE(0,0,ebeam,ebeam);
}


void kinematics::setTarget(TString tgt){
	tgt_name = tgt;
	if(tgt_name == "40Ca"){ 
		m_A = m_40Ca;
		m_Am1 = m_39K;
		n_A = A_40Ca;
	}else if(tgt_name == "48Ca"){ 
		m_A = m_48Ca;
		m_Am1 = m_47K;
		n_A = A_48Ca;
	}else if(tgt_name == "40Ar"){ 
		m_A = m_40Ar;
		m_Am1 = m_39Cl;
		n_A = A_40Ar;
	}else if(tgt_name == "12C"){ 
		m_A = m_12C;
		m_Am1 = m_11B;
		n_A = A_12C;
	}else if(tgt_name == "120Sn"){ 
		m_A = m_120Sn;
		m_Am1 = m_119In;
		n_A = A_120Sn;
	}else if(tgt_name == "He"){ 
		m_A = m_He;
		m_Am1 = m_3H;
		n_A = A_He;
	}else if(tgt_name == "2D"){ 
		m_A = m_2D;
		m_Am1 = m_1H;
		n_A = A_2D;
	}else if(tgt_name == "1H"){ 
		m_A = m_1H;
		m_Am1 = m_0H;
		n_A = A_1H;
	}else if(tgt_name == "120Sn"){
	  m_A = m_120Sn;
	  m_Am1 = m_119In;
	  n_A = A_120Sn;
	    
	}
	
	else cerr << "NO TARGET DEFINED" << endl;
	//cout << "MA  " << m_A << endl;
};


void kinematics::setParticle(TLorentzVector P){

	Pp.SetXYZM(P.X(),P.Y(),P.Z(),P.M());
};


void kinematics::calcInclusive(){

	q.SetPxPyPzE(0,0,0,0);

	q = Pbeam - Pp;
	theta_q = Pbeam.Vect().Angle(q.Vect())*TMath::RadToDeg();
	Q2 =  (q.Px()*q.Px()+q.Py()*q.Py()+q.Pz()*q.Pz()) - q.E()*q.E();
	W2 = m_p*m_p + 2*m_p*q.E() - Q2;	

	xB = Q2 / (2*m_p*q.E());
	
	calcY();
	
}


void kinematics::calcPmiss(TLorentzVector Pq){

	theta_pq = TMath::ACos((Pp.Vect().Dot(Pq.Vect())) / (Pp.Vect().Mag()*Pq.Vect().Mag()))*TMath::RadToDeg();
	pq = Pp.Vect().Mag() / Pq.Vect().Mag();

	Pmiss = Pp - Pq;
	
	Mmiss2 = pow(Pq.E()+2*m_p-Pp.E(),2) - Pmiss.Vect().Mag2();
	Mmiss = sqrt( Mmiss2 );
	
	Emiss0 = m_p - m_A + sqrt( pow(Pq.E()+m_A-Pp.E(),2) - Pmiss.Vect().Mag2() );
	Emiss = Pq.E() - Pp.E() + m_p - (sqrt(pow(m_Am1,2) + Pmiss.Vect().Mag2()) - m_Am1);

	theta_pmq = Pmiss.Vect().Angle(Pq.Vect())*TMath::RadToDeg();	
	pmissq_perp = TMath::Sin(theta_pmq*TMath::DegToRad())*Pmiss.P();
	pmissq_para = TMath::Cos(theta_pmq*TMath::DegToRad())*Pmiss.P();
		
	calcAlpha(Pq);

}


void kinematics::calcAlpha(TLorentzVector Pq){

	double pnq = TMath::Cos(theta_pq*TMath::DegToRad())*Pp.P();
	alpha_N = (Pp.E()-pnq) / (m_A/n_A);
	alpha_q = (Pq.E()-Pq.P()) / (m_A/n_A);

	alpha_miss = alpha_N - alpha_q;

}


void kinematics::calcY(){

	double q_mag = sqrt(q.Px()*q.Px()+q.Py()*q.Py()+q.Pz()*q.Pz());
        double W = sqrt( pow(m_A+q.E(),2) - q_mag*q_mag );
        double lambda = 0.5*(m_Am1*m_Am1 - m_p*m_p + W*W);

	Yscale = ((m_A+q.E())* sqrt(lambda*lambda-m_Am1*m_Am1*W*W) - q_mag*lambda) / (W*W);

}
