#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"

class kinematics{

public:

	kinematics(){
		setEbeam(5.98636);
		Pp.SetXYZM(-999.,-999.,-999.,0);
	}
	~kinematics(){}

	void setParticle(TLorentzVector P);
	void calcInclusive();
	void calcPmiss(TLorentzVector Pq);
	
	void setEbeam(double Ebeam=5.98636);
	void setTarget(TString tgt);

	TLorentzVector getMomentum(){return Pp;};
	TLorentzVector getPbeam(){return Pbeam;};

	TLorentzVector getQ(){return q;};
	double getThetaq(){return theta_q;};
	double getPQ(){return pq;};
	double getQ2(){return Q2;};
	double getXB(){return xB;};
	double getW2(){return W2;};
	double getY(){return Yscale;};
	
	double getTheta(){return Pp.Theta();};
	double getPhi(){return Pp.Phi();};
	
	TVector3 getPmiss(){return Pmiss.Vect();};
	double getThetapq(){return theta_pq;};
	double getThetaPmissQ(){return theta_pmq;};
	double getMmiss(){return Mmiss;};
	double getMmiss2(){return Mmiss2;};
	double getEmiss0(){return Emiss0;};
	double getEmiss(){return Emiss;};
	double getAlpha(){return alpha_miss;};
	double getPmissQ_para(){return pmissq_para;};
	double getPmissQ_perp(){return pmissq_perp;};

private:

	void calcAlpha(TLorentzVector Pq);
	void calcY();
	
	TString tgt_name;
	double ebeam;
	TLorentzVector Pp;
	TLorentzVector Pbeam;

	TLorentzVector q;
        double theta_q;
        double Q2;
        double xB;
        double W2;
        double Yscale;
	
	TLorentzVector Pmiss;
        double theta_pq;
        double pq;
	double Mmiss;
	double Mmiss2;
	double Emiss0;
	double Emiss;
        double theta_pmq;
        double pmissq_para;
        double pmissq_perp;
        double alpha_miss;
        double alpha_N;
        double alpha_q;

	const double clight = 29.9792458; //cm/ns
	const double m_e = 0.000511; //GeV, eletron mass
	const double m_p = 0.938272; //GeV, proton mass
	double m_A;
	double m_Am1;
	double n_A;
	const double m_40Ca = 37.224918; //GeV, 40Ca mass
	const double m_48Ca = 44.667492; //GeV, 48Ca mass
	const double m_12C = 11.177928; //GeV, 12C mass
	const double m_40Ar = 37.224724; //GeV, 40Ar mass
	const double m_120Sn = 111.6881872; //GeV, 120Sn mass
	const double m_2D = 1.876124; //GeV, D mass
	const double m_He = 3.728401; //GeV, 4He mass
	const double A_40Ca = 40; 
	const double A_48Ca = 48; 
	const double A_12C = 12; 
	const double A_40Ar = 40; 
	const double A_120Sn = 120; 
	const double A_He = 4; 
	const double A_2D = 2; 
	const double A_1H = 1; 
	const double m_39K = 36.294459;
	const double m_47K = 43.744506;
	const double m_11B = 10.255102;
	const double m_39Cl = 36.298470;
	const double m_119In = 110.760093;
	const double m_3H = 2.809432;
	const double m_1H = 0.939565;
	const double m_0H = 0.0;
	
};

#endif

