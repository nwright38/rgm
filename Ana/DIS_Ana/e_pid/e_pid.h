#ifndef E_PID_HH
#define E_PID_HH

#include "TF1.h"
#include "TCanvas.h"
#include "TVector3.h"

#include "../kinematics/kinematics.h"


class e_pid : public kinematics 
{

public:

	e_pid() : kinematics() {
		//setParamsRGB(5.98636); -> cuts are loaded in AnaClass
		intervalEpcal = 3.5; //sigma
		intervalMom = 3.5;
		color = 2;	
	}
	~e_pid(){}

	bool isElectronFull(clas12::region_part_ptr eHit);

	void setParamsRGB(double Ebeam);
	void fillParams();
	void setIntervalEpcal(double newInterval);
	void setIntervalMom(double newInterval);
	void setColor(int newColor);

	void drawEpcal(int sector, TCanvas * myCanvas);
	void drawMom(int sector, TCanvas * myCanvas);

private:
	
	bool SF_Epcal_Cut(int sector, double Epcal, double SF);
	bool SF_Mom_Cut(int sector, double p, double SF);
	bool SFpcal_SFecin_Cut(double p, double Epcal, double Eecin);
	bool DC_Track_Cut(int sector, clas12::region_part_ptr eHit);

	bool SF_Cut(int sector, double x, double SF, double params[6][6], double interval);
	bool DC_Cut(int sector, int region, TVector3 v, double params[3][12][2]);
	double meanSF(int sector, double x, double params[6][6]);
	double sigmaSF(int sector, double x, double params[6][6]);
	double maxSF(int sector, double x, double params[6][6], double interval);
	double minSF(int sector, double x, double params[6][6], double interval);
	double FF(double x, double sf1, double sf2, double sf3);
	double minDC(int sector, int region, TVector3 v, double params[3][12][2]);
	double maxDC(int sector, int region, TVector3 v, double params[3][12][2]);

	char paramFileNameEpcal[100];
	double paramsEpcal[6][6];
	double intervalEpcal;
	char paramFileNameMom[100];
	double paramsMom[6][6];
	double intervalMom;
	char paramFileNameDC1[100];
	char paramFileNameDC2[100];
	char paramFileNameDC3[100];
	double paramsDC[3][12][2];
	int color;
};

#endif

