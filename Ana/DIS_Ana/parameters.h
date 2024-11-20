#ifndef PARAMETERS_H
#define PARAMETERS_H


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
	p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void SetLorentzVectorTrue(TLorentzVector &p4, double p[4]){
	p4.SetXYZM(p[0],p[1],p[2],p[3]);
}


Double_t minXB = 0.2;//0.55;

const Double_t pbeam = 5.98636; //GeV
//const Double_t pbeam = 4.030;//GeV
const Double_t clight = 29.9792458; //cm/ns
const Double_t mp = 0.938272; //GeV
const Int_t nsec = 7; //CD + 6xFD

TString ntarget = "40Ca"; //default target
Double_t vtxX_Ca = 0.0115043;//real 40Ca  
Double_t vtxY_Ca = 0.0563248;//real 40Ca
Double_t vtxZ_Ca = -2.34291;
Double_t vtxXs_Ca = 0.0526075;
Double_t vtxYs_Ca = 0.071465;
Double_t vtxZs_Ca = 6.58553e-01;

const Double_t vtxX_40CaB = 0.08227;//mislabeled 48Ca
const Double_t vtxY_40CaB = 0.122;//mislabeled 48Ca
const Double_t vtxZ_40CaB = -2.34291;
const Double_t vtxXs_40CaB = 0.0526075;
const Double_t vtxYs_40CaB = 0.071465;
const Double_t vtxZs_40CaB = 6.58553e-01;

const Double_t vtxX_48Ca = 0.127641;
const Double_t vtxY_48Ca = 0.0712873;
const Double_t vtxZ_48Ca = -2.73387;
const Double_t vtxXs_48Ca = 0.0526075;
const Double_t vtxYs_48Ca = 0.071465;
const Double_t vtxZs_48Ca = 6.19187e-01;

const Double_t vtxX_12C = 0.0904852;
const Double_t vtxY_12C = 0.11577;
const Double_t vtxZ_12C = -2.34291;
const Double_t vtxXs_12C = 0.0526075;
const Double_t vtxYs_12C = 0.071465;
const Double_t vtxZs_12C = 1.4;//1.12;

const Double_t vtxX_He = 2.77913e-02;
const Double_t vtxY_He = -3.35426e-02;
const Double_t vtxZ_He = -2.758;
const Double_t vtxXs_He = 2.421e-02;
const Double_t vtxYs_He = 5.242e-02;
const Double_t vtxZs_He = 0.9;

const Double_t vtxX_D = 2.77913e-02;
const Double_t vtxY_D = -3.35426e-02;
const Double_t vtxZ_D = -2.758;
const Double_t vtxXs_D = 2.421e-02;
const Double_t vtxYs_D = 5.242e-02;
const Double_t vtxZs_D = 0.9;

//not correct values, but also not used anymore
const Double_t vtxX_120Sn = 0.127641;
const Double_t vtxY_120Sn = 0.0712873;
const Double_t vtxZ_120Sn = -2.73387;
const Double_t vtxXs_120Sn = 0.0526075;
const Double_t vtxYs_120Sn = 0.071465;
const Double_t vtxZs_120Sn = 1.4;

const Double_t vtxX_Sim = 0.;
const Double_t vtxY_Sim = 0.;
const Double_t vtxZ_Sim = -2.;
const Double_t vtxXs_Sim = 0.257;
const Double_t vtxYs_Sim = 0.253;
const Double_t vtxZs_Sim = 0.427;


void setPara(int tgtA){

	if(tgtA==40){
		ntarget = "40Ca";

	}else if(tgtA==402){
		ntarget = "40Ca";

		vtxX_Ca = vtxX_40CaB;
		vtxY_Ca = vtxY_40CaB;
		vtxZ_Ca = vtxZ_40CaB;
		vtxXs_Ca = vtxXs_40CaB;
		vtxYs_Ca = vtxYs_40CaB;
		vtxZs_Ca = vtxZs_40CaB;

	}else if(tgtA==48){
		ntarget = "48Ca";

		vtxX_Ca = vtxX_48Ca;
		vtxY_Ca = vtxY_48Ca;
		vtxZ_Ca = vtxZ_48Ca;
		vtxXs_Ca = vtxXs_48Ca;
		vtxYs_Ca = vtxYs_48Ca;
		vtxZs_Ca = vtxZs_48Ca;

	}else if(tgtA==404){
		ntarget = "40Ar";
		//wrong! 
		vtxX_Ca = vtxX_48Ca;
		vtxY_Ca = vtxY_48Ca;
		vtxZ_Ca = vtxZ_48Ca;
		vtxXs_Ca = vtxXs_48Ca;
		vtxYs_Ca = vtxYs_48Ca;
		vtxZs_Ca = vtxZs_48Ca;

	}else if(tgtA==12){
		ntarget = "12C";

		vtxX_Ca = vtxX_12C;
		vtxY_Ca = vtxY_12C;
		vtxZ_Ca = vtxZ_12C;
		vtxXs_Ca = vtxXs_12C;
		vtxYs_Ca = vtxYs_12C;
		vtxZs_Ca = vtxZs_12C;

		//use e ID cuts for 48Ca 
	}else if(tgtA==403){
		ntarget = "40Ca";

		vtxX_Ca = vtxX_Sim;
		vtxY_Ca = vtxY_Sim;
		vtxZ_Ca = vtxZ_Sim;
		vtxXs_Ca = vtxXs_Sim;
		vtxYs_Ca = vtxYs_Sim;
		vtxZs_Ca = vtxZs_Sim;

	}else if(tgtA==4){
		ntarget = "He";

		vtxX_Ca = vtxX_He;
		vtxY_Ca = vtxY_He;
		vtxZ_Ca = vtxZ_He;
		vtxXs_Ca = vtxXs_He;
		vtxYs_Ca = vtxYs_He;
		vtxZs_Ca = vtxZs_He;
	

	}else if(tgtA==2){
		ntarget = "2D";

		vtxX_Ca = vtxX_D;
		vtxY_Ca = vtxY_D;
		vtxZ_Ca = vtxZ_D;
		vtxXs_Ca = vtxXs_D;
		vtxYs_Ca = vtxYs_D;
		vtxZs_Ca = vtxZs_D;
	
	}else if(tgtA==1){
		ntarget = "1H";

		//wrong.
		vtxX_Ca = vtxX_D;
		vtxY_Ca = vtxY_D;
		vtxZ_Ca = vtxZ_D;
		vtxXs_Ca = vtxXs_D;
		vtxYs_Ca = vtxYs_D;
		vtxZs_Ca = vtxZs_D;
	
	}else if(tgtA==120){
		ntarget = "120Sn";

		vtxX_Ca = vtxX_120Sn;
		vtxY_Ca = vtxY_120Sn;
		vtxZ_Ca = vtxZ_120Sn;
		vtxXs_Ca = vtxXs_120Sn;
		vtxYs_Ca = vtxYs_120Sn;
		vtxZs_Ca = vtxZs_120Sn;
	}
}

#endif
