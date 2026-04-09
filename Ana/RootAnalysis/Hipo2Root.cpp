#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <string>
#include <iostream>
#include <vector>
#include "HipoChain.h"
#include "clas12reader.h"
#include "clas12ana.h"
#include "Corrections.h"
#include <fstream>
#include <vector>


using namespace std;

double Pival=TMath::Pi();





int main(int argc, char **argv){
    clas12root::HipoChain chain;
    clas12ana clasAna;
    const int isMC = 0; // set to 1 if running on MC



   //chain.Add("/w/hallb-scshelf2102/clas12/igorpars/RGMrep/output/hipos/*.hipo");
   // Allow either:
   //   1) quoted wildcard pattern:  ./Hipo2Root "/path/*.hipo" out.root
   //   2) shell-expanded list:      ./Hipo2Root file1.hipo file2.hipo ... out.root
   if(argc < 3){
       std::cerr << "Usage:\n"
                 << "  " << argv[0] << " \"<hipo_pattern>\" <output_root_file>\n"
                 << "  " << argv[0] << " <hipo1> <hipo2> ... <output_root_file>\n";
       return 1;
   }

   // Last argument is ALWAYS the output ROOT file
   const char* outFile = argv[argc-1];

   // If only one input argument (argc==3), treat argv[1] as a (possibly quoted) pattern
   // Otherwise, treat argv[1..argc-2] as explicit input files
   if(argc == 3){
       chain.Add(argv[1]);
   } else {
       for(int i = 1; i < argc-1; i++){
           chain.Add(argv[i]);
       }
   }
 
    TVector3 p3_prot;
    TVector3 p3_gamma;
    TLorentzVector vgs;
    TLorentzVector p4_electron;
    TLorentzVector p4_prot;
    TLorentzVector p4_miss;
    
        
    TLorentzVector p4_beam;
    TLorentzVector p4_target;
    
    p4_beam.SetXYZM(0,0,5.986,0);
    //p4_beam.SetXYZM(0,0,2.07,0);
//    p4_target.SetXYZM(0,0,0,3.727);
    p4_target.SetXYZM(0,0,0,1.875); //D2

    
    
    
    //TFile *rFile=TFile::Open("ntuple.root","Recreate");
    TFile *rFile = TFile::Open(outFile,"RECREATE");
    if(!rFile || rFile->IsZombie()){
        std::cerr << "ERROR: could not create output ROOT file: " << outFile << std::endl;
        return 2;
    }
    TTree *T=new TTree("T","e'p and photons");
    
    
    // Define variables
    
    Float_t mp  =0.938;
    
    Float_t startTime;
    Long_t runNum;
    Long_t evNum;
    Float_t beamQ;
       
    // ==== electron =====
    Int_t nml;
    Int_t Epid[1000];
    Int_t Echarge[1000];
    Float_t Ep[1000];
    Float_t Epx[1000];
    Float_t Epy[1000];
    Float_t Epz[1000];
    Float_t Etheta[1000];
    Float_t Ephi[1000];
    Float_t Etime[1000];
    Float_t Epath[1000];
    Float_t Evz[1000];
    Int_t Estat[1000];
    Float_t EdE[1000];
    Float_t Ebeta[1000];
    
    Int_t Esector[1000];
    Int_t Eregion[1000];
    
    
    
    
    
    Float_t ePCAL_time[1000];
    Float_t ePCAL_energy[1000];
    Float_t ePCAL_path[1000];
    Int_t ePCAL_sector[1000];
    Float_t ePCAL_x[1000];
    Float_t ePCAL_y[1000];
    Float_t ePCAL_z[1000];
    Float_t ePCAL_lu[1000];
    Float_t ePCAL_lv[1000];
    Float_t ePCAL_lw[1000];
    Int_t ePCAL_status[1000];
    
    Float_t eECIN_time[1000];
    Float_t eECIN_energy[1000];
    Float_t eECIN_path[1000];
    Int_t eECIN_sector[1000];
    Float_t eECIN_x[1000];
    Float_t eECIN_y[1000];
    Float_t eECIN_z[1000];
    Float_t eECIN_lu[1000];
    Float_t eECIN_lv[1000];
    Float_t eECIN_lw[1000];
    Int_t eECIN_status[1000];
    
    
    Float_t eECOUT_time[1000];
    Float_t eECOUT_energy[1000];
    Float_t eECOUT_path[1000];
    Int_t eECOUT_sector[1000];
    Float_t eECOUT_x[1000];
    Float_t eECOUT_y[1000];
    Float_t eECOUT_z[1000];
    Float_t eECOUT_lu[1000];
    Float_t eECOUT_lv[1000];
    Float_t eECOUT_lw[1000];
    Int_t eECOUT_status[1000];
    
    Int_t eCHindex[1000];
    Int_t eCHdetector[1000];
    Int_t eCHnphe[1000];
    Float_t eCHtime[1000];
    Float_t eCHpath[1000];
    Int_t eCHsector[1000];
    Float_t eCHx[1000];
    Float_t eCHy[1000];
    Float_t eCHz[1000];
    Float_t eCHdtheta[1000];
    Float_t eCHdphi[1000];
    Int_t eCHstatus[1000];

    
    Float_t EdcX[1000][3]; //layer 6, 18,36
    Float_t EdcY[1000][3]; //layer 6 18 36
    Float_t EdcZ[1000][3]; //layer 6 18 36
    
    
    /*  particles[i]->cal(DET)->getTime();
     particles[i]->cal(DET)->getEnergy();
     particles[i]->cal(DET)->getPath();
     particles[i]->cal(DET)->getSector();
     particles[i]->cal(DET)->getX();
     particles[i]->cal(DET)->getY();
     particles[i]->cal(DET)->getZ();
     
     particles[i]->cal(DET)->getLu();
     particles[i]->cal(DET)->getLv();
     particles[i]->cal(DET)->getLw();
     
     particles[i]->cal(DET)->getStatus();
     */
    
    
    // =====  charged particles =====
    Int_t nP;
    Float_t prp[1000];
    Float_t prpx[1000];
    Float_t prpy[1000];
    Float_t prpz[1000];
    Float_t prtheta[1000];
    Float_t prphi[1000];

    // --- proton-only (selected by clas12ana) extra info (raw where appropriate) ---
    Int_t   prpid[1000];
    Int_t   prcharge[1000];
    Float_t prvz[1000];
    Int_t   prstat[1000];
    Float_t prtime[1000];
    Float_t prpath[1000];
    Float_t prdE[1000];
    Float_t prbeta[1000];


    // Float_t prdcX[1000][3];
    // Float_t prdcY[1000][3];
    // Float_t prdcZ[1000][3];
    
    Int_t nmb;
    Int_t ppid[1000];
    Int_t pcharge[1000];
    Float_t pp[1000];
    Float_t ppx[1000];
    Float_t ppy[1000];
    Float_t ppz[1000];
    Float_t ptheta[1000];
    Float_t pphi[1000];
    Float_t ptime[1000];
    Float_t ppath[1000];
    Float_t pvz[1000];
    Int_t pstat[1000];
    Float_t pdE[1000];
    Float_t pbeta[1000];
    
    Int_t psector[1000];
    Int_t pregion[1000];

    Float_t pFTOF1Atime[1000];
    Int_t pFTOF1ADet[1000];
    Int_t pFTOF1ALay[1000];
    Int_t pFTOF1AComp[1000];
    Float_t pFTOF1Apath[1000];
    Int_t pFTOF1Asector[1000];
    Int_t pFTOF1Astatus[1000];
    Float_t pFTOF1Ax[1000];
    Float_t pFTOF1Ay[1000];
    Float_t pFTOF1Az[1000];

    Float_t pFTOF2time[1000];
    Int_t pFTOF2Det[1000];
    Int_t pFTOF2Lay[1000];
    Int_t pFTOF2Comp[1000];
    Float_t pFTOF2path[1000];
    Int_t pFTOF2sector[1000];
    Int_t pFTOF2status[1000];
    Float_t pFTOF2x[1000];
    Float_t pFTOF2y[1000];
    Float_t pFTOF2z[1000];
    
    Float_t pFTOF1Btime[1000];
    Int_t pFTOF1BDet[1000];
    Int_t pFTOF1BLay[1000];
    Int_t pFTOF1BComp[1000];
    Float_t pFTOF1Bpath[1000];
    Int_t pFTOF1Bsector[1000];
    Int_t pFTOF1Bstatus[1000];
    Float_t pFTOF1Bx[1000];
    Float_t pFTOF1By[1000];
    Float_t pFTOF1Bz[1000];

    Float_t pCND1time[1000];
    Int_t pCND1Det[1000];
    Int_t pCND1Lay[1000];
    Int_t pCND1Comp[1000];
    Float_t pCND1path[1000];
    Int_t pCND1sector[1000];
    Int_t pCND1status[1000];
    Float_t pCND1Energy[1000];
    Float_t pCND1x[1000];
    Float_t pCND1y[1000];
    Float_t pCND1z[1000];
    

    Float_t pCND2time[1000];
    Int_t pCND2Det[1000];
    Int_t pCND2Lay[1000];
    Int_t pCND2Comp[1000];
    Float_t pCND2path[1000];
    Int_t pCND2sector[1000];
    Int_t pCND2status[1000];
    Float_t pCND2Energy[1000];
    Float_t pCND2x[1000];
    Float_t pCND2y[1000];
    Float_t pCND2z[1000];
   
    Float_t pCND3time[1000];
    Int_t pCND3Det[1000];
    Int_t pCND3Lay[1000];
    Int_t pCND3Comp[1000];
    Float_t pCND3path[1000];
    Int_t pCND3sector[1000];
    Int_t pCND3status[1000];
    Float_t pCND3Energy[1000];
    Float_t pCND3x[1000];
    Float_t pCND3y[1000];
    Float_t pCND3z[1000];


    Float_t pCTOFtime[1000];
    Int_t pCTOFDet[1000];
    Int_t pCTOFLay[1000];
    Int_t pCTOFComp[1000];
    Float_t pCTOFpath[1000];
    Int_t pCTOFsector[1000];
    Int_t pCTOFstatus[1000];
    Float_t pCTOFEnergy[1000];
    Float_t pCTOFx[1000];
    Float_t pCTOFy[1000];
    Float_t pCTOFz[1000];

    Int_t pCVTpind[1000];
    Int_t pCVTDet[1000];
    Int_t pCVTLay[1000];
    Float_t pCVTx[1000];
    Float_t pCVTy[1000];
    Float_t pCVTz[1000];
    Int_t pCVTtrkpind[1000];
    Int_t pCVTtrkDet[1000];
    Int_t pCVTtrkNDF[1000];
    Int_t pCVTtrkSec[1000];
    Float_t pCVTtrkChi2[1000];
 
    Float_t pPCAL_time[1000];
    Float_t pPCAL_energy[1000];
    Float_t pPCAL_path[1000];
    Int_t pPCAL_sector[1000];
    Float_t pPCAL_x[1000];
    Float_t pPCAL_y[1000];
    Float_t pPCAL_z[1000];
    Float_t pPCAL_lu[1000];
    Float_t pPCAL_lv[1000];
    Float_t pPCAL_lw[1000];
    Int_t pPCAL_status[1000];
    
    Float_t pECIN_time[1000];
    Float_t pECIN_energy[1000];
    Float_t pECIN_path[1000];
    Int_t pECIN_sector[1000];
    Float_t pECIN_x[1000];
    Float_t pECIN_y[1000];
    Float_t pECIN_z[1000];
    Float_t pECIN_lu[1000];
    Float_t pECIN_lv[1000];
    Float_t pECIN_lw[1000];
    Int_t pECIN_status[1000];
    
    
    Float_t pECOUT_time[1000];
    Float_t pECOUT_energy[1000];
    Float_t pECOUT_path[1000];
    Int_t pECOUT_sector[1000];
    Float_t pECOUT_x[1000];
    Float_t pECOUT_y[1000];
    Float_t pECOUT_z[1000];
    Float_t pECOUT_lu[1000];
    Float_t pECOUT_lv[1000];
    Float_t pECOUT_lw[1000];
    Int_t pECOUT_status[1000];
    
    
Int_t pECOUT_det[1000];
Int_t pECIN_det[1000];
Int_t pPCAL_det[1000];

Int_t ptrajPind[1000][3];
Int_t ptrajDet[1000][3];
Int_t ptrajLay[1000][3];
Float_t ptrajX[1000][3];
Float_t ptrajY[1000][3];
Float_t ptrajZ[1000][3];
Float_t ptrajPath[1000][3];   

 
    Float_t pdcX[1000][3]; //layer 6, 18,36
    Float_t pdcY[1000][3]; //layer 6 18 36
    Float_t pdcZ[1000][3]; //layer 6 18 36
    
    
    
    
    std::vector<int> hitInd;
    std::vector<int> hitpInd;
    std::vector<int> hitDet;
    std::vector<int> hitSec;
    std::vector<int> hitLay;
    std::vector<int> hitComp;
    std::vector<float> hitEnergy;
    std::vector<float> hitTime;
    std::vector<float> hitPath;
    std::vector<float> hitChi2;
    std::vector<float> hitX;
    std::vector<float> hitY;
    std::vector<float> hitZ;
    std::vector<float> hitHX;
    std::vector<float> hitHY;
    std::vector<float> hitHZ;
    std::vector<int> hitStat;
    std::vector<float> hitDeDx;
    std::vector<int> hitSize;
    std::vector<int> hitLayMult;
    
    
    Float_t pmisp[1000];
    Float_t pmispx[1000];
    Float_t pmispy[1000];
    Float_t pmispz[1000];
    Float_t pmisM2[1000];
    Float_t pmistheta[1000];
    Float_t pmisphi[1000];

    bool store = false;
    //=========== Tree branches ===============
    
    Int_t helicity;
    T->Branch("startTime",&startTime,"startTime/F");
    T->Branch("runNum",&runNum,"runNum/L");
    T->Branch("evNum",&evNum,"evNum/L");
    T->Branch("helicity",&helicity,"helicity/I");
    T->Branch("beamQ",&beamQ,"beamQ/F");
    
    
    T->Branch("nml",&nml,"nml/I");
    T->Branch("Epid",&Epid,"Epid[nml]/I");
    T->Branch("Echarge",&Echarge,"Echarge[nml]/I");
    T->Branch("Ep",&Ep,"Ep[nml]/F");
    T->Branch("Epx",&Epx,"Epx[nml]/F");
    T->Branch("Epy",&Epy,"Epy[nml]/F");
    T->Branch("Epz",&Epz,"Epz[nml]/F");
    T->Branch("Etheta",&Etheta,"Etheta[nml]/F");
    T->Branch("Ephi",&Ephi,"Ephi[nml]/F");
    T->Branch("Etime",&Etime,"Etime[nml]/F");
    T->Branch("Epath",&Epath,"Epath[nml]/F");
    T->Branch("Evz",&Evz,"Evz[nml]/F");
    T->Branch("Estat",&Estat,"Estat[nml]/I");
    T->Branch("EdE",&EdE,"EdE[nml]/F");
    T->Branch("Ebeta",&Ebeta,"Ebeta[nml]/F");
    
    T->Branch("Esector",&Esector,"Esector[nml]/I");
    T->Branch("Eregion",&Eregion,"Eregion[nml]/I");
    
    
    T->Branch("ePCAL_time",&ePCAL_time,"ePCAL_time[nml]/F");
    T->Branch("ePCAL_energy",&ePCAL_energy,"ePCAL_energy[nml]/F");
    T->Branch("ePCAL_path",&ePCAL_path,"ePCAL_path[nml]/F");
    T->Branch("ePCAL_sector",&ePCAL_sector,"ePCAL_sector[nml]/I");
    T->Branch("ePCAL_x",&ePCAL_x,"ePCAL_x[nml]/F");
    T->Branch("ePCAL_y",&ePCAL_y,"ePCAL_y[nml]/F");
    T->Branch("ePCAL_z",&ePCAL_z,"ePCAL_z[nml]/F");
    T->Branch("ePCAL_lu",&ePCAL_lu,"ePCAL_lu[nml]/F");
    T->Branch("ePCAL_lv",&ePCAL_lv,"ePCAL_lv[nml]/F");
    T->Branch("ePCAL_lw",&ePCAL_lw,"ePCAL_lw[nml]/F");
    T->Branch("ePCAL_status",&ePCAL_status,"ePCAL_status[nml]/I");
    
    T->Branch("eECIN_time",&eECIN_time,"eECIN_time[nml]/F");
    T->Branch("eECIN_energy",&eECIN_energy,"eECIN_energy[nml]/F");
    T->Branch("eECIN_path",&eECIN_path,"eECIN_path[nml]/F");
    T->Branch("eECIN_sector",&eECIN_sector,"eECIN_sector[nml]/I");
    T->Branch("eECIN_x",&eECIN_x,"eECIN_x[nml]/F");
    T->Branch("eECIN_y",&eECIN_y,"eECIN_y[nml]/F");
    T->Branch("eECIN_z",&eECIN_z,"eECIN_z[nml]/F");
    T->Branch("eECIN_lu",&eECIN_lu,"eECIN_lu[nml]/F");
    T->Branch("eECIN_lv",&eECIN_lv,"eECIN_lv[nml]/F");
    T->Branch("eECIN_lw",&eECIN_lw,"eECIN_lw[nml]/F");
    T->Branch("eECIN_status",&eECIN_status,"eECIN_status[nml]/I");
    
    T->Branch("eECOUT_time",&eECOUT_time,"eECOUT_time[nml]/F");
    T->Branch("eECOUT_energy",&eECOUT_energy,"eECOUT_energy[nml]/F");
    T->Branch("eECOUT_path",&eECOUT_path,"eECOUT_path[nml]/F");
    T->Branch("eECOUT_sector",&eECOUT_sector,"eECOUT_sector[nml]/I");
    T->Branch("eECOUT_x",&eECOUT_x,"eECOUT_x[nml]/F");
    T->Branch("eECOUT_y",&eECOUT_y,"eECOUT_y[nml]/F");
    T->Branch("eECOUT_z",&eECOUT_z,"eECOUT_z[nml]/F");
    T->Branch("eECOUT_lu",&eECOUT_lu,"eECOUT_lu[nml]/F");
    T->Branch("eECOUT_lv",&eECOUT_lv,"eECOUT_lv[nml]/F");
    T->Branch("eECOUT_lw",&eECOUT_lw,"eECOUT_lw[nml]/F");
    T->Branch("eECOUT_status",&eECOUT_status,"eECOUT_status[nml]/I");
    
    
    T->Branch("eCHindex",&eCHindex,"eCHindex[nml]/I");
    T->Branch("eCHdetector",&eCHdetector,"eCHdetector[nml]/I");
    T->Branch("eCHnphe",&eCHnphe,"eCHnphe[nml]/I");
    T->Branch("eCHtime",&eCHtime,"eCHtime[nml]/F");
    T->Branch("eCHpath",&eCHpath,"eCHpath[nml]/F");
    T->Branch("eCHsector",&eCHsector,"eCHsector[nml]/I");
    T->Branch("eCHx",&eCHx,"eCHx[nml]/F");
    T->Branch("eCHy",&eCHy,"eCHy[nml]/F");
    T->Branch("eCHz",&eCHz,"eCHz[nml]/F");
    T->Branch("eCHdtheta",&eCHdtheta,"eCHdtheta[nml]/F");
    T->Branch("eCHdphi",&eCHdphi,"eCHdphi[nml]/F");
    T->Branch("eCHstatus",&eCHstatus,"eCHstatus[nml]/I");

    
    T->Branch("EdcX",&EdcX,"EdcX[nml][3]/F"); //layer 6, 18,36
    T->Branch("EdcY",&EdcY,"EdcY[nml][3]/F"); //layer 6 18 36
    T->Branch("EdcZ",&EdcZ,"EdcZ[nml][3]/F"); //layer 6 18 36
    
    
    ////charged particles + branch ====================================
    T->Branch("nmb",&nmb,"nmb/I");
    T->Branch("ppid",&ppid,"ppid[nmb]/I");
    T->Branch("pcharge",&pcharge,"pcharge[nmb]/I");
    T->Branch("pp",&pp,"pp[nmb]/F");
    T->Branch("ppx",&ppx,"ppx[nmb]/F");
    T->Branch("ppy",&ppy,"ppy[nmb]/F");
    T->Branch("ppz",&ppz,"ppz[nmb]/F");
    T->Branch("ptheta",&ptheta,"ptheta[nmb]/F");
    T->Branch("pphi",&pphi,"pphi[nmb]/F");
    T->Branch("ptime",&ptime,"ptime[nmb]/F");
    T->Branch("ppath",&ppath,"ppath[nmb]/F");
    T->Branch("pvz",&pvz,"pvz[nmb]/F");
    T->Branch("pstat",&pstat,"pstat[nmb]/I");
    T->Branch("pdE",&pdE,"pdE[nmb]/F");
    T->Branch("pbeta",&pbeta,"pbeta[nmb]/F");
   
    T->Branch("pFTOF1Atime",&pFTOF1Atime,"pFTOF1Atime[nmb]/F");
    T->Branch("pFTOF1ADet",&pFTOF1ADet,"pFTOF1ADet[nmb]/I");
    T->Branch("pFTOF1ALay",&pFTOF1ALay,"pFTOF1ALay[nmb]/I");
    T->Branch("pFTOF1AComp",&pFTOF1AComp,"pFTOF1AComp[nmb]/I");
    T->Branch("pFTOF1Apath",&pFTOF1Apath,"pFTOF1Apath[nmb]/F");
    T->Branch("pFTOF1Asector",&pFTOF1Asector,"pFTOF1Asector[nmb]/I");
    T->Branch("pFTOF1Astatus",&pFTOF1Astatus,"pFTOF1Astatus[nmb]/I");
    T->Branch("pFTOF1Ax",&pFTOF1Ax,"pFTOF1Ax[nmb]/F");
    T->Branch("pFTOF1Ay",&pFTOF1Ay,"pFTOF1Ay[nmb]/F");
    T->Branch("pFTOF1Az",&pFTOF1Az,"pFTOF1Az[nmb]/F");


    T->Branch("pFTOF2time",&pFTOF2time,"pFTOF2time[nmb]/F");
    T->Branch("pFTOF2Det",&pFTOF2Det,"pFTOF2Det[nmb]/I");
    T->Branch("pFTOF2Lay",&pFTOF2Lay,"pFTOF2Lay[nmb]/I");
    T->Branch("pFTOF2Comp",&pFTOF2Comp,"pFTOF2Comp[nmb]/I");
    T->Branch("pFTOF2path",&pFTOF2path,"pFTOF2path[nmb]/F");
    T->Branch("pFTOF2sector",&pFTOF2sector,"pFTOF2sector[nmb]/I");
    T->Branch("pFTOF2status",&pFTOF2status,"pFTOF2status[nmb]/I");
    T->Branch("pFTOF2x",&pFTOF2x,"pFTOF2x[nmb]/F");
    T->Branch("pFTOF2y",&pFTOF2y,"pFTOF2y[nmb]/F");
    T->Branch("pFTOF2z",&pFTOF2z,"pFTOF2z[nmb]/F");


    T->Branch("pFTOF1Btime",&pFTOF1Btime,"pFTOF1Btime[nmb]/F");
    T->Branch("pFTOF1BDet",&pFTOF1BDet,"pFTOF1BDet[nmb]/I");
    T->Branch("pFTOF1BLay",&pFTOF1BLay,"pFTOF1BLay[nmb]/I");
    T->Branch("pFTOF1BComp",&pFTOF1BComp,"pFTOF1BComp[nmb]/I");
    T->Branch("pFTOF1Bpath",&pFTOF1Bpath,"pFTOF1Bpath[nmb]/F");
    T->Branch("pFTOF1Bsector",&pFTOF1Bsector,"pFTOF1Bsector[nmb]/I");
    T->Branch("pFTOF1Bstatus",&pFTOF1Bstatus,"pFTOF1Bstatus[nmb]/I");
    T->Branch("pFTOF1Bx",&pFTOF1Bx,"pFTOF1Bx[nmb]/F");
    T->Branch("pFTOF1By",&pFTOF1By,"pFTOF1By[nmb]/F");
    T->Branch("pFTOF1Bz",&pFTOF1Bz,"pFTOF1Bz[nmb]/F");

    T->Branch("pCND1time",&pCND1time,"pCND1time[nmb]/F");
    T->Branch("pCND1Det",&pCND1Det,"pCND1Det[nmb]/I");
    T->Branch("pCND1Lay",&pCND1Lay,"pCND1Lay[nmb]/I");
    T->Branch("pCND1Comp",&pCND1Comp,"pCND1Comp[nmb]/I");
    T->Branch("pCND1path",&pCND1path,"pCND1path[nmb]/F");
    T->Branch("pCND1sector",&pCND1sector,"pCND1sector[nmb]/I");
    T->Branch("pCND1status",&pCND1status,"pCND1status[nmb]/I");
    T->Branch("pCND1Energy",&pCND1Energy,"pCND1Energy[nmb]/F");
    T->Branch("pCND1x",&pCND1x,"pCND1x[nmb]/F");
    T->Branch("pCND1y",&pCND1y,"pCND1y[nmb]/F");
    T->Branch("pCND1z",&pCND1z,"pCND1z[nmb]/F");
 

    T->Branch("pCND2time",&pCND2time,"pCND2time[nmb]/F");
    T->Branch("pCND2Det",&pCND2Det,"pCND2Det[nmb]/I");
    T->Branch("pCND2Lay",&pCND2Lay,"pCND2Lay[nmb]/I");
    T->Branch("pCND2Comp",&pCND2Comp,"pCND2Comp[nmb]/I");
    T->Branch("pCND2path",&pCND2path,"pCND2path[nmb]/F");
    T->Branch("pCND2sector",&pCND2sector,"pCND2sector[nmb]/I");
    T->Branch("pCND2status",&pCND2status,"pCND2status[nmb]/I");
    T->Branch("pCND2Energy",&pCND2Energy,"pCND2Energy[nmb]/F");
    T->Branch("pCND2x",&pCND2x,"pCND2x[nmb]/F");
    T->Branch("pCND2y",&pCND2y,"pCND2y[nmb]/F");
    T->Branch("pCND2z",&pCND2z,"pCND2z[nmb]/F");


    T->Branch("pCND3time",&pCND3time,"pCND3time[nmb]/F");
    T->Branch("pCND3Det",&pCND3Det,"pCND3Det[nmb]/I");
    T->Branch("pCND3Lay",&pCND3Lay,"pCND3Lay[nmb]/I");
    T->Branch("pCND3Comp",&pCND3Comp,"pCND3Comp[nmb]/I");
    T->Branch("pCND3path",&pCND3path,"pCND3path[nmb]/F");
    T->Branch("pCND3sector",&pCND3sector,"pCND3sector[nmb]/I");
    T->Branch("pCND3status",&pCND3status,"pCND3status[nmb]/I");
    T->Branch("pCND3Energy",&pCND3Energy,"pCND3Energy[nmb]/F");
    T->Branch("pCND3x",&pCND3x,"pCND3x[nmb]/F");
    T->Branch("pCND3y",&pCND3y,"pCND3y[nmb]/F");
    T->Branch("pCND3z",&pCND3z,"pCND3z[nmb]/F");

    T->Branch("pCTOFtime",&pCTOFtime,"pCTOFtime[nmb]/F");
    T->Branch("pCTOFDet",&pCTOFDet,"pCTOFDet[nmb]/I");
    T->Branch("pCTOFLay",&pCTOFLay,"pCTOFLay[nmb]/I");
    T->Branch("pCTOFComp",&pCTOFComp,"pCTOFComp[nmb]/I");
    T->Branch("pCTOFpath",&pCTOFpath,"pCTOFpath[nmb]/F");
    T->Branch("pCTOFsector",&pCTOFsector,"pCTOFsector[nmb]/I");
    T->Branch("pCTOFstatus",&pCTOFstatus,"pCTOFstatus[nmb]/I");
    T->Branch("pCTOFEnergy",&pCTOFEnergy,"pCTOFEnergy[nmb]/F");
    T->Branch("pCTOFx",&pCTOFx,"pCTOFx[nmb]/F");
    T->Branch("pCTOFy",&pCTOFy,"pCTOFy[nmb]/F");
    T->Branch("pCTOFz",&pCTOFz,"pCTOFz[nmb]/F");

    T->Branch("pCVTpind",&pCVTpind,"pCVTpind[nmb]/I"); //
    T->Branch("pCVTDet",&pCVTDet,"pCVTDet[nmb]/I"); //
    T->Branch("pCVTLay",&pCVTLay,"pCVTLay[nmb]/I"); //
    T->Branch("pCVTx",&pCVTx,"pCVTx[nmb]/F"); //
    T->Branch("pCVTy",&pCVTy,"pCVTy[nmb]/F"); //
    T->Branch("pCVTz",&pCVTz,"pCVTz[nmb]/F"); //
    T->Branch("pCVTtrkpind",&pCVTtrkpind,"pCVTtrkpind[nmb]/I"); //
    T->Branch("pCVTtrkDet",&pCVTtrkDet,"pCVTtrkDet[nmb]/I"); //
    T->Branch("pCVTtrkNDF",&pCVTtrkNDF,"pCVTtrkNDF[nmb]/I"); //
    T->Branch("pCVTtrkSec",&pCVTtrkSec,"pCVTtrkSec[nmb]/I"); //
    T->Branch("pCVTtrkChi2",&pCVTtrkChi2,"pCVTtrkChi2[nmb]/F"); //

    T->Branch("psector",&psector,"psector[nmb]/I");
    T->Branch("pregion",&pregion,"pregion[nmb]/I");
    
    
    T->Branch("pPCAL_time",&pPCAL_time,"pPCAL_time[nmb]/F");
    T->Branch("pPCAL_energy",&pPCAL_energy,"pPCAL_energy[nmb]/F");
    T->Branch("pPCAL_path",&pPCAL_path,"pPCAL_path[nmb]/F");
    T->Branch("pPCAL_sector",&pPCAL_sector,"pPCAL_sector[nmb]/I");
    T->Branch("pPCAL_x",&pPCAL_x,"pPCAL_x[nmb]/F");
    T->Branch("pPCAL_y",&pPCAL_y,"pPCAL_y[nmb]/F");
    T->Branch("pPCAL_z",&pPCAL_z,"pPCAL_z[nmb]/F");
    T->Branch("pPCAL_lu",&pPCAL_lu,"pPCAL_lu[nmb]/F");
    T->Branch("pPCAL_lv",&pPCAL_lv,"pPCAL_lv[nmb]/F");
    T->Branch("pPCAL_lw",&pPCAL_lw,"pPCAL_lw[nmb]/F");
    T->Branch("pPCAL_status",&pPCAL_status,"pPCAL_status[nmb]/I");
    
    T->Branch("pECIN_time",&pECIN_time,"pECIN_time[nmb]/F");
    T->Branch("pECIN_energy",&pECIN_energy,"pECIN_energy[nmb]/F");
    T->Branch("pECIN_path",&pECIN_path,"pECIN_path[nmb]/F");
    T->Branch("pECIN_sector",&pECIN_sector,"pECIN_sector[nmb]/I");
    T->Branch("pECIN_x",&pECIN_x,"pECIN_x[nmb]/F");
    T->Branch("pECIN_y",&pECIN_y,"pECIN_y[nmb]/F");
    T->Branch("pECIN_z",&pECIN_z,"pECIN_z[nmb]/F");
    T->Branch("pECIN_lu",&pECIN_lu,"pECIN_lu[nmb]/F");
    T->Branch("pECIN_lv",&pECIN_lv,"pECIN_lv[nmb]/F");
    T->Branch("pECIN_lw",&pECIN_lw,"pECIN_lw[nmb]/F");
    T->Branch("pECIN_status",&pECIN_status,"pECIN_status[nmb]/I");
    
    T->Branch("pECOUT_time",&pECOUT_time,"pECOUT_time[nmb]/F");
    T->Branch("pECOUT_energy",&pECOUT_energy,"pECOUT_energy[nmb]/F");
    T->Branch("pECOUT_path",&pECOUT_path,"pECOUT_path[nmb]/F");
    T->Branch("pECOUT_sector",&pECOUT_sector,"pECOUT_sector[nmb]/I");
    T->Branch("pECOUT_x",&pECOUT_x,"pECOUT_x[nmb]/F");
    T->Branch("pECOUT_y",&pECOUT_y,"pECOUT_y[nmb]/F");
    T->Branch("pECOUT_z",&pECOUT_z,"pECOUT_z[nmb]/F");
    T->Branch("pECOUT_lu",&pECOUT_lu,"pECOUT_lu[nmb]/F");
    T->Branch("pECOUT_lv",&pECOUT_lv,"pECOUT_lv[nmb]/F");
    T->Branch("pECOUT_lw",&pECOUT_lw,"pECOUT_lw[nmb]/F");
    T->Branch("pECOUT_status",&pECOUT_status,"pECOUT_status[nmb]/I");
    
    T->Branch("pdcX",&pdcX,"pdcX[nmb][3]/F"); //layer 6, 18,36
    T->Branch("pdcY",&pdcY,"pdcY[nmb][3]/F"); //layer 6 18 36
    T->Branch("pdcZ",&pdcZ,"pdcZ[nmb][3]/F"); //layer 6 18 36

T->Branch("ptrajPind",&ptrajPind,"ptrajPind[nmb][3]/I");
T->Branch("ptrajDet",&ptrajDet,"ptrajDet[nmb][3]/I");
T->Branch("ptrajLay",&ptrajLay,"ptrajLay[nmb][3]/I");
T->Branch("ptrajX",&ptrajX,"ptrajX[nmb][3]/F");
T->Branch("ptrajY",&ptrajY,"ptrajY[nmb][3]/F");
T->Branch("ptrajZ",&ptrajZ,"ptrajZ[nmb][3]/F");
T->Branch("ptrajPath",&ptrajPath,"ptrajPath[nmb][3]/F");

T->Branch("pECOUT_det",&pECOUT_det,"pECOUT_det[nmb]/I");
T->Branch("pECIN_det",&pECIN_det,"pECIN_det[nmb]/I");
T->Branch("pPCAL_det",&pPCAL_det,"pPCAL_det[nmb]/I");
    
    
    
    T->Branch("hitInd",&hitInd);
    T->Branch("hitpInd",&hitpInd);
    T->Branch("hitDet",&hitDet);
    T->Branch("hitSec",&hitSec);
    T->Branch("hitLay",&hitLay);
    T->Branch("hitComp",&hitComp);
    T->Branch("hitEnergy",&hitEnergy);
    T->Branch("hitTime",&hitTime);
    T->Branch("hitPath",&hitPath);
    T->Branch("hitChi2",&hitChi2);
    T->Branch("hitX",&hitX);
    T->Branch("hitY",&hitY);
    T->Branch("hitZ",&hitZ);
    T->Branch("hitHX",&hitHX);
    T->Branch("hitHY",&hitHY);
    T->Branch("hitHZ",&hitHZ);
    T->Branch("hitStat",&hitStat);
    T->Branch("hitDeDx",&hitDeDx);
    T->Branch("hitSize",&hitSize);
    T->Branch("hitLayMult",&hitLayMult);
    
    
    T->Branch("nP",&nP,"nP/I");
    T->Branch("prp",&prp,"prp[nP]/F");
    T->Branch("prpx",&prpx,"prpx[nP]/F");
    T->Branch("prpy",&prpy,"prpy[nP]/F");
    T->Branch("prpz",&prpz,"prpz[nP]/F");
    T->Branch("prtheta",&prtheta,"prtheta[nP]/F");
    T->Branch("prphi",&prphi,"prphi[nP]/F");
    T->Branch("prpid",&prpid,"prpid[nP]/I");
    T->Branch("prcharge",&prcharge,"prcharge[nP]/I");
    T->Branch("prvz",&prvz,"prvz[nP]/F");
    T->Branch("prstat",&prstat,"prstat[nP]/I");
    T->Branch("prtime",&prtime,"prtime[nP]/F");
    T->Branch("prpath",&prpath,"prpath[nP]/F");
    T->Branch("prdE",&prdE,"prdE[nP]/F");
    T->Branch("prbeta",&prbeta,"prbeta[nP]/F");

    // T->Branch("prsector",&prsector,"prsector[nP]/I");
    // T->Branch("prregion",&prregion,"prregion[nP]/I");

    // T->Branch("prdcX",&prdcX,"prdcX[nP][3]/F");
    // T->Branch("prdcY",&prdcY,"prdcY[nP][3]/F");
    // T->Branch("prdcZ",&prdcZ,"prdcZ[nP][3]/F");
    
    T->Branch("pmisp",&pmisp,"pmisp[nP]/F");
    T->Branch("pmispx",&pmispx,"pmispx[nP]/F");
    T->Branch("pmispy",&pmispy,"pmispy[nP]/F");
    T->Branch("pmispz",&pmispz,"pmispz[nP]/F");
    T->Branch("pmisM2",&pmisM2,"pmisM2[nP]/F");
    T->Branch("pmistheta",&pmistheta,"pmistheta[nP]/F");
    T->Branch("pmisphi",&pmisphi,"pmisphi[nP]/F");
    
    long index =0;
    int tempInt;

    Float_t Q2;
    Float_t Nu;
    Float_t xB;
    Float_t q;
    Float_t qx;
    Float_t qy;
    Float_t qz;
    Float_t qtheta;
    Float_t qphi;


    T->Branch("Q2",&Q2,"Q2/F");
    T->Branch("Nu",&Nu,"Nu/F");
    T->Branch("xB",&xB,"xB/F");
    T->Branch("q",&q,"q/F");
    T->Branch("qx",&qx,"qx/F");
    T->Branch("qy",&qy,"qy/F");
    T->Branch("qz",&qz,"qz/F");
    T->Branch("qtheta",&qtheta,"qtheta/F");
    T->Branch("qphi",&qphi,"qphi/F");

    
    //==============================================================
    //loop over files
    //================================================================
    for(int ifile=0; ifile<chain.GetNFiles();++ifile){
        auto c12 = std::make_unique<clas12::clas12reader>(chain.GetFileName(ifile).Data());
        c12->addExactPid(11,1);	//exactly 1 electron
	// c12->addExactPid(2212,1);  // proton
	// c12->addExactPid(211,1);   // pi+
	// c12->addExactPid(-211,1);  // pi-
//        c12->addExactPid(2212,1);	//exactly 1 proton
//        c12->addExactPid(211,1);      //exactly 1 gammas
        
        nml  = 0;
        nmb  = 0;
        
        auto idx_RECScintExtra = c12->addBank("REC::ScintExtras");
        auto iRECScintExtraDEDX =   c12->getBankOrder(idx_RECScintExtra,"dedx");
        auto iRECScintExtraSize =   c12->getBankOrder(idx_RECScintExtra,"size");
        auto iRECScintExtraLMult=   c12->getBankOrder(idx_RECScintExtra,"layermult");

        
        auto idx_RECScint = c12->addBank("REC::Scintillator");
        auto iRECScintInd = c12->getBankOrder(idx_RECScint,"index");
        auto iRECScintPind = c12->getBankOrder(idx_RECScint,"pindex");
        auto iRECScintDet   =   c12->getBankOrder(idx_RECScint,"detector");
        auto iRECScintSec   =   c12->getBankOrder(idx_RECScint,"sector");
        auto iRECScintLay   =   c12->getBankOrder(idx_RECScint,"layer");
        auto iRECScintCom   =   c12->getBankOrder(idx_RECScint,"component");
        auto iRECScintEn    =   c12->getBankOrder(idx_RECScint,"energy");
        auto iRECScintTime  =   c12->getBankOrder(idx_RECScint,"time");
        auto iRECScintPath  =   c12->getBankOrder(idx_RECScint,"path");
        auto iRECScintChi   =   c12->getBankOrder(idx_RECScint,"chi2");
        auto iRECScintX     =   c12->getBankOrder(idx_RECScint,"x");
        auto iRECScintY     =   c12->getBankOrder(idx_RECScint,"y");
        auto iRECScintZ     =   c12->getBankOrder(idx_RECScint,"z");
        auto iRECScintHX     =   c12->getBankOrder(idx_RECScint,"hx");
        auto iRECScintHY     =   c12->getBankOrder(idx_RECScint,"hy");
        auto iRECScintHZ     =   c12->getBankOrder(idx_RECScint,"hz");
        auto iRECScintStat  =   c12->getBankOrder(idx_RECScint,"status");
  


        
        while(c12->next() == true){
            //cout<<" Event ======= "<<index<<endl;
            nP=0; nmb=0; store=false;
            index++;
            if(c12->getDetParticles().empty())
                continue;
            
            
            startTime = c12->event()->getStartTime();
            helicity = c12->event()->getHelicity();
            beamQ = c12->event()->getBeamCharge();

            // Apply the same PID + fiducial/vertex cuts as example_ana.cpp
            clasAna.Run(c12);
            auto electrons = clasAna.getByPid(11);
            auto protons_sel = clasAna.getByPid(2212);
            if(electrons.size() != 1) continue;

            auto electron = electrons[0];

            // Use corrected electron 4-vector (energy-loss / angle / momentum corrections)
            TLorentzVector el_corr(0,0,0,0.000511);
            GetLorentzVector_Corrected(el_corr, electron, isMC);
            p4_electron = el_corr;

            nP=0;
            vgs = p4_beam - p4_electron;
            Q2 = -1.*vgs.Mag2();
            Nu = vgs.E();
            xB = Q2/(2*0.938*Nu);
            q = vgs.P();
            qx= vgs.Px();
            qy= vgs.Py();
            qz = vgs.Pz();
            qtheta = vgs.Theta();
            qphi = vgs.Phi();

            Epid[0]=electron->par()->getPid();
            Echarge[0]=electron->par()->getCharge();
            Ep[0]= p4_electron.P();
            Epx[0] = p4_electron.Px();
            Epy[0] = p4_electron.Py();
            Epz[0] = p4_electron.Pz();
            Etheta[0] = p4_electron.Theta();
            Ephi[0] = p4_electron.Phi();
            Etime[0] = electron->getTime();
            Epath[0] = electron->getPath();
            Evz[0] = electron->par()->getVz();
            Estat[0] = electron->par()->getStatus();
            EdE[0] = electron->getDeltaEnergy();
            Ebeta[0] = electron->par()->getBeta();
            
            Esector[0] = electron->getSector();
            Eregion[0] = electron->getRegion();
            
            
            ePCAL_time[0] = electron->cal(clas12::PCAL)->getTime();
            ePCAL_energy[0] = electron->cal(clas12::PCAL)->getEnergy();
            ePCAL_path[0] = electron->cal(clas12::PCAL)->getPath();
            ePCAL_sector[0] = electron->cal(clas12::PCAL)->getSector();
            ePCAL_x[0] = electron->cal(clas12::PCAL)->getX();
            ePCAL_y[0] = electron->cal(clas12::PCAL)->getY();
            ePCAL_z[0] = electron->cal(clas12::PCAL)->getZ();
            ePCAL_lu[0] = electron->cal(clas12::PCAL)->getLu();
            ePCAL_lv[0] = electron->cal(clas12::PCAL)->getLv();
            ePCAL_lw[0] = electron->cal(clas12::PCAL)->getLw();
            ePCAL_status[0] = electron->cal(clas12::PCAL)->getStatus();
            
            eECIN_time[0] = electron->cal(clas12::ECIN)->getTime();
            eECIN_energy[0] = electron->cal(clas12::ECIN)->getEnergy();
            eECIN_path[0] = electron->cal(clas12::ECIN)->getPath();
            eECIN_sector[0] = electron->cal(clas12::ECIN)->getSector();
            eECIN_x[0] = electron->cal(clas12::ECIN)->getX();
            eECIN_y[0] = electron->cal(clas12::ECIN)->getY();
            eECIN_z[0] = electron->cal(clas12::ECIN)->getZ();
            eECIN_lu[0] = electron->cal(clas12::ECIN)->getLu();
            eECIN_lv[0] = electron->cal(clas12::ECIN)->getLv();
            eECIN_lw[0] = electron->cal(clas12::ECIN)->getLw();
            eECIN_status[0] = electron->cal(clas12::ECIN)->getStatus();
            
            
            eECOUT_time[0] = electron->cal(clas12::ECOUT)->getTime();
            eECOUT_energy[0] = electron->cal(clas12::ECOUT)->getEnergy();
            eECOUT_path[0] = electron->cal(clas12::ECOUT)->getPath();
            eECOUT_sector[0] = electron->cal(clas12::ECOUT)->getSector();
            eECOUT_x[0] = electron->cal(clas12::ECOUT)->getX();
            eECOUT_y[0] = electron->cal(clas12::ECOUT)->getY();
            eECOUT_z[0] = electron->cal(clas12::ECOUT)->getZ();
            eECOUT_lu[0] = electron->cal(clas12::ECOUT)->getLu();
            eECOUT_lv[0] = electron->cal(clas12::ECOUT)->getLv();
            eECOUT_lw[0] = electron->cal(clas12::ECOUT)->getLw();
            eECOUT_status[0] = electron->cal(clas12::ECOUT)->getStatus();
            

            
            eCHindex[0]     = electron->che(clas12::HTCC)->getPindex();
            eCHdetector[0]  = electron->che(clas12::HTCC)->getDetector();
            eCHnphe[0]      = electron->che(clas12::HTCC)->getNphe();
            eCHtime[0]      = electron->che(clas12::HTCC)->getTime();
            eCHpath[0]      = electron->che(clas12::HTCC)->getPath();
            eCHsector[0]    = electron->che(clas12::HTCC)->getSector();
            eCHx[0]         = electron->che(clas12::HTCC)->getX();
            eCHy[0]         = electron->che(clas12::HTCC)->getY();
            eCHz[0]         = electron->che(clas12::HTCC)->getZ();
            eCHdtheta[0]    = electron->che(clas12::HTCC)->getDtheta();
            eCHdphi[0]      = electron->che(clas12::HTCC)->getDPhi();
            eCHstatus[0]    = electron->che(clas12::HTCC)->getStatus();
            
            
            EdcX[0][0] = electron->traj(6,6)->getX();
            EdcY[0][0] = electron->traj(6,6)->getY();
            EdcZ[0][0] = electron->traj(6,6)->getZ();
            
            EdcX[0][1] = electron->traj(6,18)->getX();
            EdcY[0][1] = electron->traj(6,18)->getY();
            EdcZ[0][1] = electron->traj(6,18)->getZ();
            
            EdcX[0][2] = electron->traj(6,36)->getX();
            EdcY[0][2] = electron->traj(6,36)->getY();
            EdcZ[0][2] = electron->traj(6,36)->getZ();
            
        
            
            
            nml=1; // meanwhile analyze only for one electron
            
            //======= proton ========
            //nmb = c12->getByCharge(1).size();
            nmb = c12->getDetParticles().size();
            for(auto ii=0;ii<nmb;ii++){
                auto prot = c12->getDetParticles()[ii];
                
                
                p3_prot.SetXYZ(prot->par()->getPx(),prot->par()->getPy(),prot->par()->getPz());
                
                ppid[ii]=prot->par()->getPid();
                pcharge[ii]=prot->par()->getCharge();
                pp[ii]= p3_prot.Mag();;
                ppx[ii] = p3_prot.X();
                ppy[ii] = p3_prot.Y();
                ppz[ii] = p3_prot.Z();
                ptheta[ii] = p3_prot.Theta();
                pphi[ii] = p3_prot.Phi();
                ptime[ii] = prot->getTime();
                ppath[ii] = prot->getPath();
                pvz[ii] = prot->par()->getVz();
                pstat[ii] = prot->par()->getStatus();
                pdE[ii] = prot->getDeltaEnergy();
                pbeta[ii] = prot->par()->getBeta();
                
                // if(ppid[ii] == 2212) {
                //     TLorentzVector p4_corr(0,0,0,mp);      // mp already defined above
                //     GetLorentzVector_Corrected(p4_corr, prot, isMC);
                //     prp[nP] = p4_corr.P();
                //     prpx[nP] = p4_corr.X();
                //     prpy[nP] = p4_corr.Y();
                //     prpz[nP] = p4_corr.Z();
                //     prtheta[nP] = p4_corr.Theta();
                //     prphi[nP] = p4_corr.Phi();
                //     nP++;
                // }
                
                // if(ii==0){
                //     p4_prot.SetXYZM(ppx[0],ppy[0],ppz[0],0.938);
                // }
                
                
                psector[ii] = prot->getSector();
                pregion[ii] = prot->getRegion();

                pFTOF1Atime[ii] = prot->sci(clas12::FTOF1A)->getTime();
                pFTOF1ADet[ii]  = prot->sci(clas12::FTOF1A)->getDetector();
                pFTOF1ALay[ii]  = prot->sci(clas12::FTOF1A)->getLayer();
                pFTOF1AComp[ii]  = prot->sci(clas12::FTOF1A)->getComponent();
                pFTOF1Apath[ii]  = prot->sci(clas12::FTOF1A)->getPath();
                pFTOF1Asector[ii] = prot->sci(clas12::FTOF1A)->getSector();
                pFTOF1Astatus[ii] = prot->sci(clas12::FTOF1A)->getStatus();
                pFTOF1Ax[ii] = prot->sci(clas12::FTOF1A)->getX();
                pFTOF1Ay[ii] = prot->sci(clas12::FTOF1A)->getY();
                pFTOF1Az[ii] = prot->sci(clas12::FTOF1A)->getZ();

                pFTOF2time[ii] = prot->sci(clas12::FTOF2)->getTime();
                pFTOF2Det[ii]  = prot->sci(clas12::FTOF2)->getDetector();
                pFTOF2Lay[ii]  = prot->sci(clas12::FTOF2)->getLayer();
                pFTOF2Comp[ii]  = prot->sci(clas12::FTOF2)->getComponent();
                pFTOF2path[ii]  = prot->sci(clas12::FTOF2)->getPath();
                pFTOF2sector[ii] = prot->sci(clas12::FTOF2)->getSector();
                pFTOF2status[ii] = prot->sci(clas12::FTOF2)->getStatus();
                pFTOF2x[ii] = prot->sci(clas12::FTOF2)->getX();
                pFTOF2y[ii] = prot->sci(clas12::FTOF2)->getY();
                pFTOF2z[ii] = prot->sci(clas12::FTOF2)->getZ();

                pFTOF1Btime[ii] = prot->sci(clas12::FTOF1B)->getTime();
                pFTOF1BDet[ii]  = prot->sci(clas12::FTOF1B)->getDetector();
                pFTOF1BLay[ii]  = prot->sci(clas12::FTOF1B)->getLayer();
                pFTOF1BComp[ii]  = prot->sci(clas12::FTOF1B)->getComponent();
                pFTOF1Bpath[ii]  = prot->sci(clas12::FTOF1B)->getPath();
                pFTOF1Bsector[ii] = prot->sci(clas12::FTOF1B)->getSector();
                pFTOF1Bstatus[ii] = prot->sci(clas12::FTOF1B)->getStatus();
                pFTOF1Bx[ii] = prot->sci(clas12::FTOF1B)->getX();
                pFTOF1By[ii] = prot->sci(clas12::FTOF1B)->getY();
                pFTOF1Bz[ii] = prot->sci(clas12::FTOF1B)->getZ();


                pCND1time[ii] = prot->sci(clas12::CND1)->getTime();
                pCND1Det[ii]  = prot->sci(clas12::CND1)->getDetector();
                pCND1Lay[ii]  = prot->sci(clas12::CND1)->getLayer();
                pCND1Comp[ii]  = prot->sci(clas12::CND1)->getComponent();
                pCND1path[ii]  = prot->sci(clas12::CND1)->getPath();
                pCND1sector[ii] = prot->sci(clas12::CND1)->getSector();
                pCND1status[ii] = prot->sci(clas12::CND1)->getStatus();
                pCND1Energy[ii] = prot->sci(clas12::CND1)->getEnergy();
                pCND1x[ii] = prot->sci(clas12::CND1)->getX();
                pCND1y[ii] = prot->sci(clas12::CND1)->getY();
                pCND1z[ii] = prot->sci(clas12::CND1)->getZ();

                pCND2time[ii] = prot->sci(clas12::CND2)->getTime();
                pCND2Det[ii]  = prot->sci(clas12::CND2)->getDetector();
                pCND2Lay[ii]  = prot->sci(clas12::CND2)->getLayer();
                pCND2Comp[ii]  = prot->sci(clas12::CND2)->getComponent();
                pCND2path[ii]  = prot->sci(clas12::CND2)->getPath();
                pCND2sector[ii] = prot->sci(clas12::CND2)->getSector();
                pCND2status[ii] = prot->sci(clas12::CND2)->getStatus();
                pCND2Energy[ii] = prot->sci(clas12::CND2)->getEnergy();
                pCND2x[ii] = prot->sci(clas12::CND2)->getX();
                pCND2y[ii] = prot->sci(clas12::CND2)->getY();
                pCND2z[ii] = prot->sci(clas12::CND2)->getZ();


                pCND3time[ii] = prot->sci(clas12::CND3)->getTime();
                pCND3Det[ii]  = prot->sci(clas12::CND3)->getDetector();
                pCND3Lay[ii]  = prot->sci(clas12::CND3)->getLayer();
                pCND3Comp[ii]  = prot->sci(clas12::CND3)->getComponent();
                pCND3path[ii]  = prot->sci(clas12::CND3)->getPath();
                pCND3sector[ii] = prot->sci(clas12::CND3)->getSector();
                pCND3status[ii] = prot->sci(clas12::CND3)->getStatus();
                pCND3Energy[ii] = prot->sci(clas12::CND3)->getEnergy();
                pCND3x[ii] = prot->sci(clas12::CND3)->getX();
                pCND3y[ii] = prot->sci(clas12::CND3)->getY();
                pCND3z[ii] = prot->sci(clas12::CND3)->getZ();

                pCTOFtime[ii] = prot->sci(clas12::CTOF)->getTime();
                pCTOFDet[ii]  = prot->sci(clas12::CTOF)->getDetector();
                pCTOFLay[ii]  = prot->sci(clas12::CTOF)->getLayer();
                pCTOFComp[ii]  = prot->sci(clas12::CTOF)->getComponent();
                pCTOFpath[ii]  = prot->sci(clas12::CTOF)->getPath();
                pCTOFsector[ii] = prot->sci(clas12::CTOF)->getSector();
                pCTOFstatus[ii] = prot->sci(clas12::CTOF)->getStatus();
                pCTOFEnergy[ii] = prot->sci(clas12::CTOF)->getEnergy();
                pCTOFx[ii] = prot->sci(clas12::CTOF)->getX();
                pCTOFy[ii] = prot->sci(clas12::CTOF)->getY();
                pCTOFz[ii] = prot->sci(clas12::CTOF)->getZ();

		pCVTpind[ii] = prot->traj(5,12)->getPindex();
                pCVTDet[ii] = prot->traj(5,12)->getDetector();
                pCVTLay[ii] = prot->traj(5,12)->getLayer();
                pCVTx[ii] = prot->traj(5,12)->getX();
                pCVTy[ii] = prot->traj(5,12)->getY();
                pCVTz[ii] = prot->traj(5,12)->getZ();
                pCVTtrkpind[ii] = prot->trk(5)->getPindex();
                pCVTtrkDet[ii] = prot->trk(5)->getDetector();
                pCVTtrkNDF[ii] = prot->trk(5)->getNDF();
                pCVTtrkSec[ii] = prot->trk(5)->getSector();
                pCVTtrkChi2[ii] = prot->trk(5)->getChi2();
                
                pPCAL_time[ii] = prot->cal(clas12::PCAL)->getTime();
                pPCAL_energy[ii] = prot->cal(clas12::PCAL)->getEnergy();
		pPCAL_det[ii] = prot->cal(clas12::PCAL)->getDetector();
                pPCAL_path[ii] = prot->cal(clas12::PCAL)->getPath();
                pPCAL_sector[ii] = prot->cal(clas12::PCAL)->getSector();
                pPCAL_x[ii] = prot->cal(clas12::PCAL)->getX();
                pPCAL_y[ii] = prot->cal(clas12::PCAL)->getY();
                pPCAL_z[ii] = prot->cal(clas12::PCAL)->getZ();
                pPCAL_lu[ii] = prot->cal(clas12::PCAL)->getLu();
                pPCAL_lv[ii] = prot->cal(clas12::PCAL)->getLv();
                pPCAL_lw[ii] = prot->cal(clas12::PCAL)->getLw();
                pPCAL_status[ii] = prot->cal(clas12::PCAL)->getStatus();
                
                pECIN_time[ii] = prot->cal(clas12::ECIN)->getTime();
                pECIN_energy[ii] = prot->cal(clas12::ECIN)->getEnergy();
		pECIN_det[ii] = prot->cal(clas12::ECIN)->getDetector();
                pECIN_path[ii] = prot->cal(clas12::ECIN)->getPath();
                pECIN_sector[ii] = prot->cal(clas12::ECIN)->getSector();
                pECIN_x[ii] = prot->cal(clas12::ECIN)->getX();
                pECIN_y[ii] = prot->cal(clas12::ECIN)->getY();
                pECIN_z[ii] = prot->cal(clas12::ECIN)->getZ();
                pECIN_lu[ii] = prot->cal(clas12::ECIN)->getLu();
                pECIN_lv[ii] = prot->cal(clas12::ECIN)->getLv();
                pECIN_lw[ii] = prot->cal(clas12::ECIN)->getLw();
                pECIN_status[ii] = prot->cal(clas12::ECIN)->getStatus();
                
                
                pECOUT_time[ii] = prot->cal(clas12::ECOUT)->getTime();
                pECOUT_energy[ii] = prot->cal(clas12::ECOUT)->getEnergy();
		pECOUT_det[ii] = prot->cal(clas12::ECOUT)->getDetector();
                pECOUT_path[ii] = prot->cal(clas12::ECOUT)->getPath();
                pECOUT_sector[ii] = prot->cal(clas12::ECOUT)->getSector();
                pECOUT_x[ii] = prot->cal(clas12::ECOUT)->getX();
                pECOUT_y[ii] = prot->cal(clas12::ECOUT)->getY();
                pECOUT_z[ii] = prot->cal(clas12::ECOUT)->getZ();
                pECOUT_lu[ii] = prot->cal(clas12::ECOUT)->getLu();
                pECOUT_lv[ii] = prot->cal(clas12::ECOUT)->getLv();
                pECOUT_lw[ii] = prot->cal(clas12::ECOUT)->getLw();
                pECOUT_status[ii] = prot->cal(clas12::ECOUT)->getStatus();
                
                
                
                pdcX[ii][0] = prot->traj(6,6)->getX();
                pdcY[ii][0] = prot->traj(6,6)->getY();
                pdcZ[ii][0] = prot->traj(6,6)->getZ();
                
                pdcX[ii][1] = prot->traj(6,18)->getX();
                pdcY[ii][1] = prot->traj(6,18)->getY();
                pdcZ[ii][1] = prot->traj(6,18)->getZ();
                
                pdcX[ii][2] = prot->traj(6,36)->getX();
                pdcY[ii][2] = prot->traj(6,36)->getY();
                pdcZ[ii][2] = prot->traj(6,36)->getZ();
               


		ptrajPind[ii][0] = prot->traj(7,1)->getPindex();
		ptrajDet[ii][0] = prot->traj(7,1)->getDetector();
		ptrajLay[ii][0] = prot->traj(7,1)->getLayer();
		ptrajX[ii][0] = prot->traj(7,1)->getX();
		ptrajY[ii][0] = prot->traj(7,1)->getY();
		ptrajZ[ii][0] = prot->traj(7,1)->getZ();
		ptrajPath[ii][0] = prot->traj(7,1)->getPath();

                ptrajPind[ii][1] = prot->traj(7,4)->getPindex();
                ptrajDet[ii][1] = prot->traj(7,4)->getDetector();
                ptrajLay[ii][1] = prot->traj(7,4)->getLayer();
                ptrajX[ii][1] = prot->traj(7,4)->getX();
                ptrajY[ii][1] = prot->traj(7,4)->getY();
                ptrajZ[ii][1] = prot->traj(7,4)->getZ();
                ptrajPath[ii][1] = prot->traj(7,4)->getPath();

                ptrajPind[ii][2] = prot->traj(7,7)->getPindex();
                ptrajDet[ii][2] = prot->traj(7,7)->getDetector();
                ptrajLay[ii][2] = prot->traj(7,7)->getLayer();
                ptrajX[ii][2] = prot->traj(7,7)->getX();
                ptrajY[ii][2] = prot->traj(7,7)->getY();
                ptrajZ[ii][2] = prot->traj(7,7)->getZ();
                ptrajPath[ii][2] = prot->traj(7,7)->getPath(); 

                
                
            }// loop over proton


            // ---- corrected + selected protons (fiducial/vertex/PID via clas12ana) ----
            // Fill the dedicated proton branches (pr*) ONLY from clasAna-selected protons
            nP = 0;
            for(const auto& psel : protons_sel){
                if(nP >= 1000) break; // safety

                // corrected proton kinematics (energy-loss / angle / momentum corrections)
                TLorentzVector p4_corr(0,0,0,mp);
                GetLorentzVector_Corrected(p4_corr, psel, isMC);

                prp[nP]      = p4_corr.P();
                prpx[nP]     = p4_corr.Px();
                prpy[nP]     = p4_corr.Py();
                prpz[nP]     = p4_corr.Pz();
                prtheta[nP]  = p4_corr.Theta();
                prphi[nP]    = p4_corr.Phi();

                // raw / auxiliary info for the SAME selected proton
                prpid[nP]    = psel->par()->getPid();
                prcharge[nP] = psel->par()->getCharge();
                prvz[nP]     = psel->par()->getVz();
                prstat[nP]   = psel->par()->getStatus();
                prtime[nP]   = psel->getTime();
                prpath[nP]   = psel->getPath();
                prdE[nP]     = psel->getDeltaEnergy();
                prbeta[nP]   = psel->par()->getBeta();

                

                nP++;
            }
            
            hitInd.clear();
            hitpInd.clear();
            hitDet.clear();
            hitSec.clear();
            hitLay.clear();
            hitComp.clear();
            hitEnergy.clear();
            hitTime.clear();
            hitPath.clear();
            hitChi2.clear();
            hitX.clear();
            hitY.clear();
            hitZ.clear();
            hitHX.clear();
            hitHY.clear();
            hitHZ.clear();
            hitStat.clear();
  
            
            for(auto ipa=0;ipa<c12->getBank(idx_RECScint)->getRows();ipa++){
                
                auto tind   =   c12->getBank(idx_RECScint)->getInt(iRECScintInd,ipa);
                auto tpind  =   c12->getBank(idx_RECScint)->getInt(iRECScintPind,ipa);
                auto tdet   =   c12->getBank(idx_RECScint)->getInt(iRECScintDet,ipa);
                auto tsec   =   c12->getBank(idx_RECScint)->getInt(iRECScintSec,ipa);
                auto tlay   =   c12->getBank(idx_RECScint)->getInt(iRECScintLay,ipa);
                auto tcom   =   c12->getBank(idx_RECScint)->getInt(iRECScintCom,ipa);
                auto tene   =   c12->getBank(idx_RECScint)->getFloat(iRECScintEn,ipa);
                auto ttim   =   c12->getBank(idx_RECScint)->getFloat(iRECScintTime,ipa);
                auto tpat   =   c12->getBank(idx_RECScint)->getFloat(iRECScintPath,ipa);
                auto tchi   =   c12->getBank(idx_RECScint)->getFloat(iRECScintChi,ipa);
                auto tx     =   c12->getBank(idx_RECScint)->getFloat(iRECScintX,ipa);
                auto ty     =   c12->getBank(idx_RECScint)->getFloat(iRECScintY,ipa);
                auto tz     =   c12->getBank(idx_RECScint)->getFloat(iRECScintZ,ipa);
                auto thx     =   c12->getBank(idx_RECScint)->getFloat(iRECScintHX,ipa);
                auto thy     =   c12->getBank(idx_RECScint)->getFloat(iRECScintHY,ipa);
                auto thz     =   c12->getBank(idx_RECScint)->getFloat(iRECScintHZ,ipa);
                auto tsta   =   c12->getBank(idx_RECScint)->getInt(iRECScintStat,ipa);

                hitInd.push_back(tind);
                hitpInd.push_back(tpind);
                hitDet.push_back(tdet);
                hitSec.push_back(tsec);
                hitLay.push_back(tlay);
                hitComp.push_back(tcom);
                hitEnergy.push_back(tene);
                hitTime.push_back(ttim);
                hitPath.push_back(tpat);
                hitChi2.push_back(tchi);
                hitX.push_back(tx);
                hitY.push_back(ty);
                hitZ.push_back(tz);
                hitHX.push_back(thx);
                hitHY.push_back(thy);
                hitHZ.push_back(thz);
                hitStat.push_back(tsta);
            }
            
            hitDeDx.clear();
            hitSize.clear();
            hitLayMult.clear();
            
            for(auto ipa=0;ipa<c12->getBank(idx_RECScintExtra)->getRows();ipa++){
                
                auto tdedx   =   c12->getBank(idx_RECScintExtra)->getFloat(iRECScintExtraDEDX,ipa);
                auto tsize  =   c12->getBank(idx_RECScintExtra)->getInt(iRECScintExtraSize,ipa);
                auto tlmul   =   c12->getBank(idx_RECScintExtra)->getInt(iRECScintExtraLMult,ipa);
                
                hitDeDx.push_back(tdedx);
                hitSize.push_back(tsize);
                hitLayMult.push_back(tlmul);

                
            }
            store = false;
            if(nml==1 && nP>0){
                for(Int_t jj=0;jj<nP;jj++){
                    
                    p4_prot.SetXYZM(prpx[jj],prpy[jj],prpz[jj],0.938);
        
                    
                    p4_miss = p4_beam + p4_target - p4_electron - p4_prot;
                    pmisp[jj] = p4_miss.P();
                    pmispx[jj] = p4_miss.Px();
                    pmispy[jj] = p4_miss.Py();
                    pmispz[jj] = p4_miss.Pz();
                    pmisM2[jj] = p4_miss.M2();
                    pmistheta[jj] = p4_miss.Theta();
                    pmisphi[jj] = p4_miss.Phi();
                    

		    double angleTemp = TMath::ACos((qx*prpx[jj] + qy*prpy[jj] + qz*prpz[jj])/(prp[jj]*q))*180./TMath::Pi();
		    double LeadFraq = prp[jj]/q;
                    if(TMath::Sqrt(pmisM2[jj])<1.6 && pmistheta[jj]*180/TMath::Pi() > 40 && pmistheta[jj]*180/TMath::Pi()<145) store=true;
                }
                
            }
            

            




            
            if(nml==1 && Q2>1. && xB>1.2 && nP == 1 && store == true)
                T->Fill();
            
        }
        
       // if(ifile>20) break;
    }
  rFile->Write();
  rFile->Close();
}



