#include <cmath>
#include <iostream>

#include <TBranch.h>
#include <TFile.h>
#include <TNamed.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

double getBeamEnergyFromMetadata(TFile *file, double defaultBeam)
{
  TNamed *meta = dynamic_cast<TNamed*>(file->Get("skimmer_metadata"));
  if(!meta){
    return defaultBeam;
  }

  TString text = meta->GetTitle();
  Ssiz_t start = text.Index("beam_energy=");
  if(start == kNPOS){
    return defaultBeam;
  }

  start += TString("beam_energy=").Length();
  Ssiz_t end = text.Index("\n", start);
  TString value = (end == kNPOS) ? text(start, text.Length() - start) : text(start, end - start);
  return value.Atof();
}

void appendRotatedPmissBranches(const char *filename,
                                const char *treeName = "deeppim",
                                double beamEnergyOverride = -1.)
{
  TFile *file = TFile::Open(filename, "UPDATE");
  if(!file || file->IsZombie()){
    std::cerr << "Could not open " << filename << " in UPDATE mode.\n";
    return;
  }

  TTree *tree = dynamic_cast<TTree*>(file->Get(treeName));
  if(!tree){
    std::cerr << "Could not find tree " << treeName << " in " << filename << ".\n";
    file->Close();
    return;
  }

  if(tree->GetBranch("pxMiss_rot") || tree->GetBranch("pyMiss_rot") || tree->GetBranch("pzMiss_rot")){
    std::cerr << "Rotated missing-momentum branches already exist. Leaving file unchanged.\n";
    file->Close();
    return;
  }

  double p1P = 0.;
  double p1Px = 0.;
  double p1Py = 0.;
  double p1Pz = 0.;
  double p2P = 0.;
  double p2Px = 0.;
  double p2Py = 0.;
  double p2Pz = 0.;
  double ePx = 0.;
  double ePy = 0.;
  double ePz = 0.;
  double pxMiss = 0.;
  double pyMiss = 0.;
  double pzMiss = 0.;

  tree->SetBranchAddress("p1P", &p1P);
  tree->SetBranchAddress("p1Px", &p1Px);
  tree->SetBranchAddress("p1Py", &p1Py);
  tree->SetBranchAddress("p1Pz", &p1Pz);
  tree->SetBranchAddress("p2P", &p2P);
  tree->SetBranchAddress("p2Px", &p2Px);
  tree->SetBranchAddress("p2Py", &p2Py);
  tree->SetBranchAddress("p2Pz", &p2Pz);
  tree->SetBranchAddress("ePx", &ePx);
  tree->SetBranchAddress("ePy", &ePy);
  tree->SetBranchAddress("ePz", &ePz);
  tree->SetBranchAddress("pxMiss", &pxMiss);
  tree->SetBranchAddress("pyMiss", &pyMiss);
  tree->SetBranchAddress("pzMiss", &pzMiss);

  double pxMiss_rot = -9.;
  double pyMiss_rot = -9.;
  double pzMiss_rot = -9.;
  double pMiss_rot = -9.;
  double thetaMiss_rot = -9.;
  double phiMiss_rot = -9.;
  int leadProtonIndex_rot = -1;

  TBranch *b_pxMiss_rot = tree->Branch("pxMiss_rot", &pxMiss_rot, "pxMiss_rot/D");
  TBranch *b_pyMiss_rot = tree->Branch("pyMiss_rot", &pyMiss_rot, "pyMiss_rot/D");
  TBranch *b_pzMiss_rot = tree->Branch("pzMiss_rot", &pzMiss_rot, "pzMiss_rot/D");
  TBranch *b_pMiss_rot = tree->Branch("pMiss_rot", &pMiss_rot, "pMiss_rot/D");
  TBranch *b_thetaMiss_rot = tree->Branch("thetaMiss_rot", &thetaMiss_rot, "thetaMiss_rot/D");
  TBranch *b_phiMiss_rot = tree->Branch("phiMiss_rot", &phiMiss_rot, "phiMiss_rot/D");
  TBranch *b_leadProtonIndex_rot = tree->Branch("leadProtonIndex_rot", &leadProtonIndex_rot, "leadProtonIndex_rot/I");

  double beamEnergy = (beamEnergyOverride > 0.) ? beamEnergyOverride : getBeamEnergyFromMetadata(file, 5.98636);
  TVector3 beam(0., 0., beamEnergy);

  Long64_t nEntries = tree->GetEntries();
  Long64_t nFilled = 0;
  Long64_t nBadBasis = 0;

  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);

    TVector3 miss(pxMiss, pyMiss, pzMiss);
    TVector3 q = beam - TVector3(ePx, ePy, ePz);

    pxMiss_rot = -9.;
    pyMiss_rot = -9.;
    pzMiss_rot = -9.;
    pMiss_rot = -9.;
    thetaMiss_rot = -9.;
    phiMiss_rot = -9.;
    leadProtonIndex_rot = -1;

    TVector3 p1(p1Px, p1Py, p1Pz);
    TVector3 p2(p2Px, p2Py, p2Pz);
    TVector3 lead = (p1P >= p2P) ? p1 : p2;
    leadProtonIndex_rot = (p1P >= p2P) ? 1 : 2;
    TVector3 leadMinusQ = lead - q;

    TVector3 yRaw = leadMinusQ.Cross(q);
    if(miss.Mag() > 0. && q.Mag() > 0. && leadMinusQ.Mag() > 0. && yRaw.Mag() > 0.){
      TVector3 vz = leadMinusQ.Unit();
      TVector3 vy = yRaw.Unit();
      TVector3 vx = vz.Cross(vy).Unit();
      TVector3 missRot(miss.Dot(vx), miss.Dot(vy), miss.Dot(vz));

      pxMiss_rot = missRot.X();
      pyMiss_rot = missRot.Y();
      pzMiss_rot = missRot.Z();
      pMiss_rot = missRot.Mag();
      thetaMiss_rot = missRot.Theta() * 180. / M_PI;
      phiMiss_rot = missRot.Phi() * 180. / M_PI;
      nFilled++;
    }
    else{
      nBadBasis++;
    }

    b_pxMiss_rot->Fill();
    b_pyMiss_rot->Fill();
    b_pzMiss_rot->Fill();
    b_pMiss_rot->Fill();
    b_thetaMiss_rot->Fill();
    b_phiMiss_rot->Fill();
    b_leadProtonIndex_rot->Fill();
  }

  tree->Write("", TObject::kOverwrite);
  file->Close();

  std::cout << "Appended rotated pmiss branches to " << filename
            << " using Ebeam=" << beamEnergy << " GeV. Filled " << nFilled
            << " entries";
  if(nBadBasis > 0){
    std::cout << " (" << nBadBasis << " entries had an undefined basis and were set to -9)";
  }
  std::cout << ".\n";
}
