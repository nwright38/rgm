#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TTree.h>
#include <TVector3.h>

#include "Corrections.h"
#include "HipoChain.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

const double mP = 0.9382720813;
const double mD = 1.8756;
const double mPi = 0.13957039;
const double me = 0.000511;
const double rad2deg = 180. / M_PI;

struct TopologyConfig {
  string name;
  int id = 0;
  int minProtons = 2;
  int minPiMinus = 1;
  bool useP2 = true;
  bool usePiMinus = true;
  double targetMass = mD;
  string targetName = "deuterium";
};

struct ParticleBranches {
  Double_t p = -9.;
  Double_t theta = -9.;
  Double_t phi = -9.;
  Double_t px = -9.;
  Double_t py = -9.;
  Double_t pz = -9.;
  Double_t energy = -9.;
  Double_t vz = -99.;
  Int_t status = -9;
  Int_t region = -9;

  void reset()
  {
    p = theta = phi = px = py = pz = energy = -9.;
    vz = -99.;
    status = -9;
    region = -9;
  }

  void book(TTree *tree, const string &prefix)
  {
    tree->Branch((prefix + "P").c_str(), &p, (prefix + "P/D").c_str());
    tree->Branch((prefix + "Theta").c_str(), &theta, (prefix + "Theta/D").c_str());
    tree->Branch((prefix + "Phi").c_str(), &phi, (prefix + "Phi/D").c_str());
    tree->Branch((prefix + "Px").c_str(), &px, (prefix + "Px/D").c_str());
    tree->Branch((prefix + "Py").c_str(), &py, (prefix + "Py/D").c_str());
    tree->Branch((prefix + "Pz").c_str(), &pz, (prefix + "Pz/D").c_str());
    tree->Branch((prefix + "E").c_str(), &energy, (prefix + "E/D").c_str());
    tree->Branch((prefix + "Vz").c_str(), &vz, (prefix + "Vz/D").c_str());
    tree->Branch((prefix + "Status").c_str(), &status, (prefix + "Status/I").c_str());
    tree->Branch((prefix + "Region").c_str(), &region, (prefix + "Region/I").c_str());
  }

  void fill(TLorentzVector const &p4, clas12::region_part_ptr part)
  {
    p = p4.P();
    theta = p4.Theta() * rad2deg;
    phi = p4.Phi() * rad2deg;
    px = p4.Px();
    py = p4.Py();
    pz = p4.Pz();
    energy = p4.E();
    vz = part->par()->getVz();
    status = part->getStatus();
    region = static_cast<int>(part->getRegion());
  }
};

struct CorrectedParticle {
  int index = -1;
  clas12::region_part_ptr part;
  TLorentzVector p4;
};

void Usage()
{
  cerr << "Usage: ./skim_deeppim <MC=1,Data=0> [Ebeam(GeV)] <output.root> <input.hipo> [more hipo...] [options]\n"
       << "Options:\n"
       << "  --topology deeppim|eppim|epp|ep   default: deeppim\n"
       << "  --target deuterium|hydrogen       default: topology-dependent\n"
       << "  --beam Ebeam                      default: 5.98636 GeV\n"
       << "  --max-rows N                      optional quick-test row limit\n";
}

void setCorrectedP4(TLorentzVector &p4, clas12::region_part_ptr part, double mass, bool isMC)
{
  p4.SetXYZM(0., 0., 0., mass);
  GetLorentzVector_Corrected(p4, part, isMC);
}

TopologyConfig makeTopology(string name)
{
  transform(name.begin(), name.end(), name.begin(), [](unsigned char c){ return tolower(c); });

  TopologyConfig cfg;
  cfg.name = name;

  if(name == "ep" || name == "h_elastic" || name == "hydrogen_elastic"){
    cfg.id = 1;
    cfg.name = "ep";
    cfg.minProtons = 1;
    cfg.minPiMinus = 0;
    cfg.useP2 = false;
    cfg.usePiMinus = false;
    cfg.targetMass = mP;
    cfg.targetName = "hydrogen";
  }
  else if(name == "epp"){
    cfg.id = 2;
    cfg.minProtons = 2;
    cfg.minPiMinus = 0;
    cfg.useP2 = true;
    cfg.usePiMinus = false;
    cfg.targetMass = mD;
    cfg.targetName = "deuterium";
  }
  else if(name == "eppim"){
    cfg.id = 3;
    cfg.minProtons = 1;
    cfg.minPiMinus = 1;
    cfg.useP2 = false;
    cfg.usePiMinus = true;
    cfg.targetMass = mD;
    cfg.targetName = "deuterium";
  }
  else if(name == "deeppim" || name == "epppim"){
    cfg.id = 4;
    cfg.name = "deeppim";
    cfg.minProtons = 2;
    cfg.minPiMinus = 1;
    cfg.useP2 = true;
    cfg.usePiMinus = true;
    cfg.targetMass = mD;
    cfg.targetName = "deuterium";
  }
  else{
    cerr << "Unknown topology '" << name << "'.\n";
    Usage();
    exit(-1);
  }

  return cfg;
}

bool isNumber(const string &text)
{
  char *end = nullptr;
  strtod(text.c_str(), &end);
  return end != text.c_str() && end != nullptr && *end == '\0';
}

void applyTargetOverride(TopologyConfig &cfg, string target)
{
  transform(target.begin(), target.end(), target.begin(), [](unsigned char c){ return tolower(c); });
  if(target == "hydrogen" || target == "h" || target == "proton"){
    cfg.targetMass = mP;
    cfg.targetName = "hydrogen";
  }
  else if(target == "deuterium" || target == "d" || target == "deuteron"){
    cfg.targetMass = mD;
    cfg.targetName = "deuterium";
  }
  else{
    cerr << "Unknown target '" << target << "'.\n";
    Usage();
    exit(-1);
  }
}

void fillRotatedMissingBranches(TVector3 const &vec,
                                TVector3 const &zAxis,
                                TVector3 const &planeAxis,
                                Double_t &pxMissRot,
                                Double_t &pyMissRot,
                                Double_t &pzMissRot,
                                Double_t &pMissRot,
                                Double_t &thetaMissRot,
                                Double_t &phiMissRot)
{
  pxMissRot = pyMissRot = pzMissRot = -9.;
  pMissRot = thetaMissRot = phiMissRot = -9.;

  TVector3 yRaw = planeAxis.Cross(zAxis);
  if(vec.Mag() <= 0. || zAxis.Mag() <= 0. || planeAxis.Mag() <= 0. || yRaw.Mag() <= 0.){
    return;
  }

  TVector3 vz = zAxis.Unit();
  TVector3 vy = yRaw.Unit();
  TVector3 vx = vz.Cross(vy).Unit();
  TVector3 missRot(vec.Dot(vx), vec.Dot(vy), vec.Dot(vz));

  pxMissRot = missRot.X();
  pyMissRot = missRot.Y();
  pzMissRot = missRot.Z();
  pMissRot = missRot.Mag();
  thetaMissRot = missRot.Theta() * rad2deg;
  phiMissRot = missRot.Phi() * rad2deg;
}

int main(int argc, char **argv)
{
  if(argc < 4){
    Usage();
    return -1;
  }

  bool isMC = (atoi(argv[1]) == 1);
  double Ebeam = 5.98636;
  long long maxRows = -1;
  string topologyName = "deeppim";
  string targetOverride;
  string outputName = argv[2];
  vector<string> inputFiles;

  int firstInputArg = 3;
  if(argc >= 5 && isNumber(argv[2])){
    Ebeam = atof(argv[2]);
    outputName = argv[3];
    firstInputArg = 4;
  }

  for(int k = firstInputArg; k < argc; k++){
    string arg = argv[k];
    if(arg == "--beam"){
      if(k + 1 >= argc){ Usage(); return -1; }
      Ebeam = atof(argv[++k]);
    }
    else if(arg == "--topology"){
      if(k + 1 >= argc){ Usage(); return -1; }
      topologyName = argv[++k];
    }
    else if(arg == "--target"){
      if(k + 1 >= argc){ Usage(); return -1; }
      targetOverride = argv[++k];
    }
    else if(arg == "--max-rows"){
      if(k + 1 >= argc){ Usage(); return -1; }
      maxRows = atoll(argv[++k]);
    }
    else{
      inputFiles.push_back(arg);
    }
  }

  if(inputFiles.empty()){
    Usage();
    return -1;
  }

  TopologyConfig topology = makeTopology(topologyName);
  if(!targetOverride.empty()){
    applyTargetOverride(topology, targetOverride);
  }

  TFile *outFile = new TFile(outputName.c_str(), "RECREATE");
  TTree *tree = new TTree("deeppim", "Configurable e-scattering missing-system skim");

  Int_t b_run = -9;
  Long64_t b_event = -9;
  Bool_t b_isMC = isMC;
  Float_t b_weight = 1.f;
  Int_t b_topologyId = topology.id;
  Double_t b_targetMass = topology.targetMass;
  Int_t b_nElectrons = 0;
  Int_t b_nProtons = 0;
  Int_t b_nPiMinus = 0;
  Int_t b_p1Index = -1;
  Int_t b_p2Index = -1;
  Int_t b_pimIndex = -1;
  Int_t b_leadProtonIndex_rot = -1;

  Double_t b_Q2 = -9.;
  Double_t b_xB = -9.;
  Double_t b_omega = -9.;
  Double_t b_qP = -9.;
  Double_t b_qTheta = -9.;
  Double_t b_qPhi = -9.;
  Double_t b_qPx = -9.;
  Double_t b_qPy = -9.;
  Double_t b_qPz = -9.;

  Double_t b_mMiss = -9.;
  Double_t b_mMiss2 = -9.;
  Double_t b_pMiss = -9.;
  Double_t b_thetaMiss = -9.;
  Double_t b_phiMiss = -9.;
  Double_t b_pxMiss = -9.;
  Double_t b_pyMiss = -9.;
  Double_t b_pzMiss = -9.;
  Double_t b_pxMiss_rot = -9.;
  Double_t b_pyMiss_rot = -9.;
  Double_t b_pzMiss_rot = -9.;
  Double_t b_pMiss_rot = -9.;
  Double_t b_thetaMiss_rot = -9.;
  Double_t b_phiMiss_rot = -9.;

  ParticleBranches b_e;
  ParticleBranches b_pim;
  ParticleBranches b_p1;
  ParticleBranches b_p2;

  tree->Branch("run", &b_run, "run/I");
  tree->Branch("event", &b_event, "event/L");
  tree->Branch("isMC", &b_isMC, "isMC/O");
  tree->Branch("weight", &b_weight, "weight/F");
  tree->Branch("topologyId", &b_topologyId, "topologyId/I");
  tree->Branch("targetMass", &b_targetMass, "targetMass/D");
  tree->Branch("nElectrons", &b_nElectrons, "nElectrons/I");
  tree->Branch("nProtons", &b_nProtons, "nProtons/I");
  tree->Branch("nPiMinus", &b_nPiMinus, "nPiMinus/I");
  tree->Branch("p1Index", &b_p1Index, "p1Index/I");
  tree->Branch("p2Index", &b_p2Index, "p2Index/I");
  tree->Branch("pimIndex", &b_pimIndex, "pimIndex/I");
  tree->Branch("leadProtonIndex_rot", &b_leadProtonIndex_rot, "leadProtonIndex_rot/I");

  tree->Branch("Q2", &b_Q2, "Q2/D");
  tree->Branch("xB", &b_xB, "xB/D");
  tree->Branch("omega", &b_omega, "omega/D");
  tree->Branch("qP", &b_qP, "qP/D");
  tree->Branch("qTheta", &b_qTheta, "qTheta/D");
  tree->Branch("qPhi", &b_qPhi, "qPhi/D");
  tree->Branch("qPx", &b_qPx, "qPx/D");
  tree->Branch("qPy", &b_qPy, "qPy/D");
  tree->Branch("qPz", &b_qPz, "qPz/D");

  tree->Branch("mMiss", &b_mMiss, "mMiss/D");
  tree->Branch("mMiss2", &b_mMiss2, "mMiss2/D");
  tree->Branch("pMiss", &b_pMiss, "pMiss/D");
  tree->Branch("thetaMiss", &b_thetaMiss, "thetaMiss/D");
  tree->Branch("phiMiss", &b_phiMiss, "phiMiss/D");
  tree->Branch("pxMiss", &b_pxMiss, "pxMiss/D");
  tree->Branch("pyMiss", &b_pyMiss, "pyMiss/D");
  tree->Branch("pzMiss", &b_pzMiss, "pzMiss/D");
  tree->Branch("pxMiss_rot", &b_pxMiss_rot, "pxMiss_rot/D");
  tree->Branch("pyMiss_rot", &b_pyMiss_rot, "pyMiss_rot/D");
  tree->Branch("pzMiss_rot", &b_pzMiss_rot, "pzMiss_rot/D");
  tree->Branch("pMiss_rot", &b_pMiss_rot, "pMiss_rot/D");
  tree->Branch("thetaMiss_rot", &b_thetaMiss_rot, "thetaMiss_rot/D");
  tree->Branch("phiMiss_rot", &b_phiMiss_rot, "phiMiss_rot/D");

  b_e.book(tree, "e");
  b_pim.book(tree, "pim");
  b_p1.book(tree, "p1");
  b_p2.book(tree, "p2");

  clas12root::HipoChain chain;
  for(const auto &fname : inputFiles){
    cout << "Input file " << fname << endl;
    chain.Add(fname.c_str());
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  clas12ana clasAna;

  TLorentzVector beam(0., 0., Ebeam, Ebeam);
  TLorentzVector target(0., 0., 0., topology.targetMass);
  TLorentzVector el;

  long long nRead = 0;
  long long nOneElectron = 0;
  long long nEnoughProtons = 0;
  long long nEnoughPiMinus = 0;
  long long nGoodElectronKinematics = 0;
  long long nFDLeadKinematics = 0;
  long long nWritten = 0;
  

  while(chain.Next()){
    if(nRead % 100000 == 0){
      cout << "Processing event " << nRead << "\t" << nWritten << " rows saved." << "\t" << nFDLeadKinematics << " FD e-." << endl;
    }
    nRead++;
    if(maxRows > 0 && nWritten >= maxRows){ break; }
    if(nFDLeadKinematics > 10000) break;

    clasAna.Run(c12);
    auto electrons = clasAna.getByPid(11);
    auto protons = clasAna.getByPid(2212);
    auto piminus = clasAna.getByPid(-211);

    b_nElectrons = static_cast<int>(electrons.size());
    b_nProtons = static_cast<int>(protons.size());
    b_nPiMinus = static_cast<int>(piminus.size());

    if(electrons.size() != 1){ continue; }
    nOneElectron++;
    if(static_cast<int>(protons.size()) < topology.minProtons){ continue; }
    nEnoughProtons++;
    if(static_cast<int>(piminus.size()) < topology.minPiMinus){ continue; }
    nEnoughPiMinus++;

    setCorrectedP4(el, electrons[0], me, isMC);
    if(el.P() <= 0.){ continue; }

    TLorentzVector q = beam - el;
    b_Q2 = -q.M2();
    b_omega = q.E();
    if(b_omega <= 0.){ continue; }
    b_xB = b_Q2 / (2. * mP * b_omega);
    b_qP = q.P();
    b_qTheta = q.Theta() * rad2deg;
    b_qPhi = q.Phi() * rad2deg;
    b_qPx = q.Px();
    b_qPy = q.Py();
    b_qPz = q.Pz();

    if(b_xB < .7) continue;
   

    nGoodElectronKinematics++;



    b_run = c12->runconfig()->getRun();
    b_event = c12->runconfig()->getEvent();
    b_weight = isMC ? c12->mcevent()->getWeight() : 1.f;
    b_topologyId = topology.id;
    b_targetMass = topology.targetMass;
    b_e.fill(el, electrons[0]);

    vector<CorrectedParticle> corrProtons;
    corrProtons.reserve(protons.size());
    for(int i = 0; i < static_cast<int>(protons.size()); i++){
      CorrectedParticle cp;
      cp.index = i;
      cp.part = protons[i];
      setCorrectedP4(cp.p4, protons[i], mP, isMC);
      corrProtons.push_back(cp);
    }

    vector<CorrectedParticle> corrPiMinus;
    corrPiMinus.reserve(piminus.size());
    for(int i = 0; i < static_cast<int>(piminus.size()); i++){
      CorrectedParticle cp;
      cp.index = i;
      cp.part = piminus[i];
      setCorrectedP4(cp.p4, piminus[i], mPi, isMC);
      corrPiMinus.push_back(cp);
    }

    auto fillRow = [&](CorrectedParticle const &lead,
                       CorrectedParticle const *recoil,
                       CorrectedParticle const *pim)
    {
      b_p1.reset();
      b_p2.reset();
      b_pim.reset();
      b_p1Index = lead.index;
      b_p2Index = recoil ? recoil->index : -1;
      b_pimIndex = pim ? pim->index : -1;
      b_leadProtonIndex_rot = lead.index;

      b_p1.fill(lead.p4, lead.part);
      TLorentzVector miss = beam + target - el - lead.p4;

      if(recoil){
        b_p2.fill(recoil->p4, recoil->part);
        miss -= recoil->p4;
      }
      if(pim){
        b_pim.fill(pim->p4, pim->part);
        miss -= pim->p4;
      }

      TVector3 pMiss = lead.p4.Vect() - q.Vect();
      b_mMiss = miss.M();
      b_mMiss2 = miss.M2();
      b_pMiss = pMiss.Mag();
      b_thetaMiss = pMiss.Theta() * rad2deg;
      b_phiMiss = pMiss.Phi() * rad2deg;
      b_pxMiss = pMiss.X();
      b_pyMiss = pMiss.Y();
      b_pzMiss = pMiss.Z();

      TVector3 rotZ = (topology.name == "ep") ? q.Vect() : pMiss;
      TVector3 rotPlane = (topology.name == "ep") ? pMiss : q.Vect();
      fillRotatedMissingBranches(pMiss, rotZ, rotPlane,
                                 b_pxMiss_rot, b_pyMiss_rot, b_pzMiss_rot,
                                 b_pMiss_rot, b_thetaMiss_rot, b_phiMiss_rot);

      tree->Fill();

      if(lead.p4.Theta() * rad2deg < 37.){
        nFDLeadKinematics++;
      }
      nWritten++;
    };

    if(topology.useP2){
      for(int i = 0; i < static_cast<int>(corrProtons.size()); i++){
        for(int j = i + 1; j < static_cast<int>(corrProtons.size()); j++){
          CorrectedParticle const *lead = &corrProtons[i];
          CorrectedParticle const *recoil = &corrProtons[j];
          if(recoil->p4.P() > lead->p4.P()){
            swap(lead, recoil);
          }

          if(topology.usePiMinus){
            for(auto const &pim : corrPiMinus){
              fillRow(*lead, recoil, &pim);
            }
          }
          else{
            fillRow(*lead, recoil, nullptr);
          }
        }
      }
    }
    else{
      for(auto const &lead : corrProtons){
        if(topology.usePiMinus){
          for(auto const &pim : corrPiMinus){
            fillRow(lead, nullptr, &pim);
          }
        }
        else{
          fillRow(lead, nullptr, nullptr);
        }
      }
    }
  }

  outFile->cd();
  ostringstream meta;
  time_t now = time(nullptr);
  string dateString = ctime(&now);
  while(!dateString.empty() && dateString.back() == '\n'){ dateString.pop_back(); }
  meta << "date=" << dateString << "\n"
       << "topology=" << topology.name << "\n"
       << "topology_id=" << topology.id << "\n"
       << "target=" << topology.targetName << "\n"
       << "target_mass=" << topology.targetMass << "\n"
       << "cuts=exactly_one_electron,nProtons>=" << topology.minProtons
       << ",nPiMinus>=" << topology.minPiMinus << "\n"
       << "corrections=GetLorentzVector_Corrected from Corrections.h, matching the SRC_General correction prescription\n"
       << "pMiss_convention=p1_lead-q; mMiss is from the full missing four-vector for the selected topology\n"
       << "rotation=ep z_axis=q; otherwise z_axis=p1_lead-q; y_axis=plane_axis x z_axis; x_axis=z_axis x y_axis\n"
       << "beam_energy=" << Ebeam << "\n"
       << "QADB=OFF\n"
       << "isMC_arg=" << isMC << "\n"
       << "events_read=" << nRead << "\n"
       << "cutflow_one_electron=" << nOneElectron << "\n"
       << "cutflow_enough_protons=" << nEnoughProtons << "\n"
       << "cutflow_enough_pi_minus=" << nEnoughPiMinus << "\n"
       << "cutflow_good_electron_kinematics=" << nGoodElectronKinematics << "\n"
       << "rows_written=" << nWritten;
  TNamed metadata("skimmer_metadata", meta.str().c_str());
  metadata.Write();
  tree->Write();
  outFile->Close();

  cout << "Done. Processed " << nRead << " events. Saved " << nWritten << " rows.\n"
       << "Topology " << topology.name << " on " << topology.targetName << ". Cutflow: one electron "
       << nOneElectron << ", enough protons " << nEnoughProtons
       << ", enough pi- " << nEnoughPiMinus
       << ", electron kinematics " << nGoodElectronKinematics << ".\n";
  return 0;
}
