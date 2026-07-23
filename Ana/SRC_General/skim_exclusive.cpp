#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <set>
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
const double mN = 0.9395654133;
const double mD = 1.8756;
const double mPi = 0.13957039;
const double mK = 0.493677;
const double me = 0.000511;
const double rad2deg = 180. / M_PI;

// ------------------------------------------------------------------
// Manual-edit block for new exclusive topologies.
// ------------------------------------------------------------------
namespace ExclusiveConfig {
  // Examples:
  //   H(e,e'p)         : {2212}, targetMass = mP
  //   d(e,e'pp pi-)   : {2212, 2212, -211}, targetMass = mD
  //   d(e,e'p pi-)    : {2212, -211}, targetMass = mD
  const vector<int> topologyPids = {2212, 2212, -211};
  const double targetMass = mD;
  const string topologyName = "eppim_pp";
  const bool requireExactlyOneElectron = true;
  const bool applyCorrections = true;

  TLorentzVector missingP4(TLorentzVector const &beam,
                           TLorentzVector const &target,
                           TLorentzVector const &electron,
                           vector<TLorentzVector> const &selected)
  {
    TLorentzVector miss = beam + target - electron;
    for(auto const &p4 : selected){
      miss -= p4;
    }
    return miss;
  }

  // Default SRC-style convention. Edit here if a topology needs something else.
  TVector3 pMiss3(TLorentzVector const &q, vector<TLorentzVector> const &selected)
  {
    if(selected.empty()){
      return TVector3(-9., -9., -9.);
    }
    return selected[0].Vect() - q.Vect();
  }

  // Default: rotate pMiss into z, with q setting the reaction plane.
  TVector3 rotationZAxis(TLorentzVector const &q, TVector3 const &pMiss, vector<TLorentzVector> const &selected)
  {
    (void)q;
    (void)selected;
    return pMiss;
  }

  TVector3 rotationPlaneAxis(TLorentzVector const &q, TVector3 const &pMiss, vector<TLorentzVector> const &selected)
  {
    (void)pMiss;
    (void)selected;
    return q.Vect();
  }
}

struct ParticleInfo {
  int index = -1;
  int pid = 0;
  int charge = 0;
  int status = 0;
  int region = 0;
  double beta = -9.;
  double vz = -99.;
  TLorentzVector p4;
};

void Usage()
{
  cerr << "Usage: ./skim_exclusive <MC=1,Data=0> [Ebeam(GeV)] <output.root> <input.hipo> [more hipo...] [--beam Ebeam] [--max-rows N]\n";
}

bool isNumber(const string &text)
{
  char *end = nullptr;
  strtod(text.c_str(), &end);
  return end != text.c_str() && end != nullptr && *end == '\0';
}

double massForPid(int pid)
{
  switch(pid){
    case 11:
    case -11:
      return me;
    case 2212:
      return mP;
    case 2112:
      return mN;
    case 45:
      return mD;
    case 211:
    case -211:
      return mPi;
    case 321:
    case -321:
      return mK;
    case 22:
    case 0:
      return 0.;
    default:
      return 0.;
  }
}

TLorentzVector correctedP4(clas12::region_part_ptr part, bool isMC)
{
  int pid = part->par()->getPid();
  TLorentzVector p4(0., 0., 0., massForPid(pid));
  GetLorentzVector_ReconVector(p4, part);

  if(!ExclusiveConfig::applyCorrections){
    return p4;
  }

  if(part->par()->getCharge() != 0 && (part->getRegion() == FD || part->getRegion() == CD)){
    GetLorentzVector_Corrected(p4, part, isMC);
  }
  return p4;
}

void fillRotated(TVector3 const &vec,
                 TVector3 const &zAxis,
                 TVector3 const &planeAxis,
                 double &xRot,
                 double &yRot,
                 double &zRot,
                 double &pRot,
                 double &thetaRot,
                 double &phiRot)
{
  xRot = yRot = zRot = -9.;
  pRot = thetaRot = phiRot = -9.;

  TVector3 yRaw = planeAxis.Cross(zAxis);
  if(vec.Mag() <= 0. || zAxis.Mag() <= 0. || planeAxis.Mag() <= 0. || yRaw.Mag() <= 0.){
    return;
  }

  TVector3 vz = zAxis.Unit();
  TVector3 vy = yRaw.Unit();
  TVector3 vx = vz.Cross(vy).Unit();
  TVector3 rotated(vec.Dot(vx), vec.Dot(vy), vec.Dot(vz));

  xRot = rotated.X();
  yRot = rotated.Y();
  zRot = rotated.Z();
  pRot = rotated.Mag();
  thetaRot = rotated.Theta() * rad2deg;
  phiRot = rotated.Phi() * rad2deg;
}

void findCombinations(vector<ParticleInfo> const &particles,
                      vector<int> const &topologyPids,
                      int slot,
                      int electronIndex,
                      vector<int> &selected,
                      vector<vector<int>> &out)
{
  if(slot == static_cast<int>(topologyPids.size())){
    out.push_back(selected);
    return;
  }

  int requiredPid = topologyPids[slot];
  int minIndex = -1;
  if(slot > 0 && topologyPids[slot - 1] == requiredPid){
    minIndex = selected.back();
  }

  for(int i = 0; i < static_cast<int>(particles.size()); i++){
    if(i == electronIndex){ continue; }
    if(particles[i].pid != requiredPid){ continue; }
    if(i <= minIndex){ continue; }
    if(find(selected.begin(), selected.end(), i) != selected.end()){ continue; }
    selected.push_back(i);
    findCombinations(particles, topologyPids, slot + 1, electronIndex, selected, out);
    selected.pop_back();
  }
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

  TFile *outFile = new TFile(outputName.c_str(), "RECREATE");
  TTree *tree = new TTree("exclusive", "Generic exclusive-topology skim");

  int run = -9;
  long long event = -9;
  bool b_isMC = isMC;
  float weight = 1.f;
  int nParticles = 0;
  int nTopologyCombos = 0;

  double Q2 = -9.;
  double xB = -9.;
  double omega = -9.;
  double mMiss = -9.;
  double mMiss2 = -9.;
  double pMiss = -9.;
  double thetaMiss = -9.;
  double phiMiss = -9.;
  double pxMiss = -9.;
  double pyMiss = -9.;
  double pzMiss = -9.;
  double pxMiss_rot = -9.;
  double pyMiss_rot = -9.;
  double pzMiss_rot = -9.;
  double pMiss_rot = -9.;
  double thetaMiss_rot = -9.;
  double phiMiss_rot = -9.;

  vector<int> pid;
  vector<int> charge;
  vector<int> status;
  vector<int> region;
  vector<double> p;
  vector<double> theta;
  vector<double> phi;
  vector<double> px;
  vector<double> py;
  vector<double> pz;
  vector<double> energy;
  vector<double> beta;
  vector<double> vz;
  vector<int> topologyParticleIndex;
  vector<int> topologyPid;

  tree->Branch("run", &run);
  tree->Branch("event", &event);
  tree->Branch("isMC", &b_isMC);
  tree->Branch("weight", &weight);
  tree->Branch("nParticles", &nParticles);
  tree->Branch("nTopologyCombos", &nTopologyCombos);
  tree->Branch("Q2", &Q2);
  tree->Branch("xB", &xB);
  tree->Branch("omega", &omega);
  tree->Branch("mMiss", &mMiss);
  tree->Branch("mMiss2", &mMiss2);
  tree->Branch("pMiss", &pMiss);
  tree->Branch("thetaMiss", &thetaMiss);
  tree->Branch("phiMiss", &phiMiss);
  tree->Branch("pxMiss", &pxMiss);
  tree->Branch("pyMiss", &pyMiss);
  tree->Branch("pzMiss", &pzMiss);
  tree->Branch("pxMiss_rot", &pxMiss_rot);
  tree->Branch("pyMiss_rot", &pyMiss_rot);
  tree->Branch("pzMiss_rot", &pzMiss_rot);
  tree->Branch("pMiss_rot", &pMiss_rot);
  tree->Branch("thetaMiss_rot", &thetaMiss_rot);
  tree->Branch("phiMiss_rot", &phiMiss_rot);
  tree->Branch("pid", &pid);
  tree->Branch("charge", &charge);
  tree->Branch("status", &status);
  tree->Branch("region", &region);
  tree->Branch("p", &p);
  tree->Branch("theta", &theta);
  tree->Branch("phi", &phi);
  tree->Branch("px", &px);
  tree->Branch("py", &py);
  tree->Branch("pz", &pz);
  tree->Branch("energy", &energy);
  tree->Branch("beta", &beta);
  tree->Branch("vz", &vz);
  tree->Branch("topologyParticleIndex", &topologyParticleIndex);
  tree->Branch("topologyPid", &topologyPid);

  clas12root::HipoChain chain;
  for(auto const &fname : inputFiles){
    cout << "Input file " << fname << endl;
    chain.Add(fname.c_str());
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12 = chain.GetC12Reader();
  auto &c12 = chain.C12ref();

  TLorentzVector beam(0., 0., Ebeam, Ebeam);
  TLorentzVector target(0., 0., 0., ExclusiveConfig::targetMass);

  long long nRead = 0;
  long long nWritten = 0;
  long long nOneElectron = 0;
  long long nWithTopology = 0;

  while(chain.Next()){
    if(nRead % 100000 == 0){
      cout << "Processing event " << nRead << "\t" << nWritten << " rows saved." << endl;
    }
    nRead++;
    if(maxRows > 0 && nWritten >= maxRows){ break; }

    auto detParticles = c12->getDetParticles();
    vector<ParticleInfo> particles;
    particles.reserve(detParticles.size());

    pid.clear();
    charge.clear();
    status.clear();
    region.clear();
    p.clear();
    theta.clear();
    phi.clear();
    px.clear();
    py.clear();
    pz.clear();
    energy.clear();
    beta.clear();
    vz.clear();

    int electronIndex = -1;
    for(int i = 0; i < static_cast<int>(detParticles.size()); i++){
      auto part = detParticles[i];
      ParticleInfo info;
      info.index = i;
      info.pid = part->par()->getPid();
      info.charge = part->par()->getCharge();
      info.status = part->getStatus();
      info.region = static_cast<int>(part->getRegion());
      info.beta = part->par()->getBeta();
      info.vz = part->par()->getVz();
      info.p4 = correctedP4(part, isMC);
      particles.push_back(info);

      if(info.pid == 11){
        if(electronIndex < 0){
          electronIndex = i;
        }
        else if(ExclusiveConfig::requireExactlyOneElectron){
          electronIndex = -2;
        }
      }

      pid.push_back(info.pid);
      charge.push_back(info.charge);
      status.push_back(info.status);
      region.push_back(info.region);
      p.push_back(info.p4.P());
      theta.push_back(info.p4.Theta() * rad2deg);
      phi.push_back(info.p4.Phi() * rad2deg);
      px.push_back(info.p4.Px());
      py.push_back(info.p4.Py());
      pz.push_back(info.p4.Pz());
      energy.push_back(info.p4.E());
      beta.push_back(info.beta);
      vz.push_back(info.vz);
    }

    nParticles = static_cast<int>(particles.size());
    if(electronIndex < 0){ continue; }
    nOneElectron++;

    TLorentzVector electron = particles[electronIndex].p4;
    TLorentzVector q = beam - electron;
    Q2 = -q.M2();
    omega = q.E();
    if(omega <= 0.){ continue; }
    xB = Q2 / (2. * mP * omega);

    vector<int> selected;
    vector<vector<int>> combinations;
    findCombinations(particles, ExclusiveConfig::topologyPids, 0, electronIndex, selected, combinations);
    nTopologyCombos = static_cast<int>(combinations.size());
    if(combinations.empty()){ continue; }
    nWithTopology++;

    for(auto const &combo : combinations){
      vector<TLorentzVector> selectedP4;
      topologyParticleIndex.clear();
      topologyPid.clear();

      for(int idx : combo){
        selectedP4.push_back(particles[idx].p4);
        topologyParticleIndex.push_back(particles[idx].index);
        topologyPid.push_back(particles[idx].pid);
      }

      TLorentzVector miss = ExclusiveConfig::missingP4(beam, target, electron, selectedP4);
      TVector3 pmiss3 = ExclusiveConfig::pMiss3(q, selectedP4);
      TVector3 zAxis = ExclusiveConfig::rotationZAxis(q, pmiss3, selectedP4);
      TVector3 planeAxis = ExclusiveConfig::rotationPlaneAxis(q, pmiss3, selectedP4);

      run = c12->runconfig()->getRun();
      event = c12->runconfig()->getEvent();
      weight = isMC ? c12->mcevent()->getWeight() : 1.f;

      mMiss = miss.M();
      mMiss2 = miss.M2();
      pMiss = pmiss3.Mag();
      thetaMiss = pmiss3.Theta() * rad2deg;
      phiMiss = pmiss3.Phi() * rad2deg;
      pxMiss = pmiss3.X();
      pyMiss = pmiss3.Y();
      pzMiss = pmiss3.Z();

      fillRotated(pmiss3, zAxis, planeAxis,
                  pxMiss_rot, pyMiss_rot, pzMiss_rot,
                  pMiss_rot, thetaMiss_rot, phiMiss_rot);

      tree->Fill();
      nWritten++;
      if(maxRows > 0 && nWritten >= maxRows){ break; }
    }
  }

  outFile->cd();
  ostringstream meta;
  time_t now = time(nullptr);
  string dateString = ctime(&now);
  while(!dateString.empty() && dateString.back() == '\n'){ dateString.pop_back(); }
  meta << "date=" << dateString << "\n"
       << "topology_name=" << ExclusiveConfig::topologyName << "\n"
       << "topology_pids=";
  for(size_t i = 0; i < ExclusiveConfig::topologyPids.size(); i++){
    if(i > 0){ meta << ","; }
    meta << ExclusiveConfig::topologyPids[i];
  }
  meta << "\n"
       << "target_mass=" << ExclusiveConfig::targetMass << "\n"
       << "beam_energy=" << Ebeam << "\n"
       << "corrections=" << (ExclusiveConfig::applyCorrections ? "on" : "off") << "\n"
       << "events_read=" << nRead << "\n"
       << "events_with_one_electron=" << nOneElectron << "\n"
       << "events_with_topology=" << nWithTopology << "\n"
       << "rows_written=" << nWritten;

  TNamed metadata("skimmer_metadata", meta.str().c_str());
  metadata.Write();
  tree->Write();
  outFile->Close();

  cout << "Done. Processed " << nRead << " events. Saved " << nWritten << " rows.\n";
  return 0;
}
