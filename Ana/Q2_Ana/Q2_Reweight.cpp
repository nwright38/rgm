// Builds a reusable Q2 reweighting map from data and simulation hipo files.
//
// The output ROOT file stores the raw data/sim Q2 spectra, normalized spectra,
// and h_Q2_reweight = normalized data / normalized simulation. Downstream
// analyses can load the file with Q2Reweight.h and multiply MC event weights
// by q2Weight.weight(Q2).
//
// Usage:
//   ./Q2_Reweight output.root simA simType [--selection ep|epp] [--max N] \
//       data1.hipo ... --sim sim1.hipo ...
//
// simType follows the existing convention in Main_Figs_Binned:
//   1 AV18, 2 AV4, 3 N2LO10, 4 N2LO12, 5 NV

#define MAIN_FIGS_BINNED_NO_MAIN
#include "Main_Figs_Binned.cpp"

#include <algorithm>
#include <stdexcept>

namespace {

enum class MatchSelection { EP, EPP };

void Usage() {
  cerr << "Usage:\n"
       << "  ./Q2_Reweight output.root simA simType [--selection ep|epp] [--max N] "
       << "data1.hipo ... --sim sim1.hipo ...\n\n"
       << "Example:\n"
       << "  ./Q2_Reweight q2_weights.root 4 1 data/*.hipo --sim sim/*.hipo\n";
}

Long64_t ParseLongLong(const string& s) {
  try {
    size_t pos = 0;
    Long64_t val = stoll(s, &pos);
    if (pos != s.size()) throw invalid_argument("extra characters");
    return val;
  } catch (const exception&) {
    throw runtime_error("Could not parse integer value: " + s);
  }
}

char* UTypeFromSimType(int simType) {
  if (simType == 1) return const_cast<char*>("AV18");
  if (simType == 2) return const_cast<char*>("AV4");
  if (simType == 3) return const_cast<char*>("N2LO10");
  if (simType == 4) return const_cast<char*>("N2LO12");
  if (simType == 5) return const_cast<char*>("NV");
  cerr << "Unknown simType " << simType << ", using AV18.\n";
  return const_cast<char*>("AV18");
}

void Normalize(TH1D* hist) {
  const double integral = hist->Integral(1, hist->GetNbinsX());
  if (integral > 0.0) hist->Scale(1.0 / integral);
}

bool PassesSelection(const EventKinematics& ek, MatchSelection selection) {
  if (selection == MatchSelection::EP) return ek.passep;
  return ek.passepp;
}

void FillFromHipo(const vector<string>& files, bool isMC, int nucleusA, char* uType,
                  MatchSelection selection, Long64_t maxEvents, TH1D* hist) {
  clas12ana clasAna;
  clasAna.printParams();

  clas12root::HipoChain chain;
  for (const auto& file : files) {
    cout << "Input " << (isMC ? "simulation" : "data") << " file " << file << endl;
    chain.Add(file);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12 = chain.GetC12Reader();
  auto& c12 = chain.C12ref();

  int Z = 2;
  int N = 2;
  if (isMC) {
    Z = nucleusA / 2;
    N = nucleusA / 2;
  }
  reweighter nominalWeight(beam_E, Z, N, kelly, uType, .15);
  CutVariation nominalCut = CutVariation::Nominal();

  int ctr = 0;
  Long64_t accepted = 0;
  while (chain.Next()) {
    if (maxEvents > 0 && ctr >= maxEvents) break;
    if (ctr % 1000 == 0) {
      cout << (isMC ? "Simulation" : "Data") << " event " << ctr << "\r" << flush;
    }

    EventKinematics ek = computeEventKinematics(c12, clasAna, false, ctr);
    bool passep = false;
    bool passepp = false;
    nominalCut.apply(ek, passep, passepp);
    ek.passep = passep;
    ek.passepp = passepp;
    if (!PassesSelection(ek, selection)) continue;

    double weight = 1.0;
    if (isMC) {
      const double originalWeight = c12->mcevent()->getWeight();
      if (selection == MatchSelection::EP) {
        weight = originalWeight * nominalWeight.get_weight_ep(c12->mcparts());
      } else {
        weight = originalWeight * nominalWeight.get_weight_epp(c12->mcparts());
      }
    }
    hist->Fill(ek.qSq, weight);
    accepted++;
  }
  cout << "\nAccepted " << accepted << " " << (isMC ? "simulation" : "data")
       << " events for Q2 matching.\n";
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 6) {
    Usage();
    return -1;
  }

  const string outFileName = argv[1];
  const int nucleusA = atoi(argv[2]);
  const int simType = atoi(argv[3]);
  char* uType = UTypeFromSimType(simType);

  MatchSelection selection = MatchSelection::EP;
  Long64_t maxEvents = -1;
  vector<string> dataFiles;
  vector<string> simFiles;
  bool readingSim = false;

  for (int i = 4; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "--sim") {
      readingSim = true;
      continue;
    }
    if (arg == "--selection") {
      if (i + 1 >= argc) {
        Usage();
        return -1;
      }
      string val = argv[++i];
      if (val == "ep") selection = MatchSelection::EP;
      else if (val == "epp") selection = MatchSelection::EPP;
      else {
        cerr << "Unknown selection '" << val << "'. Use ep or epp.\n";
        return -1;
      }
      continue;
    }
    if (arg.rfind("--selection=", 0) == 0) {
      string val = arg.substr(12);
      if (val == "ep") selection = MatchSelection::EP;
      else if (val == "epp") selection = MatchSelection::EPP;
      else {
        cerr << "Unknown selection '" << val << "'. Use ep or epp.\n";
        return -1;
      }
      continue;
    }
    if (arg == "--max") {
      if (i + 1 >= argc) {
        Usage();
        return -1;
      }
      try {
        maxEvents = ParseLongLong(argv[++i]);
      } catch (const exception& e) {
        cerr << e.what() << endl;
        return -1;
      }
      continue;
    }
    if (arg.rfind("--max=", 0) == 0) {
      try {
        maxEvents = ParseLongLong(arg.substr(6));
      } catch (const exception& e) {
        cerr << e.what() << endl;
        return -1;
      }
      continue;
    }

    if (readingSim) simFiles.push_back(arg);
    else dataFiles.push_back(arg);
  }

  if (dataFiles.empty() || simFiles.empty()) {
    cerr << "Need at least one data hipo file and one simulation hipo file.\n";
    Usage();
    return -1;
  }

  TH1::AddDirectory(kFALSE);

  TH1D hDataRaw("h_Q2_data_raw", "Data Q^{2};Q^{2};Counts", bE_Q2.size() - 1, bE_Q2.data());
  TH1D hSimRaw("h_Q2_sim_raw", "Simulation Q^{2};Q^{2};Counts", bE_Q2.size() - 1, bE_Q2.data());
  hDataRaw.Sumw2();
  hSimRaw.Sumw2();

  FillFromHipo(dataFiles, false, nucleusA, uType, selection, maxEvents, &hDataRaw);
  FillFromHipo(simFiles, true, nucleusA, uType, selection, maxEvents, &hSimRaw);

  TH1D hDataNorm(hDataRaw);
  TH1D hSimNorm(hSimRaw);
  hDataNorm.SetName("h_Q2_data_norm");
  hSimNorm.SetName("h_Q2_sim_norm");
  Normalize(&hDataNorm);
  Normalize(&hSimNorm);

  TH1D hWeights("h_Q2_reweight", "Q^{2} reweight;Q^{2};Data/Simulation",
                bE_Q2.size() - 1, bE_Q2.data());
  hWeights.Sumw2();
  for (int bin = 1; bin <= hWeights.GetNbinsX(); ++bin) {
    const double data = hDataNorm.GetBinContent(bin);
    const double sim = hSimNorm.GetBinContent(bin);
    hWeights.SetBinContent(bin, sim > 0.0 ? data / sim : 0.0);
    hWeights.SetBinError(bin, 0.0);
  }

  TH1D hSimReweighted(hSimRaw);
  hSimReweighted.SetName("h_Q2_sim_reweighted");
  hSimReweighted.SetTitle("Reweighted simulation Q^{2};Q^{2};Counts");
  for (int bin = 1; bin <= hSimReweighted.GetNbinsX(); ++bin) {
    hSimReweighted.SetBinContent(bin, hSimRaw.GetBinContent(bin) * hWeights.GetBinContent(bin));
    hSimReweighted.SetBinError(bin, hSimRaw.GetBinError(bin) * hWeights.GetBinContent(bin));
  }

  TFile outFile(outFileName.c_str(), "RECREATE");
  if (outFile.IsZombie()) {
    cerr << "Could not create output file " << outFileName << endl;
    return -1;
  }
  TDirectory* dir = outFile.mkdir("q2_reweighting");
  dir->cd();
  hDataRaw.Write();
  hSimRaw.Write();
  hDataNorm.Write();
  hSimNorm.Write();
  hWeights.Write();
  hSimReweighted.Write();
  outFile.Close();

  cout << "Wrote Q2 reweighting file " << outFileName << endl;
  cout << "Selection: " << (selection == MatchSelection::EP ? "ep" : "epp") << endl;
  cout << "Downstream histogram: q2_reweighting/h_Q2_reweight" << endl;

  return 0;
}
