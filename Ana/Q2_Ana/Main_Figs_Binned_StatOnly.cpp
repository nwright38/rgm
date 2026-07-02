// Fast nominal-only companion to Main_Figs_Binned.cpp.
//
// Produces the same diffTable / integratedTable / hists / graphs layout that
// Main_Figs_Binned.cpp followed by BuildGraphs.cpp produces, but skips the
// 100 cut/reweighting toys. This is intended for plots where only statistical
// uncertainties should be shown.
//
// Usage:
//   ./Main_Figs_Binned_StatOnly isMC A outputfile.root [--mode legacy|modern] inputfiles.hipo

#define MAIN_FIGS_BINNED_NO_MAIN
#include "Main_Figs_Binned.cpp"

#define BUILD_GRAPHS_NO_MAIN
#include "BuildGraphs.cpp"

void StatOnlyUsage() {
  std::cerr << "Usage: ./Main_Figs_Binned_StatOnly isMC A outputfile.root "
            << "[--mode legacy|modern] [--q2-reweight weights.root] inputfiles.hipo \n\n\n";
}

int main(int argc, char** argv) {
  if (argc < 4) {
    StatOnlyUsage();
    return -1;
  }

  bool legacyCompatMode = false;
  string q2ReweightFile;
  int inputStartArg = 4;
  while (inputStartArg < argc) {
    string opt = argv[inputStartArg];
    if (opt.rfind("--mode", 0) == 0) {
      string mode;
      if (opt == "--mode" && inputStartArg + 1 < argc) {
        mode = argv[inputStartArg + 1];
        inputStartArg += 2;
      } else if (opt.rfind("--mode=", 0) == 0) {
        mode = opt.substr(7);
        inputStartArg += 1;
      } else {
        StatOnlyUsage();
        return -1;
      }

      if (mode == "legacy") legacyCompatMode = true;
      else if (mode == "modern") legacyCompatMode = false;
      else {
        std::cerr << "Unknown mode '" << mode << "'. Use 'legacy' or 'modern'.\n";
        return -1;
      }
      continue;
    }
    if (opt == "--q2-reweight") {
      if (inputStartArg + 1 >= argc) {
        StatOnlyUsage();
        return -1;
      }
      q2ReweightFile = argv[inputStartArg + 1];
      inputStartArg += 2;
      continue;
    }
    if (opt.rfind("--q2-reweight=", 0) == 0) {
      q2ReweightFile = opt.substr(14);
      inputStartArg += 1;
      continue;
    }
    break;
  }

  TH1::AddDirectory(kFALSE);

  int isMC = atoi(argv[1]);
  char* uType = "AV18";
  if (isMC > 0) {
    if (isMC == 1) uType = "AV18";
    if (isMC == 2) uType = "AV4";
    if (isMC == 3) uType = "N2LO10";
    if (isMC == 4) uType = "N2LO12";
    if (isMC == 5) uType = "NV";
    isMC = 1;
  }

  int nucleus_A = atoi(argv[2]);
  TString outFile = argv[3];
  cout << "Output file " << outFile << endl;

  clas12ana clasAna;
  clasAna.printParams();

  clas12root::HipoChain chain;
  for (int k = inputStartArg; k < argc; k++) {
    cout << "Input file " << argv[k] << endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12 = chain.GetC12Reader();

  auto& c12 = chain.C12ref();

  int Z = 2, N = 2;
  if (isMC) {
    Z = nucleus_A / 2;
    N = nucleus_A / 2;
  }
  reweighter newWeight(beam_E, Z, N, kelly, uType, .15);
  Q2Reweight q2Reweight;
  if (!q2ReweightFile.empty()) {
    if (!isMC) {
      cout << "Ignoring --q2-reweight for data input." << endl;
    } else if (!q2Reweight.load(q2ReweightFile)) {
      return -1;
    } else {
      cout << "Loaded Q2 reweight file " << q2ReweightFile << endl;
    }
  }

  vector<FillTask<EventKinematics>> tasks = buildFillTasks(legacyCompatMode);

  HistStore<EventKinematics> nominalStore("nominal");
  nominalStore.book(tasks);
  vector<HistStore<EventKinematics>> noToyStores;

  CutVariation nominalCut = CutVariation::Nominal();

  int ctr = 0;
  int max_ev = 1000000;
  while (chain.Next() && ctr < max_ev) {
    if (ctr % 1000 == 0) {
      cout << "Event " << ctr << " of " << max_ev << ". " << '\t'
           << (double)ctr / max_ev * 100 << "%\r" << flush;
    }

    double wep = 1, wepp = 1, original_weight = 1;
    if (isMC) {
      original_weight = c12->mcevent()->getWeight();
      wep = original_weight * newWeight.get_weight_ep(c12->mcparts());
      wepp = original_weight * newWeight.get_weight_epp(c12->mcparts());
    }

    EventKinematics ek = computeEventKinematics(c12, clasAna, false, ctr);
    if (isMC && q2Reweight.enabled()) {
      const double q2Weight = q2Reweight.weight(ek.qSq);
      wep *= q2Weight;
      wepp *= q2Weight;
    }

    bool passep_nom, passepp_nom;
    nominalCut.apply(ek, passep_nom, passepp_nom);
    ek.passep = passep_nom;
    ek.passepp = passepp_nom;
    nominalStore.fill(tasks, ek, wep, wepp);
  }

  TFile* f = new TFile(outFile, "RECREATE");
  f->cd();

  TTree* diffTree = new TTree("diffTable", "Differential bin table");
  string d_task, d_sel;
  vector<int> d_axisBin;
  vector<double> d_axisCenter;
  int d_valueBin;
  double d_valueCenter, d_count, d_statErr, d_sysErr, d_sysErrUp, d_sysErrDown;
  diffTree->Branch("task_name", &d_task);
  diffTree->Branch("selection", &d_sel);
  diffTree->Branch("axis_bin", &d_axisBin);
  diffTree->Branch("axis_center", &d_axisCenter);
  diffTree->Branch("value_bin", &d_valueBin, "value_bin/I");
  diffTree->Branch("value_center", &d_valueCenter, "value_center/D");
  diffTree->Branch("count", &d_count, "count/D");
  diffTree->Branch("stat_error", &d_statErr, "stat_error/D");
  diffTree->Branch("sys_error", &d_sysErr, "sys_error/D");
  diffTree->Branch("sys_error_up", &d_sysErrUp, "sys_error_up/D");
  diffTree->Branch("sys_error_down", &d_sysErrDown, "sys_error_down/D");

  TTree* intTree = new TTree("integratedTable", "Integrated (collapsed-axis) table");
  string i_task, i_sel, i_pattern;
  vector<int> i_axisBin;
  vector<double> i_axisCenter;
  double i_count, i_statErr, i_sysErr, i_sysErrUp, i_sysErrDown;
  intTree->Branch("task_name", &i_task);
  intTree->Branch("selection", &i_sel);
  intTree->Branch("pattern", &i_pattern);
  intTree->Branch("axis_bin", &i_axisBin);
  intTree->Branch("axis_center", &i_axisCenter);
  intTree->Branch("count", &i_count, "count/D");
  intTree->Branch("stat_error", &i_statErr, "stat_error/D");
  intTree->Branch("sys_error", &i_sysErr, "sys_error/D");
  intTree->Branch("sys_error_up", &i_sysErrUp, "sys_error_up/D");
  intTree->Branch("sys_error_down", &i_sysErrDown, "sys_error_down/D");

  for (size_t t = 0; t < tasks.size(); t++) {
    auto diffRows = buildDiffRows(tasks[t], (int)t, nominalStore, noToyStores,
                                  CentralValueMode::NOMINAL);
    for (auto& r : diffRows) {
      d_task = r.task_name;
      d_sel = r.selection;
      d_axisBin = r.axis_bin;
      d_axisCenter = r.axis_center;
      d_valueBin = r.value_bin;
      d_valueCenter = r.value_center;
      d_count = r.count;
      d_statErr = r.stat_error;
      d_sysErr = r.sys_error;
      d_sysErrUp = r.sys_error_up;
      d_sysErrDown = r.sys_error_down;
      diffTree->Fill();
    }

    vector<pair<vector<int>, string>> collapsePatterns;
    collapsePatterns.push_back({{}, "collapse:none"});
    for (size_t a = 0; a < tasks[t].binAxes.size(); a++) {
      collapsePatterns.push_back({{(int)a}, "collapse:" + tasks[t].binAxes[a]->name});
    }
    if (tasks[t].binAxes.size() > 1) {
      vector<int> allAxes;
      for (size_t a = 0; a < tasks[t].binAxes.size(); a++) allAxes.push_back((int)a);
      collapsePatterns.push_back({allAxes, "collapse:all"});
    }

    for (auto& pat : collapsePatterns) {
      auto intRows = buildIntegratedRows(tasks[t], (int)t, nominalStore, noToyStores,
                                          pat.first, pat.second);
      for (auto& r : intRows) {
        i_task = r.task_name;
        i_sel = r.selection;
        i_pattern = r.pattern;
        i_axisBin = r.axis_bin;
        i_axisCenter = r.axis_center;
        i_count = r.count;
        i_statErr = r.stat_error;
        i_sysErr = r.sys_error;
        i_sysErrUp = r.sys_error_up;
        i_sysErrDown = r.sys_error_down;
        intTree->Fill();
      }
    }
  }

  auto findTaskIdx = [&](const string& name) -> int {
    for (size_t t = 0; t < tasks.size(); t++) if (tasks[t].name == name) return (int)t;
    return -1;
  };
  vector<pair<string, string>> ratioSpecs = {
      {"pMiss_epp", "pMiss_ep"},
      {"kMiss_epp", "kMiss_ep"},
      {"Q2_epp_SRC_pmiss", "Q2_ep_SRC_pmiss"},
      {"Q2_epp_SRC_kmiss", "Q2_ep_SRC_kmiss"},
  };
  for (auto& spec : ratioSpecs) {
    int numIdx = findTaskIdx(spec.first);
    int denIdx = findTaskIdx(spec.second);
    if (numIdx < 0 || denIdx < 0) continue;
    string ratioName = spec.first + "_over_" + spec.second;
    auto rows = buildRatioDiffRows(tasks[numIdx], numIdx, denIdx, ratioName,
                                   nominalStore, noToyStores,
                                   CentralValueMode::NOMINAL);
    for (auto& r : rows) {
      d_task = r.task_name;
      d_sel = r.selection;
      d_axisBin = r.axis_bin;
      d_axisCenter = r.axis_center;
      d_valueBin = r.value_bin;
      d_valueCenter = r.value_center;
      d_count = r.count;
      d_statErr = r.stat_error;
      d_sysErr = r.sys_error;
      d_sysErrUp = r.sys_error_up;
      d_sysErrDown = r.sys_error_down;
      diffTree->Fill();
    }
  }

  intTree->Write();
  diffTree->Write();

  TDirectory* histDir = f->mkdir("hists");
  histDir->cd();
  nominalStore.writeAll();

  f->cd();
  TDirectory* graphsDir = f->mkdir("graphs");
  buildDiffGraphs(diffTree, graphsDir);
  buildIntegratedGraphs(intTree, graphsDir);

  f->cd();
  f->Write("", TObject::kOverwrite);
  f->Close();
  return 0;
}
