// Rehydrates the flat diffTable / integratedTable TTrees written by
// Main_Figs_Binned.cpp into TGraphAsymmErrors objects, so the python plotting
// layer can read graphs by name exactly like the old pipeline did.
//
// Deliberately plain/procedural (loops + std::map, no templates, no
// std::function): the input tables are already a fully generic flat
// format, so this step doesn't need the genericity machinery in
// BinnedHistStore.h -- it just groups rows and builds graphs.
//
// Usage: ./BuildGraphs <file.root>
// Opens <file.root> in UPDATE mode and writes everything into a new
// "graphs" TDirectory inside it (parallel to the "hists" directory
// Main_Figs_Binned.cpp already writes).
//
// Error convention: each graph's y-errors are
//   total_up   = sqrt(stat_error^2 + sys_error_up^2)
//   total_down = sqrt(stat_error^2 + sys_error_down^2)
// If the tree only has the legacy symmetric sys_error branch, BuildGraphs
// falls back to sys_error_up = sys_error_down = sys_error.

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphAsymmErrors.h>
#include <TDirectory.h>

using namespace std;

// Joins a vector<int> like {2,1} into "2_1", or "none" if empty -- used to
// build a unique, readable suffix for each selector-bin combination.
string joinBins(const vector<int>& bins) {
  if (bins.empty()) return "none";
  string s;
  for (size_t i = 0; i < bins.size(); i++) {
    if (i > 0) s += "_";
    s += to_string(bins[i]);
  }
  return s;
}

// One point collected from a diffTable row (already grouped by task/selection/axis_bin).
struct DiffPoint {
  int value_bin;
  double value_center, count, stat_error, sys_error_up, sys_error_down;
};

// One point collected from an integratedTable row (already grouped by task/selection/pattern).
struct IntegratedPoint {
  vector<int> axis_bin;
  vector<double> axis_center;
  double count, stat_error, sys_error_up, sys_error_down;
};

void buildDiffGraphs(TTree* diffTree, TDirectory* outDir) {
  string task, sel;
  vector<int>* axisBin = nullptr;
  int valueBin = 0;
  double valueCenter = 0.0, count = 0.0, statErr = 0.0, sysErr = 0.0;
  double sysErrUp = 0.0, sysErrDown = 0.0;
  const bool hasAsymm = diffTree->GetBranch("sys_error_up") && diffTree->GetBranch("sys_error_down");

  diffTree->SetBranchAddress("task_name", &task);
  diffTree->SetBranchAddress("selection", &sel);
  diffTree->SetBranchAddress("axis_bin", &axisBin);
  diffTree->SetBranchAddress("value_bin", &valueBin);
  diffTree->SetBranchAddress("value_center", &valueCenter);
  diffTree->SetBranchAddress("count", &count);
  diffTree->SetBranchAddress("stat_error", &statErr);
  diffTree->SetBranchAddress("sys_error", &sysErr);
  if (hasAsymm) {
    diffTree->SetBranchAddress("sys_error_up", &sysErrUp);
    diffTree->SetBranchAddress("sys_error_down", &sysErrDown);
  }

  // group key -> all points for that (task, selection, selector-bin) combo
  map<string, vector<DiffPoint>> groups;

  const Long64_t nEntries = diffTree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    diffTree->GetEntry(i);
    string key = task + "|" + sel + "|" + joinBins(*axisBin);
    DiffPoint p;
    p.value_bin = valueBin;
    p.value_center = valueCenter;
    p.count = count;
    p.stat_error = statErr;
    p.sys_error_up = hasAsymm ? sysErrUp : sysErr;
    p.sys_error_down = hasAsymm ? sysErrDown : sysErr;
    groups[key].push_back(p);
  }

  outDir->cd();
  for (auto& kv : groups) {
    vector<DiffPoint> pts = kv.second;
    sort(pts.begin(), pts.end(),
         [](const DiffPoint& a, const DiffPoint& b) { return a.value_bin < b.value_bin; });

    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    g->SetName(kv.first.c_str());
    for (auto& p : pts) {
      double errUp = sqrt(p.stat_error * p.stat_error + p.sys_error_up * p.sys_error_up);
      double errDown = sqrt(p.stat_error * p.stat_error + p.sys_error_down * p.sys_error_down);
      int n = g->GetN();
      g->SetPoint(n, p.value_center, p.count);
      g->SetPointError(n, 0.0, 0.0, errDown, errUp);
    }
    g->Write();
  }
  cout << "Wrote " << groups.size() << " differential graphs." << endl;
}

void buildIntegratedGraphs(TTree* intTree, TDirectory* outDir) {
  string task, sel, pattern;
  vector<int>* axisBin = nullptr;
  vector<double>* axisCenter = nullptr;
  double count = 0.0, statErr = 0.0, sysErr = 0.0;
  double sysErrUp = 0.0, sysErrDown = 0.0;
  const bool hasAsymm = intTree->GetBranch("sys_error_up") && intTree->GetBranch("sys_error_down");

  intTree->SetBranchAddress("task_name", &task);
  intTree->SetBranchAddress("selection", &sel);
  intTree->SetBranchAddress("pattern", &pattern);
  intTree->SetBranchAddress("axis_bin", &axisBin);
  intTree->SetBranchAddress("axis_center", &axisCenter);
  intTree->SetBranchAddress("count", &count);
  intTree->SetBranchAddress("stat_error", &statErr);
  intTree->SetBranchAddress("sys_error", &sysErr);
  if (hasAsymm) {
    intTree->SetBranchAddress("sys_error_up", &sysErrUp);
    intTree->SetBranchAddress("sys_error_down", &sysErrDown);
  }

  // group key -> all points for that (task, selection, pattern)
  map<string, vector<IntegratedPoint>> groups;

  const Long64_t nEntries = intTree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    intTree->GetEntry(i);
    string key = task + "|" + sel + "|" + pattern;
    IntegratedPoint p;
    p.axis_bin = *axisBin;
    p.axis_center = *axisCenter;
    p.count = count;
    p.stat_error = statErr;
    p.sys_error_up = hasAsymm ? sysErrUp : sysErr;
    p.sys_error_down = hasAsymm ? sysErrDown : sysErr;
    groups[key].push_back(p);
  }

  outDir->cd();
  int nWritten = 0, nSkipped = 0;
  for (auto& kv : groups) {
    vector<IntegratedPoint> pts = kv.second;

    // Only 0 or 1 remaining selector axes make sense as a 1D graph (x-axis
    // is that axis's bin center, or x=0 for the single "grand total" point
    // when no axis remains). Two or more remaining axes can't be reduced to
    // one x value without picking a dimension arbitrarily -- skip those for
    // now rather than silently guessing.
    size_t naxes = pts.empty() ? 0 : pts[0].axis_bin.size();
    if (naxes > 1) {
      nSkipped++;
      continue;
    }

    sort(pts.begin(), pts.end(), [](const IntegratedPoint& a, const IntegratedPoint& b) {
      int ba = a.axis_bin.empty() ? 0 : a.axis_bin[0];
      int bb = b.axis_bin.empty() ? 0 : b.axis_bin[0];
      return ba < bb;
    });

    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    g->SetName(kv.first.c_str());
    for (auto& p : pts) {
      double x = p.axis_center.empty() ? 0.0 : p.axis_center[0];
      double errUp = sqrt(p.stat_error * p.stat_error + p.sys_error_up * p.sys_error_up);
      double errDown = sqrt(p.stat_error * p.stat_error + p.sys_error_down * p.sys_error_down);
      int n = g->GetN();
      g->SetPoint(n, x, p.count);
      g->SetPointError(n, 0.0, 0.0, errDown, errUp);
    }
    g->Write();
    nWritten++;
  }
  cout << "Wrote " << nWritten << " integrated graphs (" << nSkipped
       << " skipped: more than one remaining selector axis)." << endl;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: ./BuildGraphs <file.root>\n";
    return -1;
  }

  TFile* f = TFile::Open(argv[1], "UPDATE");
  if (!f || f->IsZombie()) {
    cerr << "Could not open " << argv[1] << " in UPDATE mode.\n";
    return -1;
  }

  TTree* diffTree = (TTree*)f->Get("diffTable");
  TTree* intTree = (TTree*)f->Get("integratedTable");
  if (!diffTree || !intTree) {
    cerr << "diffTable/integratedTable not found in " << argv[1]
         << " -- is this a Main_Figs_Binned output file?\n";
    return -1;
  }

  TDirectory* graphsDir = f->GetDirectory("graphs");
  if (!graphsDir) graphsDir = f->mkdir("graphs");

  buildDiffGraphs(diffTree, graphsDir);
  buildIntegratedGraphs(intTree, graphsDir);

  f->cd();
  f->Write("", TObject::kOverwrite);
  f->Close();
  return 0;
}
