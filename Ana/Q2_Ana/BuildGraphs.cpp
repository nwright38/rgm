// Rehydrates the flat diffTable / integratedTable TTrees written by
// Main_Figs_Binned.cpp into TGraphErrors objects, so the python plotting
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
// Error convention: each graph's y-error is sqrt(stat_error^2 +
// sys_error^2), matching the old pipeline's getGraphWithError() -- stat and
// sys are combined in quadrature into one number per point, not kept as
// separate graphs.

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
#include <TGraphErrors.h>
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
  double value_center, count, stat_error, sys_error;
};

// One point collected from an integratedTable row (already grouped by task/selection/pattern).
struct IntegratedPoint {
  vector<int> axis_bin;
  vector<double> axis_center;
  double count, stat_error, sys_error;
};

void buildDiffGraphs(TTree* diffTree, TDirectory* outDir) {
  TTreeReader reader(diffTree);
  TTreeReaderValue<string> r_task(reader, "task_name");
  TTreeReaderValue<string> r_sel(reader, "selection");
  TTreeReaderValue<vector<int>> r_axisBin(reader, "axis_bin");
  TTreeReaderValue<int> r_valueBin(reader, "value_bin");
  TTreeReaderValue<double> r_valueCenter(reader, "value_center");
  TTreeReaderValue<double> r_count(reader, "count");
  TTreeReaderValue<double> r_statErr(reader, "stat_error");
  TTreeReaderValue<double> r_sysErr(reader, "sys_error");

  // group key -> all points for that (task, selection, selector-bin) combo
  map<string, vector<DiffPoint>> groups;

  while (reader.Next()) {
    string key = *r_task + "|" + *r_sel + "|" + joinBins(*r_axisBin);
    DiffPoint p;
    p.value_bin = *r_valueBin;
    p.value_center = *r_valueCenter;
    p.count = *r_count;
    p.stat_error = *r_statErr;
    p.sys_error = *r_sysErr;
    groups[key].push_back(p);
  }

  outDir->cd();
  for (auto& kv : groups) {
    vector<DiffPoint> pts = kv.second;
    sort(pts.begin(), pts.end(),
         [](const DiffPoint& a, const DiffPoint& b) { return a.value_bin < b.value_bin; });

    TGraphErrors* g = new TGraphErrors();
    g->SetName(kv.first.c_str());
    for (auto& p : pts) {
      double err = sqrt(p.stat_error * p.stat_error + p.sys_error * p.sys_error);
      int n = g->GetN();
      g->SetPoint(n, p.value_center, p.count);
      g->SetPointError(n, 0.0, err);
    }
    g->Write();
  }
  cout << "Wrote " << groups.size() << " differential graphs." << endl;
}

void buildIntegratedGraphs(TTree* intTree, TDirectory* outDir) {
  TTreeReader reader(intTree);
  TTreeReaderValue<string> r_task(reader, "task_name");
  TTreeReaderValue<string> r_sel(reader, "selection");
  TTreeReaderValue<string> r_pattern(reader, "pattern");
  TTreeReaderValue<vector<int>> r_axisBin(reader, "axis_bin");
  TTreeReaderValue<vector<double>> r_axisCenter(reader, "axis_center");
  TTreeReaderValue<double> r_count(reader, "count");
  TTreeReaderValue<double> r_statErr(reader, "stat_error");
  TTreeReaderValue<double> r_sysErr(reader, "sys_error");

  // group key -> all points for that (task, selection, pattern)
  map<string, vector<IntegratedPoint>> groups;

  while (reader.Next()) {
    string key = *r_task + "|" + *r_sel + "|" + *r_pattern;
    IntegratedPoint p;
    p.axis_bin = *r_axisBin;
    p.axis_center = *r_axisCenter;
    p.count = *r_count;
    p.stat_error = *r_statErr;
    p.sys_error = *r_sysErr;
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

    TGraphErrors* g = new TGraphErrors();
    g->SetName(kv.first.c_str());
    for (auto& p : pts) {
      double x = p.axis_center.empty() ? 0.0 : p.axis_center[0];
      double err = sqrt(p.stat_error * p.stat_error + p.sys_error * p.sys_error);
      int n = g->GetN();
      g->SetPoint(n, x, p.count);
      g->SetPointError(n, 0.0, err);
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
