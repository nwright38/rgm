// ExtractFitQuantities.C
//
// ROOT macro (run interpreted, not compiled into the build) that computes
// the Q^2-dependent quantities that can't come from simple bin sums and
// therefore were deliberately left out of diffTable/integratedTable when
// Main_Figs_Binned.cpp was written:
//   - sigma_pcmx, sigma_pcmy: Gaussian-fit width of the C.M. momentum
//     components, vs Q^2.
//   - std dev / mean of E1miss, E2miss vs Q^2, in each pMiss bin and each
//     kMiss bin.
//
// Same methodology as the old Main_Figs_Sys_Err.cpp: for each of the 100
// systematic-variation ("toy") histograms, extract the value (a Gaussian
// fit sigma, or a plain TH1::GetStdDev()/GetMean()), then the point plotted
// is the MEAN of those 100 toy values and its error is their STDDEV --
// the nominal histogram itself is not used as the central value here,
// matching the original getSigma()/getStdDev()/getMean() + getMeanStddev()
// pattern exactly. (If you'd rather center on the nominal fit and treat the
// toy spread as a separate systematic on top of the nominal's own stat
// error -- closer to how diffTable/integratedTable handle other
// quantities -- that's a deliberate methodology change from the original
// and should be a separate decision, not something silently swapped in
// here.)
//
// Reads the hists/nominal and hists/toy_NNN histograms Main_Figs_Binned.cpp
// already wrote (nominal is opened but currently unused, kept available for
// the alternative methodology above). Writes one TGraphErrors per quantity
// into the same "graphs" TDirectory BuildGraphs.cpp writes into, using the
// SAME naming convention as graph_names.py's diff_graph_name(): one entry
// per quantity is named "<quantity>|<selection>|<bin-suffix>", so the
// existing python graph_io/graph_names code can read these with no changes.
//
// Usage:
//   root -l -b -q 'ExtractFitQuantities.C("Data_He.root")'

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>

using namespace std;

const int kNToys = 100;
vector<double> bE_Q2 = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
const int kNQ2Bins = (int)bE_Q2.size() - 1;

double Q2center(int j) { return 0.5 * (bE_Q2[j] + bE_Q2[j + 1]); }

// Same normalized-Gaussian form as the original G() in Main_Figs_Sys_Err.cpp.
double G(double x, double N, double mu, double sigma) {
  double z = (x - mu) / sigma;
  return (N / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * z * z);
}

// Same as the original getSigma(): Gaussian-fits h over [min,max] and
// returns the fit sigma, falling back to the histogram's own GetStdDev()
// if there are too few entries to fit reliably.
double fitSigma(TH1D* h, double min, double max) {
  if (!h || h->GetEntries() == 0) return 0.0;
  if (h->GetEntries() < 20) return h->GetStdDev();

  TF1* gFit = new TF1("GausFit",
                       [](double* x, double* p) { return G(x[0], p[0], p[1], p[2]); },
                       min, max, 3);
  gFit->SetParameter(0, h->GetMaximum() / G(0, 1, 0, 0.1));
  gFit->SetParameter(1, (max + min) / 2);
  gFit->SetParLimits(1, min, max);
  gFit->SetParameter(2, (max - min) / 4);
  gFit->SetParLimits(2, 0.0, max - min);

  TFitResultPtr fr = h->Fit(gFit, "SrBeqn", "", min, max);
  double sigma = (fr.Get() != nullptr) ? fr->Parameter(2) : h->GetStdDev();
  delete gFit;
  return sigma;
}

double plainStdDev(TH1D* h) {
  if (!h || h->GetEntries() == 0) return 0.0;
  return h->GetStdDev();
}

double plainMean(TH1D* h) {
  if (!h || h->GetEntries() == 0) return 0.0;
  return h->GetMean();
}

void meanStddev(const vector<double>& v, double& mean, double& stddev) {
  mean = 0.0;
  stddev = 0.0;
  if (v.empty()) return;
  for (double x : v) mean += x;
  mean /= v.size();
  double s2 = 0.0;
  for (double x : v) s2 += (x - mean) * (x - mean);
  stddev = sqrt(s2 / v.size());
}

string toyLabel(int i) {
  char buf[32];
  snprintf(buf, sizeof(buf), "toy_%03d", i);
  return string(buf);
}

// Matches the histogram-naming convention HistStore::bookTask() uses in
// BinnedHistStore.h: "<label>_<task><axisName><axisBin>_Q2<q2bin>" when
// there is a selector axis besides Q2, or "<label>_<task>_Q2<q2bin>" when
// Q2 is the only selector axis (e.g. the pcm tasks).
string histName(const string& label, const string& sourceTask,
                const string& axisName, int axisBin, int q2Bin) {
  string n = label + "_" + sourceTask;
  if (!axisName.empty()) n += "_" + axisName + to_string(axisBin);
  n += "_Q2" + to_string(q2Bin);
  return n;
}

TH1D* getHist(TDirectory* histsDir, const string& label, const string& name) {
  TDirectory* d = histsDir->GetDirectory(label.c_str());
  if (!d) return nullptr;
  return (TH1D*)d->Get(name.c_str());
}

// Matches how Main_Figs_Sys_Err.cpp defines Q2_mean[]:
// Q2_mean[j] = dataGroup.h_Q2_epp_SRC_Q2[j]->GetMean().
double q2MeanFromNominal(TDirectory* histsDir, int q2Bin) {
  string hname = histName("nominal", "Q2_epp_SRC_Q2", "", 0, q2Bin);
  TH1D* h = getHist(histsDir, "nominal", hname);
  if (!h || h->GetEntries() == 0) return Q2center(q2Bin);
  return h->GetMean();
}

enum Method { STDDEV, MEAN, FIT_SIGMA };

// One row per quantity to extract. nAxisBins=1 with axisName="" means "no
// selector besides Q2" (the pcm tasks); nAxisBins=4 means one graph per
// pMiss or kMiss bin (axisBin = 0..3).
struct FitSpec {
  string quantityName;  // becomes the "<quantity>" part of the graph name
  string selection;     // "ep" or "epp" -- just a label here, the source
                        // histogram already reflects the right selection
  string sourceTask;    // FillTask name from Main_Figs_Binned.cpp
  string axisName;      // "" , "pMiss", or "kMiss"
  int nAxisBins;
  Method method;
  double fitMin, fitMax;  // only used for FIT_SIGMA
};

vector<FitSpec> buildSpecs() {
  vector<FitSpec> specs;

  specs.push_back({"sigma_pcmx", "epp", "pcmx_epp_SRC_Q2", "", 1, FIT_SIGMA, -0.2, 0.2});
  specs.push_back({"sigma_pcmy", "epp", "pcmy_epp_SRC_Q2", "", 1, FIT_SIGMA, -0.2, 0.2});

  struct EmissTask { string suffix; string selection; string task; };
  vector<EmissTask> emissTasks = {
      {"E1miss_ep_pmiss",   "ep",  "E1miss_ep_SRC_pmiss_Q2"},
      {"E1miss_epp_pmiss",  "epp", "E1miss_epp_SRC_pmiss_Q2"},
      {"E2miss_epp_pmiss",  "epp", "E2miss_epp_SRC_pmiss_Q2"},
      {"E1miss_ep_kmiss",   "ep",  "E1miss_ep_SRC_kmiss_Q2"},
      {"E1miss_epp_kmiss",  "epp", "E1miss_epp_SRC_kmiss_Q2"},
      {"E2miss_epp_kmiss",  "epp", "E2miss_epp_SRC_kmiss_Q2"},
  };
  for (auto& t : emissTasks) {
    string axisName = (t.suffix.find("pmiss") != string::npos) ? "pMiss" : "kMiss";
    specs.push_back({"stddev_" + t.suffix, t.selection, t.task, axisName, 4, STDDEV, 0, 0});
    specs.push_back({"mean_" + t.suffix,   t.selection, t.task, axisName, 4, MEAN,   0, 0});
  }

  return specs;
}

string joinBins(const vector<int>& bins) {
  if (bins.empty()) return "none";
  string s;
  for (size_t i = 0; i < bins.size(); i++) {
    if (i > 0) s += "_";
    s += to_string(bins[i]);
  }
  return s;
}

double extract(const FitSpec& spec, TH1D* h) {
  if (spec.method == STDDEV) return plainStdDev(h);
  if (spec.method == MEAN) return plainMean(h);
  return fitSigma(h, spec.fitMin, spec.fitMax);
}

void ExtractFitQuantities(const char* filename) {
  TFile* f = TFile::Open(filename, "UPDATE");
  if (!f || f->IsZombie()) {
    cerr << "Could not open " << filename << " in UPDATE mode.\n";
    return;
  }

  TDirectory* histsDir = f->GetDirectory("hists");
  if (!histsDir) {
    cerr << "No 'hists' directory in " << filename
         << " -- is this a Main_Figs_Binned output file?\n";
    return;
  }

  TDirectory* graphsDir = f->GetDirectory("graphs");
  if (!graphsDir) graphsDir = f->mkdir("graphs");

  vector<double> q2x(kNQ2Bins, 0.0);
  for (int j = 0; j < kNQ2Bins; j++) q2x[j] = q2MeanFromNominal(histsDir, j);

  vector<FitSpec> specs = buildSpecs();
  int nWritten = 0;
  int nSkipped = 0;
  int nMissing = 0;

  for (auto& spec : specs) {
    for (int axisBin = 0; axisBin < spec.nAxisBins; axisBin++) {
      TGraphErrors* g = new TGraphErrors();
      vector<int> binKey;
      if (!spec.axisName.empty()) binKey.push_back(axisBin);
      string gname = spec.quantityName + "|" + spec.selection + "|" + joinBins(binKey);
      g->SetName(gname.c_str());

      for (int j = 0; j < kNQ2Bins; j++) {
        vector<double> toyVals;
        toyVals.reserve(kNToys);
        int missingThisPoint = 0;
        string firstMissing;
        for (int i = 0; i < kNToys; i++) {
          string label = toyLabel(i);
          string hname = histName(label, spec.sourceTask, spec.axisName, axisBin, j);
          TH1D* h = getHist(histsDir, label, hname);
          if (!h) {
            missingThisPoint++;
            if (firstMissing.empty()) firstMissing = label + "/" + hname;
            continue;
          }
          toyVals.push_back(extract(spec, h));
        }

        if (missingThisPoint > 0) {
          nMissing += missingThisPoint;
          cerr << "[ExtractFitQuantities] Missing " << missingThisPoint
               << " toy hist(s) for graph '" << gname << "', Q2 bin " << j
               << ". Example missing path: hists/" << firstMissing << "\n";
        }

        // Preserve the original methodology: central value is the mean over
        // toy-derived quantities and error is the toy spread. If no toy hists
        // were found for this point, skip it rather than silently writing zero.
        if (toyVals.empty()) {
          continue;
        }
        double mean, stddev;
        meanStddev(toyVals, mean, stddev);
        int n = g->GetN();
        g->SetPoint(n, q2x[j], mean);
        g->SetPointError(n, 0.0, stddev);
      }

      if (g->GetN() == 0) {
        nSkipped++;
        delete g;
        continue;
      }

      graphsDir->cd();
      g->Write("", TObject::kOverwrite);
      nWritten++;
    }
  }

  cout << "Wrote " << nWritten << " fit-derived graphs into " << filename << ":/graphs"
       << " (" << nSkipped << " skipped with no valid points, "
       << nMissing << " missing toy histogram lookups).\n";

  f->cd();
  f->Write("", TObject::kOverwrite);
  f->Close();
}
