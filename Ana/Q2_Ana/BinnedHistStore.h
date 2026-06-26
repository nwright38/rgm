#ifndef BINNEDHISTSTORE_HH
#define BINNEDHISTSTORE_HH

// Generic N-dimensional binning + histogram bookkeeping, with no physics in it.
// Analysis-specific code (kinematics, cuts, the list of plots to make) lives in
// Main_Figs_Binned.cpp and is passed in as data (FillTask) or via the EventT
// template parameter -- this header never names a physics variable directly.

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "TH1D.h"
#include "TDirectory.h"

// One binning axis: a single source of truth for a set of bin edges, so the
// same vector<double> is never retyped/duplicated in multiple places.
struct Axis {
  std::string name;
  std::vector<double> edges;

  int nbins() const { return (int)edges.size() - 1; }

  // Same convention as the original binX(): -1 if out of range.
  int findBin(double x) const {
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
      if (x >= edges[i] && x < edges[i + 1]) return (int)i;
    }
    return -1;
  }

  double center(int bin) const {
    return 0.5 * (edges.at(bin) + edges.at(bin + 1));
  }
};

enum class Selection { EP, EPP };

enum class CentralValueMode { NOMINAL, TOY_MEAN };

inline const char* selectionName(Selection s) {
  return (s == Selection::EP) ? "ep" : "epp";
}

// Describes one quantity to histogram: which selector axes pick the
// histogram (e.g. {pMissAxis, Q2Axis}), what value goes into that histogram
// (e.g. E1miss), what binning the histogrammed value itself uses, and which
// selection (ep/epp) gates whether an event is filled at all.
template <typename EventT>
struct FillTask {
  std::string name;
  Selection selection = Selection::EP;

  // Returns whether this event is filled at all under this task (e.g. ek.passep).
  std::function<bool(const EventT&)> passFn;

  // Selector axes: which bin of each axis routes the event to which histogram.
  std::vector<const Axis*> binAxes;
  std::vector<std::function<double(const EventT&)>> axisValueFns; // one per binAxes entry

  // The quantity actually filled into the chosen histogram.
  std::function<double(const EventT&)> valueFn;

  // Histogram binning for valueFn. If valueEdges is non-empty it takes
  // priority (variable-width bins); otherwise nbinsHist/xmin/xmax is used.
  int nbinsHist = 1;
  double xmin = 0.0, xmax = 1.0;
  std::vector<double> valueEdges;
};

// Holds one full set of booked/filled histograms -- either the nominal set
// or a single systematic-variation ("toy") set. Construct one HistStore per
// toy plus one for nominal; book() must be called with the same tasks vector
// (by content) every time so bin keys line up across stores.
template <typename EventT>
class HistStore {
 public:
  explicit HistStore(std::string label) : label_(std::move(label)) {}

  void book(const std::vector<FillTask<EventT>>& tasks) {
    hists_.clear();
    for (size_t t = 0; t < tasks.size(); ++t) bookTask(tasks[t], (int)t);
  }

  void fill(const std::vector<FillTask<EventT>>& tasks, const EventT& ev,
            double wep, double wepp) {
    for (size_t t = 0; t < tasks.size(); ++t) {
      const FillTask<EventT>& task = tasks[t];
      if (task.passFn && !task.passFn(ev)) continue;

      std::vector<int> binIdx;
      binIdx.reserve(task.binAxes.size());
      bool ok = true;
      for (size_t a = 0; a < task.binAxes.size(); ++a) {
        double x = task.axisValueFns[a](ev);
        int b = task.binAxes[a]->findBin(x);
        if (b < 0) { ok = false; break; }
        binIdx.push_back(b);
      }
      if (!ok) continue;

      TH1D* h = get((int)t, binIdx);
      if (!h) continue;  // book() not called with a matching task list
      double w = (task.selection == Selection::EP) ? wep : wepp;
      h->Fill(task.valueFn(ev), w);
    }
  }

  TH1D* get(int taskIdx, const std::vector<int>& binIdx) const {
    auto it = hists_.find(makeKey(taskIdx, binIdx));
    return it == hists_.end() ? nullptr : it->second;
  }

  // All booked (taskIdx ++ binIdx) -> histogram entries, for output/integration code.
  const std::map<std::vector<int>, TH1D*>& raw() const { return hists_; }

  const std::string& label() const { return label_; }

  // Writes every booked histogram into a subdirectory of the current
  // directory named after this store's label (e.g. "nominal", "toy_007").
  void writeAll() const {
    TDirectory* parent = gDirectory;
    TDirectory* dir = parent->mkdir(label_.c_str());
    if (!dir) dir = parent->GetDirectory(label_.c_str());
    dir->cd();
    for (auto& kv : hists_) kv.second->Write();
    parent->cd();
  }

 private:
  std::string label_;
  std::map<std::vector<int>, TH1D*> hists_;

  static std::vector<int> makeKey(int taskIdx, const std::vector<int>& binIdx) {
    std::vector<int> key;
    key.reserve(binIdx.size() + 1);
    key.push_back(taskIdx);
    key.insert(key.end(), binIdx.begin(), binIdx.end());
    return key;
  }

  void bookTask(const FillTask<EventT>& task, int taskIdx) {
    std::vector<int> nb;
    nb.reserve(task.binAxes.size());
    for (auto* ax : task.binAxes) nb.push_back(ax->nbins());

    std::vector<int> idx(nb.size(), 0);
    bool done = false;
    do {
      std::string hn = label_ + "_" + task.name;
      for (size_t a = 0; a < idx.size(); ++a) {
        hn += "_" + task.binAxes[a]->name + std::to_string(idx[a]);
      }
      TH1D* h = makeHist(hn, task);
      h->Sumw2();
      hists_[makeKey(taskIdx, idx)] = h;

      if (nb.empty()) break;

      // Odometer-style increment over all selector-axis bin combinations.
      int pos = (int)idx.size() - 1;
      while (pos >= 0) {
        idx[pos]++;
        if (idx[pos] < nb[pos]) break;
        idx[pos] = 0;
        pos--;
      }
      done = (pos < 0);
    } while (!done);
  }

  static TH1D* makeHist(const std::string& name, const FillTask<EventT>& task) {
    if (!task.valueEdges.empty()) {
      return new TH1D(name.c_str(), name.c_str(),
                       (int)task.valueEdges.size() - 1, task.valueEdges.data());
    }
    return new TH1D(name.c_str(), name.c_str(), task.nbinsHist, task.xmin, task.xmax);
  }
};

// Mean/stddev of a list of toy values -- same definition as the original
// getMeanStddev(), reused for both per-bin and integrated systematic errors.
inline void meanStddev(const std::vector<double>& v, double& mean, double& stddev) {
  mean = 0.0;
  stddev = 0.0;
  if (v.empty()) return;
  mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
  double s2 = 0.0;
  for (double x : v) s2 += (x - mean) * (x - mean);
  stddev = std::sqrt(s2 / v.size());
}

// Percentile helper with linear interpolation between adjacent ordered values.
// p in [0,1], e.g. 0.16 / 0.50 / 0.84.
inline double percentile(std::vector<double> v, double p) {
  if (v.empty()) return 0.0;
  if (p <= 0.0) p = 0.0;
  if (p >= 1.0) p = 1.0;
  std::sort(v.begin(), v.end());
  double x = p * (double)(v.size() - 1);
  size_t i0 = (size_t)std::floor(x);
  size_t i1 = (size_t)std::ceil(x);
  if (i0 == i1) return v[i0];
  double t = x - (double)i0;
  return (1.0 - t) * v[i0] + t * v[i1];
}

// One row of the flat differential output table: which task, which selector
// bin, which value-histogram bin, and its nominal count / stat error / sys error.
struct DiffRow {
  std::string task_name;
  std::string selection;
  std::vector<int> axis_bin;      // selector-axis bin indices, in task.binAxes order
  std::vector<double> axis_center;
  int value_bin = 0;              // bin index within the histogrammed quantity
  double value_center = 0.0;
  double count = 0.0;
  double stat_error = 0.0;
  double sys_error = 0.0;
  double sys_error_up = 0.0;
  double sys_error_down = 0.0;
};

// Builds the differential table for one task: every (selector-bin, value-bin)
// combination, with sys_error taken from the spread of the same bin's content
// across all toy stores.
template <typename EventT>
std::vector<DiffRow> buildDiffRows(const FillTask<EventT>& task, int taskIdx,
                                    const HistStore<EventT>& nominal,
                                    const std::vector<HistStore<EventT>>& toys,
                                    CentralValueMode centralMode = CentralValueMode::NOMINAL) {
  std::vector<DiffRow> rows;
  for (auto& kv : nominal.raw()) {
    const std::vector<int>& key = kv.first;
    if (key.front() != taskIdx) continue;
    std::vector<int> binIdx(key.begin() + 1, key.end());
    TH1D* hNom = kv.second;

    std::vector<double> centers;
    for (size_t a = 0; a < task.binAxes.size(); ++a)
      centers.push_back(task.binAxes[a]->center(binIdx[a]));

    for (int vb = 1; vb <= hNom->GetNbinsX(); ++vb) {
      DiffRow row;
      row.task_name = task.name;
      row.selection = selectionName(task.selection);
      row.axis_bin = binIdx;
      row.axis_center = centers;
      row.value_bin = vb - 1;
      row.value_center = hNom->GetXaxis()->GetBinCenter(vb);
      row.count = hNom->GetBinContent(vb);
      row.stat_error = hNom->GetBinError(vb);

      std::vector<double> toyVals;
      toyVals.reserve(toys.size());
      for (auto& toy : toys) {
        TH1D* hToy = toy.get(taskIdx, binIdx);
        toyVals.push_back(hToy ? hToy->GetBinContent(vb) : 0.0);
      }
      double m, s;
      meanStddev(toyVals, m, s);
      if (centralMode == CentralValueMode::TOY_MEAN) row.count = m;
      row.sys_error = s;
      row.sys_error_up = s;
      row.sys_error_down = s;

      rows.push_back(std::move(row));
    }
  }
  return rows;
}

// One row of the integrated output table: same task, but with one or more
// selector axes (and optionally the value axis itself) summed over.
struct IntegratedRow {
  std::string task_name;
  std::string selection;
  std::string pattern;            // human-readable description, e.g. "collapse:pMiss"
  std::vector<int> axis_bin;      // remaining selector-axis bin indices (after collapse)
  std::vector<double> axis_center;
  double count = 0.0;
  double stat_error = 0.0;
  double sys_error = 0.0;
  double sys_error_up = 0.0;
  double sys_error_down = 0.0;
};

// Sums every booked histogram's full integral (Integral() over all value
// bins) for a task, collapsing the selector axes listed in collapseAxisPos
// (positions into task.binAxes) and keeping the rest as the row's bin key.
// sys_error is computed by summing each TOY's integrals first (over the same
// collapsed selector bins), then taking the stddev across toys -- this is
// what correctly captures correlation between bins from a shared toy cut
// variation, instead of naively combining the per-bin sys errors in quadrature.
template <typename EventT>
std::vector<IntegratedRow> buildIntegratedRows(
    const FillTask<EventT>& task, int taskIdx,
    const HistStore<EventT>& nominal, const std::vector<HistStore<EventT>>& toys,
    const std::vector<int>& collapseAxisPos, const std::string& patternLabel) {
  // Group nominal bin keys by their "kept" (non-collapsed) selector indices.
  std::map<std::vector<int>, std::vector<int>> groups;  // keptKey -> list of full binIdx (encoded as flat ints via lookup below)
  std::map<std::vector<int>, std::vector<std::vector<int>>> groupMembers;

  for (auto& kv : nominal.raw()) {
    const std::vector<int>& key = kv.first;
    if (key.front() != taskIdx) continue;
    std::vector<int> binIdx(key.begin() + 1, key.end());

    std::vector<int> kept;
    for (size_t a = 0; a < binIdx.size(); ++a) {
      bool collapse = false;
      for (int c : collapseAxisPos) if (c == (int)a) collapse = true;
      if (!collapse) kept.push_back(binIdx[a]);
    }
    groupMembers[kept].push_back(binIdx);
  }

  std::vector<IntegratedRow> rows;
  for (auto& g : groupMembers) {
    const std::vector<int>& kept = g.first;
    const std::vector<std::vector<int>>& members = g.second;

    IntegratedRow row;
    row.task_name = task.name;
    row.selection = selectionName(task.selection);
    row.pattern = patternLabel;

    // Map kept indices back to their original axis identity for centers.
    size_t keptPos = 0;
    for (size_t a = 0; a < task.binAxes.size(); ++a) {
      bool collapse = false;
      for (int c : collapseAxisPos) if (c == (int)a) collapse = true;
      if (!collapse) {
        row.axis_bin.push_back(kept[keptPos]);
        row.axis_center.push_back(task.binAxes[a]->center(kept[keptPos]));
        keptPos++;
      }
    }

    double nomSum = 0.0, statSq = 0.0;
    for (auto& binIdx : members) {
      TH1D* h = nominal.get(taskIdx, binIdx);
      if (!h) continue;
      double err = 0.0;
      nomSum += h->IntegralAndError(1, h->GetNbinsX(), err);
      statSq += err * err;
    }
    row.count = nomSum;
    row.stat_error = std::sqrt(statSq);

    std::vector<double> toySums;
    toySums.reserve(toys.size());
    for (auto& toy : toys) {
      double s = 0.0;
      for (auto& binIdx : members) {
        TH1D* h = toy.get(taskIdx, binIdx);
        if (h) s += h->Integral();
      }
      toySums.push_back(s);
    }
    double m, s;
    meanStddev(toySums, m, s);
    row.sys_error = s;
    row.sys_error_up = s;
    row.sys_error_down = s;

    rows.push_back(std::move(row));
  }
  return rows;
}

// Ratio of two tasks that share the same selector-axis structure (e.g.
// epp/ep yield ratio). Computed via TH1::Divide() on cloned histograms (never
// mutates the originals, unlike calling Divide() in place), for the nominal
// store and every toy store, then read out exactly like buildDiffRows.
//
// Systematic spread for ratios is evaluated from toy percentiles:
//   p16 = percentile(toyVals, 0.16), p50 = percentile(toyVals, 0.50),
//   p84 = percentile(toyVals, 0.84)
// and stored both as a legacy symmetric additive uncertainty
//   sigma_add = 0.5 * (p84 - p16)
// and as asymmetric distances from the chosen central value:
//   sys_down = max(0, central - p16)
//   sys_up   = max(0, p84 - central)
//
// If centralMode == TOY_MEAN, row.count uses p50 (toy median central value).
template <typename EventT>
std::vector<DiffRow> buildRatioDiffRows(const FillTask<EventT>& numTask, int numIdx,
                                         int denIdx, const std::string& ratioName,
                                         const HistStore<EventT>& nominal,
                                         const std::vector<HistStore<EventT>>& toys,
                                         CentralValueMode centralMode = CentralValueMode::NOMINAL) {
  std::vector<DiffRow> rows;
  for (auto& kv : nominal.raw()) {
    const std::vector<int>& key = kv.first;
    if (key.front() != numIdx) continue;
    std::vector<int> binIdx(key.begin() + 1, key.end());
    TH1D* hNumNom = kv.second;
    TH1D* hDenNom = nominal.get(denIdx, binIdx);
    if (!hDenNom) continue;

    TH1D* ratioNom = (TH1D*)hNumNom->Clone((ratioName + "_ratio_nominal_tmp").c_str());
    ratioNom->Divide(hDenNom);

    std::vector<double> centers;
    for (size_t a = 0; a < numTask.binAxes.size(); ++a)
      centers.push_back(numTask.binAxes[a]->center(binIdx[a]));

    for (int vb = 1; vb <= ratioNom->GetNbinsX(); ++vb) {
      DiffRow row;
      row.task_name = ratioName;
      row.selection = "ratio";
      row.axis_bin = binIdx;
      row.axis_center = centers;
      row.value_bin = vb - 1;
      row.value_center = ratioNom->GetXaxis()->GetBinCenter(vb);
      row.count = ratioNom->GetBinContent(vb);
      row.stat_error = ratioNom->GetBinError(vb);

      std::vector<double> toyVals;
      toyVals.reserve(toys.size());
      for (auto& toy : toys) {
        TH1D* hN = toy.get(numIdx, binIdx);
        TH1D* hD = toy.get(denIdx, binIdx);
        if (hN && hD) {
          TH1D* r = (TH1D*)hN->Clone("ratio_toy_tmp");
          r->Divide(hD);
          double v = r->GetBinContent(vb);
          if (std::isfinite(v)) toyVals.push_back(v);
          delete r;
        }
      }
      if (!toyVals.empty()) {
        double p16 = percentile(toyVals, 0.16);
        double p50 = percentile(toyVals, 0.50);
        double p84 = percentile(toyVals, 0.84);
        if (centralMode == CentralValueMode::TOY_MEAN) row.count = p50;
        row.sys_error_down = std::max(0.0, row.count - p16);
        row.sys_error_up = std::max(0.0, p84 - row.count);
        row.sys_error = 0.5 * (row.sys_error_up + row.sys_error_down);
        if (row.sys_error < 0.0) row.sys_error = 0.0;
      } else {
        row.sys_error = 0.0;
        row.sys_error_up = 0.0;
        row.sys_error_down = 0.0;
      }

      rows.push_back(std::move(row));
    }
    delete ratioNom;
  }
  return rows;
}

#endif
