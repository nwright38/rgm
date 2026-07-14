#include "SigmaCMPlotOutput.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TObject.h>
#include <TAxis.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace sigmacm {
namespace {

constexpr int kScanBins = 100;
constexpr double kSigmaMin = 0.08;
constexpr double kSigmaMax = 0.25;

struct Axis {
  const char* shortName;
  const char* histName;
  const char* fitName;
  const char* q2HistPrefix;
  const char* chi2Name;
  const char* scaleName;
  const char* sigIntName;
  const char* sigQ2Name;
  int bins;
  double min;
  double max;
  double fitMin;
  double fitMax;
};

std::array<Axis, 4> axes(const Config& cfg) {
  return {{{"pcmx", "pcmx_epp", "pcmx_epp_fit", "h_pcmx_epp_SRC_Q2",
            "g_chi2_pcmx_epp", "g_scale_pcmx_epp", "sigmacmx_int", "sigmacmx_Q2",
            cfg.binsX, cfg.histXMin, cfg.histXMax, -cfg.cutRangeXY, cfg.cutRangeXY},
           {"pcmy", "pcmy_epp", "pcmy_epp_fit", "h_pcmy_epp_SRC_Q2",
            "g_chi2_pcmy_epp", "g_scale_pcmy_epp", "sigmacmy_int", "sigmacmy_Q2",
            cfg.binsY, cfg.histYMin, cfg.histYMax, -cfg.cutRangeXY, cfg.cutRangeXY},
           {"pcmz", "pcmz_epp", "pcmz_epp_fit", "h_pcmz_epp_SRC_Q2",
            "g_chi2_pcmz_epp", "g_scale_pcmz_epp", "sigmacmz_int", "sigmacmz_Q2",
            cfg.binsZ, cfg.histZMin, cfg.histZMax, cfg.fitZMin, cfg.fitZMax},
           {"pcmT", "pcmT_epp", "pcmT_epp_fit", "h_pcmT_epp_SRC_Q2",
            "g_chi2_pcmT_epp", "g_scale_pcmT_epp", "sigmacmT_int", "sigmacmT_Q2",
            cfg.binsX, 0.0, std::sqrt(2.0) * cfg.cutRangeXY, 0.0, std::sqrt(2.0) * cfg.cutRangeXY}}};
}

double value(const Event& e, int axis) {
  if (axis == 0) return e.pcmX;
  if (axis == 1) return e.pcmY;
  if (axis == 2) return e.pcmZ;
  return std::sqrt(e.pcmX * e.pcmX + e.pcmY * e.pcmY);
}

double sigmaFor(const Result& r, int axis) {
  if (axis == 0) return r.sigmaX;
  if (axis == 1) return r.sigmaY;
  if (axis == 2) return r.sigmaZ;
  return std::sqrt(0.5 * (r.sigmaX * r.sigmaX + r.sigmaY * r.sigmaY));
}

double sigmaErrFor(const Result& r, int axis, bool high) {
  if (axis == 0) return high ? r.sigmaXErrHigh : r.sigmaXErrLow;
  if (axis == 1) return high ? r.sigmaYErrHigh : r.sigmaYErrLow;
  if (axis == 2) return high ? r.sigmaZErrHigh : r.sigmaZErrLow;
  const double ex = high ? r.sigmaXErrHigh : r.sigmaXErrLow;
  const double ey = high ? r.sigmaYErrHigh : r.sigmaYErrLow;
  return 0.5 * std::sqrt(ex * ex + ey * ey);
}

std::unique_ptr<TH1D> makeHist(const std::string& name, const Axis& axis) {
  auto h = std::make_unique<TH1D>(name.c_str(), name.c_str(), axis.bins, axis.min, axis.max);
  h->Sumw2();
  return h;
}

void fillData(TH1D& h, const std::vector<Event>& events, const Config& cfg, int axis) {
  for (const auto& e : events) {
    if (passesCuts(e, cfg)) h.Fill(value(e, axis));
  }
}

void fillMC(TH1D& h, const std::vector<Event>& events, double sigmaGen,
            const Result& r, const Config& cfg, int axis,
            double sigmaOverride, int overrideAxis) {
  double sx = r.sigmaX;
  double sy = r.sigmaY;
  double sz = r.sigmaZ;
  if (overrideAxis == 0) sx = sigmaOverride;
  if (overrideAxis == 1) sy = sigmaOverride;
  if (overrideAxis == 2) sz = sigmaOverride;
  if (overrideAxis == 3) {
    sx = sigmaOverride;
    sy = sigmaOverride;
  }
  for (const auto& e : events) {
    if (!passesCuts(e, cfg)) continue;
    const double w = eventWeight(e, sx, sy, sz, sigmaGen, cfg.auxWeightBranch);
    h.Fill(value(e, axis), w);
  }
}

double chi2WithScale(const TH1D& data, const TH1D& mc, const Axis& axis, double scale) {
  double chi2 = 0.0;
  const int first = data.GetXaxis()->FindBin(axis.fitMin);
  const int last = data.GetXaxis()->FindBin(axis.fitMax);
  for (int b = first; b <= last; ++b) {
    const double d = data.GetBinContent(b);
    const double e = mc.GetBinContent(b);
    const double de = data.GetBinError(b);
    const double ee = mc.GetBinError(b);
    const double var = de * de + scale * scale * ee * ee;
    if (var <= 0.0) continue;
    chi2 += (d - scale * e) * (d - scale * e) / var;
  }
  return chi2;
}

std::pair<double, double> bestScaleChi2(const TH1D& data, const TH1D& mc, const Axis& axis) {
  const double mcIntegral = mc.Integral();
  const double intScale = mcIntegral > 0.0 ? data.Integral() / mcIntegral : 1.0;
  const double start = std::max(1.0e-9, 0.1 * intScale);
  const double stop = std::max(start, 3.0 * intScale);
  double bestScale = intScale;
  double bestChi2 = 1.0e100;
  for (int i = 0; i < 100; ++i) {
    const double scale = start + (stop - start) * i / 99.0;
    const double chi2 = chi2WithScale(data, mc, axis, scale);
    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestScale = scale;
    }
  }
  return {bestScale, bestChi2};
}

std::pair<double, double> bestJointScaleChi2(const std::array<TH1D*, 3>& data,
                                             const std::array<TH1D*, 3>& mc,
                                             const std::array<Axis, 4>& axisDefs,
                                             const Config& cfg) {
  if (!cfg.sharedScale) {
    double chi2 = 0.0;
    double scale = 0.0;
    for (int a = 0; a < 3; ++a) {
      auto [axisScale, axisChi2] = bestScaleChi2(*data[a], *mc[a], axisDefs[a]);
      scale += axisScale;
      chi2 += axisChi2;
    }
    return {scale / 3.0, chi2};
  }

  double dataIntegral = 0.0;
  double mcIntegral = 0.0;
  for (int a = 0; a < 3; ++a) {
    dataIntegral += data[a]->Integral();
    mcIntegral += mc[a]->Integral();
  }
  const double intScale = mcIntegral > 0.0 ? dataIntegral / mcIntegral : 1.0;
  const double start = std::max(1.0e-9, 0.1 * intScale);
  const double stop = std::max(start, 3.0 * intScale);
  double bestScale = intScale;
  double bestChi2 = 1.0e100;
  for (int i = 0; i < 120; ++i) {
    const double scale = start + (stop - start) * i / 119.0;
    double chi2 = 0.0;
    for (int a = 0; a < 3; ++a) {
      chi2 += chi2WithScale(*data[a], *mc[a], axisDefs[a], scale);
    }
    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestScale = scale;
    }
  }
  return {bestScale, bestChi2};
}

int resultIndexFor(const std::vector<Result>& results, int q2Bin) {
  for (size_t i = 0; i < results.size(); ++i) {
    const int idx = results[i].config.integratedQ2 ? -1 : results[i].config.q2BinIndex;
    if (idx == q2Bin) return static_cast<int>(i);
  }
  return -1;
}

void writeObject(TObject* obj) {
  obj->Write(obj->GetName(), TObject::kOverwrite);
}

void writeOverlayCanvas(const std::string& name, TH1D& data, TH1D& fit) {
  TCanvas c(name.c_str(), name.c_str(), 1000, 760);
  data.SetLineColor(kBlack);
  data.SetMarkerColor(kBlack);
  data.SetMarkerStyle(20);
  fit.SetLineColor(kRed + 1);
  fit.SetLineWidth(2);
  data.Draw("E");
  fit.Draw("HIST SAME");
  c.Write(name.c_str(), TObject::kOverwrite);
}

void writeGraphCanvas(const std::string& name, TGraph& graph) {
  TCanvas c(name.c_str(), name.c_str(), 1000, 760);
  graph.Draw("AL");
  c.Write(name.c_str(), TObject::kOverwrite);
}

void fillSigmaPoint(TGraphAsymmErrors& g, const Result& r, int axis, double x,
                    double xLow, double xHigh) {
  const double s = sigmaFor(r, axis);
  const int p = g.GetN();
  g.SetPoint(p, x, s);
  g.SetPointError(p, xLow, xHigh, sigmaErrFor(r, axis, false), sigmaErrFor(r, axis, true));
}

void writeJointChi2Scans(const std::string& suffix, const std::vector<Event>& dataEvents,
                         const std::vector<Event>& mcEvents, double sigmaGen,
                         const Result& r, const Config& cfg,
                         const std::array<Axis, 4>& axisDefs) {
  std::array<std::unique_ptr<TH1D>, 3> dataOwned;
  std::array<TH1D*, 3> data{};
  for (int a = 0; a < 3; ++a) {
    dataOwned[a] = makeHist("joint_data_tmp_" + std::to_string(a), axisDefs[a]);
    fillData(*dataOwned[a], dataEvents, cfg, a);
    data[a] = dataOwned[a].get();
  }

  const char* sigmaNames[] = {"sigmaX", "sigmaY", "sigmaZ"};
  for (int varied = 0; varied < 3; ++varied) {
    TGraph chi2;
    chi2.SetName(("g_chi2_joint_" + std::string(sigmaNames[varied]) + suffix).c_str());
    TGraph scale;
    scale.SetName(("g_scale_joint_" + std::string(sigmaNames[varied]) + suffix).c_str());
    for (int s = 0; s < kScanBins; ++s) {
      const double sigma = kSigmaMin + (kSigmaMax - kSigmaMin) * s / kScanBins;
      std::array<std::unique_ptr<TH1D>, 3> mcOwned;
      std::array<TH1D*, 3> mc{};
      for (int a = 0; a < 3; ++a) {
        mcOwned[a] = makeHist("joint_scan_tmp_" + std::to_string(a), axisDefs[a]);
        fillMC(*mcOwned[a], mcEvents, sigmaGen, r, cfg, a, sigma, varied);
        mc[a] = mcOwned[a].get();
      }
      auto [bestScale, bestChi2] = bestJointScaleChi2(data, mc, axisDefs, cfg);
      chi2.SetPoint(chi2.GetN(), sigma, bestChi2);
      scale.SetPoint(scale.GetN(), sigma, bestScale);
    }
    writeObject(&chi2);
    writeObject(&scale);
    writeGraphCanvas("c_" + std::string(chi2.GetName()), chi2);
    writeGraphCanvas("c_" + std::string(scale.GetName()), scale);
  }
}

}  // namespace

void writePlottingRootObjects(const std::string& path,
                              const std::vector<Event>& dataEvents,
                              const std::vector<Event>& mcEvents,
                              double sigmaGen,
                              const std::vector<Result>& results) {
  if (results.empty()) return;
  const int integratedIndex = resultIndexFor(results, -1);
  if (integratedIndex < 0) return;
  const Result& integrated = results[static_cast<size_t>(integratedIndex)];
  const Config& baseCfg = integrated.config;
  const auto axisDefs = axes(baseCfg);

  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "UPDATE"));
  if (!file || file->IsZombie()) throw std::runtime_error("Could not update ROOT file '" + path + "'");
  TNamed note("sigmacm_plotting_surface",
              "Standard ROOT hists/graphs/canvases written by Sigma_CM skim fitter");
  note.Write(note.GetName(), TObject::kOverwrite);

  std::array<std::unique_ptr<TGraphAsymmErrors>, 4> q2Graphs;
  for (int a = 0; a < 4; ++a) {
    q2Graphs[a] = std::make_unique<TGraphAsymmErrors>();
    q2Graphs[a]->SetName(axisDefs[a].sigQ2Name);
  }

  const int nQ2 = static_cast<int>(baseCfg.q2Edges.size()) - 1;
  for (int q2 = -1; q2 < nQ2; ++q2) {
    const int rIndex = resultIndexFor(results, q2);
    if (rIndex < 0) continue;
    const Result& r = results[static_cast<size_t>(rIndex)];
    Config cfg = r.config;
    const bool integratedMode = q2 < 0;
    const std::string suffix = integratedMode ? "" : ("_SRC_Q2_" + std::to_string(q2));
    const std::string jointSuffix = integratedMode ? "_epp" : ("_epp_SRC_Q2_" + std::to_string(q2));

    if (!integratedMode) {
      auto hq2 = makeHist("h_Q2_" + std::to_string(q2), {"Q2", "", "", "", "", "", "", "",
                          80, cfg.q2Edges[q2], cfg.q2Edges[q2 + 1], cfg.q2Edges[q2], cfg.q2Edges[q2 + 1]});
      for (const auto& e : dataEvents) {
        if (passesCuts(e, cfg)) hq2->Fill(e.Q2);
      }
      writeObject(hq2.get());
      const double center = 0.5 * (cfg.q2Edges[q2] + cfg.q2Edges[q2 + 1]);
      for (int a = 0; a < 4; ++a) {
        fillSigmaPoint(*q2Graphs[a], r, a, center, center - cfg.q2Edges[q2], cfg.q2Edges[q2 + 1] - center);
      }
    }

    writeJointChi2Scans(jointSuffix, dataEvents, mcEvents, sigmaGen, r, cfg, axisDefs);

    for (int a = 0; a < 4; ++a) {
      const Axis& axis = axisDefs[a];
      const std::string dataName = integratedMode ? axis.histName : axis.q2HistPrefix + suffix.substr(7);
      auto hData = makeHist(dataName, axis);
      fillData(*hData, dataEvents, cfg, a);
      writeObject(hData.get());

      auto hBest = makeHist(integratedMode ? axis.fitName : dataName + "_fit", axis);
      fillMC(*hBest, mcEvents, sigmaGen, r, cfg, a, 0.0, -1);
      auto [bestScale, bestChi2] = bestScaleChi2(*hData, *hBest, axis);
      hBest->Scale(bestScale);
      writeObject(hBest.get());

      TGraph chi2;
      TGraph scale;
      chi2.SetName(integratedMode ? axis.chi2Name
                                  : ("g_chi2_" + std::string(axis.shortName) + "_epp_SRC_Q2_" + std::to_string(q2)).c_str());
      scale.SetName(integratedMode ? axis.scaleName
                                   : ("g_scale_" + std::string(axis.shortName) + "_epp_SRC_Q2_" + std::to_string(q2)).c_str());
      for (int s = 0; s < kScanBins; ++s) {
        const double sigma = kSigmaMin + (kSigmaMax - kSigmaMin) * s / kScanBins;
        auto hScan = makeHist("scan_tmp", axis);
        fillMC(*hScan, mcEvents, sigmaGen, r, cfg, a, sigma, a);
        auto [scanScale, scanChi2] = bestScaleChi2(*hData, *hScan, axis);
        chi2.SetPoint(chi2.GetN(), sigma, scanChi2);
        scale.SetPoint(scale.GetN(), sigma, scanScale);
      }
      writeObject(&chi2);
      writeObject(&scale);

      if (integratedMode) {
        TGraphAsymmErrors g;
        g.SetName(axis.sigIntName);
        fillSigmaPoint(g, r, a, 1.0, 0.0, 2.0);
        writeObject(&g);
        writeGraphCanvas("c_chi2_" + std::string(axis.shortName) + "_epp", chi2);
        writeGraphCanvas("c_scale_" + std::string(axis.shortName) + "_epp", scale);
        writeOverlayCanvas("c_overlay_" + std::string(axis.shortName) + "_epp", *hData, *hBest);
      } else {
        writeGraphCanvas("c_chi2_" + std::string(axis.shortName) + "_Q2_" + std::to_string(q2), chi2);
        writeGraphCanvas("c_scale_" + std::string(axis.shortName) + "_Q2_" + std::to_string(q2), scale);
        writeOverlayCanvas("c_overlay_" + std::string(axis.shortName) + "_Q2_" + std::to_string(q2),
                           *hData, *hBest);
      }

      (void)bestChi2;
    }
  }

  for (int a = 0; a < 4; ++a) {
    writeObject(q2Graphs[a].get());
    writeGraphCanvas("c_" + std::string(axisDefs[a].sigQ2Name), *q2Graphs[a]);
  }
  file->Close();
}

}  // namespace sigmacm
