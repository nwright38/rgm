#include "SigmaCMExtractor.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace sigmacm {
namespace {

struct AxisDef {
  int bins = 0;
  double min = 0.0;
  double max = 0.0;
  double fitMin = 0.0;
  double fitMax = 0.0;
};

struct Hist {
  std::vector<double> sum;
  std::vector<double> sumw2;
};

std::array<AxisDef, 3> axes(const Config& cfg) {
  return {{{cfg.binsX, cfg.histXMin, cfg.histXMax, -cfg.cutRangeXY, cfg.cutRangeXY},
           {cfg.binsY, cfg.histYMin, cfg.histYMax, -cfg.cutRangeXY, cfg.cutRangeXY},
           {cfg.binsZ, cfg.histZMin, cfg.histZMax, cfg.fitZMin, cfg.fitZMax}}};
}

double component(const Event& e, int axis) {
  if (axis == 0) return e.pcmX;
  if (axis == 1) return e.pcmY;
  return e.pcmZ;
}

int binFor(double x, const AxisDef& a) {
  if (x < a.min || x >= a.max) return -1;
  const int b = static_cast<int>((x - a.min) / (a.max - a.min) * a.bins);
  return (b >= 0 && b < a.bins) ? b : -1;
}

double binCenter(int bin, const AxisDef& a) {
  return a.min + (bin + 0.5) * (a.max - a.min) / a.bins;
}

void fill(Hist& h, const AxisDef& a, double value, double weight) {
  const int b = binFor(value, a);
  if (b < 0) return;
  h.sum[b] += weight;
  h.sumw2[b] += weight * weight;
}

struct Prepared {
  std::vector<const Event*> data;
  std::vector<const Event*> mc;
};

Prepared prepare(const std::vector<Event>& dataEvents, const std::vector<Event>& mcEvents, const Config& cfg) {
  Prepared p;
  for (const auto& e : dataEvents) {
    if (passesCuts(e, cfg)) p.data.push_back(&e);
  }
  for (const auto& e : mcEvents) {
    if (passesCuts(e, cfg)) p.mc.push_back(&e);
  }
  return p;
}

class Chi2Model {
 public:
  Chi2Model(const Prepared& prepared, double sigmaGen, const Config& cfg, bool fixedAxis = false,
            int axis = -1, double fixedSigma = 0.0)
      : prepared_(prepared), sigmaGen_(sigmaGen), cfg_(cfg), axes_(axes(cfg)),
        fixedAxis_(fixedAxis), fixedAxisIndex_(axis), fixedSigma_(fixedSigma) {
    for (int a = 0; a < 3; ++a) {
      data_[a].sum.assign(axes_[a].bins, 0.0);
      data_[a].sumw2.assign(axes_[a].bins, 0.0);
    }
    for (const Event* e : prepared_.data) {
      const double w = cfg_.useBootstrapWeights ? e->bootstrapWeight : 1.0;
      for (int a = 0; a < 3; ++a) fill(data_[a], axes_[a], component(*e, a), w);
    }
  }

  unsigned int nPars() const { return (fixedAxis_ ? 2 : 3) + (cfg_.sharedScale ? 1 : 3); }

  double operator()(const double* pars) const {
    std::array<double, 3> sigma{};
    int pi = 0;
    for (int a = 0; a < 3; ++a) sigma[a] = (fixedAxis_ && a == fixedAxisIndex_) ? fixedSigma_ : pars[pi++];
    std::array<double, 3> scale{};
    if (cfg_.sharedScale) {
      scale.fill(pars[pi]);
    } else {
      for (int a = 0; a < 3; ++a) scale[a] = pars[pi + a];
    }
    if (sigma[0] <= 0.02 || sigma[1] <= 0.02 || sigma[2] <= 0.02 ||
        scale[0] <= 0.0 || scale[1] <= 0.0 || scale[2] <= 0.0) {
      return 1.0e30;
    }
    auto mc = buildMC(sigma);
    return chi2For(mc, scale, nullptr);
  }

  std::array<Hist, 3> buildMC(const std::array<double, 3>& sigma) const {
    std::array<Hist, 3> mc;
    for (int a = 0; a < 3; ++a) {
      mc[a].sum.assign(axes_[a].bins, 0.0);
      mc[a].sumw2.assign(axes_[a].bins, 0.0);
    }
    for (const Event* e : prepared_.mc) {
      const double w = eventWeight(*e, sigma[0], sigma[1], sigma[2], sigmaGen_, cfg_.auxWeightBranch);
      for (int a = 0; a < 3; ++a) fill(mc[a], axes_[a], component(*e, a), w);
    }
    return mc;
  }

  double chi2For(const std::array<Hist, 3>& mc, const std::array<double, 3>& scale,
                 std::array<double, 3>* perProjection) const {
    double total = 0.0;
    if (perProjection) perProjection->fill(0.0);
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < axes_[a].bins; ++b) {
        const double c = binCenter(b, axes_[a]);
        if (c < axes_[a].fitMin || c > axes_[a].fitMax) continue;
        const double d = data_[a].sum[b];
        const double e = mc[a].sum[b];
        const double var = data_[a].sumw2[b] + scale[a] * scale[a] * mc[a].sumw2[b];
        if (var <= 0.0) continue;
        const double contribution = (d - scale[a] * e) * (d - scale[a] * e) / var;
        total += contribution;
        if (perProjection) (*perProjection)[a] += contribution;
      }
    }
    return total;
  }

  int nFitBins() const {
    int n = 0;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < axes_[a].bins; ++b) {
        const double c = binCenter(b, axes_[a]);
        if (c >= axes_[a].fitMin && c <= axes_[a].fitMax) ++n;
      }
    }
    return n;
  }

  int nScalePars() const { return cfg_.sharedScale ? 1 : 3; }

  double effectiveMC(const std::array<double, 3>& sigma) const {
    double sw = 0.0;
    double sw2 = 0.0;
    for (const Event* e : prepared_.mc) {
      const double w = eventWeight(*e, sigma[0], sigma[1], sigma[2], sigmaGen_, cfg_.auxWeightBranch);
      sw += w;
      sw2 += w * w;
    }
    return sw2 > 0.0 ? sw * sw / sw2 : 0.0;
  }

 private:
  const Prepared& prepared_;
  double sigmaGen_ = 0.0;
  Config cfg_;
  std::array<AxisDef, 3> axes_;
  bool fixedAxis_ = false;
  int fixedAxisIndex_ = -1;
  double fixedSigma_ = 0.0;
  std::array<Hist, 3> data_;
};

std::unique_ptr<ROOT::Math::Minimizer> makeMinimizer() {
  std::unique_ptr<ROOT::Math::Minimizer> minimizer(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (!minimizer) throw std::runtime_error("Could not create ROOT Minuit2 minimizer");
  minimizer->SetMaxFunctionCalls(20000);
  minimizer->SetMaxIterations(20000);
  minimizer->SetTolerance(1e-3);
  minimizer->SetPrintLevel(0);
  return minimizer;
}

Result minimizeModel(const Chi2Model& model, const Config& cfg, bool fixedAxis = false,
                     int fixedAxisIndex = -1, double fixedSigma = 0.0,
                     bool computeErrors = true) {
  ROOT::Math::Functor f(model, model.nPars());
  auto minimizer = makeMinimizer();
  minimizer->SetFunction(f);

  int pi = 0;
  for (int a = 0; a < 3; ++a) {
    if (fixedAxis && a == fixedAxisIndex) continue;
    const char* name = (a == 0) ? "sigmaX" : (a == 1) ? "sigmaY" : "sigmaZ";
    minimizer->SetLimitedVariable(pi++, name, 0.16, 0.003, 0.04, 0.40);
  }
  const int firstScaleIndex = pi;
  if (cfg.sharedScale) {
    minimizer->SetLimitedVariable(pi++, "scale", 1.0, 0.01, 1e-6, 1e6);
  } else {
    minimizer->SetLimitedVariable(pi++, "scaleX", 1.0, 0.01, 1e-6, 1e6);
    minimizer->SetLimitedVariable(pi++, "scaleY", 1.0, 0.01, 1e-6, 1e6);
    minimizer->SetLimitedVariable(pi++, "scaleZ", 1.0, 0.01, 1e-6, 1e6);
  }
  const bool ok = minimizer->Minimize();
  if (computeErrors) minimizer->Hesse();

  const double* xs = minimizer->X();
  const double* errs = minimizer->Errors();
  Result r;
  r.config = cfg;
  r.converged = ok;
  r.status = ok ? "OK" : "Minuit2 did not report convergence";
  pi = 0;
  std::array<double, 3> sigma{};
  for (int a = 0; a < 3; ++a) {
    if (fixedAxis && a == fixedAxisIndex) sigma[a] = fixedSigma;
    else sigma[a] = xs[pi++];
  }
  const int scaleIndex = firstScaleIndex;
  std::array<double, 3> scale{};
  std::array<double, 3> scaleErr{};
  if (cfg.sharedScale) {
    scale.fill(xs[firstScaleIndex]);
    scaleErr.fill(errs ? errs[firstScaleIndex] : 0.0);
  } else {
    scale = {{xs[firstScaleIndex], xs[firstScaleIndex + 1], xs[firstScaleIndex + 2]}};
    scaleErr = {{errs ? errs[firstScaleIndex] : 0.0, errs ? errs[firstScaleIndex + 1] : 0.0,
                 errs ? errs[firstScaleIndex + 2] : 0.0}};
  }
  r.sigmaX = sigma[0];
  r.sigmaY = sigma[1];
  r.sigmaZ = sigma[2];
  r.scale = (scale[0] + scale[1] + scale[2]) / 3.0;
  r.scaleErr = std::sqrt(scaleErr[0] * scaleErr[0] + scaleErr[1] * scaleErr[1] +
                         scaleErr[2] * scaleErr[2]) / 3.0;
  auto mc = model.buildMC(sigma);
  r.chi2 = model.chi2For(mc, scale, &r.projectionChi2);
  r.ndf = model.nFitBins() - static_cast<int>(model.nPars());
  r.effectiveMCEntries = model.effectiveMC(sigma);

  auto covariance = [&](int i, int j) {
    if (!computeErrors) return 0.0;
    return minimizer->CovMatrix(i, j);
  };
  auto corr = [&](int i, int j) {
    const double vii = covariance(i, i);
    const double vjj = covariance(j, j);
    if (vii <= 0.0 || vjj <= 0.0) return 0.0;
    return covariance(i, j) / std::sqrt(vii * vjj);
  };

  int sxIndex = fixedAxis && fixedAxisIndex == 0 ? -1 : 0;
  int syIndex = fixedAxis ? (fixedAxisIndex == 0 ? 0 : (fixedAxisIndex == 1 ? -1 : 1)) : 1;
  int szIndex = fixedAxis ? (fixedAxisIndex == 2 ? -1 : 2) : 2;
  auto errFor = [&](int idx) { return (computeErrors && idx >= 0 && errs) ? errs[idx] : 0.0; };
  r.sigmaXErrLow = r.sigmaXErrHigh = errFor(sxIndex);
  r.sigmaYErrLow = r.sigmaYErrHigh = errFor(syIndex);
  r.sigmaZErrLow = r.sigmaZErrHigh = errFor(szIndex);
  r.rhoScaleSigmaX = (cfg.sharedScale && sxIndex >= 0) ? corr(scaleIndex, sxIndex) : 0.0;
  r.rhoScaleSigmaY = (cfg.sharedScale && syIndex >= 0) ? corr(scaleIndex, syIndex) : 0.0;
  r.rhoScaleSigmaZ = (cfg.sharedScale && szIndex >= 0) ? corr(scaleIndex, szIndex) : 0.0;
  return r;
}

}  // namespace

double gaussianRatio(double p, double sigmaTarget, double sigmaGen) {
  if (sigmaTarget <= 0.0 || sigmaGen <= 0.0) {
    throw std::runtime_error("Gaussian widths must be positive");
  }
  return (sigmaGen / sigmaTarget) *
         std::exp(-0.5 * p * p * (1.0 / (sigmaTarget * sigmaTarget) -
                                   1.0 / (sigmaGen * sigmaGen)));
}

double eventWeight(const Event& event, double sigmaX, double sigmaY, double sigmaZ,
                   double sigmaGen, const std::string& auxWeightBranch) {
  double aux = 1.0;
  if (!auxWeightBranch.empty()) {
    auto it = event.auxWeights.find(auxWeightBranch);
    if (it == event.auxWeights.end()) {
      throw std::runtime_error("Requested auxiliary MC weight branch '" + auxWeightBranch +
                               "' was not loaded");
    }
    aux = it->second;
  }
  return event.genWeight * aux * gaussianRatio(event.pcmXTruth, sigmaX, sigmaGen) *
         gaussianRatio(event.pcmYTruth, sigmaY, sigmaGen) *
         gaussianRatio(event.pcmZTruth, sigmaZ, sigmaGen);
}

bool passesCuts(const Event& e, const Config& cfg) {
  if (!(e.xB > cfg.xBLower)) return false;
  if (!(e.Q2 > cfg.q2Lower && e.Q2 < cfg.q2Upper)) return false;
  if (!(e.mMiss > cfg.mMissLower && e.mMiss < cfg.mMissUpper)) return false;
  if (!(e.kMiss > cfg.kMissLower)) return false;
  if (!(e.pLead > cfg.pLeadLower)) return false;
  if (!(e.hasRecoil && e.pRec > cfg.pRecLower)) return false;
  if (cfg.requirePcmLtPrel && !(e.pCM < e.pRel)) return false;
  if (cfg.leadMode == LeadMode::FD) {
    if (!(e.leadRegion == cfg.fdLeadRegionValue && e.thetaLead < cfg.thetaFDUpper)) return false;
  } else if (cfg.leadMode == LeadMode::CD) {
    if (!(e.leadRegion == cfg.cdLeadRegionValue && e.thetaLead > cfg.thetaCDLower)) return false;
  } else {
    const bool passFD = e.leadRegion == cfg.fdLeadRegionValue && e.thetaLead < cfg.thetaFDUpper;
    const bool passCD = e.leadRegion == cfg.cdLeadRegionValue && e.thetaLead > cfg.thetaCDLower;
    if (!(passFD || passCD)) return false;
  }
  if (!cfg.integratedQ2) {
    if (cfg.q2BinIndex < 0 || cfg.q2BinIndex + 1 >= static_cast<int>(cfg.q2Edges.size())) {
      throw std::runtime_error("Invalid q2BinIndex in Config");
    }
    if (!(e.Q2 >= cfg.q2Edges[cfg.q2BinIndex] && e.Q2 < cfg.q2Edges[cfg.q2BinIndex + 1])) {
      return false;
    }
  }
  return true;
}

Result extract(const std::vector<Event>& dataEvents, const std::vector<Event>& mcEvents,
               double sigmaGen, const Config& cfg) {
  if (sigmaGen <= 0.0) throw std::runtime_error("sigma_gen metadata must be positive");
  Config local = cfg;
  const Prepared prepared = prepare(dataEvents, mcEvents, local);
  if (prepared.data.empty()) throw std::runtime_error("No data events pass Config cuts");
  if (prepared.mc.empty()) throw std::runtime_error("No MC events pass Config cuts");
  Chi2Model model(prepared, sigmaGen, local);
  Result r = minimizeModel(model, local);
  r.nEventsData = static_cast<long long>(prepared.data.size());
  return r;
}

Result extractProfileScan(const std::vector<Event>& dataEvents, const std::vector<Event>& mcEvents,
                          double sigmaGen, const Config& cfg, int axis,
                          double scanMin, double scanMax, int nPoints) {
  if (sigmaGen <= 0.0) throw std::runtime_error("sigma_gen metadata must be positive");
  if (axis < 0 || axis > 2) throw std::runtime_error("Profile scan axis must be 0, 1, or 2");
  Config local = cfg;
  const Prepared prepared = prepare(dataEvents, mcEvents, local);
  if (prepared.data.empty()) throw std::runtime_error("No data events pass Config cuts");
  if (prepared.mc.empty()) throw std::runtime_error("No MC events pass Config cuts");
  Chi2Model nominal(prepared, sigmaGen, local);
  Result base = minimizeModel(nominal, local, false, -1, 0.0, false);
  base.nEventsData = static_cast<long long>(prepared.data.size());
  base.profileAxis = axis;
  base.profile.clear();
  for (int i = 0; i < nPoints; ++i) {
    const double s = scanMin + (scanMax - scanMin) * i / std::max(1, nPoints - 1);
    try {
      Chi2Model fixed(prepared, sigmaGen, local, true, axis, s);
      Result rr = minimizeModel(fixed, local, true, axis, s, false);
      base.profile.push_back({s, rr.chi2, rr.sigmaX, rr.sigmaY, rr.sigmaZ, rr.scale});
    } catch (const std::exception& e) {
      base.status += "; profile point failed at sigma=" + std::to_string(s) + ": " + e.what();
    }
  }
  if (base.profile.empty()) throw std::runtime_error("All profile scan points failed");
  return base;
}

}  // namespace sigmacm
