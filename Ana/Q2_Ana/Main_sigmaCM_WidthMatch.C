// Main_sigmaCM_WidthMatch.C
//
// Standalone stage-2 extractor for the Main_sigmaCM stage-1 ROOT output.
// It extracts input sigma_CM by fitting the reco-level data and each
// reco-level reweighted simulation template with the same width proxy, then
// interpolating the sigma_reco(sigma_CM input) map.
//
// Usage:
//   root -l -b -q 'Ana/Q2_Ana/Main_sigmaCM_WidthMatch.C("stage1.root","sigmaCM_widthmatch.root")'

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

const int linbin = 100;
const double min_sigma = 0.1;
const double max_sigma = 0.5;

const vector<double> bE_Q2 = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
const int bQ2 = (int)bE_Q2.size() - 1;

double sq(double x) { return x * x; }

double sCM(int j) {
  return (double(j) / double(linbin)) * (max_sigma - min_sigma) + min_sigma;
}

double G(double x, double N, double mu, double sigma) {
  return (N / (sigma * sqrt(2.0 * M_PI))) * exp(-0.5 * sq((x - mu) / sigma));
}

double Rayleigh(double x, double N, double sigma) {
  if (x < 0.0 || sigma <= 0.0) return 0.0;
  return N * (x / sq(sigma)) * exp(-0.5 * sq(x / sigma));
}

struct FitWidth {
  bool ok = false;
  double sigma = numeric_limits<double>::quiet_NaN();
  double err = numeric_limits<double>::quiet_NaN();
  double chi2ndf = numeric_limits<double>::quiet_NaN();
};

struct MatchResult {
  bool ok = false;
  int bestTemplate = -1;
  double sigmaCM = numeric_limits<double>::quiet_NaN();
  double errLow = numeric_limits<double>::quiet_NaN();
  double errHigh = numeric_limits<double>::quiet_NaN();
  double dataSigma = numeric_limits<double>::quiet_NaN();
  double dataSigmaErr = numeric_limits<double>::quiet_NaN();
};

struct Component {
  string label;
  string dataIntegrated;
  string simIntegratedPrefix;
  string dataQ2Prefix;
  string simQ2Prefix;
  string graphIntName;
  string graphQ2Name;
  double fitMin;
  double fitMax;
  bool transverse;
};

TH1D* getHist(TFile* f, const string& name, bool required = true) {
  TH1D* h = dynamic_cast<TH1D*>(f->Get(name.c_str()));
  if (!h && required) {
    cerr << "ERROR: missing histogram " << name << endl;
  }
  return h;
}

FitWidth fitWidth(TH1D* h, const string& fitName, double min, double max, bool transverse) {
  FitWidth out;
  if (!h || h->Integral() <= 0.0 || h->GetEntries() < 10) return out;

  // pcmT is a transverse magnitude and is not truly Gaussian. By default it
  // is still fit with G() as a consistent data/sim width proxy. If visual QA
  // shows those T fits are poor, flip this fallback on; it is applied
  // identically to data and simulation.
  const bool useTransverseFallback = false;

  TF1* f = nullptr;
  if (transverse && useTransverseFallback) {
    f = new TF1(fitName.c_str(), [](double* x, double* p) { return Rayleigh(x[0], p[0], p[1]); }, min, max, 2);
    f->SetParameter(0, h->Integral("width"));
    f->SetParameter(1, std::max(0.05, h->GetMean() / sqrt(M_PI / 2.0)));
    f->SetParLimits(1, 0.001, max - min);
  } else {
    f = new TF1(fitName.c_str(), [](double* x, double* p) { return G(x[0], p[0], p[1], p[2]); }, min, max, 3);
    const double seedSigma = std::max(0.05, (max - min) / 4.0);
    f->SetParameter(0, h->GetMaximum() / G(0.0, 1.0, 0.0, 0.1));
    f->SetParameter(1, transverse ? h->GetMean() : 0.0);
    f->SetParLimits(1, min, max);
    f->SetParameter(2, seedSigma);
    f->SetParLimits(2, 0.001, max - min);
  }

  TFitResultPtr fr = h->Fit(f, "SrBeqn", "", min, max);
  if (fr.Get() != nullptr && int(fr) == 0) {
    const int sigmaPar = (transverse && useTransverseFallback) ? 1 : 2;
    out.ok = true;
    out.sigma = std::abs(fr->Parameter(sigmaPar));
    out.err = std::abs(fr->ParError(sigmaPar));
    const int ndf = fr->Ndf();
    if (ndf > 0) out.chi2ndf = fr->Chi2() / double(ndf);
  }
  delete f;
  return out;
}

bool isMonotone(const vector<double>& y) {
  int direction = 0;
  for (size_t i = 1; i < y.size(); ++i) {
    if (!isfinite(y[i - 1]) || !isfinite(y[i])) continue;
    const double d = y[i] - y[i - 1];
    if (std::abs(d) < 1e-9) continue;
    if (direction == 0) direction = d > 0.0 ? 1 : -1;
    if (direction * d < -1e-9) return false;
  }
  return true;
}

double invertLinear(const vector<double>& x, const vector<double>& y, double target) {
  vector<pair<double, double>> pts;
  for (size_t i = 0; i < x.size(); ++i) {
    if (isfinite(x[i]) && isfinite(y[i])) pts.push_back({x[i], y[i]});
  }
  if (pts.size() < 2 || !isfinite(target)) return numeric_limits<double>::quiet_NaN();

  const bool increasing = pts.back().second >= pts.front().second;
  if (increasing) {
    if (target <= pts.front().second) return pts.front().first;
    if (target >= pts.back().second) return pts.back().first;
  } else {
    if (target >= pts.front().second) return pts.front().first;
    if (target <= pts.back().second) return pts.back().first;
  }

  for (size_t i = 1; i < pts.size(); ++i) {
    const double y0 = pts[i - 1].second;
    const double y1 = pts[i].second;
    const bool crossed = increasing ? (target >= y0 && target <= y1) : (target <= y0 && target >= y1);
    if (!crossed) continue;
    if (std::abs(y1 - y0) < 1e-12) return 0.5 * (pts[i - 1].first + pts[i].first);
    const double t = (target - y0) / (y1 - y0);
    return pts[i - 1].first + t * (pts[i].first - pts[i - 1].first);
  }
  return numeric_limits<double>::quiet_NaN();
}

int nearestTemplateBin(double sigmaCM) {
  int best = 0;
  double bestDist = numeric_limits<double>::max();
  for (int j = 0; j < linbin; ++j) {
    const double dist = std::abs(sCM(j) - sigmaCM);
    if (dist < bestDist) {
      bestDist = dist;
      best = j;
    }
  }
  return best;
}

double fitScale(TH1D* data, TH1D* sim, double min, double max) {
  if (!data || !sim || sim->Integral() <= 0.0) return 1.0;
  double num = 0.0;
  double den = 0.0;
  for (int b = data->FindBin(min); b <= data->FindBin(max); ++b) {
    const double d = data->GetBinContent(b);
    const double s = sim->GetBinContent(b);
    const double e2 = sq(data->GetBinError(b)) + sq(sim->GetBinError(b));
    if (e2 <= 0.0) continue;
    num += d * s / e2;
    den += s * s / e2;
  }
  return den > 0.0 ? num / den : 1.0;
}

MatchResult matchOne(TFile* outFile,
                     const Component& comp,
                     TH1D* data,
                     const vector<TH1D*>& sim,
                     const string& tag) {
  MatchResult result;
  if (!data) return result;

  vector<double> x;
  vector<double> y;
  vector<double> ey;
  auto* gMap = new TGraphErrors();
  gMap->SetName(Form("g_sigma_reco_vs_sCM_%s_%s", comp.label.c_str(), tag.c_str()));
  gMap->SetTitle(Form("%s %s;input #sigma_{CM};reco fit #sigma", comp.label.c_str(), tag.c_str()));

  for (int j = 0; j < linbin; ++j) {
    FitWidth fw = fitWidth(sim[j], Form("fit_%s_%s_sim_%d", comp.label.c_str(), tag.c_str(), j),
                           comp.fitMin, comp.fitMax, comp.transverse);
    if (!fw.ok) continue;
    x.push_back(sCM(j));
    y.push_back(fw.sigma);
    ey.push_back(fw.err);
    gMap->SetPoint(gMap->GetN(), sCM(j), fw.sigma);
    gMap->SetPointError(gMap->GetN() - 1, 0.0, fw.err);
  }
  outFile->cd();
  gMap->Write();

  if (!isMonotone(y)) {
    cerr << "WARNING: sigma_reco(sCM) is not monotone for " << comp.label << " " << tag
         << "; linear inversion will use the first crossing." << endl;
  }

  FitWidth dataFit = fitWidth(data, Form("fit_%s_%s_data", comp.label.c_str(), tag.c_str()),
                              comp.fitMin, comp.fitMax, comp.transverse);
  if (!dataFit.ok) {
    cerr << "WARNING: data fit failed for " << comp.label << " " << tag << endl;
    return result;
  }

  const double center = invertLinear(x, y, dataFit.sigma);
  const double lo = invertLinear(x, y, dataFit.sigma - dataFit.err);
  const double hi = invertLinear(x, y, dataFit.sigma + dataFit.err);
  if (!isfinite(center) || !isfinite(lo) || !isfinite(hi)) return result;

  result.ok = true;
  result.sigmaCM = center;
  result.errLow = std::abs(center - lo);
  result.errHigh = std::abs(hi - center);
  result.dataSigma = dataFit.sigma;
  result.dataSigmaErr = dataFit.err;
  result.bestTemplate = nearestTemplateBin(center);
  return result;
}

TH1D* matchedTemplateClone(const Component& comp,
                           TH1D* data,
                           TH1D* sim,
                           const string& name) {
  if (!sim) return nullptr;
  TH1D* clone = (TH1D*)sim->Clone(name.c_str());
  clone->SetDirectory(nullptr);
  clone->Scale(fitScale(data, sim, comp.fitMin, comp.fitMax));
  clone->SetLineColor(kRed + 1);
  clone->SetLineWidth(3);
  return clone;
}

void drawOverlay(TFile* outFile,
                 TCanvas* canvas,
                 const char* cname,
                 const char* ctitle,
                 TH1D* data,
                 TH1D* sim) {
  if (!data || !sim) return;
  canvas->Clear();
  canvas->Divide(1, 1);
  canvas->cd(1);
  data->SetLineColor(kBlack);
  data->SetMarkerStyle(20);
  data->Draw("E");
  sim->Draw("HIST SAME");
  data->Draw("E SAME");

  outFile->cd();
  TCanvas* cOut = new TCanvas(cname, ctitle, canvas->GetWw(), canvas->GetWh());
  canvas->DrawClonePad();
  cOut->Write();
  delete cOut;
}

int Main_sigmaCM_WidthMatch(const char* inFileName = nullptr, const char* outFileName = nullptr) {
  if (!inFileName || !outFileName) {
    cerr << "Usage:\n"
         << "  root -l -b -q 'Ana/Q2_Ana/Main_sigmaCM_WidthMatch.C+(\"stage1.root\",\"sigmaCM_widthmatch.root\")'\n";
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile* inFile = TFile::Open(inFileName, "READ");
  if (!inFile || inFile->IsZombie()) {
    cerr << "ERROR: could not open input file " << inFileName << endl;
    return 1;
  }

  TFile* outFile = TFile::Open(outFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    cerr << "ERROR: could not open output file " << outFileName << endl;
    return 1;
  }

  vector<TH1D*> hQ2(bQ2, nullptr);
  for (int i = 0; i < bQ2; ++i) {
    hQ2[i] = getHist(inFile, Form("h_Q2_%d", i));
  }

  vector<Component> comps = {
    {"pcmx", "pcmx_epp", "h_pcmx_epp_simSCM", "h_pcmx_epp_SRC_Q2", "h_pcmx_epp_SRC_simSCM_Q2",
     "sigmacmx_int", "sigmacmx_Q2", -0.2, 0.2, false},
    {"pcmy", "pcmy_epp", "h_pcmy_epp_simSCM", "h_pcmy_epp_SRC_Q2", "h_pcmy_epp_SRC_simSCM_Q2",
     "sigmacmy_int", "sigmacmy_Q2", -0.2, 0.2, false},
    {"pcmz", "pcmz_epp", "h_pcmz_epp_simSCM", "h_pcmz_epp_SRC_Q2", "h_pcmz_epp_SRC_simSCM_Q2",
     "sigmacmz_int", "sigmacmz_Q2", -0.2, 0.2, false},
    {"pcmT", "pcmT_epp", "h_pcmT_epp_simSCM", "h_pcmT_epp_SRC_Q2", "h_pcmT_epp_SRC_simSCM_Q2",
     "sigmacmT_int", "sigmacmT_Q2", 0.0, 0.5, true},
  };

  const int pixelx = 1980;
  const int pixely = 1530;
  TCanvas* myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);

  for (const Component& comp : comps) {
    TH1D* hDataInt = getHist(inFile, comp.dataIntegrated);
    vector<TH1D*> hSimInt(linbin, nullptr);
    for (int j = 0; j < linbin; ++j) {
      hSimInt[j] = getHist(inFile, Form("%s_%d", comp.simIntegratedPrefix.c_str(), j));
    }

    MatchResult intResult = matchOne(outFile, comp, hDataInt, hSimInt, "int");
    auto* gInt = new TGraphAsymmErrors();
    gInt->SetName(comp.graphIntName.c_str());
    if (intResult.ok) {
      gInt->SetPoint(0, 1.0, intResult.sigmaCM);
      gInt->SetPointError(0, 0.0, 2.0, intResult.errLow, intResult.errHigh);
    }
    outFile->cd();
    if (hDataInt) hDataInt->Write();
    gInt->Write();

    TH1D* hFitInt = nullptr;
    if (intResult.ok) {
      hFitInt = matchedTemplateClone(comp, hDataInt, hSimInt[intResult.bestTemplate],
                                     Form("%s_epp_fit", comp.label.c_str()));
      if (hFitInt) hFitInt->Write();
      drawOverlay(outFile, myCanvas, Form("c_overlay_%s_epp", comp.label.c_str()),
                  Form("data+matched sim overlay - %s integrated", comp.label.c_str()),
                  hDataInt, hFitInt);
    }

    auto* gQ2 = new TGraphAsymmErrors();
    gQ2->SetName(comp.graphQ2Name.c_str());
    vector<MatchResult> q2Results(bQ2);

    for (int i = 0; i < bQ2; ++i) {
      TH1D* hDataQ2 = getHist(inFile, Form("%s_%d", comp.dataQ2Prefix.c_str(), i));
      vector<TH1D*> hSimQ2(linbin, nullptr);
      for (int j = 0; j < linbin; ++j) {
        hSimQ2[j] = getHist(inFile, Form("%s_%d_%d", comp.simQ2Prefix.c_str(), j, i));
      }

      q2Results[i] = matchOne(outFile, comp, hDataQ2, hSimQ2, Form("Q2_%d", i));
      if (!q2Results[i].ok || !hQ2[i]) continue;

      const double q2Center = hQ2[i]->GetMean();
      const double q2Low = q2Center - bE_Q2[i];
      const double q2High = bE_Q2[i + 1] - q2Center;
      const int p = gQ2->GetN();
      gQ2->SetPoint(p, q2Center, q2Results[i].sigmaCM);
      gQ2->SetPointError(p, q2Low, q2High, q2Results[i].errLow, q2Results[i].errHigh);

      TH1D* hFitQ2 = matchedTemplateClone(comp, hDataQ2, hSimQ2[q2Results[i].bestTemplate],
                                          Form("%s_epp_SRC_Q2_%d_fit", comp.label.c_str(), i));
      if (hFitQ2) {
        outFile->cd();
        hFitQ2->Write();
        drawOverlay(outFile, myCanvas, Form("c_overlay_%s_Q2_%d", comp.label.c_str(), i),
                    Form("data+matched sim %s Q2bin=%d", comp.label.c_str(), i),
                    hDataQ2, hFitQ2);
      }
    }

    outFile->cd();
    gQ2->Write();
    myCanvas->Clear();
    myCanvas->Divide(1, 1);
    myCanvas->cd(1);
    gQ2->Draw("AP");
    TCanvas* cQ2 = new TCanvas(Form("c_%s", comp.graphQ2Name.c_str()),
                               Form("%s vs Q2", comp.graphQ2Name.c_str()), pixelx, pixely);
    myCanvas->DrawClonePad();
    cQ2->Write();
    delete cQ2;
  }

  outFile->Close();
  inFile->Close();
  cout << "Wrote width-matched sigma_CM output to " << outFileName << endl;
  return 0;
}
