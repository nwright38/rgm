// simple.C
//
// Q2-binned 1D shape overlays from a TTree.
//
// Usage examples:
//   root -l -b -q 'scratch/simple.C'
//   root -l -b -q 'scratch/simple.C("~/data/RGM_DATA/c12_src_skim.root")'
//   root -l -b -q 'scratch/simple.C("input.root","srcTree","my_q2_overlays.pdf","xB>1.2")'
//
// Edit variablesToPlot() below to add/remove variables.  For each variable, the
// PDF gets one e'p page and one e'pp page.  Every Q2-bin histogram is scaled to
// unit area after filling.  The tree is read once and all histograms are filled
// in that single event loop.

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeFormula.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;

int ep_bins = 40;
int epp_bins = 20;

namespace {

const vector<double> bE_Q2 = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
const vector<double> bE_pmiss = {0.4,0.45,0.5,0.55,0.6,0.65,0.75,0.85,0.95,1.05,1.2};

struct PlotVar {
  string name;     // short name used in ROOT object names
  string expr;     // TTreeFormula expression, e.g. "pMiss" or "leadTheta*180./TMath::Pi()"
  string xTitle;   // axis title
  double xmin;
  double xmax;
  bool drawEp;
  bool drawEpp;
};

vector<PlotVar> variablesToPlot() {
  return {
      {"pMiss", "pMiss", "p_{miss} [GeV]", 0.3, 1.1, true, true},
      {"kMiss", "kMiss", "k_{miss} [GeV]", 0.2, 1.2, true, true},
      {"xB", "xB", "x_{B}", 0.8, 2.2, true, true},
      {"leadP", "leadP", "p_{lead} [GeV]", 0.8, 4.0, true, true},
      {"leadTheta", "leadTheta*180./TMath::Pi()", "#theta_{lead} [deg]", 0.0, 45.0, true, true},

      // Pair/recoil quantities: e'pp only.
      {"pCMx", "pCMx", "p_{CM,x} [GeV]", -0.8, 0.8, false, true},
      {"pCMy", "pCMy", "p_{CM,y} [GeV]", -0.8, 0.8, false, true},
      {"pCMz", "pCMz", "p_{CM,z} [GeV]", -0.8, 0.8, false, true},
      {"pCM", "pCM", "|p_{CM}| [GeV]", 0.0, 1.2, false, true},
      {"pRel", "pRel", "|p_{rel}| [GeV]", 0.0, 1.2, false, true},
      {"recP", "recP", "p_{rec} [GeV]", 0.28, 1.2, false, true},
      // Examples if/when you want recoil plots:
      // {"recP", "recP", "p_{rec} [GeV]", 0.0, 1.2, false, true},
      // {"recTheta", "recTheta*180./TMath::Pi()", "#theta_{rec} [deg]", 0.0, 180.0, false, true},
  };
}

struct Channel {
  string key;
  string label;
  string cutExpr;
  string weightExpr;
};

class Formula {
 public:
  Formula(TTree* tree, const string& name, const string& expr)
      : expr_(expr), formula_(nullptr) {
    if (!expr_.empty()) formula_ = new TTreeFormula(name.c_str(), expr_.c_str(), tree);
  }

  ~Formula() { delete formula_; }

  bool ok() const { return expr_.empty() || (formula_ && formula_->GetNdim() > 0); }
  const string& expr() const { return expr_; }

  double eval(double emptyValue = 1.0) const {
    if (!formula_) return emptyValue;
    return formula_->EvalInstance();
  }

 private:
  string expr_;
  TTreeFormula* formula_;
};

int findQ2Bin(double q2) {
  for (int q = 0; q < (int)bE_Q2.size() - 1; q++) {
    if (q2 >= bE_Q2[q] && q2 < bE_Q2[q + 1]) return q;
  }
  return -1;
}

bool passes(const Formula& f) {
  return f.expr().empty() || f.eval(1.0) != 0.0;
}

double weight(const Formula& f) {
  return f.eval(1.0);
}

int binsForChannel(const PlotVar& v, const Channel& channel) {
  if (channel.key == "ep" && v.drawEp) return ep_bins;
  if (channel.key == "epp" && v.drawEpp) return epp_bins;
  return epp_bins;
}

void styleHist(TH1D* h, int q2Bin) {
  const int colors[] = {
      kBlack, kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1, kCyan + 2};
  const int styles[] = {1, 1, 1, 1, 1, 1, 1};

  h->SetDirectory(nullptr);
  h->SetStats(0);
  h->SetLineColor(colors[q2Bin % 7]);
  h->SetMarkerColor(colors[q2Bin % 7]);
  h->SetLineStyle(styles[q2Bin % 7]);
  h->SetLineWidth(3);
}

double unitAreaScale(TH1D* h) {
  if (!h) return 0.0;
  const double area = h->Integral("width");
  if (area > 0.0) return 1.0 / area;

  const double counts = h->Integral();
  if (counts > 0.0) return 1.0 / counts;
  return 0.0;
}

int findPlotVarIndex(const vector<PlotVar>& vars, const string& name) {
  for (int i = 0; i < (int)vars.size(); i++) {
    if (vars[i].name == name) return i;
  }
  return -1;
}

void drawMissingPage(TCanvas* c, const string& outPdf, const string& message) {
  c->Clear();
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.035);
  t.DrawLatex(0.10, 0.58, message.c_str());
  c->Print(outPdf.c_str(), "pdf");
}

void drawOverlayPage(TCanvas* c,
                     const string& outPath,
                     const PlotVar& v,
                     const Channel& channel,
                     const vector<vector<vector<TH1D*>>>& hists,
                     int channelIndex,
                     int varIndex,
                     const string& baseCut) {
  double ymax = 0.0;
  for (int q = 0; q < (int)bE_Q2.size() - 1; q++) {
    TH1D* h = hists[channelIndex][varIndex][q];
    const double scale = unitAreaScale(h);
    if (scale > 0.0) h->Scale(scale);
    ymax = std::max(ymax, h->GetMaximum());
  }

  c->Clear();
  if (ymax <= 0.0) {
    drawMissingPage(c, outPath, "No entries found for " + channel.label + " " + v.name + ".");
    return;
  }

  TLegend* leg = new TLegend(0.58, 0.52, 0.93, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.027);

  for (int q = 0; q < (int)bE_Q2.size() - 1; q++) {
    TH1D* h = hists[channelIndex][varIndex][q];
    h->GetXaxis()->SetTitle(v.xTitle.c_str());
    h->GetYaxis()->SetTitle("Unit-area density");
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.038);
    h->GetYaxis()->SetLabelSize(0.038);
    h->GetYaxis()->SetRangeUser(0.0, 1.25 * ymax);

    const string drawOpt = (q == 0) ? "PE" : "PE SAME";
    h->Draw(drawOpt.c_str());
    leg->AddEntry(h, Form("%.2f < Q^{2} < %.2f", bE_Q2[q], bE_Q2[q + 1]), "l");
  }

  TLatex title;
  title.SetNDC();
  title.SetTextFont(62);
  title.SetTextSize(0.038);
  title.DrawLatex(0.12, 0.94,
                  Form("%s: %s by Q^{2} bin", channel.label.c_str(), v.xTitle.c_str()));

  TLatex note;
  note.SetNDC();
  note.SetTextSize(0.025);
  note.DrawLatex(0.12, 0.89, "Each curve normalized to unit area");
  if (!baseCut.empty()) note.DrawLatex(0.12, 0.85, Form("Common cut: %s", baseCut.c_str()));
  if (!channel.cutExpr.empty()) note.DrawLatex(0.12, 0.81, Form("%s cut: %s", channel.key.c_str(), channel.cutExpr.c_str()));
  if (!channel.weightExpr.empty()) note.DrawLatex(0.12, 0.77, Form("%s weight: %s", channel.key.c_str(), channel.weightExpr.c_str()));

  leg->Draw();
  c->Print(outPath.c_str(), "pdf");
}

void drawRatioOverlayPage(TCanvas* c,
                          const string& outPath,
                          const vector<TH1D*>& hEppPmiss,
                          const vector<TH1D*>& hEpPmiss,
                          const string& baseCut,
                          const Channel& epChannel,
                          const Channel& eppChannel) {
  c->Clear();

  vector<TH1D*> ratios;
  ratios.reserve(bE_Q2.size() - 1);
  double ymax = 0.0;

  for (int q = 0; q < (int)bE_Q2.size() - 1; q++) {
    if (!hEppPmiss[q] || !hEpPmiss[q] ||
        hEppPmiss[q]->Integral() <= 0.0 || hEpPmiss[q]->Integral() <= 0.0) {
      ratios.push_back(nullptr);
      continue;
    }

    TH1D* hRatio = (TH1D*)hEppPmiss[q]->Clone(Form("h_ratio_epp_over_ep_pMiss_Q2_%d", q));
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hEppPmiss[q], hEpPmiss[q], 1.0, 1.0);
    styleHist(hRatio, q);
    hRatio->SetMarkerStyle(20 + q);
    hRatio->SetMarkerSize(0.85);
    hRatio->SetLineWidth(2);
    hRatio->GetXaxis()->SetTitle("p_{miss} [GeV]");
    hRatio->GetYaxis()->SetTitle("e'pp / e'p");
    hRatio->GetXaxis()->SetTitleSize(0.045);
    hRatio->GetYaxis()->SetTitleSize(0.045);
    hRatio->GetXaxis()->SetLabelSize(0.038);
    hRatio->GetYaxis()->SetLabelSize(0.038);

    for (int b = 1; b <= hRatio->GetNbinsX(); b++) {
      ymax = std::max(ymax, hRatio->GetBinContent(b) + hRatio->GetBinError(b));
    }
    ratios.push_back(hRatio);
  }

  if (ymax <= 0.0) {
    drawMissingPage(c, outPath, "No valid Q^{2} bins found for e'pp/e'p vs p_{miss}.");
    return;
  }

  TLegend* leg = new TLegend(0.58, 0.52, 0.93, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.027);

  bool drewFirst = false;
  TH1D* axisHist = nullptr;
  for (int q = 0; q < (int)ratios.size(); q++) {
    TH1D* hRatio = ratios[q];
    if (!hRatio) continue;
    hRatio->GetYaxis()->SetRangeUser(0.0, 0.5);
    hRatio->Draw(drewFirst ? "PE SAME" : "PE");
    drewFirst = true;
    if (!axisHist) axisHist = hRatio;
    leg->AddEntry(hRatio, Form("%.2f < Q^{2} < %.2f", bE_Q2[q], bE_Q2[q + 1]), "lep");
  }

  TLine* zero = new TLine(axisHist->GetXaxis()->GetXmin(), 0.0,
                          axisHist->GetXaxis()->GetXmax(), 0.0);
  zero->SetLineColor(kGray + 2);
  zero->SetLineStyle(2);
  zero->Draw("SAME");

  TLatex title;
  title.SetNDC();
  title.SetTextFont(62);
  title.SetTextSize(0.038);
  title.DrawLatex(0.12, 0.94, "e'pp / e'p vs p_{miss} by Q^{2} bin");

  TLatex note;
  note.SetNDC();
  note.SetTextSize(0.025);
  note.DrawLatex(0.12, 0.89, "One ratio curve per Q^{2} bin");
  if (!baseCut.empty()) note.DrawLatex(0.12, 0.85, Form("Common cut: %s", baseCut.c_str()));
  if (!epChannel.weightExpr.empty()) note.DrawLatex(0.12, 0.81, Form("e'p weight: %s", epChannel.weightExpr.c_str()));
  if (!eppChannel.cutExpr.empty()) note.DrawLatex(0.12, 0.77, Form("e'pp cut: %s", eppChannel.cutExpr.c_str()));
  if (!eppChannel.weightExpr.empty()) note.DrawLatex(0.12, 0.73, Form("e'pp weight: %s", eppChannel.weightExpr.c_str()));

  leg->Draw();
  c->Print(outPath.c_str(), "pdf");
}

}  // namespace

void simple(const char* inputFile = "~/data/RGM_DATA/c12_src_skim.root",
            const char* treeName = "srcTree",
            const char* outPdf = "myPlots/pdf/scratch/q2_bin_overlays.pdf",
            const char* baseCut = "",
            const char* epCut = "",
            const char* eppCut = "pCM > 0",
            const char* epWeightExpr = "weight_ep",
            const char* eppWeightExpr = "weight_epp") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(inputFile, "READ");
  if (!f || f->IsZombie()) {
    cout << "Could not open input file: " << inputFile << endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(f->Get(treeName));
  if (!tree) {
    cout << "Could not find TTree '" << treeName << "' in " << inputFile << endl;
    f->Close();
    return;
  }

  const string outPath = outPdf;
  const string outDir = gSystem->DirName(outPdf);
  if (!outDir.empty() && outDir != ".") gSystem->mkdir(outDir.c_str(), kTRUE);

  TCanvas* c = new TCanvas("c_q2_bin_overlays", "Q2-binned overlays", 900, 700);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);

  const vector<PlotVar> vars = variablesToPlot();
  const vector<Channel> channels = {
      {"ep", "e'p", epCut, epWeightExpr},
      {"epp", "e'pp", eppCut, eppWeightExpr},
  };
  const int pMissVarIndex = findPlotVarIndex(vars, "pMiss");

  vector<TH1D*> hRatioEpPmiss;
  vector<TH1D*> hRatioEppPmiss;
  if (pMissVarIndex >= 0) {
    const PlotVar& pMissVar = vars[pMissVarIndex];
    hRatioEpPmiss.reserve(bE_Q2.size() - 1);
    hRatioEppPmiss.reserve(bE_Q2.size() - 1);
    for (int q = 0; q < (int)bE_Q2.size() - 1; q++) {
      TH1D* hEp = new TH1D(Form("h_ratio_den_ep_pMiss_Q2_%d", q), "",
                           (int)bE_pmiss.size() - 1, bE_pmiss.data());
      TH1D* hEpp = new TH1D(Form("h_ratio_num_epp_pMiss_Q2_%d", q), "",
                            (int)bE_pmiss.size() - 1, bE_pmiss.data());
      hEp->SetDirectory(nullptr);
      hEpp->SetDirectory(nullptr);
      hEp->Sumw2();
      hEpp->Sumw2();
      hRatioEpPmiss.push_back(hEp);
      hRatioEppPmiss.push_back(hEpp);
    }
  }

  vector<vector<vector<TH1D*>>> hists(
      channels.size(), vector<vector<TH1D*>>(vars.size()));
  for (int cidx = 0; cidx < (int)channels.size(); cidx++) {
    for (int vidx = 0; vidx < (int)vars.size(); vidx++) {
      hists[cidx][vidx].reserve(bE_Q2.size() - 1);
      for (int q = 0; q < (int)bE_Q2.size() - 1; q++) {
        const PlotVar& v = vars[vidx];
        const string hname = Form("h_%s_%s_Q2_%d",
                                  channels[cidx].key.c_str(), v.name.c_str(), q);
        TH1D* h = new TH1D(hname.c_str(), "", binsForChannel(v, channels[cidx]), v.xmin, v.xmax);
        h->Sumw2();
        styleHist(h, q);
        hists[cidx][vidx].push_back(h);
      }
    }
  }

  Formula q2Formula(tree, "q2_formula", "Q2");
  Formula commonCutFormula(tree, "common_cut_formula", baseCut);
  Formula epCutFormula(tree, "ep_cut_formula", epCut);
  Formula eppCutFormula(tree, "epp_cut_formula", eppCut);
  Formula epWeightFormula(tree, "ep_weight_formula", epWeightExpr);
  Formula eppWeightFormula(tree, "epp_weight_formula", eppWeightExpr);

  vector<Formula*> varFormulas;
  varFormulas.reserve(vars.size());
  for (int vidx = 0; vidx < (int)vars.size(); vidx++) {
    varFormulas.push_back(new Formula(tree, Form("var_formula_%d", vidx), vars[vidx].expr));
  }

  bool formulasOK = q2Formula.ok() && commonCutFormula.ok() && epCutFormula.ok() &&
                    eppCutFormula.ok() && epWeightFormula.ok() && eppWeightFormula.ok();
  for (const auto* vf : varFormulas) formulasOK = formulasOK && vf->ok();

  if (!formulasOK) {
    cout << "At least one formula failed to compile. Check variable/cut/weight expressions." << endl;
    for (auto* vf : varFormulas) delete vf;
    f->Close();
    return;
  }

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nEntries; entry++) {
    tree->GetEntry(entry);

    const int q2Bin = findQ2Bin(q2Formula.eval(0.0));
    if (q2Bin < 0) continue;
    if (!passes(commonCutFormula)) continue;

    const bool passEp = passes(epCutFormula);
    const bool passEpp = passes(eppCutFormula);
    const double wEp = passEp ? weight(epWeightFormula) : 0.0;
    const double wEpp = passEpp ? weight(eppWeightFormula) : 0.0;
    if (wEp == 0.0 && wEpp == 0.0) continue;

    if (pMissVarIndex >= 0) {
      const double pMissValue = varFormulas[pMissVarIndex]->eval(0.0);
      if (passEp && wEp != 0.0) hRatioEpPmiss[q2Bin]->Fill(pMissValue, wEp);
      if (passEpp && wEpp != 0.0) hRatioEppPmiss[q2Bin]->Fill(pMissValue, wEpp);
    }

    for (int vidx = 0; vidx < (int)vars.size(); vidx++) {
      const double value = varFormulas[vidx]->eval(0.0);
      if (vars[vidx].drawEp && passEp && wEp != 0.0) hists[0][vidx][q2Bin]->Fill(value, wEp);
      if (vars[vidx].drawEpp && passEpp && wEpp != 0.0) hists[1][vidx][q2Bin]->Fill(value, wEpp);
    }
  }

  c->Print((outPath + "[").c_str(), "pdf");

  for (int vidx = 0; vidx < (int)vars.size(); vidx++) {
    if (vars[vidx].drawEp) drawOverlayPage(c, outPath, vars[vidx], channels[0], hists, 0, vidx, baseCut);
    if (vars[vidx].drawEpp) drawOverlayPage(c, outPath, vars[vidx], channels[1], hists, 1, vidx, baseCut);
  }
  if (pMissVarIndex >= 0) {
    drawRatioOverlayPage(c, outPath, hRatioEppPmiss, hRatioEpPmiss, baseCut, channels[0], channels[1]);
  }

  c->Print((outPath + "]").c_str(), "pdf");
  cout << "Wrote " << outPath << endl;

  for (auto* vf : varFormulas) delete vf;
  f->Close();
}
