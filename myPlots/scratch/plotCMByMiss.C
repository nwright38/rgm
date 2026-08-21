// STEP 2 OF 2: read the ROOT histograms made by fillCMByMiss.C, fit them,
// and create a multipage PDF.
//
// Default workflow (after running fillCMByMiss.C):
//   root -l -b -q 'myPlots/scratch/plotCMByMiss.C()'
//     Input:  cm_by_miss.root
//     Output: cm_by_miss.pdf
//
// Custom ROOT input and PDF output:
//   root -l -b -q 'myPlots/scratch/plotCMByMiss.C("my_hists.root","my_plots.pdf")'
//
// Several files (histograms are added before fitting):
//   root -l -b -q 'myPlots/scratch/plotCMByMiss.C("run1.root,run2.root","cm.pdf")'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

struct Bin {
  std::string variable;
  int index;
  double low;
  double high;
};

struct FitPoint {
  double center;
  double halfWidth;
  double sigma;
  double sigmaError;
  bool valid;
};

std::vector<std::string> splitFiles(const char *text) {
  std::vector<std::string> files;
  std::stringstream stream(text ? text : "");
  std::string item;
  while (std::getline(stream, item, ','))
    if (!item.empty()) files.push_back(item);
  return files;
}

std::vector<Bin> readBins(TFile *file) {
  std::vector<Bin> bins;
  TTree *tree = dynamic_cast<TTree *>(file->Get("bin_info"));
  if (!tree) return bins;
  char variable[16] = "";
  int bin = 0;
  double low = 0., high = 0.;
  tree->SetBranchAddress("variable", variable);
  tree->SetBranchAddress("bin", &bin);
  tree->SetBranchAddress("low", &low);
  tree->SetBranchAddress("high", &high);
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    bins.push_back({variable, bin, low, high});
  }
  return bins;
}

TH1D *combinedHistogram(const std::vector<TFile *> &files,
                        const std::string &name) {
  TH1D *sum = nullptr;
  for (TFile *file : files) {
    TH1D *source = dynamic_cast<TH1D *>(file->Get(name.c_str()));
    if (!source) continue;
    if (!sum) {
      sum = dynamic_cast<TH1D *>(source->Clone((name + "_sum").c_str()));
      sum->SetDirectory(nullptr);
    } else {
      sum->Add(source);
    }
  }
  return sum;
}

FitPoint drawAndFit(TH1D *hist, const Bin &bin, const char component) {
  FitPoint point = {(bin.low + bin.high) / 2., (bin.high - bin.low) / 2.,
                    0., 0., false};
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.050);
  if (!hist || hist->GetEffectiveEntries() < 5.) {
    text.DrawLatex(0.18, 0.80, "Insufficient entries");
    return point;
  }
  hist->SetStats(false);
  hist->SetLineColor(component == 'x' ? kBlue + 1 : kRed + 1);
  hist->SetMarkerColor(hist->GetLineColor());
  hist->SetMarkerStyle(20);
  hist->Draw("E");
  TF1 *fit = new TF1(Form("fit_%s_%c_%02d", bin.variable.c_str(), component,
                          bin.index),
                     "gaus", -0.15, 0.15);
  fit->SetLineColor(kBlack);
  TFitResultPtr result = hist->Fit(fit, "RQSN");
  if (result.Get() && static_cast<int>(result) == 0 &&
      std::isfinite(result->Parameter(2))) {
    point.sigma = std::abs(result->Parameter(2));
    point.sigmaError = result->ParError(2);
    point.valid = true;
    fit->Draw("SAME");
    text.DrawLatex(0.17, 0.82,
                   Form("#sigma_{CM,%c} = %.4f #pm %.4f GeV/c", component,
                        point.sigma, point.sigmaError));
  } else {
    text.DrawLatex(0.18, 0.82, "Gaussian fit failed");
  }
  return point;
}

TGraphErrors *makeGraph(const std::vector<FitPoint> &points, int color,
                        int marker, const char *name) {
  TGraphErrors *graph = new TGraphErrors();
  graph->SetName(name);
  for (const FitPoint &point : points) {
    if (!point.valid) continue;
    const int n = graph->GetN();
    graph->SetPoint(n, point.center, point.sigma);
    graph->SetPointError(n, point.halfWidth, point.sigmaError);
  }
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  graph->SetLineWidth(2);

  return graph;
}

}  // namespace

void plotCMByMiss(const char *inputFileNames = "cm_by_miss_he.root",
                  const char *outputPdf = "cm_by_miss_he.pdf") {
  gStyle->SetOptStat(0);
  std::vector<std::unique_ptr<TFile>> ownedFiles;
  std::vector<TFile *> files;
  for (const std::string &name : splitFiles(inputFileNames)) {
    std::unique_ptr<TFile> file(TFile::Open(name.c_str(), "READ"));
    if (!file || file->IsZombie()) {
      std::cerr << "[plotCMByMiss] Cannot open " << name << '\n';
      continue;
    }
    files.push_back(file.get());
    ownedFiles.push_back(std::move(file));
  }
  if (files.empty()) return;
  std::vector<Bin> bins = readBins(files.front());
  if (bins.empty()) {
    std::cerr << "[plotCMByMiss] bin_info is missing.\n";
    return;
  }

  TCanvas canvas("cm_canvas", "CM widths", 1200, 850);
  canvas.Print(Form("%s[", outputPdf));
  std::vector<FitPoint> pX, pY, kX, kY;
  for (const std::string variable : {"pmiss", "kmiss"}) {
    std::vector<Bin> selected;
    for (const Bin &bin : bins)
      if (bin.variable == variable) selected.push_back(bin);
    std::sort(selected.begin(), selected.end(),
              [](const Bin &a, const Bin &b) { return a.index < b.index; });

    for (int first = 0; first < static_cast<int>(selected.size()); first += 4) {
      canvas.Clear();
      canvas.Divide(2, 2);
      std::vector<TH1D *> pageHistograms;
      for (int offset = 0; offset < 4 && first + offset < static_cast<int>(selected.size());
           ++offset) {
        const Bin &bin = selected[first + offset];
        canvas.cd(offset + 1);
        TH1D *hx = combinedHistogram(
            files, Form("h_pCMx_%s_bin%02d", variable.c_str(), bin.index));
        pageHistograms.push_back(hx);
        FitPoint px = drawAndFit(hx, bin, 'x');
        (variable == "pmiss" ? pX : kX).push_back(px);
      }
      canvas.Print(outputPdf);
      for (TH1D *hist : pageHistograms) delete hist;
    }
    for (int first = 0; first < static_cast<int>(selected.size()); first += 4) {
      canvas.Clear();
      canvas.Divide(2, 2);
      std::vector<TH1D *> pageHistograms;
      for (int offset = 0; offset < 4 && first + offset < static_cast<int>(selected.size());
           ++offset) {
        const Bin &bin = selected[first + offset];
        canvas.cd(offset + 1);
        TH1D *hy = combinedHistogram(
            files, Form("h_pCMy_%s_bin%02d", variable.c_str(), bin.index));
        pageHistograms.push_back(hy);
        FitPoint py = drawAndFit(hy, bin, 'y');
        (variable == "pmiss" ? pY : kY).push_back(py);
      }
      canvas.Print(outputPdf);
      for (TH1D *hist : pageHistograms) delete hist;
    }
  }

  const std::vector<std::string> variables = {"pmiss", "kmiss"};
  const std::vector<std::vector<FitPoint> *> xs = {&pX, &kX};
  const std::vector<std::vector<FitPoint> *> ys = {&pY, &kY};
  for (int i = 0; i < 2; ++i) {
    canvas.Clear();
    canvas.cd();
    TGraphErrors *gx = makeGraph(*xs[i], kBlue + 1, 20,
                                 Form("g_sigma_x_%s", variables[i].c_str()));
    TGraphErrors *gy = makeGraph(*ys[i], kRed + 1, 21,
                                 Form("g_sigma_y_%s", variables[i].c_str()));
    TMultiGraph *multi = new TMultiGraph(
        Form("mg_%s", variables[i].c_str()),
        Form("#sigma_{CM} vs %s;%s [GeV/c];#sigma_{CM} [GeV/c]",
             variables[i].c_str(),
             variables[i] == "pmiss" ? "p_{miss}" : "k_{miss}"));
    multi->Add(gx, "P");
    multi->Add(gy, "P");
    multi->Draw("A");
    multi->GetYaxis()->SetRangeUser(.1,.3);

    TLegend *legend = new TLegend(0.62, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(gx, "p_{CM,x}", "lp");
    legend->AddEntry(gy, "p_{CM,y}", "lp");
    legend->Draw();
    canvas.Print(outputPdf);
  }
  canvas.Print(Form("%s]", outputPdf));
  std::cout << "[plotCMByMiss] Wrote " << outputPdf << '\n';
}
