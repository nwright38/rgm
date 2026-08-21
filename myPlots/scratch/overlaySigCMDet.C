// Overlay data and simulation histograms from two extrSigCMDet.C outputs.
//
// Usage:
//   root -l -b -q \
//     'overlaySigCMDet.C("extrSigCMDet_data.root",
//                        "extrSigCMDet_sim.root",
//                        "extrSigCMDet_overlays.root",true,
//                        "extrSigCMDet_overlays.pdf",
//                        "lead_cd_rec_fd,lead_cd_rec_cd")'
//
// normalize=true scales each histogram to unit area.  Set it to false to
// compare the weighted yields stored in the two input files.  The PDF contains
// one overlay per page.  skipSimulationFor is a comma- or semicolon-separated
// list of detector combinations for which only data should be plotted.

#include <cmath>
#include <cctype>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <vector>

namespace {

struct ExtractedValues {
  bool found = false;
  double sigmaX = 0.;
  double sigmaXError = 0.;
  double sigmaY = 0.;
  double sigmaYError = 0.;
};

ExtractedValues getExtractedValues(TFile *file,
                                   const char *wantedCombination) {
  ExtractedValues values;
  TTree *tree = dynamic_cast<TTree *>(file->Get("fit_results"));
  if (!tree) return values;

  char combination[64] = "";
  double sigmaX = 0.;
  double sigmaXError = 0.;
  double sigmaY = 0.;
  double sigmaYError = 0.;
  tree->SetBranchAddress("combination", combination);
  tree->SetBranchAddress("sigma_cm_x", &sigmaX);
  tree->SetBranchAddress("sigma_cm_x_error", &sigmaXError);
  tree->SetBranchAddress("sigma_cm_y", &sigmaY);
  tree->SetBranchAddress("sigma_cm_y_error", &sigmaYError);
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    if (TString(combination) != wantedCombination) continue;
    values.found = true;
    values.sigmaX = sigmaX;
    values.sigmaXError = sigmaXError;
    values.sigmaY = sigmaY;
    values.sigmaYError = sigmaYError;
    break;
  }
  tree->ResetBranchAddresses();
  return values;
}

TH1D *getHistogram(TFile *file, const char *combination,
                   const char *variable) {
  return dynamic_cast<TH1D *>(file->Get(
      Form("%s/h_%s_%s", combination, variable, combination)));
}

TF1 *getFit(TFile *file, const char *combination, const char *variable) {
  return dynamic_cast<TF1 *>(file->Get(
      Form("%s/fit_%s_%s", combination, variable, combination)));
}

bool compatibleBinning(const TH1D *a, const TH1D *b) {
  if (!a || !b || a->GetNbinsX() != b->GetNbinsX()) return false;
  for (int bin = 1; bin <= a->GetNbinsX() + 1; ++bin) {
    if (std::abs(a->GetXaxis()->GetBinLowEdge(bin) -
                 b->GetXaxis()->GetBinLowEdge(bin)) > 1.e-12)
      return false;
  }
  return true;
}

void style(TH1D *hist, int color, int marker) {
  hist->SetDirectory(nullptr);
  hist->SetStats(false);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(0.75);
  hist->SetLineWidth(2);
}

void scaleHistogramAndFit(TH1D *hist, TF1 *fit, bool normalize) {
  if (!normalize) return;
  const double integral = hist->Integral();
  if (integral == 0.) return;
  const double factor = 1. / integral;
  hist->Scale(factor);
  // Parameter zero is the Gaussian amplitude in extrSigCMDet.C.
  if (fit) fit->SetParameter(0, factor * fit->GetParameter(0));
}

bool shouldSkipSimulation(const char *combination,
                          const char *skipSimulationFor) {
  if (!skipSimulationFor || skipSimulationFor[0] == '\0') return false;

  std::string token;
  const char *cursor = skipSimulationFor;
  while (*cursor) {
    if (*cursor == ',' || *cursor == ';') {
      if (!token.empty() && token == combination) return true;
      token.clear();
    } else if (!std::isspace(static_cast<unsigned char>(*cursor))) {
      token.push_back(*cursor);
    }
    ++cursor;
  }
  return !token.empty() && token == combination;
}

}  // namespace

void overlaySigCMDet(
    const char *dataFileName = "extrSigCMDet_he4_dat.root",
    const char *simulationFileName = "extrSigCMDet_he4_sim.root",
    const char *outputFileName = "extrSigCMDet_overlays.root",
    bool normalize = true,
  const char *outputPdfName = "extrSigCMDet_overlays.pdf",
  const char *skipSimulationFor = "lead_cd_rec_fd,lead_cd_rec_cd") {

  TFile *dataFile = TFile::Open(dataFileName, "READ");
  TFile *simulationFile = TFile::Open(simulationFileName, "READ");
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "[overlaySigCMDet] Could not open data file "
              << dataFileName << '\n';
    if (simulationFile) simulationFile->Close();
    return;
  }
  if (!simulationFile || simulationFile->IsZombie()) {
    std::cerr << "[overlaySigCMDet] Could not open simulation file "
              << simulationFileName << '\n';
    dataFile->Close();
    return;
  }

  TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "[overlaySigCMDet] Could not create " << outputFileName
              << '\n';
    dataFile->Close();
    simulationFile->Close();
    return;
  }

  const char *combinations[] = {
      "lead_fd_rec_fd", "lead_fd_rec_cd",
      "lead_cd_rec_fd", "lead_cd_rec_cd"};
  const char *variables[] = {"pCMx", "pCMy", "cos2chi"};

  std::vector<TCanvas *> overlayCanvases;
  for (const char *combination : combinations) {
    TDirectory *directory = outputFile->mkdir(combination);
    directory->cd();
    const ExtractedValues dataValues =
        getExtractedValues(dataFile, combination);
    const ExtractedValues simulationValues =
        getExtractedValues(simulationFile, combination);

    for (int variableIndex = 0; variableIndex < 3; ++variableIndex) {
      const char *variable = variables[variableIndex];
      const bool skipSimulation =
          shouldSkipSimulation(combination, skipSimulationFor);
      TH1D *sourceData = getHistogram(dataFile, combination, variable);
      TH1D *sourceSimulation =
          skipSimulation ? nullptr
                         : getHistogram(simulationFile, combination, variable);
      if (!skipSimulation && !sourceSimulation) {
        std::cerr << "[overlaySigCMDet] Missing " << combination << "/h_"
                  << variable << "_" << combination << " in "
                  << simulationFileName
                  << "; continuing with data-only plot" << '\n';
      }
      if (!sourceData) {
        std::cerr << "[overlaySigCMDet] Missing " << combination << "/h_"
                  << variable << "_" << combination << " in "
                  << dataFileName
                  << '\n';
        continue;
      }
      if (sourceSimulation && !compatibleBinning(sourceData, sourceSimulation)) {
        std::cerr << "[overlaySigCMDet] Incompatible binning for "
                  << combination << ", " << variable << '\n';
        continue;
      }

      TH1D *dataHist = dynamic_cast<TH1D *>(sourceData->Clone(
          Form("h_%s_data_%s", variable, combination)));
      TH1D *simulationHist = nullptr;
      if (sourceSimulation) {
        simulationHist = dynamic_cast<TH1D *>(sourceSimulation->Clone(
            Form("h_%s_simulation_%s", variable, combination)));
      }
      style(dataHist, kBlack, 20);
      if (simulationHist) style(simulationHist, kRed + 1, 24);

      TF1 *dataFit = nullptr;
      TF1 *simulationFit = nullptr;
      if (variableIndex < 2) {
        TF1 *sourceDataFit = getFit(dataFile, combination, variable);
        TF1 *sourceSimulationFit =
            sourceSimulation ? getFit(simulationFile, combination, variable)
                             : nullptr;
        if (sourceDataFit) {
          dataFit = dynamic_cast<TF1 *>(sourceDataFit->Clone(
              Form("fit_%s_data_%s", variable, combination)));
          dataFit->SetLineColor(kBlack);
          dataFit->SetLineWidth(2);
        }
        if (sourceSimulationFit) {
          simulationFit = dynamic_cast<TF1 *>(sourceSimulationFit->Clone(
              Form("fit_%s_simulation_%s", variable, combination)));
          simulationFit->SetLineColor(kRed + 1);
          simulationFit->SetLineStyle(2);
          simulationFit->SetLineWidth(2);
        }
      }

      scaleHistogramAndFit(dataHist, dataFit, normalize);
      if (simulationHist)
        scaleHistogramAndFit(simulationHist, simulationFit, normalize);
      dataHist->GetYaxis()->SetTitle(
          normalize ? "Normalized counts" : "Weighted counts");
      double maximum = dataHist->GetMaximum();
      if (simulationHist)
        maximum = TMath::Max(maximum, simulationHist->GetMaximum());
      dataHist->SetMaximum(maximum > 0. ? 1.25 * maximum : 1.);

      TCanvas *canvas = new TCanvas(
          Form("c_overlay_%s_%s", variable, combination),
          Form("Data/simulation overlay: %s, %s", combination, variable),
          900, 700);
      canvas->SetLeftMargin(0.13);
      canvas->SetBottomMargin(0.12);
      canvas->cd();
      dataHist->Draw("E");
      if (simulationHist) simulationHist->Draw("E same");
      if (dataFit) dataFit->Draw("same");
      if (simulationFit) simulationFit->Draw("same");

      TLine *dataMean = nullptr;
      TLine *simulationMean = nullptr;
      if (variableIndex == 2) {
        dataMean = new TLine(dataHist->GetMean(), 0., dataHist->GetMean(),
                             dataHist->GetMaximum());
        if (simulationHist) {
          simulationMean =
              new TLine(simulationHist->GetMean(), 0.,
                        simulationHist->GetMean(),
                        simulationHist->GetMaximum());
        }
        dataMean->SetLineColor(kBlack);
        dataMean->SetLineStyle(2);
        dataMean->SetLineWidth(2);
        if (simulationMean) {
          simulationMean->SetLineColor(kRed + 1);
          simulationMean->SetLineStyle(2);
          simulationMean->SetLineWidth(2);
        }
        dataMean->Draw("same");
        if (simulationMean) simulationMean->Draw("same");
      }

      TLegend *legend = new TLegend(0.57, 0.70, 0.89, 0.89);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->AddEntry(dataHist, "Data", "lep");
      if (simulationHist)
        legend->AddEntry(simulationHist, "Simulation", "lep");
      legend->Draw();

      TLatex valuesText;
      valuesText.SetNDC();
      valuesText.SetTextSize(0.035);
      if (variableIndex == 0) {
        if (dataValues.found) {
          valuesText.SetTextColor(kBlack);
          valuesText.DrawLatex(
              0.15, 0.86,
              Form("Data #sigma_{CM,x} = %.4f #pm %.4f",
                   dataValues.sigmaX, dataValues.sigmaXError));
        }
        if (simulationHist && simulationValues.found) {
          valuesText.SetTextColor(kRed + 1);
          valuesText.DrawLatex(
              0.15, 0.81,
              Form("Simulation #sigma_{CM,x} = %.4f #pm %.4f",
                   simulationValues.sigmaX,
                   simulationValues.sigmaXError));
        }
      } else if (variableIndex == 1) {
        if (dataValues.found) {
          valuesText.SetTextColor(kBlack);
          valuesText.DrawLatex(
              0.15, 0.86,
              Form("Data #sigma_{CM,y} = %.4f #pm %.4f",
                   dataValues.sigmaY, dataValues.sigmaYError));
        }
        if (simulationHist && simulationValues.found) {
          valuesText.SetTextColor(kRed + 1);
          valuesText.DrawLatex(
              0.15, 0.81,
              Form("Simulation #sigma_{CM,y} = %.4f #pm %.4f",
                   simulationValues.sigmaY,
                   simulationValues.sigmaYError));
        }
      } else {
        if (dataValues.found) {
          valuesText.SetTextColor(kBlack);
          valuesText.DrawLatex(
              0.15, 0.86,
              Form("Data #sigma_{CM,x/y} = %.4f / %.4f",
                   dataValues.sigmaX, dataValues.sigmaY));
        }
        if (simulationHist && simulationValues.found) {
          valuesText.SetTextColor(kRed + 1);
          valuesText.DrawLatex(
              0.15, 0.81,
              Form("Simulation #sigma_{CM,x/y} = %.4f / %.4f",
                   simulationValues.sigmaX, simulationValues.sigmaY));
        }
      }

      directory->cd();
      dataHist->Write();
      if (simulationHist) simulationHist->Write();
      if (dataFit) dataFit->Write();
      if (simulationFit) simulationFit->Write();
      if (dataMean)
        dataMean->Write(Form("mean_marker_cos2chi_data_%s", combination));
      if (simulationMean)
        simulationMean->Write(
            Form("mean_marker_cos2chi_simulation_%s", combination));
      canvas->Write();
      overlayCanvases.push_back(canvas);
    }
  }

  if (outputPdfName && outputPdfName[0] != '\0' &&
      !overlayCanvases.empty()) {
    overlayCanvases.front()->Print(Form("%s[", outputPdfName));
    for (TCanvas *canvas : overlayCanvases)
      canvas->Print(outputPdfName, "pdf");
    overlayCanvases.back()->Print(Form("%s]", outputPdfName));
  }

  outputFile->Close();
  dataFile->Close();
  simulationFile->Close();

  std::cout << "[overlaySigCMDet] Wrote " << overlayCanvases.size()
            << " overlay canvases to " << outputFileName;
  if (outputPdfName && outputPdfName[0] != '\0')
    std::cout << " and " << outputPdfName;
  std::cout << '\n';
}
