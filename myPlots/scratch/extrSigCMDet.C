// Extract sigma_CM and <cos(2 chi)> for all lead/recoil FD/CD combinations.
//
// The four detector combinations are
//   leadTheta < 37 deg, recTheta < 37 deg
//   leadTheta < 37 deg, recTheta > 45 deg
//   leadTheta > 45 deg, recTheta < 37 deg
//   leadTheta > 45 deg, recTheta > 45 deg
//
// Usage:
//   root -l -b -q 'extrSigCMDet.C()'
//   root -l -b -q 'extrSigCMDet.C("input.root","output.root")'
//
// The output contains one directory per detector combination.  Each directory
// has the three histograms, their fit/mean functions, and a canvas showing the
// overlays.  The "fit_results" TTree provides the extracted numbers in a form
// that is easy to inspect or use in another macro.  Because the FD/FD sample
// is small and non-Gaussian, its sigma_cm_x/y values are histogram RMS values;
// its histogram means are saved explicitly as mean_cm_x/y.

#include <cmath>
#include <cstring>
#include <iostream>

#include <TCanvas.h>
#include <TCut.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

namespace {

struct DetectorRegion {
  const char *key;
  const char *label;
  const char *cut;
};

struct Combination {
  const DetectorRegion *lead;
  const DetectorRegion *recoil;
};

bool validInput(TFile *file, TTree *tree, const char *fileName,
                const char *treeName) {
  if (!file || file->IsZombie()) {
    std::cerr << "[extrSigCMDet] Could not open input file " << fileName << '\n';
    return false;
  }
  if (!tree) {
    std::cerr << "[extrSigCMDet] Could not find TTree \"" << treeName
              << "\" in " << fileName << '\n';
    return false;
  }

  const char *requiredBranches[] = {
      "leadTheta", "recTheta", "pCMx", "pCMy", "chi", "weight_epp"};
  for (const char *branch : requiredBranches) {
    if (!tree->GetBranch(branch) && !tree->GetLeaf(branch)) {
      std::cerr << "[extrSigCMDet] Required branch \"" << branch
                << "\" is missing from " << treeName << '\n';
      return false;
    }
  }
  return true;
}

void styleHistogram(TH1D *hist, int color) {
  hist->Sumw2();
  hist->SetStats(false);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.7);
  hist->GetYaxis()->SetTitleOffset(1.35);
  (void)color;
}

}  // namespace

void extrSigCMDet(
    const char *dataFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *outputFileName = "extrSigCMDet.root",
    const char *treeName = "srcTree",
    const char *extraCut = "pCM > 0 && (weight_epp < 200)",
    const char *weightExpression = "(weight_epp)") {

  TFile *dataFile = TFile::Open(dataFileName, "READ");
  TTree *dataTree =
      dataFile && !dataFile->IsZombie()
          ? dynamic_cast<TTree *>(dataFile->Get(treeName))
          : nullptr;
  if (!validInput(dataFile, dataTree, dataFileName, treeName)) {
    if (dataFile) dataFile->Close();
    return;
  }

  TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "[extrSigCMDet] Could not create output file "
              << outputFileName << '\n';
    dataFile->Close();
    return;
  }

  // The gap from 37 to 45 degrees is deliberately excluded.
  const DetectorRegion fd = {
      "fd", "FD (#theta < 37#circ)",
      "%sTheta*180./TMath::Pi() < 37."};
  const DetectorRegion cd = {
      "cd", "CD (#theta > 45#circ)",
      "%sTheta*180./TMath::Pi() > 45."};
  const Combination combinations[] = {
      {&fd, &fd}, {&fd, &cd}, {&cd, &fd}, {&cd, &cd}};

  char combinationName[64] = "";
  char extractionMethod[32] = "";
  char selection[512] = "";
  double sigmaCMx = 0.;
  double sigmaCMxError = 0.;
  double sigmaCMy = 0.;
  double sigmaCMyError = 0.;
  double meanCMx = 0.;
  double meanCMxError = 0.;
  double meanCMy = 0.;
  double meanCMyError = 0.;
  double meanCos2Chi = 0.;
  double meanCos2ChiError = 0.;
  double sumWeights = 0.;
  Long64_t entries = 0;
  int fitStatusX = -1;
  int fitStatusY = -1;

  outputFile->cd();
  TTree *results = new TTree(
      "fit_results",
      "sigma_CM Gaussian fits and histogram mean of cos(2 chi)");
  results->Branch("combination", combinationName, "combination/C");
  results->Branch("extraction_method", extractionMethod,
                  "extraction_method/C");
  results->Branch("selection", selection, "selection/C");
  results->Branch("sigma_cm_x", &sigmaCMx, "sigma_cm_x/D");
  results->Branch("sigma_cm_x_error", &sigmaCMxError,
                  "sigma_cm_x_error/D");
  results->Branch("sigma_cm_y", &sigmaCMy, "sigma_cm_y/D");
  results->Branch("sigma_cm_y_error", &sigmaCMyError,
                  "sigma_cm_y_error/D");
  results->Branch("mean_cm_x", &meanCMx, "mean_cm_x/D");
  results->Branch("mean_cm_x_error", &meanCMxError, "mean_cm_x_error/D");
  results->Branch("mean_cm_y", &meanCMy, "mean_cm_y/D");
  results->Branch("mean_cm_y_error", &meanCMyError, "mean_cm_y_error/D");
  results->Branch("mean_cos2chi", &meanCos2Chi, "mean_cos2chi/D");
  results->Branch("mean_cos2chi_error", &meanCos2ChiError,
                  "mean_cos2chi_error/D");
  results->Branch("entries", &entries, "entries/L");
  results->Branch("sum_weights", &sumWeights, "sum_weights/D");
  results->Branch("fit_status_x", &fitStatusX, "fit_status_x/I");
  results->Branch("fit_status_y", &fitStatusY, "fit_status_y/I");

  const double componentMin = -0.8;
  const double componentMax = 0.8;
  const double fitMin = -0.2;
  const double fitMax = 0.2;

  for (const Combination &combination : combinations) {
    snprintf(combinationName, sizeof(combinationName), "lead_%s_rec_%s",
             combination.lead->key, combination.recoil->key);
    const bool useHistogramMoments =
        !strcmp(combination.lead->key, "fd") &&
        !strcmp(combination.recoil->key, "fd");
    snprintf(extractionMethod, sizeof(extractionMethod), "%s",
             useHistogramMoments ? "histogram_rms" : "gaussian_fit");

    const TString leadCut =
        Form(combination.lead->cut, "lead");
    const TString recoilCut =
        Form(combination.recoil->cut, "rec");
    const TString booleanCut =
        Form("(%s) && (%s) && (%s)", extraCut, leadCut.Data(),
             recoilCut.Data());
    const TString weightedCut =
        Form("(%s)*(%s)", weightExpression, booleanCut.Data());
    snprintf(selection, sizeof(selection), "%s", booleanCut.Data());

    TDirectory *directory = outputFile->mkdir(combinationName);
    directory->cd();

    TH1D *hX = new TH1D(
        Form("h_pCMx_%s", combinationName),
        Form("%s, %s;p_{CM,x} [GeV/c];Weighted counts",
             combination.lead->label, combination.recoil->label),
        50, componentMin, componentMax);
    TH1D *hY = new TH1D(
        Form("h_pCMy_%s", combinationName),
        Form("%s, %s;p_{CM,y} [GeV/c];Weighted counts",
             combination.lead->label, combination.recoil->label),
        50, componentMin, componentMax);
    TH1D *hCos = new TH1D(
        Form("h_cos2chi_%s", combinationName),
        Form("%s, %s;cos(2#chi);Weighted counts",
             combination.lead->label, combination.recoil->label),
        40, -1., 1.);
    styleHistogram(hX, kBlue + 1);
    styleHistogram(hY, kRed + 1);
    styleHistogram(hCos, kGreen + 2);

    dataTree->Project(hX->GetName(), "pCMx", weightedCut);
    dataTree->Project(hY->GetName(), "pCMy", weightedCut);
    dataTree->Project(hCos->GetName(), "TMath::Cos(2*chi)", weightedCut);

    TF1 *fitX = new TF1(Form("fit_pCMx_%s", combinationName), "gaus",
                        fitMin, fitMax);
    TF1 *fitY = new TF1(Form("fit_pCMy_%s", combinationName), "gaus",
                        fitMin, fitMax);
    fitX->SetLineColor(kBlue + 1);
    fitY->SetLineColor(kRed + 1);
    fitX->SetLineWidth(2);
    fitY->SetLineWidth(2);

    fitStatusX = -1;
    fitStatusY = -1;
    meanCMx = hX->GetMean();
    meanCMxError = hX->GetMeanError();
    meanCMy = hY->GetMean();
    meanCMyError = hY->GetMeanError();
    sigmaCMx = sigmaCMxError = TMath::QuietNaN();
    sigmaCMy = sigmaCMyError = TMath::QuietNaN();
    if (useHistogramMoments) {
      sigmaCMx = hX->GetRMS();
      sigmaCMxError = hX->GetRMSError();
      sigmaCMy = hY->GetRMS();
      sigmaCMyError = hY->GetRMSError();
    } else if (hX->GetEffectiveEntries() >= 10.) {
      TFitResultPtr fitResult = hX->Fit(fitX, "RQS");
      fitStatusX = static_cast<int>(fitResult);
      if (fitResult.Get()) {
        sigmaCMx = std::abs(fitResult->Parameter(2));
        sigmaCMxError = fitResult->ParError(2);
      }
    }
    if (!useHistogramMoments && hY->GetEffectiveEntries() >= 10.) {
      TFitResultPtr fitResult = hY->Fit(fitY, "RQS");
      fitStatusY = static_cast<int>(fitResult);
      if (fitResult.Get()) {
        sigmaCMy = std::abs(fitResult->Parameter(2));
        sigmaCMyError = fitResult->ParError(2);
      }
    }

    meanCos2Chi = hCos->GetMean();
    meanCos2ChiError = hCos->GetMeanError();
    entries = static_cast<Long64_t>(hCos->GetEntries());
    sumWeights = hCos->Integral(0, hCos->GetNbinsX() + 1);

    // A constant function at the histogram mean makes the extracted average
    // visible on the saved cos(2 chi) plot.  It is not a fit to bin contents.
    TF1 *meanLine = new TF1(Form("mean_cos2chi_%s", combinationName),
                            "[0]", -1., 1.);
    meanLine->SetParameter(0, meanCos2Chi);
    meanLine->FixParameter(0, meanCos2Chi);
    meanLine->SetLineColor(kGreen + 2);
    meanLine->SetLineStyle(2);
    meanLine->SetLineWidth(2);

    TCanvas *canvas = new TCanvas(
        Form("c_%s", combinationName), combinationName, 1500, 500);
    canvas->Divide(3, 1);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.045);

    canvas->cd(1);
    hX->Draw("E");
    if (!useHistogramMoments && fitStatusX >= 0) fitX->Draw("same");
    text.DrawLatex(
        0.16, 0.84,
        Form("%s_{CM,x} = %.4f #pm %.4f",
             useHistogramMoments ? "RMS" : "#sigma", sigmaCMx,
             sigmaCMxError));
    text.DrawLatex(0.16, 0.76,
                   Form("Mean_{CM,x} = %.4f #pm %.4f", meanCMx,
                        meanCMxError));

    canvas->cd(2);
    hY->Draw("E");
    if (!useHistogramMoments && fitStatusY >= 0) fitY->Draw("same");
    text.DrawLatex(
        0.16, 0.84,
        Form("%s_{CM,y} = %.4f #pm %.4f",
             useHistogramMoments ? "RMS" : "#sigma", sigmaCMy,
             sigmaCMyError));
    text.DrawLatex(0.16, 0.76,
                   Form("Mean_{CM,y} = %.4f #pm %.4f", meanCMy,
                        meanCMyError));

    canvas->cd(3);
    hCos->Draw("E");
    const double ymax = hCos->GetMaximum() > 0. ? 1.05 * hCos->GetMaximum() : 1.;
    TLine *verticalMean =
        new TLine(meanCos2Chi, 0., meanCos2Chi, ymax);
    verticalMean->SetLineColor(kGreen + 2);
    verticalMean->SetLineStyle(2);
    verticalMean->SetLineWidth(2);
    verticalMean->Draw("same");
    text.DrawLatex(0.16, 0.84,
                   Form("#LTcos(2#chi)#GT = %.4f #pm %.4f",
                        meanCos2Chi, meanCos2ChiError));

    hX->Write();
    hY->Write();
    hCos->Write();
    if (!useHistogramMoments) {
      fitX->Write();
      fitY->Write();
    }
    meanLine->Write();
    verticalMean->Write(Form("mean_marker_cos2chi_%s", combinationName));
    canvas->Write();
    results->Fill();

    std::cout << combinationName << ": sigma_cm_x = " << sigmaCMx
              << " +/- " << sigmaCMxError
              << ", sigma_cm_y = " << sigmaCMy
              << " +/- " << sigmaCMyError
              << ", mean_cm_x = " << meanCMx
              << " +/- " << meanCMxError
              << ", mean_cm_y = " << meanCMy
              << " +/- " << meanCMyError
              << ", <cos(2 chi)> = " << meanCos2Chi
              << " +/- " << meanCos2ChiError << '\n';
  }

  outputFile->cd();
  results->Write();
  outputFile->Write();
  outputFile->Close();
  dataFile->Close();

  std::cout << "[extrSigCMDet] Wrote histograms, overlays, and fit_results to "
            << outputFileName << '\n';
}
