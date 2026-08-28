// overlayPCMComponentsByDet.C
//
// Build overlay pages for lab-frame (pcmx_lab, pcmy_lab, pcmz_lab) and
// CM-frame (pCMx, pCMy, pCMz) components. Two overlay organizations are
// produced for each frame:
//   1) fixed detector combination, overlay x/y/z
//   2) fixed coordinate (x, y, z), overlay detector combinations
//
// Gaussian fits are performed for all histograms to extract sigma values, but
// the output PDF contains overlay pages only.
//
// Usage example:
//   root -l -b -q 'myPlots/scratch/overlayPCMComponentsByDet.C()'
//   root -l -b -q 'myPlots/scratch/overlayPCMComponentsByDet.C("~/data/RGM_DATA/c12_src_skim.root","srcTree","overlay_pcm_components.pdf",true,"pCM > 0 && pMiss < 1. && recP < 1.","1",-0.15,0.15)'

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

struct DetectorCombination {
  const char *key;
  const char *label;
  const char *cut;
};

struct ComponentSpec {
  const char *key;
  const char *expr;
  const char *axisTitle;
};

struct FrameSpec {
  const char *key;
  const char *label;
  ComponentSpec components[3];
};

struct FitSummary {
  bool fitted = false;
  bool usedMoments = false;
  int fitStatus = -999;
  double nEff = 0.;
  double mean = TMath::QuietNaN();
  double meanErr = TMath::QuietNaN();
  double sigma = TMath::QuietNaN();
  double sigmaErr = TMath::QuietNaN();
};

const int kAxisColors[3] = {kBlack, kRed + 1, kBlue + 1};
const int kComboColors[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};

bool hasBranchOrLeaf(TTree *tree, const char *name) {
  return tree && (tree->GetBranch(name) || tree->GetLeaf(name));
}

void styleAxisOverlay(TH1D *hist, int axisIdx) {
  static const int kMarkers[3] = {20, 21, 22};
  hist->SetStats(false);
  hist->SetLineColor(kAxisColors[axisIdx]);
  hist->SetMarkerColor(kAxisColors[axisIdx]);
  hist->SetMarkerStyle(kMarkers[axisIdx]);
  hist->SetMarkerSize(0.8);
  hist->SetLineWidth(2);
}

void styleComboOverlay(TH1D *hist, int comboIdx) {
  static const int kMarkers[4] = {20, 21, 22, 23};
  hist->SetStats(false);
  hist->SetLineColor(kComboColors[comboIdx]);
  hist->SetMarkerColor(kComboColors[comboIdx]);
  hist->SetMarkerStyle(kMarkers[comboIdx]);
  hist->SetMarkerSize(0.8);
  hist->SetLineWidth(2);
}

FitSummary fitHistogram(TH1D *hist, TF1 *fit, double fitMin, double fitMax) {
  FitSummary out;
  if (!hist) return out;
  out.nEff = hist->GetEffectiveEntries();
  out.mean = hist->GetMean();
  out.meanErr = hist->GetMeanError();

  fit->SetRange(fitMin, fitMax);
  fit->SetParameter(0, hist->GetMaximum());
  fit->SetParameter(1, hist->GetMean());
  fit->SetParameter(2, TMath::Max(0.05, 0.5 * hist->GetRMS()));

  if (out.nEff >= 10.) {
    TFitResultPtr fr = hist->Fit(fit, "RQS0");
    out.fitStatus = static_cast<int>(fr);
    if (fr.Get()) {
      out.fitted = true;
      out.sigma = std::abs(fr->Parameter(2));
      out.sigmaErr = fr->ParError(2);
      out.mean = fr->Parameter(1);
      out.meanErr = fr->ParError(1);
      return out;
    }
  }

  out.usedMoments = true;
  out.fitStatus = -1;
  out.sigma = hist->GetRMS();
  out.sigmaErr = hist->GetRMSError();
  return out;
}

TH1D *makeOverlayClone(const TH1D *source, const char *newName,
                       bool normalizeToUnity) {
  if (!source) return nullptr;
  TH1D *h = dynamic_cast<TH1D *>(source->Clone(newName));
  h->SetDirectory(nullptr);
  if (normalizeToUnity) {
    const double integral = h->Integral();
    if (integral > 0.) h->Scale(1. / integral);
  }
  return h;
}

void drawSigmaLine(double x, double y, int color, const char *label,
                   const FitSummary &s) {
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.036);
  t.SetTextColor(color);
  if (std::isfinite(s.sigma) && std::isfinite(s.sigmaErr)) {
    t.DrawLatex(x, y, Form("%s: #sigma_{CM}=%.3f #pm %.3f", label, s.sigma,
                           s.sigmaErr));
  } else {
    t.DrawLatex(x, y, Form("%s: #sigma_{CM}=nan", label));
  }
}

}  // namespace

void overlayPCMComponentsByDet(
    const char *inputFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *treeName = "srcTree",
    const char *outputPdfName = "overlay_pcm_components.pdf",
    bool normalizeOverlays = true,
    const char *baseCut = "weight_epp < 200.",
    const char *weightExpression = "weight_epp",
    double fitMin = -0.2,
    double fitMax = 0.2,
    int bins = 50,
    double xmin = -0.8,
    double xmax = 0.8) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "[overlayPCMComponentsByDet] Could not open input file "
              << inputFileName << '\n';
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(inputFile->Get(treeName));
  if (!tree) {
    std::cerr << "[overlayPCMComponentsByDet] Could not find tree \""
              << treeName << "\" in " << inputFileName << '\n';
    inputFile->Close();
    return;
  }

  const char *requiredBranches[] = {
      "leadTheta", "recTheta", "pcmx_lab", "pcmy_lab", "pcmz_lab",
      "pCMx", "pCMy", "pCMz"};
  for (const char *branch : requiredBranches) {
    if (!hasBranchOrLeaf(tree, branch)) {
      std::cerr << "[overlayPCMComponentsByDet] Missing required branch \""
                << branch << "\" in " << treeName << '\n';
      inputFile->Close();
      return;
    }
  }

  const DetectorCombination combos[4] = {
      {"lead_fd_rec_fd", "Lead FD Rec FD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_fd_rec_cd", "Lead FD Rec CD",
       "(leadTheta*180./TMath::Pi()<37.) && (recTheta*180./TMath::Pi()>45.)"},
      {"lead_cd_rec_fd", "Lead CD Rec FD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()<37.)"},
      {"lead_cd_rec_cd", "Lead CD Rec CD",
       "(leadTheta*180./TMath::Pi()>45.) && (recTheta*180./TMath::Pi()>45.)"}};

  const FrameSpec frames[2] = {
      {"lab", "Lab frame",
       {{"x", "pcmx_lab", "p_{CM,x}^{lab} [GeV/c]"},
        {"y", "pcmy_lab", "p_{CM,y}^{lab} [GeV/c]"},
        {"z", "pcmz_lab", "p_{CM,z}^{lab} [GeV/c]"}}},
      {"cm", "CM frame",
       {{"x", "pCMx", "p_{CM,x} [GeV/c]"},
        {"y", "pCMy", "p_{CM,y} [GeV/c]"},
        {"z", "pCMz", "p_{CM,z} [GeV/c]"}}}};

  TH1D *hRaw[2][4][3] = {{{nullptr}}};
  TF1 *fGaus[2][4][3] = {{{nullptr}}};
  FitSummary fitSummary[2][4][3];
  const int kFirstIncludedCombo = 1;  // Omit FD+FD everywhere.

  for (int iFrame = 0; iFrame < 2; ++iFrame) {
    for (int iCombo = kFirstIncludedCombo; iCombo < 4; ++iCombo) {
      const TString comboCut = combos[iCombo].cut;
      const TString fullWeight =
          Form("(%s)*((%s)&&(%s))", weightExpression, baseCut,
               comboCut.Data());
      for (int iAxis = 0; iAxis < 3; ++iAxis) {
        const ComponentSpec &comp = frames[iFrame].components[iAxis];
        const TString histName =
            Form("h_%s_%s_%s", frames[iFrame].key, comp.key, combos[iCombo].key);

        hRaw[iFrame][iCombo][iAxis] =
            new TH1D(histName, Form("%s;%s;Weighted counts", combos[iCombo].label,
                                   comp.axisTitle),
                     bins, xmin, xmax);
        const Long64_t nProjected = tree->Project(histName, comp.expr, fullWeight);
        hRaw[iFrame][iCombo][iAxis]->SetDirectory(nullptr);
        if (nProjected <= 0) {
          std::cerr << "[overlayPCMComponentsByDet] Warning: zero projected entries for "
                    << histName << " with cut " << fullWeight << '\n';
        }

        fGaus[iFrame][iCombo][iAxis] =
            new TF1(Form("fit_%s", histName.Data()), "gaus", fitMin, fitMax);
        fGaus[iFrame][iCombo][iAxis]->SetLineColor(kMagenta + 2);
        fGaus[iFrame][iCombo][iAxis]->SetLineWidth(2);

        fitSummary[iFrame][iCombo][iAxis] =
            fitHistogram(hRaw[iFrame][iCombo][iAxis],
                         fGaus[iFrame][iCombo][iAxis], fitMin, fitMax);
      }
    }
  }

  TCanvas *canvas =
      new TCanvas("c_overlay_pcm_components", "PCM component overlays", 1400, 900);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TString pdfName(outputPdfName);
  canvas->Print(pdfName + "[", "pdf");

  // 1) Same detector config: overlay x/y/z (no fit drawn).
  for (int iFrame = 0; iFrame < 2; ++iFrame) {
    for (int iCombo = kFirstIncludedCombo; iCombo < 4; ++iCombo) {
      canvas->Clear();
      canvas->cd(1);
      gPad->SetLeftMargin(0.12);
      gPad->SetBottomMargin(0.12);

      TH1D *hOverlay[3] = {nullptr, nullptr, nullptr};
      double yMax = 0.;
      for (int iAxis = 0; iAxis < 3; ++iAxis) {
        hOverlay[iAxis] = makeOverlayClone(
            hRaw[iFrame][iCombo][iAxis],
            Form("h_overlay_combo_%s_%s_%d", frames[iFrame].key,
                 combos[iCombo].key, iAxis),
            normalizeOverlays);
        styleAxisOverlay(hOverlay[iAxis], iAxis);
        yMax = std::max(yMax, hOverlay[iAxis]->GetMaximum());
      }

      if (yMax <= 0.) yMax = 1.;
      hOverlay[0]->SetMaximum(1.25 * yMax);
      hOverlay[0]->SetTitle(Form("%s | %s", frames[iFrame].label, combos[iCombo].label));
      hOverlay[0]->GetYaxis()->SetTitle(normalizeOverlays ? "Normalized counts"
                                                          : "Weighted counts");
      hOverlay[0]->Draw("E");
      hOverlay[1]->Draw("E same");
      hOverlay[2]->Draw("E same");

      // Print fitted sigma values without drawing fit curves on overlays.
      const double yStart = 0.86;
      const double yStep = 0.06;
      for (int iAxis = 0; iAxis < 3; ++iAxis) {
        drawSigmaLine(0.56, yStart - iAxis * yStep, kAxisColors[iAxis],
                      frames[iFrame].components[iAxis].axisTitle,
                      fitSummary[iFrame][iCombo][iAxis]);
      }

      canvas->Print(pdfName, "pdf");

      for (int iAxis = 0; iAxis < 3; ++iAxis) delete hOverlay[iAxis];

    }
  }

  // 2) Same coordinate: overlay detector combinations (no fit drawn).
  for (int iFrame = 0; iFrame < 2; ++iFrame) {
    for (int iAxis = 0; iAxis < 3; ++iAxis) {
      canvas->Clear();
      canvas->cd(1);
      gPad->SetLeftMargin(0.12);
      gPad->SetBottomMargin(0.12);

      TH1D *hOverlay[4] = {nullptr, nullptr, nullptr, nullptr};
      double yMax = 0.;
      for (int iCombo = kFirstIncludedCombo; iCombo < 4; ++iCombo) {
        hOverlay[iCombo] = makeOverlayClone(
            hRaw[iFrame][iCombo][iAxis],
            Form("h_overlay_axis_%s_%s_%d", frames[iFrame].key,
                 frames[iFrame].components[iAxis].key, iCombo),
            normalizeOverlays);
        styleComboOverlay(hOverlay[iCombo], iCombo);
        yMax = std::max(yMax, hOverlay[iCombo]->GetMaximum());
      }

      if (yMax <= 0.) yMax = 1.;
      hOverlay[kFirstIncludedCombo]->SetMaximum(1.25 * yMax);
      hOverlay[kFirstIncludedCombo]->SetTitle(Form("%s | %s", frames[iFrame].label,
                                                   frames[iFrame].components[iAxis].axisTitle));
      hOverlay[kFirstIncludedCombo]->GetYaxis()->SetTitle(
          normalizeOverlays ? "Normalized counts" : "Weighted counts");
      hOverlay[kFirstIncludedCombo]->Draw("E");
      for (int iCombo = kFirstIncludedCombo + 1; iCombo < 4; ++iCombo) {
        hOverlay[iCombo]->Draw("E same");
      }

      // Print fitted sigma values without drawing fit curves on overlays.
      const double yStart = 0.86;
      const double yStep = 0.06;
      for (int iCombo = kFirstIncludedCombo; iCombo < 4; ++iCombo) {
        drawSigmaLine(0.54, yStart - (iCombo - kFirstIncludedCombo) * yStep,
                      kComboColors[iCombo], combos[iCombo].label,
                      fitSummary[iFrame][iCombo][iAxis]);
      }

      canvas->Print(pdfName, "pdf");

      for (int iCombo = kFirstIncludedCombo; iCombo < 4; ++iCombo) {
        delete hOverlay[iCombo];
      }

    }
  }

  canvas->Print(pdfName + "]", "pdf");

  inputFile->Close();
}