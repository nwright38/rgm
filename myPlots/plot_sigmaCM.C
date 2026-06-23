// plot_sigmaCM.C
//
// Reads the output ROOT file from the sigma_CM fitting code and produces
// a multi-page PDF with one page per sigma_CM component.
//
// Usage (from ROOT prompt):
//   root -l -q 'plot_sigmaCM.C("fitOutput.root","sigmaCM_plots.pdf")'
//
// Or compiled:
//   root -l -q 'plot_sigmaCM.C+("fitOutput.root","sigmaCM_plots.pdf")'
//
// Contents expected in the input file:
//   sigmacmx_int  -- TGraphAsymmErrors, integrated sigma_CMx (single point)
//   sigmacmy_int  -- TGraphAsymmErrors, integrated sigma_CMy (single point)
//   sigmacmz_int  -- TGraphAsymmErrors, integrated sigma_CMz (single point)
//   sigmacmT_int  -- TGraphAsymmErrors, integrated sigma_CMT (single point)
//   sigmacmx_Q2   -- TGraphAsymmErrors, sigma_CMx vs Q2 bins
//   sigmacmy_Q2   -- TGraphAsymmErrors, sigma_CMy vs Q2 bins
//   sigmacmz_Q2   -- TGraphAsymmErrors, sigma_CMz vs Q2 bins
//   sigmacmT_Q2   -- TGraphAsymmErrors, sigma_CMT vs Q2 bins
//
// Note: the _int graphs store absolute sigma positions as the Y errors
//       (not deltas); the _Q2 graphs store proper deltas. Both are handled.

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <cstdio>

// Q2 bin edges (must match what was used to produce the file)
static const std::vector<double> bE_Q2 = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};
static const int bQ2 = (int)bE_Q2.size() - 1;

static const double min_sigma = 0.050;
static const double max_sigma = 0.250;

// -----------------------------------------------------------------------
// Draw one page: sigma_CM component vs Q2
// -----------------------------------------------------------------------
void drawPage(TCanvas *c,
              TGraphAsymmErrors *gQ2,
              TGraphAsymmErrors *gInt,
              const char *label,   // e.g. "CMx", "CMT"
              const char *pdfFile,
              bool isLast = false)
{
    c->Clear();
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13);

    // --- Q2-differential graph (main panel) ---
    gQ2->SetMarkerStyle(20);
    gQ2->SetMarkerSize(1.2);
    gQ2->SetMarkerColor(kBlue+1);
    gQ2->SetLineColor(kBlue+1);
    gQ2->SetLineWidth(2);
    gQ2->SetFillColorAlpha(kBlue+1, 0.20);

    // Axis labels
    gQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2}/c^{2})");
    gQ2->GetYaxis()->SetTitle(Form("#sigma_{%s} (GeV/c)", label));
    gQ2->GetXaxis()->SetTitleSize(0.05);
    gQ2->GetYaxis()->SetTitleSize(0.05);
    gQ2->GetXaxis()->SetLabelSize(0.04);
    gQ2->GetYaxis()->SetLabelSize(0.04);
    gQ2->GetYaxis()->SetRangeUser(0.0, max_sigma * 1.2);
    gQ2->GetXaxis()->SetRangeUser(bE_Q2.front() * 0.9, bE_Q2.back() * 1.05);

    gQ2->Draw("A P Z");   // Z = draw error bar end-caps

    // --- Integrated value as a horizontal band ---
    if (gInt && gInt->GetN() > 0) {
        double x_dummy, sigma_int;
        gInt->GetPoint(0, x_dummy, sigma_int);

        // These are already error deltas, not absolute lower/upper sigma positions.
        double err_lo = gInt->GetErrorYlow(0);
        double err_hi = gInt->GetErrorYhigh(0);

        // Shaded band
        double xlo = bE_Q2.front() * 0.9;
        double xhi = bE_Q2.back() * 1.05;

        TBox *band = new TBox(xlo, sigma_int - err_lo,
                              xhi, sigma_int + err_hi);
        band->SetFillColorAlpha(kRed-4, 0.25);
        band->SetLineColor(kRed+1);
        band->SetLineWidth(1);
        band->Draw("SAME");

        TLine *lmid = new TLine(xlo, sigma_int, xhi, sigma_int);
        lmid->SetLineColor(kRed+1);
        lmid->SetLineWidth(2);
        lmid->SetLineStyle(2);
        lmid->Draw("SAME");

        // Redraw points on top
        gQ2->Draw("P Z SAME");

        // Legend
        TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.038);
        leg->AddEntry(gQ2,  Form("#sigma_{%s} vs Q^{2}", label), "lpe");
        leg->AddEntry(lmid, Form("Integrated #sigma_{%s} = %.3f^{+%.3f}_{-%.3f} GeV/c",
                                 label, sigma_int, err_hi, err_lo), "l");
        leg->Draw();

        // Numeric label in upper-left
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.038);
        tex.SetTextColor(kRed+1);
        tex.DrawLatex(0.16, 0.85,
                      Form("Integrated: #sigma_{%s} = %.3f^{+%.3f}_{-%.3f} GeV/c",
                           label, sigma_int, err_hi, err_lo));
    }

    // Page title
    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.048);
    title.SetTextFont(62);
    title.DrawLatex(0.14, 0.92, Form("Centre-of-Mass Momentum Width: #sigma_{%s}", label));

    c->Print(pdfFile, "pdf");
}

// -----------------------------------------------------------------------
// MAIN macro function
// -----------------------------------------------------------------------
void plot_sigmaCM(const char *inFileName  = "/work/clas12/users/nwright/rgm_andrew/build/Ana/Q2_Ana/test.root",
                  const char *outPDFFile  = "sigmaCM_plots.pdf")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetGridColor(kGray);

    // ------------------------------------------------------------------
    // Open input file
    // ------------------------------------------------------------------
    TFile *f = TFile::Open(inFileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << inFileName << "\n";
        return;
    }

    // ------------------------------------------------------------------
    // Load integrated TGraphAsymmErrors
    // ------------------------------------------------------------------
    TGraphAsymmErrors *g_int_x = (TGraphAsymmErrors*)f->Get("sigmacmx_int");
    TGraphAsymmErrors *g_int_y = (TGraphAsymmErrors*)f->Get("sigmacmy_int");
    TGraphAsymmErrors *g_int_z = (TGraphAsymmErrors*)f->Get("sigmacmz_int");
    TGraphAsymmErrors *g_int_T = (TGraphAsymmErrors*)f->Get("sigmacmT_int");

    // ------------------------------------------------------------------
    // Load Q2-binned results (errors stored as deltas, same as sigmacmT_Q2)
    // ------------------------------------------------------------------
    TGraphAsymmErrors *g_Q2_x = (TGraphAsymmErrors*)f->Get("sigmacmx_Q2");
    TGraphAsymmErrors *g_Q2_y = (TGraphAsymmErrors*)f->Get("sigmacmy_Q2");
    TGraphAsymmErrors *g_Q2_z = (TGraphAsymmErrors*)f->Get("sigmacmz_Q2");
    TGraphAsymmErrors *g_Q2_T = (TGraphAsymmErrors*)f->Get("sigmacmT_Q2");

    if (!g_Q2_x) std::cerr << "WARNING: sigmacmx_Q2 not found in file\n";
    if (!g_Q2_y) std::cerr << "WARNING: sigmacmy_Q2 not found in file\n";
    if (!g_Q2_z) std::cerr << "WARNING: sigmacmz_Q2 not found in file\n";
    if (!g_Q2_T) std::cerr << "WARNING: sigmacmT_Q2 not found in file\n";

    // ------------------------------------------------------------------
    // Canvas & PDF output
    // ------------------------------------------------------------------
    int pixelx = 900, pixely = 650;
    TCanvas *c = new TCanvas("cSigmaCM", "sigma_CM vs Q2", pixelx, pixely);
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13);

    // Open PDF
    char openTag[256];
    snprintf(openTag, sizeof(openTag), "%s[", outPDFFile);
    c->Print(openTag);

    // ------------------------------------------------------------------
    // Page 1: sigma_CMx vs Q2
    // ------------------------------------------------------------------
    if (g_Q2_x) {
        g_Q2_x->SetTitle("#sigma_{CMx} vs Q^{2};Q^{2} (GeV^{2}/c^{2});#sigma_{CMx} (GeV/c)");
        drawPage(c, g_Q2_x, g_int_x, "CMx", outPDFFile);
    }

    // ------------------------------------------------------------------
    // Page 2: sigma_CMy vs Q2
    // ------------------------------------------------------------------
    if (g_Q2_y) {
        g_Q2_y->SetTitle("#sigma_{CMy} vs Q^{2};Q^{2} (GeV^{2}/c^{2});#sigma_{CMy} (GeV/c)");
        drawPage(c, g_Q2_y, g_int_y, "CMy", outPDFFile);
    }

    // ------------------------------------------------------------------
    // Page 3: sigma_CMz vs Q2
    // ------------------------------------------------------------------
    if (g_Q2_z) {
        g_Q2_z->SetTitle("#sigma_{CMz} vs Q^{2};Q^{2} (GeV^{2}/c^{2});#sigma_{CMz} (GeV/c)");
        drawPage(c, g_Q2_z, g_int_z, "CMz", outPDFFile);
    }

    // ------------------------------------------------------------------
    // Page 4: sigma_CMT vs Q2
    // ------------------------------------------------------------------
    if (g_Q2_T) {
        g_Q2_T->SetTitle("#sigma_{CMT} vs Q^{2};Q^{2} (GeV^{2}/c^{2});#sigma_{CMT} (GeV/c)");
        drawPage(c, g_Q2_T, g_int_T, "CMT", outPDFFile, true);
    }

    // Close PDF
    char closeTag[256];
    snprintf(closeTag, sizeof(closeTag), "%s]", outPDFFile);
    c->Print(closeTag);

    std::cout << "\nDone. Output written to: " << outPDFFile << "\n";

    f->Close();
}