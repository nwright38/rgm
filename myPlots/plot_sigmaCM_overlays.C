// plot_sigmaCM_overlays.C
//
// Usage:
//   root -l -q 'plot_sigmaCM_overlays.C("output.root","sigmaCM_overlays.pdf")'
//
// Makes:
//   Page 1: integrated over Q2, 2x2 panels for pcmx, pcmy, pcmz, pcmT
//   Pages 2...: one 2x2 page per Q2 bin, with all sigma_CM sim histograms
//               drawn with HIST and data overlaid with points/errors.

#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"

using std::string;
using std::vector;

TH1D* getHistClone(TFile* f, const string& name, const string& cloneName)
{
    TH1D* h = nullptr;
    f->GetObject(name.c_str(), h);

    if (!h) {
        std::cerr << "Missing histogram: " << name << std::endl;
        return nullptr;
    }

    TH1D* hc = (TH1D*)h->Clone(cloneName.c_str());
    hc->SetDirectory(nullptr);
    return hc;
}

void drawSimOverlayWithData(TFile* f,
                            const string& dataName,
                            const string& simPrefix,
                            int nSigma,
                            bool normalizeSimToData,
                            const string& title,
                            const string& xTitle)
{
    TH1D* hData = getHistClone(f, dataName, dataName + "_clone");
    if (!hData) return;

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.7);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);
    hData->SetLineWidth(2);

    vector<TH1D*> sims;
    double ymax = hData->GetMaximum();

    for (int j = 0; j < nSigma; j++) {
        string simName = simPrefix + "_" + std::to_string(j);
        TH1D* hSim = getHistClone(f, simName, simName + "_clone");
        if (!hSim) continue;

        if (normalizeSimToData && hSim->Integral() > 0 && hData->Integral() > 0) {
            hSim->Scale(hData->Integral() / hSim->Integral());
        }

        int color = TColor::GetColorPalette(
            int((double)j / std::max(1, nSigma - 1) * 255.0)
        );

        hSim->SetLineColor(color);
        hSim->SetLineWidth(1);
        hSim->SetFillStyle(0);

        ymax = std::max(ymax, hSim->GetMaximum());
        sims.push_back(hSim);
    }

    bool drewSomething = false;

    for (size_t i = 0; i < sims.size(); i++) {
        sims[i]->SetTitle(title.c_str());
        sims[i]->GetXaxis()->SetTitle(xTitle.c_str());
        sims[i]->GetYaxis()->SetTitle("Counts");
        sims[i]->SetMaximum(1.25 * ymax);

        if (i == 0) {
            sims[i]->Draw("HIST");
            drewSomething = true;
        } else {
            sims[i]->Draw("HIST SAME");
        }
    }

    if (drewSomething) {
        hData->Draw("E SAME");
    } else {
        hData->SetTitle(title.c_str());
        hData->GetXaxis()->SetTitle(xTitle.c_str());
        hData->GetYaxis()->SetTitle("Counts");
        hData->SetMaximum(1.25 * ymax);
        hData->Draw("E");
    }

    TLegend* leg = new TLegend(0.58, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    if (!sims.empty()) {
        leg->AddEntry(sims.front(), "sim #sigma_{CM} scan", "l");
    }
    leg->AddEntry(hData, "data", "lep");
    leg->Draw();
}

void plot_sigmaCM_overlays(const char* inputFile = "/work/clas12/users/nwright/rgm_andrew/build/Ana/Q2_Ana/he_sigmaCM.root",
                           const char* outputPdf = "sigmaCM_overlays.pdf",
                           int nSigma = 100,
                           int nQ2 = 7,
                           bool normalizeSimToData = true)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open file: " << inputFile << std::endl;
        return;
    }

    vector<string> comp = {"pcmx", "pcmy", "pcmz", "pcmT"};
    vector<string> pretty = {
        "p_{CM}^{x}",
        "p_{CM}^{y}",
        "p_{CM}^{z}",
        "p_{CM}^{T}"
    };
    vector<string> xTitles = {
        "p_{CM}^{x} [GeV/c]",
        "p_{CM}^{y} [GeV/c]",
        "p_{CM}^{z} [GeV/c]",
        "p_{CM}^{T} [GeV/c]"
    };

    double q2Edges[8] = {1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0};

    TCanvas* c = new TCanvas("c", "sigma_CM overlays", 1200, 900);

    // Open multipage PDF
    c->Print((string(outputPdf) + "[").c_str());

    // ------------------------------------------------------------------
    // Page 1: integrated over Q2
    // ------------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int k = 0; k < 4; k++) {
        c->cd(k + 1);

        string dataName =  comp[k] + "_epp";
        string simPrefix = "h_" + comp[k] + "_epp_simSCM";

        drawSimOverlayWithData(f,
                               dataName,
                               simPrefix,
                               nSigma,
                               normalizeSimToData,
                               pretty[k] + ", integrated over Q^{2}",
                               xTitles[k]);
    }

    c->Print(outputPdf);

    // ------------------------------------------------------------------
    // Pages 2...: one 2x2 canvas per Q2 bin
    // ------------------------------------------------------------------
    for (int iq = 0; iq < nQ2; iq++) {
        c->Clear();
        c->Divide(2, 2);

        for (int k = 0; k < 4; k++) {
            c->cd(k + 1);

            string dataName =
                "h_" + comp[k] + "_epp_SRC_Q2_" + std::to_string(iq);

            string simPrefix =
                "h_" + comp[k] + "_epp_SRC_simSCM_Q2";

            // The sim histograms are named:
            // h_pcmx_epp_SRC_simSCM_Q2_<sigmaIndex>_<q2Index>
            // so here we need a custom prefix that includes the sigma index later.
            //
            // To reuse drawSimOverlayWithData, make a temporary wrapper pattern:
            // simPrefix + "_<j>_<iq>"
            //
            // The function expects simPrefix + "_<j>", so pass:
            // h_pcmx_epp_SRC_simSCM_Q2 and manually draw below instead.
            TH1D* hData = getHistClone(f, dataName, dataName + "_clone");
            if (!hData) continue;

            hData->SetMarkerStyle(20);
            hData->SetMarkerSize(0.7);
            hData->SetLineColor(kBlack);
            hData->SetMarkerColor(kBlack);
            hData->SetLineWidth(2);

            vector<TH1D*> sims;
            double ymax = hData->GetMaximum();

            for (int j = 0; j < nSigma; j++) {
                string simName =
                    "h_" + comp[k] + "_epp_SRC_simSCM_Q2_" +
                    std::to_string(j) + "_" + std::to_string(iq);

                TH1D* hSim = getHistClone(f, simName, simName + "_clone");
                if (!hSim) continue;

                if (normalizeSimToData && hSim->Integral() > 0 && hData->Integral() > 0) {
                    hSim->Scale(hData->Integral() / hSim->Integral());
                }

                int color = TColor::GetColorPalette(
                    int((double)j / std::max(1, nSigma - 1) * 255.0)
                );

                hSim->SetLineColor(color);
                hSim->SetLineWidth(1);
                hSim->SetFillStyle(0);

                ymax = std::max(ymax, hSim->GetMaximum());
                sims.push_back(hSim);
            }

            string title = pretty[k] + Form(", %.2f < Q^{2} < %.2f",
                                            q2Edges[iq], q2Edges[iq + 1]);

            bool drewSomething = false;

            for (size_t is = 0; is < sims.size(); is++) {
                sims[is]->SetTitle(title.c_str());
                sims[is]->GetXaxis()->SetTitle(xTitles[k].c_str());
                sims[is]->GetYaxis()->SetTitle("Counts");
                sims[is]->SetMaximum(1.25 * ymax);

                if (is == 0) {
                    sims[is]->Draw("HIST");
                    drewSomething = true;
                } else {
                    sims[is]->Draw("HIST SAME");
                }
            }

            if (drewSomething) {
                hData->Draw("E SAME");
            } else {
                hData->SetTitle(title.c_str());
                hData->GetXaxis()->SetTitle(xTitles[k].c_str());
                hData->GetYaxis()->SetTitle("Counts");
                hData->SetMaximum(1.25 * ymax);
                hData->Draw("E");
            }

            TLegend* leg = new TLegend(0.58, 0.70, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            if (!sims.empty()) {
                leg->AddEntry(sims.front(), "sim #sigma_{CM} scan", "l");
            }
            leg->AddEntry(hData, "data", "lep");
            leg->Draw();
        }

        c->Print(outputPdf);
    }

    // Close multipage PDF
    c->Print((string(outputPdf) + "]").c_str());

    f->Close();

    std::cout << "Wrote " << outputPdf << std::endl;
}