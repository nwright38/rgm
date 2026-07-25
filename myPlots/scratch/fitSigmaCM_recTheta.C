#include "TCanvas.h"
#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

void fitSigmaCM_recTheta(const char *pdfName = "fitSigmaCM_recTheta.pdf"){

    TFile *f_in = new TFile("~/data/RGM_DATA/c12_src_skim.root","READONLY");
    TTree *t_in = (TTree*)f_in->Get("srcTree");

    const int n_components = 2;
    const char *components[n_components] = {"pCMx", "pCMy"};
    const int colors[n_components] = {kBlue+1, kRed+1};

    const int n_bins = 5;
    const double sin_rt_2_edges[n_bins+1] = {0.6,0.8,.85,0.9,.95,1.0};

    //const int n_bins = 5;

    // const double sin_rt_2_edges[n_bins+1] = {0,0.4,0.6,0.8,0.9,1.0};

    TCut baseCut = "pCM > 0";
    // baseCut = baseCut && "leadTheta*180./TMath::Pi()< 37.";

    TCut weight = "weight_epp";
    baseCut = baseCut*weight;

    TH1D *componentHists[n_components][n_bins];
    TF1 *componentFits[n_components][n_bins];
    TGraph *sigmaGraphs[n_components];

    for(int component = 0; component < n_components; component++){
        sigmaGraphs[component] = new TGraph(n_bins);
        sigmaGraphs[component]->SetName(Form("sigma_%s_vs_sin2theta",components[component]));
        sigmaGraphs[component]->SetTitle(Form("#sigma_{%s}^{2}",components[component]));
        sigmaGraphs[component]->SetMarkerStyle(20 + component);
        sigmaGraphs[component]->SetMarkerColor(colors[component]);
        sigmaGraphs[component]->SetLineColor(colors[component]);

        for(int i = 0; i < n_bins; i++){
            componentHists[component][i] = new TH1D(
                Form("%s_binHist_%d",components[component],i),
                Form("%s: %.2f < sin^{2}(#theta_{recoil}) < %.2f;%s;Counts",
                     components[component],sin_rt_2_edges[i],sin_rt_2_edges[i+1],
                     components[component]),
                50,-0.8,0.8);
            componentFits[component][i] = new TF1(
                Form("%s_fitFunc_%d",components[component],i),"gaus",-0.2,0.2);
            componentFits[component][i]->SetLineColor(colors[component]);

            TCut binCut = Form(
                "TMath::Power(sin(pMissTheta),2) > %f && "
                "TMath::Power(sin(pMissTheta),2) < %f",
                sin_rt_2_edges[i],sin_rt_2_edges[i+1]);

            t_in->Draw(Form("%s>>%s_binHist_%d",
                            components[component],components[component],i),
                       baseCut && binCut,"goff");
            componentHists[component][i]->Fit(componentFits[component][i],"RQ");

            const double binCenter =
                0.5*(sin_rt_2_edges[i] + sin_rt_2_edges[i+1]);
            sigmaGraphs[component]->SetPoint(
                i,binCenter,componentFits[component][i]->GetParameter(2)*componentFits[component][i]->GetParameter(2));
        }
    }

    TCanvas *fitCanvas = new TCanvas("fitCanvas","pCM component fits",1200,700);
    for(int component = 0; component < n_components; component++){
        fitCanvas->Clear();
        fitCanvas->Divide(3,2);
        for(int i = 0; i < n_bins; i++){
            fitCanvas->cd(i + 1);
            componentHists[component][i]->Draw();
            componentFits[component][i]->Draw("same");
        }
        fitCanvas->Print(
            component == 0 ? Form("%s(",pdfName) : pdfName,"pdf");
    }

    TCanvas *sigmaCanvas =
        new TCanvas("sigmaCanvas","#sigma_{CM} vs sin^{2}(#theta_{recoil})",800,650);
    sigmaGraphs[0]->SetTitle(
        ";sin^{2}(#theta_{recoil});#sigma_{CM}^{2} [GeV/c]^{2}");

    double sigmaMin = sigmaGraphs[0]->GetY()[0];
    double sigmaMax = sigmaMin;
    for(int component = 0; component < n_components; component++){
        for(int i = 0; i < n_bins; i++){
            sigmaMin = TMath::Min(sigmaMin,sigmaGraphs[component]->GetY()[i]);
            sigmaMax = TMath::Max(sigmaMax,sigmaGraphs[component]->GetY()[i]);
        }
    }
    const double sigmaPadding =
        0.1*TMath::Max(sigmaMax - sigmaMin,TMath::Max(TMath::Abs(sigmaMax),1.e-3));
    sigmaGraphs[0]->SetMinimum(sigmaMin - sigmaPadding);
    sigmaGraphs[0]->SetMaximum(sigmaMax + sigmaPadding);

    sigmaGraphs[0]->Draw("APL");
    sigmaGraphs[1]->Draw("PL same");

    TLegend *legend = new TLegend(0.68,0.72,0.88,0.88);
    legend->AddEntry(sigmaGraphs[0],"pCMx","lp");
    legend->AddEntry(sigmaGraphs[1],"pCMy","lp");
    legend->Draw();

    sigmaCanvas->Print(Form("%s)",pdfName),"pdf");
}
