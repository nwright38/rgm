#include "TCanvas.h"
#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

void fitSigmaCM_recTheta(
    const char *pdfName = "fitSigmaCM_recTheta.pdf",
    const char *dataFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *simFileName = "~/data/RGM_DATA/c12_sim_skim.root"){

    const int n_samples = 2;
    const char *sampleNames[n_samples] = {"Data", "Simulation"};
    const char *inputFileNames[n_samples] = {dataFileName, simFileName};
    TFile *inputFiles[n_samples];
    TTree *inputTrees[n_samples];
    for(int sample = 0; sample < n_samples; sample++){
        inputFiles[sample] = new TFile(inputFileNames[sample],"READONLY");
        inputTrees[sample] = (TTree*)inputFiles[sample]->Get("srcTree");
    }

    const int n_components = 2;
    const char *components[n_components] = {"pCMx", "pCMy"};
    const int colors[n_components] = {kBlue+1, kRed+1};

    const int n_bins = 6;
    //const double bin_edges[n_bins+1] = {0,150,160,180};
    const double bin_edges[n_bins+1] = {0,0.4, 0.6, 0.8, 0.9, 1.0,1.5};

    // Plot labels: update these when changing the binned or sigma variable.
    const char *binVariableTitle = "Var";
    const char *sigmaVariableTitle = "#sigma_{CM}[(GeV/c)]";

    //const int n_bins = 5;

    // const double bin_edges[n_bins+1] = {0,0.4,0.6,0.8,0.9,1.0};

    double angleAllow = 20.;
    TCut avoidGaps = Form("((recPhi*180./TMath::Pi() + 108)<%f ) && ((recPhi*180./TMath::Pi() -15 )< %f) && ((recPhi*180./TMath::Pi() -137) < %f)", angleAllow, angleAllow, angleAllow);

    TCut baseCut = "pCM > 0";
    // baseCut = baseCut && "leadTheta*180./TMath::Pi()< 37.";
    // baseCut = baseCut && avoidGaps;

    TCut weight = "weight_epp";
    baseCut = baseCut*weight;

    TH1D *componentHists[n_samples][n_components][n_bins];
    TF1 *componentFits[n_samples][n_components][n_bins];
    TGraphErrors *sigmaGraphs[n_samples][n_components];

    for(int sample = 0; sample < n_samples; sample++){
      for(int component = 0; component < n_components; component++){
        sigmaGraphs[sample][component] = new TGraphErrors(n_bins);
        sigmaGraphs[sample][component]->SetName(
            Form("sigma_%s_%s_vs_binVariable",components[component],sampleNames[sample]));
        sigmaGraphs[sample][component]->SetTitle(Form("#sigma_{%s}",components[component]));
        sigmaGraphs[sample][component]->SetMarkerStyle(20 + component);
        sigmaGraphs[sample][component]->SetMarkerColor(colors[component]);
        sigmaGraphs[sample][component]->SetLineColor(colors[component]);

        for(int i = 0; i < n_bins; i++){
            componentHists[sample][component][i] = new TH1D(
                Form("%s_%s_binHist_%d",sampleNames[sample],components[component],i),
                Form("%s %s: %.2f < %s < %.2f;%s;Counts",
                     sampleNames[sample],components[component],bin_edges[i],
                     binVariableTitle,bin_edges[i+1],
                     components[component]),
                50,-0.8,0.8);
            componentFits[sample][component][i] = new TF1(
                Form("%s_%s_fitFunc_%d",sampleNames[sample],components[component],i),
                "gaus",-0.25,0.2);
            componentFits[sample][component][i]->SetLineColor(colors[component]);

            TCut binCut = Form(
                "pMiss > %f && "
                "pMiss < %f",
                bin_edges[i],bin_edges[i+1]);

            inputTrees[sample]->Draw(
                Form("%s>>%s_%s_binHist_%d",
                     components[component],sampleNames[sample],components[component],i),
                baseCut*binCut,"goff");
            componentHists[sample][component][i]->Fit(
                componentFits[sample][component][i],"RQ");

            const double binCenter =
                0.5*(bin_edges[i] + bin_edges[i+1]);
            const double sigma = componentFits[sample][component][i]->GetParameter(2);
            const double sigmaError =
                componentFits[sample][component][i]->GetParError(2);

            const double mean = componentFits[sample][component][i]->GetParameter(1);
            const double meanError =
                componentFits[sample][component][i]->GetParError(1);
            sigmaGraphs[sample][component]->SetPoint(i,binCenter,sigma);
            sigmaGraphs[sample][component]->SetPointError(i,0.,sigmaError);
        }
      }
    }

    TCanvas *fitCanvas = new TCanvas("fitCanvas","pCM component fits",1200,700);
    bool firstPdfPage = true;
    for(int sample = 0; sample < n_samples; sample++){
      for(int component = 0; component < n_components; component++){
        fitCanvas->Clear();
        fitCanvas->Divide(3,2);
        for(int i = 0; i < n_bins; i++){
            fitCanvas->cd(i + 1);
            componentHists[sample][component][i]->Draw();
            componentFits[sample][component][i]->Draw("same");
        }
        fitCanvas->Print(firstPdfPage ? Form("%s(",pdfName) : pdfName,"pdf");
        firstPdfPage = false;
      }
    }

    TCanvas *sigmaCanvas =
        new TCanvas("sigmaCanvas","#sigma_{CM} dependence",800,650);
    for(int sample = 0; sample < n_samples; sample++){
        sigmaCanvas->Clear();
        sigmaGraphs[sample][0]->SetTitle(
            Form("%s;%s;%s",sampleNames[sample],binVariableTitle,sigmaVariableTitle));

        double sigmaMin = sigmaGraphs[sample][0]->GetY()[0];
        double sigmaMax = sigmaMin;
        for(int component = 0; component < n_components; component++){
            for(int i = 0; i < n_bins; i++){
                const double y = sigmaGraphs[sample][component]->GetY()[i];
                const double yError = sigmaGraphs[sample][component]->GetErrorY(i);
                sigmaMin = TMath::Min(sigmaMin,y - yError);
                sigmaMax = TMath::Max(sigmaMax,y + yError);
            }
        }
        const double sigmaPadding =
            0.1*TMath::Max(sigmaMax - sigmaMin,
                           TMath::Max(TMath::Abs(sigmaMax),1.e-3));
        sigmaGraphs[sample][0]->SetMinimum(sigmaMin - sigmaPadding);
        sigmaGraphs[sample][0]->SetMaximum(sigmaMax + sigmaPadding);

        sigmaGraphs[sample][0]->Draw("APL");
        sigmaGraphs[sample][1]->Draw("PL same");

        TLegend *legend = new TLegend(0.68,0.72,0.88,0.88);
        legend->AddEntry(sigmaGraphs[sample][0],"pCMx","lp");
        legend->AddEntry(sigmaGraphs[sample][1],"pCMy","lp");
        legend->Draw();
        sigmaCanvas->Print(pdfName,"pdf");
    }

    // Final page: compare sigma_CMx between data and simulation.
    sigmaCanvas->Clear();
    sigmaGraphs[0][1]->SetTitle(
        Form("#sigma_{pCMy}: Data vs Simulation;%s;%s",
             binVariableTitle,sigmaVariableTitle));
    sigmaGraphs[0][1]->SetMarkerStyle(20);
    sigmaGraphs[0][1]->SetMarkerColor(kBlack);
    sigmaGraphs[0][1]->SetLineColor(kBlack);
    sigmaGraphs[1][1]->SetMarkerStyle(24);
    sigmaGraphs[1][1]->SetMarkerColor(kBlue+1);
    sigmaGraphs[1][1]->SetLineColor(kBlue+1);

    double comparisonMin = sigmaGraphs[0][1]->GetY()[0];
    double comparisonMax = comparisonMin;
    for(int sample = 0; sample < n_samples; sample++){
        for(int i = 0; i < n_bins; i++){
            const double y = sigmaGraphs[sample][1]->GetY()[i];
            const double yError = sigmaGraphs[sample][1]->GetErrorY(i);
            comparisonMin = TMath::Min(comparisonMin,y - yError);
            comparisonMax = TMath::Max(comparisonMax,y + yError);
        }
    }
    const double comparisonPadding =
        0.1*TMath::Max(comparisonMax - comparisonMin,
                       TMath::Max(TMath::Abs(comparisonMax),1.e-3));
    sigmaGraphs[0][1]->SetMinimum(comparisonMin - comparisonPadding);
    sigmaGraphs[0][1]->SetMaximum(comparisonMax + comparisonPadding);
    sigmaGraphs[0][1]->Draw("APL");
    sigmaGraphs[1][1]->Draw("PL same");

    TLegend *comparisonLegend = new TLegend(0.68,0.72,0.88,0.88);
    comparisonLegend->AddEntry(sigmaGraphs[0][1],"Data","lp");
    comparisonLegend->AddEntry(sigmaGraphs[1][1],"Simulation","lp");
    comparisonLegend->Draw();
    sigmaCanvas->Print(Form("%s)",pdfName),"pdf");
}
