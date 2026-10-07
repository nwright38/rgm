


void quickDraw(){


    TTree *tDat = (TTree*) TFile::Open("~/data/RGM_DATA/c12_src_skim_loosest_train.root")->Get("srcTree");
    TTree *tSim = (TTree*) TFile::Open("~/data/RGM_DATA/c12_sim_skim_loosest.root")->Get("srcTree");


    const int nXbBins = 18;
    double XbEdges[nXbBins+1];
    double xB_min = 1.;
    double xB_max = 1.3;

    for(int i = 0; i <= nXbBins; i++){
        XbEdges[i] = xB_min + (xB_max-xB_min)/nXbBins*i;
    }

    TCut cut = "weight_ep*(weight_ep < 150)*(Q2 > 1.5 && leadP > 1. && kMiss > 0.3 && leadP*180./TMath::Pi() > 45. && leadP/qP < 0.95)";
    // int colors[7] = {kRed,kBlue,kGreen+2,kOrange+7,kViolet+7,kCyan+2,kMagenta};

    const int nTrees = 2;
    TTree *t[nTrees] = {tDat,tSim};
    TH1D *h[nTrees][nXbBins];
    double centers[nTrees][nXbBins] = {0};
    double sigmas[nTrees][nXbBins] = {0};

    int nHistBins = 200;
    TF1 *f = new TF1("f","gaus",.8,1.2);
    for(int j = 0; j < nTrees; j++) {
      for (int i = 0; i < nXbBins; i++) {
        h[j][i] = new TH1D(Form("h%s_%d", j == 0 ? "Dat" : "Sim", i), Form("Xb > %.02f", XbEdges[i]), nHistBins, 0, 2);
        t[j]->Draw(Form("mMiss>>h%s_%d", j == 0 ? "Dat" : "Sim", i), cut && TCut(Form("xB > %f", XbEdges[i])), "goff");
        h[j][i]->Fit(f, "RQ");
        cout << "Tree " << j << " Xb bin " << i << ": xB > " << XbEdges[i] << '\t' << "Center = " << f->GetParameter(1) << '\t' << "Sigma = " << f->GetParameter(2) << endl;
        centers[j][i] = f->GetParameter(1);
        sigmas[j][i] = f->GetParameter(2);
      }
    }

    TGraphErrors *g[nTrees];
    for(int j = 0; j < nTrees; j++){
        g[j] = new TGraphErrors(nXbBins);
        for(int i = 0; i < nXbBins; i++){
            g[j]->SetPoint(i, 0.5*(XbEdges[i]+XbEdges[i+1]), centers[j][i]);
            g[j]->SetPointError(i, 0.5*(XbEdges[i+1]-XbEdges[i]), sigmas[j][i]);
        }
    }

    TString pdfFile = "quickDraw.pdf";
    TCanvas *c = new TCanvas("c","c",800,600);
    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    c->cd();
    c->Print(pdfFile+"[","pdf");
    
    for(int j = 0; j < nTrees; j++){
        c->Divide(3,3);
        c2->Divide(3,3);
        for(int i = 0; i < nXbBins; i++){
            if( i < 9) {
                c->cd(i+1);
            } else {
                c2->cd(i-9+1);
            }
            h[j][i]->Draw();
            TLine l(centers[j][i],0,centers[j][i],h[j][i]->GetMaximum());
            l.SetLineColor(kRed);
            l.SetLineWidth(2);
            l.Draw();
            TLine l2(centers[j][i]-sigmas[j][i],0,centers[j][i]-sigmas[j][i],h[j][i]->GetMaximum());
            l2.SetLineColor(kBlue);
            l2.SetLineWidth(2);
            l2.Draw();
            TLine l3(centers[j][i]+sigmas[j][i],0,centers[j][i]+sigmas[j][i],h[j][i]->GetMaximum());
            l3.SetLineColor(kBlue);
            l3.SetLineWidth(2);
            l3.Draw();
        }
        c->Print(pdfFile,"pdf");
        c->Clear();
        c2->Print(pdfFile,"pdf");
        c2->Clear();
    }

    c->cd();
    g[0]->SetLineColor(kRed);
    g[0]->SetMarkerColor(kRed);
    g[0]->SetMarkerStyle(20);
    g[0]->SetTitle("Data vs Sim: mMiss vs xB");
    g[0]->GetXaxis()->SetTitle("xB");
    g[0]->GetYaxis()->SetTitle("mMiss [GeV/c^{2}]");
    g[0]->Draw("AP");
    g[1]->SetLineColor(kBlue);
    g[1]->SetMarkerColor(kBlue);
    g[1]->SetMarkerStyle(21);
    g[1]->Draw("P SAME");

    g[0]->GetYaxis()->SetRangeUser(0.7,1.3);

    c->Print(pdfFile,"pdf");

    c->Clear();
    c->Print(pdfFile+"]","pdf");


}

