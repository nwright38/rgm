void fitSigmaCM_recTheta(){


    TFile *f_in = new TFile("~/data/RGM_DATA/c12_src_skim.root","READONLY");
    TTree *t_in = (TTree*)f_in->Get("srcTree");

    const int n_bins = 5;
    const double sin_rt_2_edges[n_bins+1] = {0,0.4,0.6,0.8,0.9,1.0};

    TCut currCut = "pCM > 0";

    TH1D *binHists[n_bins];
    TF1 *fitFuncs[n_bins];
    for(int i = 0; i < n_bins; i++){
        binHists[i] = new TH1D(Form("binHist_%d",i),Form("binHist_%d",i),50,-0.8,0.8);
        fitFuncs[i] = new TF1(Form("fitFunc_%d",i),"gaus",-0.2,0.2);

        TCut binCut = Form("TMath::Power(sin(recTheta),2) > %f && TMath::Power(sin(recTheta),2) < %f",sin_rt_2_edges[i],sin_rt_2_edges[i+1]);
        currCut = currCut && binCut;

        t_in->Draw(Form("pCMy>>binHist_%d",i),currCut);
        binHists[i]->Fit(Form("fitFunc_%d",i),"R");

    }

    
    
}