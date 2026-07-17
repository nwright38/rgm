void drawSigCMres(){

    double FIT_RANGE_SIG[2] = {-.2,.2};
    double FIT_RANGE_RES[2] = {-.07,.07};

    TCut cdCut = "leadTheta*180./TMath::Pi() > 45.";
    TCut fdCut = "leadTheta*180./TMath::Pi() < 37.";
    TCut weight = "weight_epp";

    TCut currCut = weight*cdCut;

    TFile *simFile = new TFile("~/data/RGM_DATA/c12_sim_skim.root");
    TFile *dataFile = new TFile("~/data/RGM_DATA/c12_sim_skim.root");

    TTree *simTree = (TTree*)simFile->Get("srcTree"); 
    TTree *dataTree = (TTree*)dataFile->Get("srcTree");

    TH1D *res_x = new TH1D("res_x", "res_x", 50, -0.3, 0.3);
    TH1D *res_y = new TH1D("res_y", "res_y", 50, -0.3, 0.3);

    TF1 *fit_res_x = new TF1("fit_res_x", "gaus", FIT_RANGE_RES[0], FIT_RANGE_RES[1]);
    TF1 *fit_res_y = new TF1("fit_res_y", "gaus", FIT_RANGE_RES[0], FIT_RANGE_RES[1]);

   

    TH1D *pcm_x = new TH1D("pcm_x", "pcm_x", 50, -0.7, 0.7);
    TH1D *pcm_y = new TH1D("pcm_y", "pcm_y", 50, -0.7, 0.7);

    TF1 *fit_x = new TF1("fit_x", "gaus", FIT_RANGE_SIG[0], FIT_RANGE_SIG[1]);
    TF1 *fit_y = new TF1("fit_y", "gaus", FIT_RANGE_SIG[0], FIT_RANGE_SIG[1]);

    simTree->Project("res_x", "pCMx-pCMx_truth", currCut);
    simTree->Project("res_y", "pCMy-pCMy_truth", currCut);

    res_x->Fit("fit_res_x", "RQ");
    res_y->Fit("fit_res_y", "RQ");

    dataTree->Project("pcm_x", "pCMx", currCut);
    dataTree->Project("pcm_y", "pCMy", currCut);

    pcm_x->Fit("fit_x", "RQ");
    pcm_y->Fit("fit_y", "RQ");

    double res_x_sigma = fit_res_x->GetParameter(2);
    double res_y_sigma = fit_res_y->GetParameter(2);

    double sig_x_sigma = fit_x->GetParameter(2);
    double sig_y_sigma = fit_y->GetParameter(2);

    double res_x_sigma_err = fit_res_x->GetParError(2);
    double res_y_sigma_err = fit_res_y->GetParError(2);

    double sig_x_sigma_err = fit_x->GetParError(2);
    double sig_y_sigma_err = fit_y->GetParError(2);

    cout << "RES X SIGMA = " << res_x_sigma << " +/- " << res_x_sigma_err << endl;
    cout << "RES Y SIGMA = " << res_y_sigma << " +/- " << res_y_sigma_err << endl;
    cout << "SIG X SIGMA = " << sig_x_sigma << " +/- " << sig_x_sigma_err << endl;
    cout << "SIG Y SIGMA = " << sig_y_sigma << " +/- " << sig_y_sigma_err << endl;

    double extr_x_sigma = sqrt(sig_x_sigma*sig_x_sigma - res_x_sigma*res_x_sigma);
    double extr_y_sigma = sqrt(sig_y_sigma*sig_y_sigma - res_y_sigma*res_y_sigma);

    double extr_x_sigma_err = sqrt(pow(sig_x_sigma*sig_x_sigma_err/sig_x_sigma,2) + pow(res_x_sigma*res_x_sigma_err/res_x_sigma,2));
    double extr_y_sigma_err = sqrt(pow(sig_y_sigma*sig_y_sigma_err/sig_y_sigma,2) + pow(res_y_sigma*res_y_sigma_err/res_y_sigma,2));
    cout << "\n" << endl;
    cout << "EXTR X SIGMA = " << extr_x_sigma << " +/- " << extr_x_sigma_err << endl;
    cout << "EXTR Y SIGMA = " << extr_y_sigma << " +/- " << extr_y_sigma_err << endl;


    // double corr_x = 0.220071 - 0.15;
    // double corr_y = 0.216312 - 0.15;

    // cout << "\n" << endl;
    // cout << "EXTR X SIGMA (corr) = " << sig_x_sigma - corr_x << endl;
    // cout << "EXTR Y SIGMA (corr) = " << sig_y_sigma - corr_y << endl;
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(2,2);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.06);
    c->cd(1);   
    res_x->Draw();
    fit_res_x->Draw("same");
    text.DrawLatex(0.18, 0.82, Form("#sigma_{res,x} = %.3f #pm %.3f", res_x_sigma, res_x_sigma_err));
    c->cd(2);
    res_y->Draw();
    fit_res_y->Draw("same");
    text.DrawLatex(0.18, 0.82, Form("#sigma_{res,y} = %.3f #pm %.3f", res_y_sigma, res_y_sigma_err));
    c->cd(3);
    pcm_x->Draw();
    fit_x->Draw("same");
    text.DrawLatex(0.18, 0.82, Form("#sigma_{fit,x} = %.3f #pm %.3f", sig_x_sigma, sig_x_sigma_err));
    text.DrawLatex(0.18, 0.74, Form("#sigma_{extr,x} = %.3f #pm %.3f", extr_x_sigma, extr_x_sigma_err));
    c->cd(4);
    pcm_y->Draw();
    fit_y->Draw("same");
    text.DrawLatex(0.18, 0.82, Form("#sigma_{fit,y} = %.3f #pm %.3f", sig_y_sigma, sig_y_sigma_err));
    text.DrawLatex(0.18, 0.74, Form("#sigma_{extr,y} = %.3f #pm %.3f", extr_y_sigma, extr_y_sigma_err));
    c->SaveAs("sigCMres.pdf");
    cout << "saved sigCMres.pdf" << endl;

    return;



}
