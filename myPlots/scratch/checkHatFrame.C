// checkHatFrame.C
//
// Rebuild the (xHat, yHat, zHat) frame used in simpleSRCSkim_archive.cpp
// from saved srcTree branches and validate orthonormality/closure.
//
// Usage:
//   root -l -b -q 'myPlots/scratch/checkHatFrame.C()'
//   root -l -b -q 'myPlots/scratch/checkHatFrame.C("~/data/RGM_DATA/c12_src_skim.root")'
//   root -l -b -q 'myPlots/scratch/checkHatFrame.C("in.root","srcTree","hat_checks.root","hat_checks.pdf")'

#include <cmath>
#include <iostream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

void checkHatFrame(
    const char *inputFileName = "~/data/RGM_DATA/c12_src_skim.root",
    const char *treeName = "srcTree",
    const char *outputRootName = "hat_frame_checks.root",
    const char *outputPdfName = "hat_frame_checks.pdf") {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);

  TFile *inFile = TFile::Open(inputFileName, "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "Could not open input file: " << inputFileName << "\n";
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(inFile->Get(treeName));
  if (!tree) {
    std::cerr << "Could not find tree \"" << treeName << "\" in "
              << inputFileName << "\n";
    inFile->Close();
    return;
  }

  Float_t pMiss = -9.f, pMissTheta = -9.f, pMissPhi = -9.f;
  Float_t qP = -9.f, qTheta = -9.f, qPhi = -9.f;
  Float_t recP = -9.f, recTheta = -9.f, recPhi = -9.f;

  tree->SetBranchAddress("pMiss", &pMiss);
  tree->SetBranchAddress("pMissTheta", &pMissTheta);
  tree->SetBranchAddress("pMissPhi", &pMissPhi);
  tree->SetBranchAddress("qP", &qP);
  tree->SetBranchAddress("qTheta", &qTheta);
  tree->SetBranchAddress("qPhi", &qPhi);
  tree->SetBranchAddress("recP", &recP);
  tree->SetBranchAddress("recTheta", &recTheta);
  tree->SetBranchAddress("recPhi", &recPhi);

  TH1D *h1 = new TH1D("h1", "zHat dot yHat;zHat #upoint yHat;Counts", 200, -1e-12, 1e-12);
  TH1D *h2 = new TH1D("h2", "xHat dot yHat;xHat #upoint yHat;Counts", 200, -1e-12, 1e-12);
  TH1D *h3 = new TH1D("h3", "|yHat|;|yHat|;Counts", 200, 0.999999999, 1.000000001);
  TH1D *h4 = new TH1D("h4", "|pCM|^{2} closure;|pCM|^{2} - (pCMx^{2}+pCMy^{2}+pCMz^{2}) [(GeV/c)^{2}];Counts", 200, -1e-12, 1e-12);
  TH1D *h5 = new TH1D("h5", "pMiss dot xHat;pMiss #upoint xHat [GeV/c];Counts", 200, -1e-12, 1e-12);
  TH1D *h6 = new TH1D("h6", "pMiss dot yHat;pMiss #upoint yHat [GeV/c];Counts", 200, -1e-12, 1e-12);

  Long64_t nEntries = tree->GetEntries();
  Long64_t nUsed = 0;

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    if (pMiss <= 0.f || qP <= 0.f || recP <= 0.f) continue;

    TVector3 pmiss;
    TVector3 q;
    TVector3 prec;

    pmiss.SetMagThetaPhi(pMiss, pMissTheta, pMissPhi);
    q.SetMagThetaPhi(qP, qTheta, qPhi);
    prec.SetMagThetaPhi(recP, recTheta, recPhi);

    if (pmiss.Mag2() == 0. || q.Mag2() == 0.) continue;

    const TVector3 zhat = pmiss.Unit();
    TVector3 yhat = pmiss.Cross(q);
    if (yhat.Mag2() == 0.) continue;
    yhat = yhat.Unit();

    TVector3 xhat = zhat.Cross(yhat);
    if (xhat.Mag2() == 0.) continue;
    xhat = xhat.Unit();

    const TVector3 pcm = pmiss + prec;
    const double pcmx = pcm.Dot(xhat);
    const double pcmy = pcm.Dot(yhat);
    const double pcmz = pcm.Dot(zhat);

    h1->Fill(zhat.Dot(yhat));
    h2->Fill(xhat.Dot(yhat));
    h3->Fill(yhat.Mag());
    h4->Fill(pcm.Mag2() - (pcmx * pcmx + pcmy * pcmy + pcmz * pcmz));
    h5->Fill(pmiss.Dot(xhat));
    h6->Fill(pmiss.Dot(yhat));

    ++nUsed;
  }

  TFile *outFile = TFile::Open(outputRootName, "RECREATE");
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  outFile->Close();

  TCanvas *c = new TCanvas("c_hat_checks", "Hat frame checks", 1200, 800);
  c->Divide(3, 2);
  c->cd(1); h1->Draw();
  c->cd(2); h2->Draw();
  c->cd(3); h3->Draw();
  c->cd(4); h4->Draw();
  c->cd(5); h5->Draw();
  c->cd(6); h6->Draw();
  c->SaveAs(outputPdfName);

  std::cout << "Processed entries: " << nEntries << "\n";
  std::cout << "Used entries:      " << nUsed << "\n";
  std::cout << "Wrote: " << outputRootName << " and " << outputPdfName << "\n";

  inFile->Close();
}
