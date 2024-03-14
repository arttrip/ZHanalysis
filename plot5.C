#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void plot5(){
 
  TFile *fDY=new TFile("histos_a20.root");
  TH2F *h= (TH2F*)fDY->Get("h_fjet1_XbbXccXqq_pt");
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle(" fjet1 p_{T}(GeV) ");
  h->GetYaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD on fjet1");
  h->SetTitle("fjet1: ParticleNetMD_XbbXccXqq_vsQCD  vs Pt ");
  h->Draw("COLZ");
}