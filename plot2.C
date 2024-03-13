#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void plot2(){
 
  TFile *fDY=new TFile("histos_a20.root");
  TH2F *h= (TH2F*)fDY->Get("h_fj1_vs_fj2_btagXbb");
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle(" ParticleNetMD_Xbb_vsQCD on fjet1  ");
  //h->GetYaxis()->SetRangeUser(0,150);
  h->GetYaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD on fjet2");
  h->SetTitle("ParticleNetMD_Xbb_vsQCD : on fjet1 vs on fjet2 (ZH m20 signal) ");
  h->Draw("COLZ");
 
  //TLine *line = new TLine(0,0,500,500);
  //line->SetLineColor(kRed);
  //line->Draw();
}