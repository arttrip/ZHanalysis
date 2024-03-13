#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void plot3(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  //TFile *fsig2=new TFile("histos_ana20.root");
  TFile *fDY=new TFile("histos_anadyzvv.root");
  // Get Histograms from signal file:
  
  TH1F *hsig_mn_eta = new TH1F("","",100,0,8);
  hsig_mn_eta=(TH1F*)fsig->Get("fjet_mult_after2");
 // hsig_mn_eta->Rebin(2);
 
  TH1F *hbkg_mn_eta=new TH1F("","",100,80,100);
  hbkg_mn_eta= (TH1F*)fDY->Get("fjet_mult_after2");
//  hbkg_mn_eta->Rebin(2);

 // TH1F *hsig2_mn_eta=new TH1F("","",100,80,100);
 // hsig2_mn_eta= (TH1F*)fsig2->Get("jet_cc_mult3");
  //hsig2_mn_eta->Rebin(4);
  // Draw the histos
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  hbkg_mn_eta->SetLineColor(kBlack);
  hbkg_mn_eta->GetXaxis()->SetTitle("no of fjets");
 // hbkg_mn_eta->GetXaxis()->SetRangeUser(0,150);
  hbkg_mn_eta->GetYaxis()->SetTitle("Entries");
  hbkg_mn_eta->SetTitle("fjet mult (step 3 of the cut flow)");
  hbkg_mn_eta->Draw("hist");
  hbkg_mn_eta->Scale(1./hbkg_mn_eta->Integral());// Normalize distribution to 1
  
 // hsig_mn_eta->SetLineColor(kRed);
 // hsig_mn_eta->Draw("histsames");
 // hsig_mn_eta->Scale(1./hsig_mn_eta->Integral());

  hsig_mn_eta->SetLineColor(kRed);
  hsig_mn_eta->Draw("histsames");
  hsig_mn_eta->Scale(1./hsig_mn_eta->Integral()); // Normalize distribution to 1
  
 
  
 
  TLegend *leg = new TLegend(0.65,0.7,.9,.5);
 
  leg->AddEntry(hsig_mn_eta,"ZH(m20) signal ","l");
  
  leg->AddEntry(hbkg_mn_eta,"z->vv","l");
  // leg->AddEntry(hsig2_mn_eta," step 3","l");
  leg->SetTextSize(0.04); 
  leg->Draw();
		       
}
