#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void plot3(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  
 
  TH1F *hsig =(TH1F*)fsig->Get("sub_mult1");
 
  TH1F *hbkg= (TH1F*)fsig->Get("sub_mult1_mat");

  // Draw the histos
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  c1->SetLogy();
  hbkg->SetLineColor(kBlack);
  hbkg->GetXaxis()->SetTitle("no of sub-jets");
  hbkg->GetYaxis()->SetTitle("Entries");
  hbkg->SetTitle("sub-jet mult in fjet1");
  hbkg->Draw("hist");
  hbkg->Scale(1./hbkg->Integral());// Normalize distribution to 1


  hsig->SetLineColor(kRed);
  hsig->Draw("histsames");
  hsig->Scale(1./hsig->Integral()); // Normalize distribution to 1
  
 
  
 
  TLegend *leg = new TLegend(0.65,0.7,.9,.5);
 
  leg->AddEntry(hsig,"fjet1 ","l");
  
  leg->AddEntry(hbkg,"fjet1 matched","l");
  leg->SetTextSize(0.04); 
  leg->Draw();
		       
}
