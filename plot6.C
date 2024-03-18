#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//pt of fjets 
void plot6(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  
  // Get Histograms from signal file:
  
  TH1F *hsig =(TH1F*)fsig->Get("h_fjet2_pt");
  hsig->Rebin(4);
  TH1F *hbkg= (TH1F*)fsig->Get("h_fj2_pt_mat");
  hbkg->Rebin(4);
  // Draw the histos
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  hbkg->SetLineColor(kBlack);
  hbkg->GetXaxis()->SetTitle("p_{T}(GeV)");
  hbkg->GetYaxis()->SetTitle("Entries");
  hbkg->SetTitle("fjet1 pt");
  hbkg->Draw("hist");
  hbkg->Scale(1./hbkg->Integral());// Normalize distribution to 1


  hsig->SetLineColor(kRed);
  hsig->Draw("histsames");
  hsig->Scale(1./hsig->Integral()); // Normalize distribution to 1
  
 
  
 
  TLegend *leg = new TLegend(0.65,0.7,.9,.5);
 
  leg->AddEntry(hsig,"fjet2","l");
  
  leg->AddEntry(hbkg,"fjet2 matched","l");
  leg->SetTextSize(0.04); 
  leg->Draw();
		       
}
