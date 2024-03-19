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
  
  TH1F *hsig =(TH1F*)fsig->Get("h_fjet1_pt");
  hsig->Rebin(4);
  TH1F *hbkg= (TH1F*)fsig->Get("h_fj1_pt_mat");
  hbkg->Rebin(4);
  // Draw the histos
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  c1->SetLogy();
  hbkg->SetLineColor(kBlack);
  hbkg->GetXaxis()->SetTitle("p_{T} (GeV)");
  hbkg->GetYaxis()->SetTitle("Entries");
  hbkg->SetTitle("fjet1 pt");
  hbkg->Draw("hist");
  hbkg->Scale(1./hbkg->Integral());// Normalize distribution to 1


  hsig->SetLineColor(kRed);
  hsig->Draw("histsames");
  hsig->Scale(1./hsig->Integral()); // Normalize distribution to 1
  
 
  
 
  TLegend *leg = new TLegend(0.65,0.7,.9,.5);
 
  leg->AddEntry(hsig,"fjet1","l");
  
  leg->AddEntry(hbkg,"fjet1 matched","l");
  leg->SetTextSize(0.04); 
  leg->Draw();

  TH1F *hsig2 =(TH1F*)fsig->Get("h_fjet2_pt");
  hsig2->Rebin(4);
  TH1F *hbkg2= (TH1F*)fsig->Get("h_fj2_pt_mat");
  hbkg2->Rebin(4);
  // Draw the histos
  TCanvas *c2 = new TCanvas ("c2","c2",900,800);
  c2->SetLogy();
  hbkg2->SetLineColor(kBlack);
  hbkg2->GetXaxis()->SetTitle("p_{T} (GeV)");
  hbkg2->GetYaxis()->SetTitle("Entries");
  hbkg2->SetTitle("fjet2 pt");
  hbkg2->Draw("hist");
  hbkg2->Scale(1./hbkg2->Integral());// Normalize distribution to 1


  hsig2->SetLineColor(kRed);
  hsig2->Draw("histsames");
  hsig2->Scale(1./hsig2->Integral()); // Normalize distribution to 1
  
 
  
 
  TLegend *leg2 = new TLegend(0.65,0.7,.9,.5);
 
  leg2->AddEntry(hsig2,"fjet2","l");
  
  leg2->AddEntry(hbkg2,"fjet2 matched","l");
  leg2->SetTextSize(0.04); 
  leg2->Draw();
		       
}
