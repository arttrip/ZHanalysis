#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void plot1(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  TFile *fDY=new TFile("histos_anadyzvv.root");
  
  TH1F *hsig = (TH1F*)fsig->Get("h_fjet2_btagXbbXccXqq");
  //hsig->Rebin(2);
  TH1F *hbkg =(TH1F*)fDY->Get("h_fjet2_btagXbbXccXqq");
  //hbkg->Rebin(2);
  std::cout << "signal :entries underflow: " << hsig->Integral(hsig->GetBin(1), hsig->GetBin(20)) << std::endl;	  
  std::cout << "z->vv :: entries underflow: " << hbkg->Integral(hbkg->GetBin(1), hbkg->GetBin(20)) << std::endl;	
  std::cout << "signal :entries overflow: " << hsig->Integral(hsig->GetBin(121), hsig->GetBin(140)) << std::endl;	  
  std::cout << "z->vv :: entries overflow: " << hbkg->Integral(hbkg->GetBin(121), hbkg->GetBin(140)) << std::endl;
  std::cout << "signal :entries in interval(0,1): " << hsig->Integral(hsig->GetBin(21), hsig->GetBin(120)) << std::endl;	
  std::cout << "z->vv :entries in interval(0,1): " << hbkg->Integral(hbkg->GetBin(21), hbkg->GetBin(120)) << std::endl;
  std::cout << "z->vv :entries in bin 130(1.1): " << hbkg->GetBinContent(130) << std::endl;
  // Draw the histos
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  c1->SetLogy();

  hsig->SetLineColor(kRed);
  hsig->SetFillColor(kRed);
  hsig->SetFillStyle(3001);
  hsig->GetXaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD");
  hsig->GetYaxis()->SetTitle("Entries");
  hsig->SetTitle("ParticleNetMD_XbbXccXqqvsQCD on fjet1  ");
  hsig->Draw("hist");
  hsig->Scale(1./hsig->Integral());// Normalize distribution to 1
  
  hbkg->SetLineColor(kBlack);
  hbkg->SetFillColor(kBlack);
  hbkg->SetFillStyle(3001);
  hbkg->Draw("histsames");
  hbkg->Scale(1./hbkg->Integral());
  
  TLegend *leg = new TLegend(0.65,0.7,.9,.5);
  leg->AddEntry(hsig,"ZH(m20) signal","l");
  leg->AddEntry(hbkg,"z->vv","l");
  leg->SetTextSize(0.04); 
  leg->Draw();
	
}