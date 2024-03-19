#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void plot1(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  TFile *fDY=new TFile("histos_zvv.root");
  
  TH1F *hsig = (TH1F*)fsig->Get("h_fjet1_btagXbbXccXqq");
  //hsig->Rebin(2);
  TH1F *hbkg =(TH1F*)fDY->Get("h_fjet1_btagXbbXccXqq");
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

 

   
  
  TH1F *hsig2 = (TH1F*)fsig->Get("h_fjet2_btagXbbXccXqq");
  //hsig->Rebin(2);
  TH1F *hbkg2 =(TH1F*)fDY->Get("h_fjet2_btagXbbXccXqq");
  //hbkg->Rebin(2);
 
  // Draw the histos
  TCanvas *c2 = new TCanvas ("c2","c2",900,800);
  c2->SetLogy();

  hsig2->SetLineColor(kRed);
  hsig2->SetFillColor(kRed);
  hsig2->SetFillStyle(3001);
  hsig2->GetXaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD");
  hsig2->GetYaxis()->SetTitle("Entries");
  hsig2->SetTitle("ParticleNetMD_XbbXccXqqvsQCD on fjet2  ");
  hsig2->Draw("hist");
  hsig2->Scale(1./hsig2->Integral());// Normalize distribution to 1
  
  hbkg2->SetLineColor(kBlack);
  hbkg2->SetFillColor(kBlack);
  hbkg2->SetFillStyle(3001);
  hbkg2->Draw("histsames");
  hbkg2->Scale(1./hbkg2->Integral());
  
  TLegend *leg2 = new TLegend(0.65,0.7,.9,.5);
  leg2->AddEntry(hsig2,"ZH(m20) signal","l");
  leg2->AddEntry(hbkg2,"z->vv","l");
  leg2->SetTextSize(0.04); 
  leg2->Draw();

   
  
  TH1F *hsig3 = (TH1F*)fsig->Get("h_fjet1_btagXbb");
  //hsig->Rebin(2);
  TH1F *hbkg3 =(TH1F*)fDY->Get("h_fjet1_btagXbb");
  //hbkg->Rebin(2);
 
  // Draw the histos
  TCanvas *c3 = new TCanvas ("c3","c3",900,800);
  c3->SetLogy();

  hsig3->SetLineColor(kRed);
  hsig3->SetFillColor(kRed);
  hsig3->SetFillStyle(3001);
  hsig3->GetXaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD");
  hsig3->GetYaxis()->SetTitle("Entries");
  hsig3->GetYaxis()->SetRange(0,1);
  hsig3->SetTitle("ParticleNetMD_XbbvsQCD on fjet1  ");
  hsig3->Draw("hist");
  hsig3->Scale(1./hsig3->Integral());// Normalize distribution to 1
  
  hbkg3->SetLineColor(kBlack);
  hbkg3->SetFillColor(kBlack);
  hbkg3->SetFillStyle(3001);
  hbkg3->Draw("histsames");
  hbkg3->Scale(1./hbkg3->Integral());
  
  TLegend *leg3 = new TLegend(0.65,0.7,.9,.5);
  leg3->AddEntry(hsig3,"ZH(m20) signal","l");
  leg3->AddEntry(hbkg3,"z->vv","l");
  leg3->SetTextSize(0.04); 
  leg3->Draw();


   TH1F *hsig4 = (TH1F*)fsig->Get("h_fjet2_btagXbb");
  //hsig->Rebin(2);
  TH1F *hbkg4 =(TH1F*)fDY->Get("h_fjet2_btagXbb");
  //hbkg->Rebin(2);
 
  // Draw the histos
  TCanvas *c4 = new TCanvas ("c4","c4",900,800);
  c4->SetLogy();

  hsig4->SetLineColor(kRed);
  hsig4->SetFillColor(kRed);
  hsig4->SetFillStyle(3001);
  hsig4->GetXaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD");
  hsig4->GetYaxis()->SetTitle("Entries");
  hsig4->GetYaxis()->SetRange(0,1);
  hsig4->SetTitle("ParticleNetMD_XbbvsQCD on fjet2  ");
  hsig4->Draw("hist");
  hsig4->Scale(1./hsig4->Integral());// Normalize distribution to 1
  
  hbkg4->SetLineColor(kBlack);
  hbkg4->SetFillColor(kBlack);
  hbkg4->SetFillStyle(3001);
  hbkg4->Draw("histsames");
  hbkg4->Scale(1./hbkg4->Integral());
  
  TLegend *leg4 = new TLegend(0.65,0.7,.9,.5);
  leg4->AddEntry(hsig3,"ZH(m20) signal","l");
  leg4->AddEntry(hbkg3,"z->vv","l");
  leg4->SetTextSize(0.04); 
  leg4->Draw();
}
void plot2(){
 
  TFile *fDY=new TFile("histos_a20.root");
  TH2F *h= (TH2F*)fDY->Get("h_fj1_vs_fj2_btagXbb");
  TCanvas *c = new TCanvas ("c","c",900,800);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle(" ParticleNetMD_Xbb_vsQCD on fjet1  ");
  h->GetYaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD on fjet2");
  h->SetTitle("ParticleNetMD_Xbb_vsQCD : on fjet1 vs on fjet2 (ZH m20 signal) ");
  h->Draw("COLZ");


 
  TH2F *h2= (TH2F*)fDY->Get("h_fj1_fj2_btagXbbXccXqq");
  TCanvas *c2 = new TCanvas ("c2","c2",900,800);
  h2->SetLineColor(kBlack);
  h2->GetXaxis()->SetTitle(" ParticleNetMD_XbbXccXqq_vsQCD on fjet1  ");
  h2->GetYaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD on fjet2");
  h2->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD : on fjet1 vs on fjet2 (ZH m20 signal) ");
  h2->Draw("COLZ");
}