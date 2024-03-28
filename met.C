#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include "tdrstyle.C"


void met(){
  // Open the input files:
  TFile *fsig=new TFile("histos_anal20.root");
  TFile *fbkg1=new TFile("histos_zvv100.root");
  TFile *fbkg2=new TFile("histos_zvv200.root");
  TFile *fbkg3=new TFile("histos_tt_had.root");
  
  // met from signal, z->vv, TTbar->had
  TH1F *hjs1 =(TH1F*)fsig->Get("met_pt");
  TH1F *hjb1= (TH1F*)fbkg1->Get("met_pt");
  TH1F *hjb2= (TH1F*)fbkg2->Get("met_pt");
  hjb1->Add(hjb2);
  TH1F *hjb3= (TH1F*)fbkg3->Get("met_pt");
  setTDRStyle();

  TCanvas *cas1 = new TCanvas ("cas1","cas1",900,800);
  cas1->SetLogy();
  hjs1->SetLineColor(kBlack);
  hjs1->SetLineWidth(2);
  hjs1->SetFillColor(kBlack);
  hjs1->SetFillStyle(3004);
  hjs1->GetXaxis()->SetTitle(" MET (GeV)");
  hjs1->GetYaxis()->SetTitle("Entries");
  hjs1->SetTitle(" total MET ");
  hjs1->Draw("hist");
  hjs1->Scale(1./hjs1->Integral());// Normalize distribution to 1

  hjb3->SetLineColor(kGreen);
  hjb3->Draw("histsames");
  hjb3->Scale(1./hjb3->Integral());
  
  hjb1->SetLineColor(kRed);
  hjb1->Draw("histsames");
  hjb1->Scale(1./hjb1->Integral()); // Normalize distribution to 1
  TLegend *le1 = new TLegend(0.55,0.7,.8,.9);
  le1->AddEntry(hjs1,"ZH (20) signal ","l");
  le1->AddEntry(hjb1,"Z #rightarrow #nu#nu","l");
  le1->AddEntry(hjb3,"TTbar  (had)","l");
  le1->SetTextSize(0.03); 
  le1->Draw("SAME");






  
 
  //met before -after trigger
  
  TH1F *h1 =(TH1F*)fsig->Get("met_vv_pt");
  h1->Rebin(4);
  TH1F *h2= (TH1F*)fsig->Get("met_vv_trig_pt");
  h2->Rebin(4);
  
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  gPad->SetGridx();
  gPad->SetGridy();
  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  Eff->BayesDivide(h2,h1);
  Eff->GetXaxis()->SetTitle("MET (GeV)");
  //Eff->GetYaxis()->SetRangeUser(0,1.2);
  Eff->GetYaxis()->SetTitle(" Efficiency ");
  Eff->SetTitle("trigger efficiency: MET from Z#rightarrow #nu#nu  after / before trigger requirement");
  Eff->SetMarkerStyle(kFullCircle);
  Eff->Draw("AP");
  
  TCanvas *ca1 = new TCanvas ("ca1","ca1",900,800);
  ca1->SetLogy();
  h1->SetLineColor(kBlack);
  h1->GetXaxis()->SetTitle("MET (GeV)");
  h1->GetYaxis()->SetTitle("Entries");
  h1->SetTitle("MET from Z#rightarrow #nu#nu ");
  h1->Draw("hist");
  h1->SetFillColor(1);
  h1->SetFillStyle(3004);
  //h2->Scale(1./h2->Integral());// Normalize distribution to 1
  h2->SetLineColor(kMagenta);
  h2->Draw("histsames");
  h2->SetLineWidth(2);
  //h1->Scale(1./h1->Integral()); // Normalize distribution to 1
  TLegend *l1 = new TLegend(0.55,0.7,.8,.9);
  l1->AddEntry(h1,"before trigger ","l");
  l1->AddEntry(h2,"after trigger","l");
  l1->SetTextSize(0.03); 
  l1->Draw("SAME");


  
 




  
  /////////////
  TH1F *hq =(TH1F*)fsig->Get("met_qq_l_pt");
	TCanvas *c3 = new TCanvas ("c3","c3",900,800);
  hq->SetLineColor(kBlack);
  hq->GetXaxis()->SetTitle("MET (GeV)");
  hq->GetYaxis()->SetTitle("Entries");
  hq->SetTitle("MET from Z#rightarrow qq ");
  hq->Draw("hist");
  hq->SetFillColor(1);
  hq->SetFillStyle(3004);	

  TH1F *hb =(TH1F*)fsig->Get("met_bb_pt");
	TCanvas *c4 = new TCanvas ("c4","c4",900,800);
  hb->SetLineColor(kBlack);
  hb->GetXaxis()->SetTitle("MET (GeV)");
  hb->GetYaxis()->SetTitle("Entries");
  hb->SetTitle("MET from Z#rightarrow bb ");
  hb->Draw("hist");
  hb->SetFillColor(1);
  hb->SetFillStyle(3004);       
}
