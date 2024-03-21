#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include "tdrstyle.C"


void met(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
 
  // Get Histograms from signal file:
  
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


  
  TH1F *h3 =(TH1F*)fsig->Get("met_pt");
  h3->Rebin(4);
  TH1F *h4= (TH1F*)fsig->Get("met_trig_pt");
  h4->Rebin(4);
  TCanvas *c2 = new TCanvas ("c2","c2",900,800);
  gPad->SetGridx();
  gPad->SetGridy();
  TGraphAsymmErrors *Eff2 = new TGraphAsymmErrors();
  Eff2->BayesDivide(h4,h3);
  Eff2->GetXaxis()->SetTitle("MET (GeV)");
  //Eff->GetYaxis()->SetRangeUser(0,1.2);
  Eff2->GetYaxis()->SetTitle(" Efficiency ");
  Eff2->SetTitle("trigger efficiency: total MET  after / before trigger requirement");
  Eff2->SetMarkerStyle(kFullCircle);
  Eff2->Draw("AP");
  TCanvas *ca2 = new TCanvas ("ca2","ca2",900,800);
  ca2->SetLogy();
  h3->SetLineColor(kBlack);
  h3->GetXaxis()->SetTitle("MET (GeV)");
  h3->GetYaxis()->SetTitle("Entries");
  h3->SetTitle("total MET");
  h3->Draw("hist");
  h3->SetFillColor(1);
  h3->SetFillStyle(3004);
  //h2->Scale(1./h2->Integral());// Normalize distribution to 1
  h4->SetLineColor(kMagenta);
  h4->Draw("histsames");
  h4->SetLineWidth(2);
  //h1->Scale(1./h1->Integral()); // Normalize distribution to 1
  TLegend *l2 = new TLegend(0.55,0.7,.8,.9);
  l2->AddEntry(h3,"before trigger ","l");
  l2->AddEntry(h4,"after trigger","l");
  l2->SetTextSize(0.03); 
  l2->Draw("SAME");

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
