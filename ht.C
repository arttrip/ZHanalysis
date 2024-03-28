#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle.C"

void ht(){
  // Open the input files:
  TFile *fsig=new TFile("histos_anal20.root");
  TFile *fbkg1=new TFile("histos_zvv100.root");
  TFile *fbkg2=new TFile("histos_zvv200.root");
  TFile *fbkg3=new TFile("histos_tt_had.root");

  
  TH1F *hjs1 =(TH1F*)fsig->Get("Ht");
  TH1F *hjb1= (TH1F*)fbkg1->Get("Ht");
  TH1F *hjb2= (TH1F*)fbkg2->Get("Ht");
  hjb1->Add(hjb2);
  TH1F *hjb3= (TH1F*)fbkg3->Get("Ht");
  setTDRStyle();
  TCanvas *cas1 = new TCanvas ("cas1","cas1",900,800);
  hjs1->Rebin(2);
  hjs1->SetLineColor(kBlack);
  hjs1->Draw("hist");
  hjs1->SetLineWidth(2);
  hjs1->SetFillColor(kBlack);
  hjs1->SetFillStyle(3004);
  hjs1->GetXaxis()->SetTitle("H_{T} (GeV)");
  hjs1->GetYaxis()->SetTitle("Entries");
  hjs1->SetTitle("H_{T} ");
  hjs1->Scale(1./hjs1->Integral());
  
  hjb1->SetLineColor(kRed);
  hjb1->Rebin(2);
  hjb1->Draw("histsames");
  
  hjb1->Scale(1./hjb1->Integral());// Normalize distribution to 1

  hjb3->SetLineColor(kGreen);
  hjb3->Rebin(2);
  hjb3->Draw("histsames");
  hjb3->Scale(1./hjb3->Integral());
  
  // Normalize distribution to 1
  TLegend *le1 = new TLegend(0.55,0.7,.8,.9);
  le1->AddEntry(hjs1,"ZH (20) signal ","l");
  le1->AddEntry(hjb1,"Z #rightarrow #nu#nu","l");
  le1->AddEntry(hjb3,"TTbar #rightarrow (had)","l");
  le1->SetTextSize(0.03); 
  le1->Draw("SAME");
}
