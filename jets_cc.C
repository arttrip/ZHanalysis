#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle.C"

void jets_cc(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  TFile *fbkg=new TFile("histos_zvv100.root");
  TFile *fbkg2=new TFile("histos_zvv200.root");
  TFile *fbkg3=new TFile("histos_tt_had.root");
  //CROSS CLEANED jet mult in step 2
  
  TH1F *hjs1 =(TH1F*)fsig->Get("jet_cc_mult2");
  TH1F *hjb1= (TH1F*)fbkg->Get("jet_cc_mult2");
  TH1F *hjb11= (TH1F*)fbkg2->Get("jet_cc_mult2");
  TH1F *hjb111= (TH1F*)fbkg3->Get("jet_cc_mult2");
  hjb1->Add(hjb11);
  setTDRStyle();
  TCanvas *cas1 = new TCanvas ("cas1","cas1",900,800);


  //cas1->SetLogy();
  hjb1->SetLineColor(kRed);
  hjb1->GetXaxis()->SetTitle("no of jets");
  hjb1->GetYaxis()->SetTitle("Entries");
  hjb1->SetTitle("cross cleaned jet mult in step 2 of cut flow ");
  hjb1->Draw("hist");
  
  hjb1->Scale(1./hjb1->Integral());

  hjb111->SetLineColor(kGreen);
  hjb111->Draw("histsames");
  hjb111->Scale(1./hjb111->Integral());

  hjs1->SetLineColor(kBlack);
  hjs1->SetFillColor(kBlack);
  hjs1->SetFillStyle(3004);
  hjs1->Draw("histsames");
  hjs1->SetLineWidth(2);
  hjs1->Scale(1./hjs1->Integral()); // Normalize distribution to 1
  TLegend *le1 = new TLegend(0.55,0.7,.8,.9);
  le1->AddEntry(hjs1,"ZH (20) signal ","l");
  le1->AddEntry(hjb1,"Z #rightarrow #nu#nu","l");
  le1->AddEntry(hjb111,"TTbar #rightarrow (had)","l");
  le1->SetTextSize(0.03); 
  le1->Draw("SAME");

  // return;
  //CROSS CLEANED (FROM FJETS) jet mult in step 3
  
  TH1F *hjs2 =(TH1F*)fsig->Get("jet_cc_mult3");
  TH1F *hjb2= (TH1F*)fbkg->Get("jet_cc_mult3");
  TH1F *hjb22= (TH1F*)fbkg2->Get("jet_cc_mult3");
  hjb2->Add(hjb22);
  TH1F *hjb222= (TH1F*)fbkg3->Get("jet_cc_mult3");
  TCanvas *cas2 = new TCanvas ("cas2","cas2",900,800);
  //cas2->SetLogy();
  hjb2->SetLineColor(kRed);
  hjb2->GetXaxis()->SetTitle("no of jets");
  hjb2->GetYaxis()->SetTitle("Entries");
  hjb2->SetTitle("cross cleaned jet mult in step 3 of cut flow ");
  hjb2->Draw("hist");
  hjb2->Scale(1./hjb2->Integral());// Normalize distribution to 1

  hjb222->SetLineColor(kGreen);
  hjb222->Draw("histsames");
  hjb222->Scale(1./hjb222->Integral());
  
  hjs2->SetLineColor(kBlack);
  hjs2->SetFillColor(kBlack);
  hjs2->SetFillStyle(3004);
  hjs2->Draw("histsames");
  hjs2->SetLineWidth(2);
  hjs2->Scale(1./hjs2->Integral()); // Normalize distribution to 1
  TLegend *le2 = new TLegend(0.55,0.7,.8,.9);
  le2->AddEntry(hjs2,"ZH (20) signal ","l");
  le2->AddEntry(hjb2,"Z #rightarrow #nu#nu","l");
  le2->AddEntry(hjb222," TTbar #rightarrow (had)","l");
  le2->SetTextSize(0.03); 
  le2->Draw("SAME");
}
