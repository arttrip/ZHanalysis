#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle.C"

void ht(){
  // Open the input files:
  TFile *fsig=new TFile("histos_a20.root");
  TFile *fbkg=new TFile("histos_zvv.root");

  //CROS CLEANED jet mult in step 2
  
  TH1F *hjs1 =(TH1F*)fsig->Get("Ht");
  TH1F *hjb1= (TH1F*)fbkg->Get("Ht");
  setTDRStyle();
  TCanvas *cas1 = new TCanvas ("cas1","cas1",900,800);
  hjb1->SetLineColor(kBlack);
  hjb1->GetXaxis()->SetTitle("H_{T} (GeV)");
  hjb1->GetYaxis()->SetTitle("Entries");
  hjb1->SetTitle("H_{T} ");
  hjb1->Draw("hist");
  hjb1->SetFillColor(1);
  hjb1->SetFillStyle(3004);
  hjb1->Scale(1./hjb1->Integral());// Normalize distribution to 1

  hjs1->SetLineColor(kRed);
  hjs1->Draw("histsames");
  hjs1->SetLineWidth(2);
  hjs1->Scale(1./hjs1->Integral()); // Normalize distribution to 1
  TLegend *le1 = new TLegend(0.55,0.7,.8,.9);
  le1->AddEntry(hjs1,"ZH (20) signal ","l");
  le1->AddEntry(hjb1,"Z #rightarrow #nu#nu","l");
  le1->SetTextSize(0.03); 
  le1->Draw("SAME");
}