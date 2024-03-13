#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>


void plot4(){
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
  Eff->SetTitle("trigger efficiency: MET from z-> vv after / before trigger requirement");
  Eff->SetMarkerStyle(kFullCircle);
  
  Eff->Draw("AP");
 
		       
}
