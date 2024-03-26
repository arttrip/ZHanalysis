#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle.C"

void plot_fjets(){
  // Open the input files:
  TFile *fsig=new TFile("histos_an20.root");
  TFile *fbkg1=new TFile("histos_zvv100.root");
  TFile *fbkg2=new TFile("histos_zvv200.root");
  

  //fjet mult in step 1
  
  TH1F *hjs1 =(TH1F*)fsig->Get("fjet_mult");
  TH1F *hjb1= (TH1F*)fbkg1->Get("fjet_mult");
  TH1F *hjb11= (TH1F*)fbkg2->Get("fjet_mult");
  hjb1->Add(hjb11);
  setTDRStyle();
  TCanvas *cas1 = new TCanvas ("cas1","cas1",900,800);


  //cas1->SetLogy();
  hjb1->SetLineColor(kRed);
  hjb1->GetXaxis()->SetTitle("no of fjets");
  hjb1->GetYaxis()->SetTitle("Entries");
  hjb1->SetTitle("fjet mult in step 1 of cut flow ");
  hjb1->Draw("hist");
  hjb1->SetFillColor(kRed);
  hjb1->SetFillStyle(3004);
  hjb1->Scale(1./hjb1->Integral());// Normalize distribution to 1

  hjs1->SetLineColor(kBlack);
  hjs1->Draw("histsames");
  hjs1->SetLineWidth(2);
  hjs1->Scale(1./hjs1->Integral()); // Normalize distribution to 1
  TLegend *le1 = new TLegend(0.55,0.7,.8,.9);
  le1->AddEntry(hjs1,"ZH (20) signal ","l");
  le1->AddEntry(hjb1,"Z #rightarrow #nu#nu","l");
  le1->SetTextSize(0.03); 
  le1->Draw("SAME");

  // return;
  //fjet mult in step 2
  
  TH1F *hjs2 =(TH1F*)fsig->Get("fjet_mult_after1");
  TH1F *hjb2= (TH1F*)fbkg1->Get("fjet_mult_after1");
  TH1F *hjb22= (TH1F*)fbkg2->Get("fjet_mult_after1");
  hjb2->Add(hjb1);
  TCanvas *cas2 = new TCanvas ("cas2","cas2",900,800);
  //cas2->SetLogy();
  hjb2->SetLineColor(kRed);
  hjb2->GetXaxis()->SetTitle("no of fjets");
  hjb2->GetYaxis()->SetTitle("Entries");
  hjb2->SetTitle("fjet mult in step 2 of cut flow ");
  hjb2->Draw("hist");
  hjb2->SetFillColor(kRed);
  hjb2->SetFillStyle(3004);
  hjb2->Scale(1./hjb2->Integral());// Normalize distribution to 1

  hjs2->SetLineColor(kBlack);
  hjs2->Draw("histsames");
  hjs2->SetLineWidth(2);
  hjs2->Scale(1./hjs2->Integral()); // Normalize distribution to 1
  TLegend *le2 = new TLegend(0.55,0.7,.8,.9);
  le2->AddEntry(hjs2,"ZH (20) signal ","l");
  le2->AddEntry(hjb2," z #rightarrow #nu#nu","l");
  le2->SetTextSize(0.03); 
  le2->Draw("SAME");


  //fjet1 subjet multiplicity
  TH1F *hsu1 =(TH1F*)fsig->Get("sub_mult1");
  TH1F *hsum1= (TH1F*)fsig->Get("sub_mult1_mat");
  TCanvas *cu1 = new TCanvas ("cu1","cu1",900,800);
  cu1->SetLogy();
  hsum1->SetLineColor(kBlack);
  hsum1->GetXaxis()->SetTitle("no of sub-jets");
  hsum1->GetYaxis()->SetTitle("Entries");
  hsum1->SetTitle("sub-jet mult in fjet1");
  hsum1->Draw("hist");
  hsum1->Scale(1./hsum1->Integral());// Normalize distribution to 1

  hsu1->SetLineColor(kRed);
  hsu1->Draw("histsames");
  hsu1->Scale(1./hsu1->Integral()); // Normalize distribution to 1
  TLegend *lg1 = new TLegend(0.55,0.7,.8,.9);
  lg1->AddEntry(hsu1,"fjet1 ","l");
  lg1->AddEntry(hsum1,"fjet1 matched","l");
  lg1->SetTextSize(0.03); 
  lg1->Draw("SAME");
  //fjet2 subjet multiplicity
  TH1F *hsu2 =(TH1F*)fsig->Get("sub_mult2");
  TH1F *hsum2= (TH1F*)fsig->Get("sub_mult2_mat");
  TCanvas *cu2 = new TCanvas ("cu2","cu2",900,800);
  cu2->SetLogy();
  hsum2->SetLineColor(kBlack);
  hsum2->GetXaxis()->SetTitle("no of sub-jets");
  hsum2->GetYaxis()->SetTitle("Entries");
  hsum2->SetTitle("sub-jet mult in fjet2");
  hsum2->Draw("hist");
  hsum2->Scale(1./hsum2->Integral());// Normalize distribution to 1

  hsu2->SetLineColor(kRed);
  hsu2->Draw("histsames");
  hsu2->Scale(1./hsu2->Integral()); // Normalize distribution to 1
  TLegend *lg2 = new TLegend(0.55,0.7,.8,.9);
  lg2->AddEntry(hsu2,"fjet2 ","l");
  lg2->AddEntry(hsum2,"fjet2 matched","l");
  lg2->SetTextSize(0.03); 
  lg2->Draw("SAME");


  //fjet1 pt
  TH1F *hs1 =(TH1F*)fsig->Get("h_fjet1_pt");
  hs1->Rebin(4);
  TH1F *hsm1= (TH1F*)fsig->Get("h_fj1_pt_mat");
  hsm1->Rebin(4);
  // Draw the histos
  TCanvas *ca1 = new TCanvas ("ca1","ca1",900,800);
  ca1->SetLogy();
  hsm1->SetLineColor(kBlack);
  hsm1->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsm1->GetYaxis()->SetTitle("Entries");
  hsm1->SetTitle("fjet1 pt");
  hsm1->Draw("hist");
  hsm1->Scale(1./hsm1->Integral());// Normalize distribution to 1


  hs1->SetLineColor(kRed);
  hs1->Draw("histsames");
  hs1->Scale(1./hs1->Integral()); // Normalize distribution to 1
  
  TLegend *l = new TLegend(0.55,0.7,.8,.9);
  l->AddEntry(hs1,"fjet1","l");
  l->AddEntry(hsm1,"fjet1 matched","l");
  l->SetTextSize(0.03); 
  l->Draw("SAME");

  //fjet2 pt
  TH1F *hs2 =(TH1F*)fsig->Get("h_fjet2_pt");
  hs2->Rebin(4);
  TH1F *hsm2= (TH1F*)fsig->Get("h_fj2_pt_mat");
  hsm2->Rebin(4);
  // Draw the histos
  TCanvas *ca2 = new TCanvas ("ca2","ca2",900,800);
  ca2->SetLogy();
  hsm2->SetLineColor(kBlack);
  hsm2->GetXaxis()->SetTitle("p_{T} (GeV)");
  hsm2->GetYaxis()->SetTitle("Entries");
  hsm2->SetTitle("fjet2 pt");
  hsm2->Draw("hist");
  //hsm2->Scale(1./hsm2->Integral());// Normalize distribution to 1

  hs2->SetLineColor(kRed);
  hs2->Draw("histsames");
  //hs2->Scale(1./hs2->Integral()); // Normalize distribution to 1

  TLegend *l2 = new TLegend(0.55,0.7,.8,.9);
  l2->AddEntry(hs2,"fjet2","l");
  l2->AddEntry(hsm2,"fjet2 matched","l");
  l2->SetTextSize(0.03); 
  l2->Draw();
  
  //ZH SIGNAL : fjet1 matched pt vs bb pair pt
  
  TH2F *hf1= (TH2F*)fsig->Get("h_fj1_pt_bb_pt");
  TCanvas *cn1 = new TCanvas ("cn1","cn1",900,800);
  hf1->GetXaxis()->SetTitle(" fjet1 p_{T} (GeV) ");
  hf1->GetYaxis()->SetTitle("bb pair p_{T} (GeV)");
  hf1->SetTitle("fjet1 matched pt vs pt(bb) ");
  hf1->SetStats(0);
  hf1->Draw("COLZ");
  TLine *line1 = new TLine(0,0,500,500);
  line1->SetLineColor(kRed);
  line1->Draw();

  //ZH signal: fjet2 matched pt vs bb pair pt
  
  TH2F *hf2= (TH2F*)fsig->Get("h_fj2_pt_bb_pt");
  TCanvas *cn2 = new TCanvas ("cn2","cn2",900,800);
  hf2->GetXaxis()->SetTitle(" fjet2 p_{T} (GeV) ");
  hf2->GetYaxis()->SetTitle("bb pair p_{T} (GeV)");
  hf2->SetStats(0);
  hf2->SetTitle("fjet2 matched pt vs pt(bb) ");
  hf2->Draw("COLZ");
  TLine *line2 = new TLine(0,0,500,500);
  line2->SetLineColor(kRed);
  line2->Draw();

  
  //ZH SIGNAL : fjet1 : pt vs sd mass
  
  TH2F *hsd1= (TH2F*)fsig->Get("fj_pt_sd_mass1");
  TCanvas *can1 = new TCanvas ("can1","can1",900,800);
  hsd1->GetXaxis()->SetTitle(" fjet1 p_{T} (GeV) ");
  hsd1->GetXaxis()->SetRangeUser(0,150);
  hsd1->GetYaxis()->SetRangeUser(0,150);
  hsd1->SetStats(0);
  hsd1->GetYaxis()->SetTitle("sd mass (GeV)");
  hsd1->SetTitle("fjet1: pt vs sd mass ");
  hsd1->Draw("COLZ");
  

  //ZH signal: fjet2 : pt vs sd mass
  TH2F *hsd2= (TH2F*)fsig->Get("fj_pt_sd_mass2");
  TCanvas *can2 = new TCanvas ("can2","can2",900,800);
  hsd2->GetXaxis()->SetTitle(" fjet2 p_{T} (GeV) ");
  hsd2->GetXaxis()->SetRangeUser(0,150);
  hsd2->GetYaxis()->SetRangeUser(0,150);
  hsd2->SetStats(0);
  hsd2->GetYaxis()->SetTitle("sd mass (GeV)");
  hsd2->SetTitle("fjet2: pt vs sd mass ");
  hsd2->Draw("COLZ");

  //ZH signal: fjet1 sd mass
  TH1F *m1 =(TH1F*)fsig->Get("fj_sd_mass1");
  TH1F *mm1= (TH1F*)fsig->Get("fj_sd_mass1_matched");

  TCanvas *ma1 = new TCanvas ("ma1","ma1",900,800);
  //ma1->SetLogy();
  m1->SetLineColor(kRed);
  m1->GetXaxis()->SetTitle("sd mass (GeV)");
  m1->GetXaxis()->SetRangeUser(0,150);
  m1->GetYaxis()->SetTitle("Entries");
  m1->SetTitle("fjet1 sd mass ");
  m1->Draw("hist");
  //m1->Scale(1./m1->Integral());// Normalize distribution to 1
  mm1->SetLineColor(kBlack);
  mm1->Draw("histsames");
  //mm1->Scale(1./mm1->Integral()); // Normalize distribution to 1
  TLegend *lm1 = new TLegend(0.55,0.7,.8,.9);
  lm1->AddEntry(m1,"fjet1","l");
  lm1->AddEntry(mm1,"fjet1 matched","l");
  lm1->SetTextSize(0.03); 
  lm1->Draw("SAME");

  //ZH signal: fjet2 sd mass
  TH1F *m2 =(TH1F*)fsig->Get("fj_sd_mass2");
  TH1F *mm2= (TH1F*)fsig->Get("fj_sd_mass2_matched");

  TCanvas *ma2 = new TCanvas ("ma2","ma2",900,800);
  //ma2->SetLogy();
  m2->SetLineColor(kRed);
  m2->GetXaxis()->SetRangeUser(0,150);
  m2->GetXaxis()->SetTitle("sd mass (GeV)");
  m2->GetYaxis()->SetTitle("Entries");
  m2->SetTitle("fjet2 sd mass ");
  m2->Draw("hist");
  //m2->Scale(1./m2->Integral());// Normalize distribution to 1
  mm2->SetLineColor(kBlack);
  mm2->Draw("histsames");
  //mm2->Scale(1./mm2->Integral()); // Normalize distribution to 1
  TLegend *lm2 = new TLegend(0.55,0.7,.8,.9);
  lm2->AddEntry(m2,"fjet2","l");
  lm2->AddEntry(mm2,"fjet2 matched","l");
  lm2->SetTextSize(0.03); 
  lm2->Draw("SAME");
  
 

  //fjet1 XbbXccXqq
  TH1F *hsig = (TH1F*)fsig->Get("h_fjet1_btagXbbXccXqq");
  TH1F *hbkg1 =(TH1F*)fbkg1->Get("h_fjet1_btagXbbXccXqq");
  TH1F *hbkg11 =(TH1F*)fbkg2->Get("h_fjet1_btagXbbXccXqq");
  hbkg1->Add(hbkg11);

  // Draw the histos
  TCanvas *c1 = new TCanvas ("c1","c1",900,800);
  c1->SetLogy();
  hsig->SetLineColor(kBlack);
  hsig->SetLineWidth(2);;
  hsig->GetXaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD");
  hsig->GetYaxis()->SetTitle("Entries");
  hsig->SetTitle("ParticleNetMD_XbbXccXqqvsQCD on fjet1  ");
  hsig->Draw("hist");
  hsig->Scale(1./hsig->Integral());// Normalize distribution to 1
  
  hbkg1->SetLineColor(kRed);
  hbkg1->SetFillColor(kRed);
  hbkg1->SetFillStyle(3001);
  hbkg1->Draw("histsames");
  hbkg1->Scale(1./hbkg1->Integral());
  
  TLegend *leg = new TLegend(0.55,0.7,.8,.9);
  leg->AddEntry(hsig,"ZH signal","l");
  leg->AddEntry(hbkg1,"Z #rightarrow #nu#nu","l");
  leg->SetTextSize(0.03); 
  leg->Draw("SAME");

 

  //fjet2 XbbXccXqq
  
  TH1F *hsig2 = (TH1F*)fsig->Get("h_fjet2_btagXbbXccXqq");
  //hsig->Rebin(2);
  TH1F *hbkg2 =(TH1F*)fbkg1->Get("h_fjet2_btagXbbXccXqq");
  //hbkg->Rebin(2);
  TH1F *hbkg22 =(TH1F*)fbkg2->Get("h_fjet1_btagXbbXccXqq");
  hbkg2->Add(hbkg22);
  // Draw the histos
  TCanvas *c2 = new TCanvas ("c2","c2",900,800);
  c2->SetLogy();

  hsig2->SetLineColor(kBlack);
  hsig2->SetLineWidth(2);
  hsig2->GetXaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD");
  hsig2->GetYaxis()->SetTitle("Entries");
  hsig2->SetTitle("ParticleNetMD_XbbXccXqqvsQCD on fjet2  ");
  hsig2->Draw("hist");
  hsig2->Scale(1./hsig2->Integral());// Normalize distribution to 1
  
  hbkg2->SetLineColor(kRed);
  hbkg2->SetFillColor(kRed);
  hbkg2->SetFillStyle(3001);
  hbkg2->Draw("histsames");
  hbkg2->Scale(1./hbkg2->Integral());
  
  TLegend *leg2 = new TLegend(0.55,0.7,.8,.9);
  leg2->AddEntry(hsig2,"ZH(20) signal","l");
  leg2->AddEntry(hbkg2,"z #rightarrow #nu#nu","l");
  leg2->SetTextSize(0.03); 
  leg2->Draw("SAME");

   
  //fjet1 Xbb
  TH1F *hsig3 = (TH1F*)fsig->Get("h_fjet1_btagXbb");
  TH1F *hbkg3 =(TH1F*)fbkg1->Get("h_fjet1_btagXbb");
  TH1F *hbkg33 =(TH1F*)fbkg2->Get("h_fjet1_btagXbbXccXqq");
  hbkg3->Add(hbkg33);
  
 
  // Draw the histos
  TCanvas *c3 = new TCanvas ("c3","c3",900,800);
  c3->SetLogy();

  hsig3->SetLineColor(kBlack);
  hsig3->SetLineWidth(2);
  hsig3->GetXaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD");
  hsig3->GetYaxis()->SetTitle("Entries");
  hsig3->GetYaxis()->SetRange(0,1);
  hsig3->SetTitle("ParticleNetMD_XbbvsQCD on fjet1  ");
  hsig3->Draw("hist");
  hsig3->Scale(1./hsig3->Integral());// Normalize distribution to 1
  
  hbkg3->SetLineColor(kRed);
  hbkg3->SetFillColor(kRed);
  hbkg3->SetFillStyle(3001);
  hbkg3->Draw("histsames");
  hbkg3->Scale(1./hbkg3->Integral());
  
  TLegend *leg3 = new TLegend(0.55,0.7,.8,.9);
  leg3->AddEntry(hsig3,"ZH(20) signal","l");
  leg3->AddEntry(hbkg3,"z #rightarrow #nu#nu","l");
  leg3->SetTextSize(0.03); 
  leg3->Draw("SAME");

 //fjet2 xbb
  TH1F *hsig4 = (TH1F*)fsig->Get("h_fjet2_btagXbb");
  TH1F *hbkg4 =(TH1F*)fbkg1->Get("h_fjet2_btagXbb");
  TH1F *hbkg44 =(TH1F*)fbkg2->Get("h_fjet1_btagXbbXccXqq");
  hbkg4->Add(hbkg44);
 
  // Draw the histos
  TCanvas *c4 = new TCanvas ("c4","c4",900,800);
  c4->SetLogy();

  hsig4->SetLineColor(kBlack);
  hsig4->SetLineWidth(2);
  hsig4->GetXaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD");
  hsig4->GetYaxis()->SetTitle("Entries");
  hsig4->GetYaxis()->SetRange(0,1);
  hsig4->SetTitle("ParticleNetMD_XbbvsQCD on fjet2  ");
  hsig4->Draw("hist");
  hsig4->Scale(1./hsig4->Integral());// Normalize distribution to 1
  
  hbkg4->SetLineColor(kRed);
  hbkg4->SetFillColor(kRed);
  hbkg4->SetFillStyle(3001);
  hbkg4->Draw("histsames");
  hbkg4->Scale(1./hbkg4->Integral());
  
  TLegend *leg4 = new TLegend(0.55,0.7,.8,.9);
  leg4->AddEntry(hsig4,"ZH(20) signal","l");
  leg4->AddEntry(hbkg4,"z #rightarrow #nu#nu","l");
  leg4->SetTextSize(0.03); 
  leg4->Draw("SAME");

 
  // ZH signal fjet1-fjet2 xbb
  TH2F *h= (TH2F*)fsig->Get("h_fj1_vs_fj2_btagXbb");
  TCanvas *c5 = new TCanvas ("c5","c5",900,800);
  h->SetStats(0);
  h->GetXaxis()->SetTitle(" ParticleNetMD_Xbb_vsQCD on fjet1  ");
  h->GetYaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD on fjet2");
  h->SetTitle("ParticleNetMD_Xbb_vsQCD : on fjet1 vs on fjet2 (ZH 20 signal) ");
  h->Draw("COLZ");


  //ZH signal fjet1-fjet2 xbbccqq
  TH2F *h2= (TH2F*)fsig->Get("h_fj1_fj2_btagXbbXccXqq");
  TCanvas *c6 = new TCanvas ("c6","c6",900,800);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle(" ParticleNetMD_XbbXccXqq_vsQCD on fjet1  ");
  h2->GetYaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD on fjet2");
  h2->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD : on fjet1 vs on fjet2 (ZH m20 signal) ");
  h2->Draw("COLZ");

  // z #rightarrow #nu#nu fjet1-fjet2 xbb
  TH2F *hb1= (TH2F*)fbkg1->Get("h_fj1_vs_fj2_btagXbb");
  TH2F *hb11= (TH2F*)fbkg2->Get("h_fj1_vs_fj2_btagXbb");
  hb1->Add(hb11);
  TCanvas *c7 = new TCanvas ("c7","c7",900,800);
  hb1->SetStats(0);
  hb1->GetXaxis()->SetTitle(" ParticleNetMD_Xbb_vsQCD on fjet1  ");
  hb1->GetYaxis()->SetTitle("ParticleNetMD_Xbb_vsQCD on fjet2");
  hb1->SetTitle("ParticleNetMD_Xbb_vsQCD : on fjet1 vs on fjet2 (z #rightarrow #nu#nu) ");
  hb1->Draw("COLZ");


  //z #rightarrow #nu#nu  fjet1-fjet2 xbbccqq
  TH2F *hb2= (TH2F*)fbkg1->Get("h_fj1_fj2_btagXbbXccXqq");
  TH2F *hb22= (TH2F*)fbkg2->Get("h_fj1_fj2_btagXbbXccXqq");
  hb2->Add(hb22);
  TCanvas *c8 = new TCanvas ("c8","c8",900,800);
  hb2->SetStats(0);
  hb2->GetXaxis()->SetTitle(" ParticleNetMD_XbbXccXqq_vsQCD on fjet1  ");
  hb2->GetYaxis()->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD on fjet2");
  hb2->SetTitle("ParticleNetMD_XbbXccXqq_vsQCD : on fjet1 vs on fjet2 (z #rightarrow #nu#nu) ");
  hb2->Draw("COLZ");
}
