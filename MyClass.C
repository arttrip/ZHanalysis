#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
//#include "tabulate.hpp"
#include "TLorentzVector.h"


#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>

#include <iterator>
#include <cstdint>
#include <array>






double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2);
bool jet_matched(TLorentzVector jet,std::vector<TLorentzVector> vecb);
bool fjet_matched(TLorentzVector jet,std::vector<TLorentzVector> vecb) ; 
double dR_j_fj_min(TLorentzVector jet,std::vector<TLorentzVector> vecfjet);
float fjet_bkg_matched_qg_pt(TLorentzVector fjet,std::vector<TLorentzVector> vecb) ;
bool fjet_bkg_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb) ;
double Ht(std::vector<TLorentzVector> vecjet,std::vector<TLorentzVector> vecfjet);
void MyClass::Loop()
{
  if (fChain == 0)
    return;
  TString file_name = "a20.root";
  bool isSignal(false);
  if(file_name.Contains("a20") || file_name.Contains("a60")) isSignal=true;

  bool verbose(true); 
 
  Long64_t nentries = fChain->GetEntriesFast();

  //----Histograms------------------------


  // control plots
  TH1F *h_en_mult = new TH1F("en_mult", "electrons", 5, 0, 5);
  TH1F *h_mn_mult = new TH1F("mn_mult", "muons", 5, 0, 5);
  TH1F *h_lepton_mult = new TH1F("lepton_mult", "leptons", 5, 0, 5);
  TH1F *h_jet_mult = new TH1F("jet_mult", "jets", 10, 0, 10);
  TH1F *h_bjet_mult = new TH1F("bjet_mult", "bJets", 10, 0, 10);
  TH1F *h_jet_cc_mult1 = new TH1F("jet_cc_mult1", "jets", 10, 0, 10);
  TH1F *h_jet_cc_mult2 = new TH1F("jet_cc_mult2", "jets", 10, 0, 10);
  TH1F *h_jet_cc_mult3 = new TH1F("jet_cc_mult3", "jets", 10, 0, 10);
  // Multiplicities after cuts
  TH1F *h_en_mult_after = new TH1F("en_mult_after", "electrons", 5, 0, 5);
  TH1F *h_mn_mult_after = new TH1F("mn_mult_after", "muons", 5, 0, 5);
  TH1F *h_lepton_mult_after = new TH1F("lepton_mult", "leptons", 5, 0, 5);
  TH1F *h_jet_mult_after = new TH1F("jet_mult_after", "Jets", 10, 0, 10);
  TH1F *h_bjet_mult_after = new TH1F("bjet_mult_after", "bJets", 10, 0, 10);
 //electrons
  TH1F *h_en1_pt = new TH1F("h_en1_pt", "en1_pt", 100, 0, 300);
  TH1F *h_en1_eta = new TH1F("h_en1_eta", "en1_eta", 100, -5, 5);
  TH1F *h_en1_phi = new TH1F("h_en1_phi", "en1_phi", 100, -TMath::Pi(), TMath::Pi());

  TH1F *h_en2_pt = new TH1F("h_en2_pt", "en2_pt", 250, 0, 300);
  TH1F *h_en2_eta = new TH1F("h_en2_eta", "en2_eta", 100, -5, 5);
  TH1F *h_en2_phi = new TH1F("h_en2_phi", "en2_phi", 100, -TMath::Pi(), TMath::Pi());
  //muons
  TH1F *h_mn1_pt = new TH1F("h_mn1_pt", "mn1_pt", 100, 0, 300);
  TH1F *h_mn1_eta = new TH1F("h_mn1_eta", "mn1_eta", 100, -5, 5);
  TH1F *h_mn1_phi = new TH1F("h_mn1_phi", "muon1 phi", 100, -TMath::Pi(), TMath::Pi());

  TH1F *h_mn2_pt = new TH1F("h_mn2_pt", "mn2_pt", 250, 0, 300);
  TH1F *h_mn2_eta = new TH1F("h_mn2_eta", "mn2_eta", 100, -5, 5);
  TH1F *h_mn2_phi = new TH1F("h_mn2_phi", "muon2 phi", 100, -TMath::Pi(), TMath::Pi());
  //dilepton_inv_m
  TH1F *h_2en_inv_m = new TH1F("h_2en_inv_m", "2en_inv_m", 200,0,200);
  TH1F *h_2mn_inv_m1 = new TH1F("h_2mn_inv_m1", "2mn_inv_m1", 200,0, 200);
  TH1F *h_2mn_inv_m2 = new TH1F("h_2mn_inv_m2", "2mn_inv_m2", 200,0, 200);
  TH1F *h_2mn_inv_m3 = new TH1F("h_2mn_inv_m3", "2mn_inv_m3", 200,0, 200);
  //dilepton_Pt
  TH1F *h_2en_pt= new TH1F("h_2en_pt", "2en_pt", 600,0,600);
  TH1F *h_2mn_pt = new TH1F("h_2mn_pt", "2mn_pt", 600,0,600);
  // jets

  TH1F *h_jet1_pt = new TH1F("h_jet1_pt", "jet1_pt", 100, 0, 300);
  TH1F *h_jet1_eta = new TH1F("h_jet1_eta", "jet1_eta", 100, -5, 5);
  TH1F *h_jet1_phi = new TH1F("h_jet1_phi", "jet1_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_jet2_pt = new TH1F("h_jet2_pt", "jet2_pt", 100, 0, 300);
  TH1F *h_jet2_eta = new TH1F("h_jet2_eta", "jet2_eta", 100, -5, 5);
  TH1F *h_jet2_phi = new TH1F("h_jet2_phi", "jet2_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_jet3_pt = new TH1F("h_jet3_pt", "jet3_pt", 250, 0, 300);
  TH1F *h_jet3_eta = new TH1F("h_jet3_eta", "jet3_eta", 100, -5, 5);
  TH1F *h_jet3_phi = new TH1F("h_jet3_phi", "jet3_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_jet4_pt = new TH1F("h_jet4_pt", "jet4_pt", 200, 0, 300);
  TH1F *h_jet4_eta = new TH1F("h_jet4_eta", "jet4_eta", 100, -5, 5);
  TH1F *h_jet4_phi = new TH1F("h_jet4_phi", "jet4_phi", 100, -TMath::Pi(), TMath::Pi());
   // bjets
  
  TH1F *h_bjet1_pt = new TH1F("h_bjet1_pt", "bjet1_pt", 100, 0, 300);
  TH1F *h_bjet1_eta = new TH1F("h_bjet1_eta", "bjet1_eta", 100, -5, 5);
  TH1F *h_bjet1_phi = new TH1F("h_bjet1_phi", "bjet1_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_bjet2_pt = new TH1F("h_bjet2_pt", "bjet2_pt", 100, 0, 300);
  TH1F *h_bjet2_eta = new TH1F("h_bjet2_eta", "bjet2_eta", 100, -5, 5);
  TH1F *h_bjet2_phi = new TH1F("h_bjet2_phi", "bjet2_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_bjet3_pt = new TH1F("h_bjet3_pt", "bjet3_pt", 250, 0, 300);
  TH1F *h_bjet3_eta = new TH1F("h_bjet3_eta", "bjet3_eta", 100, -5, 5);
  TH1F *h_bjet3_phi = new TH1F("h_bjet3_phi", "bjet3_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_bjet4_pt = new TH1F("h_bjet4_pt", "bjet4_pt", 200, 0, 300);
  TH1F *h_bjet4_eta = new TH1F("h_bjet4_eta", "bjet4_eta", 100, -5, 5);
  TH1F *h_bjet4_phi = new TH1F("h_bjet4_phi", "bjet4_phi", 100, -TMath::Pi(), TMath::Pi());
  //3-4jet_inv_m
  TH1F *h_higgs_inv_m=new TH1F("h_higgs_inv_m","higgs_inv_m",500,0,1000);
  TH1F *h_higgs_pt=new TH1F("h_higgs_pt","higgs_pt",250,0,500);
  TH1F *h_higgs_eta=new TH1F("h_higgs_eta","higgs_eta",100,-5,5);
  TH1F *h_higgs_phi=new TH1F("h_higgs_phi","higgs_phi",100, -TMath::Pi(), TMath::Pi());
  
  // fjets
  TH1F *h_fjet_mult = new TH1F("fjet_mult", "fjets", 10, 0, 10);
  TH1F *h_fjet_mult_after1 = new TH1F("fjet_mult_after1", "fJets", 10, 0, 10);
  TH1F *h_fjet_mult_after2 = new TH1F("fjet_mult_after2", "fJets", 10, 0, 10);
  TH1F *h_fjet1_pt = new TH1F("h_fjet1_pt", "fjet1_pt", 100, 0, 500);
  TH1F *h_fjet1_eta = new TH1F("h_fjet1_eta", "fjet1_eta", 100, -5, 5);
  TH1F *h_fjet1_phi = new TH1F("h_fjet1_phi", "fjet1_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_fjet1_btagXbb = new TH1F("h_fjet1_btagXbb", "fjet1_btagXbb", 50, 0, 1);
  TH1F *h_fjet1_btagXbbXccXqq = new TH1F("h_fjet1_btagXbbXccXqq", "fjet1_btagXbbXccXqq", 50, 0, 1);
  TH1F *h_fjet2_pt = new TH1F("h_fjet2_pt", "fjet2_pt", 100, 0, 500);
  TH1F *h_fjet2_eta = new TH1F("h_fjet2_eta", "fjet2_eta", 100, -5, 5);
  TH1F *h_fjet2_phi = new TH1F("h_fjet2_phi", "fjet2_phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_fjet2_btagXbb = new TH1F("h_fjet2_btagXbb", "fjet2_btagXbb", 50, 0, 1);
  TH1F *h_fjet2_btagXbbXccXqq = new TH1F("h_fjet2_btagXbbXccXqq", "fjet2_btagXbbXccXqq", 50, 0, 1);

  TH2F *h_fj1_fj2_btagXbb = new TH2F("h_fj1_vs_fj2_btagXbb", "fjet_btagXbb",50,0,1, 50, 0, 1);
  TH2F *h_fj1_fj2_btagXbbXccXqq = new TH2F("h_fj1_fj2_btagXbbXccXqq", "fjet_btagXbbXccXqq",50,0,1,50, 0, 1);



  // -- -- -- Delta R -- -- --
  TH1F *h_dR_jet_en_before = new TH1F("h_dR_jet_en_before", "dR_fjet_en_before", 100, 0, 8);
  TH1F *h_dR_jet_mn_before = new TH1F("h_dR_jet_mn_before", "dR_fjet_mn_before", 100, 0, 8);
  TH1F *h_dR_jet_en_after = new TH1F("h_dR_jet_en_after", "dR_fjet_en_after", 100, 0, 8);
  TH1F *h_dR_jet_mn_after = new TH1F("h_dR_jet_mn_after", "dR_fjet_mn_after", 100, 0, 8);
  // f dR
  TH1F *h_dR_fjet_en_before = new TH1F("h_dR_fjet_en_before", "dR_fjet_en_before", 100, 0, 8);
  TH1F *h_dR_fjet_mn_before = new TH1F("h_dR_fjet_mn_before", "dR_fjet_mn_before", 100, 0, 8);
  TH1F *h_dR_fjet_en_after = new TH1F("h_dR_fjet_en_after", "dR_fjet_en_after", 100, 0, 8);
  TH1F *h_dR_fjet_mn_after = new TH1F("h_dR_fjet_mn_after", "dR_fjet_mn_after", 100, 0, 8);
  TH1F *h_dR_j_fj = new TH1F("h_dR_j-fj", "dR_jet-fjet", 100, 0, 8);
  TH1F *h_dR_j_fj_min = new TH1F("h_dR_j-fj_min", "dR_jet-fjet_min", 100, 0,8);
  //fjet matched pt
  TH1F *h_fj1_pt_mat = new TH1F("h_fj1_pt_mat", "fj1_pt_mat", 100, 0,500);
  TH1F *h_fj2_pt_mat = new TH1F("h_fj2_pt_mat", "fj2_pt_mat", 100, 0,500);
  //pt if dr_imin<0.8
  TH1F *pt_j = new TH1F("pt_j", "pt_j", 100, 0,500);
  TH1F *pt_fj = new TH1F("pt_fj", "pt_fj", 100, 0,500);
  TH1F *pt_fj_matched = new TH1F("pt_fj_matched", "pt_fj_matched", 100, 0,500);
  TH1F *pt_j_matched = new TH1F("pt_j_matched", "pt_j_matched", 100, 0,500);
    
  //2D
  TH2F *h_pt_dR_min = new TH2F("h_dR_pt", "dR_pt",100,0,500, 100, 0,8);
  TH2F *h_fj1_pt_bb_pt = new TH2F("h_fj1_pt_bb_pt", "fj1_pt_bb_pt",100,0,500, 100, 0,500);
  TH2F *h_fj2_pt_bb_pt = new TH2F("h_fj2_pt_bb_pt", "fj2_pt_bb_pt",100,0,500, 100, 0,500);
  TH2F *fj_sd_mass1_2=new TH2F("fj_sd_mass1_2","fj_sd_mass1_2",100,0,300,100,0,300);
  TH2F *fj_sd_mass1_2_mat=new TH2F("fj_sd_mass1_2_mat","fj_sd_mass1_2",100,0,300,100,0,300);

  //fjet and fjet matched soft drop mass 

  TH1F *fj_sd_mass1=new TH1F("fj_sd_mass1","fj_sd_mass1",100,0,300);
  TH1F *fj_sd_mass2=new TH1F("fj_sd_mass2","fj_sd_mass2",100,0,300);
  TH1F *fj_sd_mass1_mat=new TH1F("fj_sd_mass1_matched","fj_sd_mass1_matched",100,0,200);
  TH1F *fj_sd_mass2_mat=new TH1F("fj_sd_mass2_matched","fj_sd_mass2_matched",100,0,200);
  TH2F *fj_pt_sd_mass1 = new TH2F("fj_pt_sd_mass1", "fj_pt_sd_mass1",100,0,300, 100, 0,300);
  TH2F *fj_pt_sd_mass2 = new TH2F("fj_pt_sd_mass2", "fj_pt_sd_mass2",100,0,300, 100, 0,300);
  //subjet count
  TH1F *subcount1 = new TH1F("sub_mult1", "sub_mult1", 5, 0, 5);
  TH1F *subcount1_mat = new TH1F("sub_mult1_mat", "sub_mult1", 5, 0, 5);
  TH1F *subcount2 = new TH1F("sub_mult2", "sub_mult1", 5, 0, 5);
  TH1F *subcount2_mat = new TH1F("sub_mult2_mat", "sub_mult1", 5, 0, 5);
  
 
 //generator level histos
 //HIGGS
  TH1F *h_ptH=new TH1F("H_HIGGS_PT","HIGGS_PT",500,0,500); 
  TH1F *h_etaH=new TH1F("H_HIGGS_ETA","HIGGS_ETA",100,-5,5); 
  TH1F *h_phiH=new TH1F("H_HIGGS_PHI","HIGGS_PHI",100,-TMath::Pi(), TMath::Pi()); 
  TH1F *h_mH=new TH1F("H_HIGGS_M","HIGGS_M",50,100,150); 
  //dileptons
  TH1F *h_pt_dlep=new TH1F("H_dilept_PT","dilept_PT",500,0,500); 
  TH1F *h_eta_dlep=new TH1F("H_dilept_ETA","dilept_ETA",100,-5,5); 
  TH1F *h_phi_dlep=new TH1F("H_dilept_PHI","dilept_PHI",100,-TMath::Pi(), TMath::Pi()); 
  TH1F *h_m_dlep=new TH1F("H_dilept_M","dilept_M",100,50,150); 
 
  //Z boson
  TH1F *h_pt_z=new TH1F("h_z_PT","z_PT",500,0,500); 
  TH1F *h_eta_z=new TH1F("H_z_eta","z_eta",100,-5,5); 
  TH1F *h_phi_z=new TH1F("H_z_phi","z_phi",100,-TMath::Pi(), TMath::Pi()); 
  TH1F *h_m_z=new TH1F("H_4b_m","4b_m",100,50,150); 
  //A
  TH1F *h_pt_a=new TH1F("h_A_PT","A_PT",500,0,500); 
  TH1F *h_pt_a1=new TH1F("h_A_PT1","A_PT",500,0,500); 
  TH1F *h_pt_a2=new TH1F("h_A_PT2","A_PT",500,0,500); 
  TH1F *h_eta_a=new TH1F("H_A_eta","A_eta",100,-5,5); 
  TH1F *h_phi_a=new TH1F("H_A_phi","A_phi",100,-TMath::Pi(), TMath::Pi()); 
  TH1F *h_m_a=new TH1F("H_A_m","A_m",500,0,500); 

  // multiplicities before acceptance cuts plots
  TH1F *h_engen_mult_bef = new TH1F("en_gen_mult_bef", "electrons_mult", 5, 0, 5);
  TH1F *h_mngen_mult_bef = new TH1F("mn_gen_mult_bef", "muons_mult", 5, 0, 5);
  TH1F *h_leptgen_mult_bef = new TH1F("lepton_gen_mult_bef", "leptons_mult", 5, 0, 5);
  TH1F *h_qg_mult_bef = new TH1F("qg_gen_mult_bef", "quark-gluon_mult", 12, 0, 12);
  TH1F *h_bq_mult_bef= new TH1F("bquark_gen_mult_bef", "bquark_mult", 10, 0, 10);

  // Multiplicities after cuts
  TH1F *h_engen_mult = new TH1F("en_gen_mult", "electrons_mult", 5, 0, 5);
  TH1F *h_mngen_mult = new TH1F("mn_gen_mult", "muons_mult", 5, 0, 5);
  TH1F *h_leptgen_mult = new TH1F("lepton_gen_mult", "leptons_mult", 5, 0, 5);
  TH1F *h_qg_mult = new TH1F("qg_gen_mult", "quark-gluon_mult", 10, 0, 10);
  TH1F *h_bq_mult= new TH1F("bquark_gen_mult", "bquark_mult", 10, 0, 10);
 
  //deltar between aas
  TH1F *h_dr_aa = new TH1F("h_dR_aa", "dR_aa", 100, 0, 8);
  // deltar between bbs of same mom
  TH1F *h_dr_bb1 = new TH1F("h_dR_bb1", "dR_bb", 100, 0, 8);
  TH1F *h_dr_bb2 = new TH1F("h_dR_bb2", "dR_bb", 100, 0, 8);
 
 
  //scalar sum of pt between bbs of same mom
  TH1F *h_pt_bb1_vec=new TH1F("h_bb1_PT_vec","bb_PT_vec",100,0,500);
  TH1F *h_pt_bb1_scal=new TH1F("h_bb1_PT_scal","bb_PT_scal",100,0,500);

  TH1F *h_pt_bb2=new TH1F("h_bb2_PT","bb_PT",100,0,500);
  TH2F *h_bb_a_vec=new TH2F("h_a_bb_vec","bb_PT",100,0,500,100,0,500);
  TH2F *h_bb_a_scal=new TH2F("h_a_bb_scal","bb_PT",100,0,500,100,0,500);
  


  TH1F *h_pt_vv_vec=new TH1F("h_vv_PT_vec","vv_PT_vec",100,0,500);
  TH1F *h_pt_vv_scal=new TH1F("h_vv_PT_scal","vv_PT_scal",100,0,500);
  //MET
  TH1F *h_met_pt=new TH1F("met_pt","met_PT",100,0,500);
  TH1F *h_met_qq_pt=new TH1F("met_qq_l_pt","met_PT",100,0,500);
  TH1F *h_met_bb_pt=new TH1F("met_bb_pt","met_PT",100,0,500);
  TH1F *h_met_vv_pt=new TH1F("met_vv_pt","met_PT",100,0,500);

  TH1F *h_met_phi=new TH1F("met_phi","met_phi",100,-TMath::Pi(), TMath::Pi()); 
  TH1F *h_Ht=new TH1F("Ht","Ht",100,0,700);
  // plots of 4 b quarks
  TH1F *h_pt_b1=new TH1F("h_b1_PT","b1_PT",100,0,250); 
  TH1F *h_eta_b1=new TH1F("H_b1_eta","b1_eta",100,-5,5); 

  TH1F *h_pt_b2=new TH1F("h_b2_PT","b2_PT",100,0,250); 
  TH1F *h_eta_b2=new TH1F("H_b2_eta","b2_eta",100,-5,5); 

  TH1F *h_pt_b3=new TH1F("h_b3_PT","b3_PT",100,0,250); 
  TH1F *h_eta_b3=new TH1F("H_b3_eta","b3_eta",100,-5,5); 

  TH1F *h_pt_b4=new TH1F("h_b4_PT","b4_PT",100,0,250); 
  TH1F *h_eta_b4=new TH1F("H_b4_eta","b4_eta",100,-5,5); 
   // plots of 4 q-g 
  TH1F *h_pt_qg1=new TH1F("h_qg1_PT","qg1_PT",100,0,250); 
  TH1F *h_eta_qg1=new TH1F("H_qg1_eta","qg1_eta",100,-5,5); 

  TH1F *h_pt_qg2=new TH1F("h_qg2_PT","qg2_PT",100,0,250); 
  TH1F *h_eta_qg2=new TH1F("H_qg2_eta","qg2_eta",100,-5,5); 

  TH1F *h_pt_qg3=new TH1F("h_qg3_PT","qg3_PT",100,0,250); 
  TH1F *h_eta_qg3=new TH1F("H_qg3_eta","qg3_eta",100,-5,5); 

  TH1F *h_pt_qg4=new TH1F("h_qg4_PT","qg4_PT",100,0,250); 
  TH1F *h_eta_qg4=new TH1F("H_qg4_eta","qg4_eta",100,-5,5); 
  //plots of lepton 1,2
  TH1F *h_pt_l1=new TH1F("h_l1_PT","l1_PT",100,0,250); 
  TH1F *h_eta_l1=new TH1F("H_l1_eta","l1_eta",100,-5,5); 

  TH1F *h_pt_l2=new TH1F("h_l2_PT","l2_PT",100,0,250); 
  TH1F *h_eta_l2=new TH1F("H_l2_eta","l2_eta",100,-5,5); 

  //------Histograms -END ----------------
  //------Histograms -END ----------------

  //------Event Counters-------------
  Int_t n_event_lepton_test = 0;
  Int_t n_event_jet_test = 0;
  Int_t n_event_bjet_test = 0;
  Int_t n_fjets=0;
  // Int_t n_event_fjet_test 
  Int_t n_event_m2el_good=0;
  Int_t n_event_m2mn_good=0;
  Int_t n_event_2ll_good = 0;
  Int_t n_3j = 0; 
  Int_t n_4j=0;
  Int_t n_3bj = 0; 
  Int_t n_4bj = 0;
  Int_t n_2el=0;
  Int_t n_2mn=0;
  Int_t INVM=0;
  Int_t METCUT=0;
  //z->:
  Int_t zqq_light=0;
  Int_t zbb=0;
  Int_t zvv=0;
  Int_t zll=0;
  
  Int_t st4_qq=0;
  Int_t st4_bb=0;
  Int_t st4_vv=0;
  Int_t st4_ll=0;



  //---------------------------------

  // ------ Start event loop---------------
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    GetEntry(jentry);

    // Create new object vectors after configuration
    // Jets
    std::vector<TLorentzVector> vec_jet;
    std::vector<TLorentzVector> vec_bjets;
    std::vector<TLorentzVector> vec_fjet;
    std::vector<TLorentzVector> vec_jet_cc;
    // Leptons
    std::vector<TLorentzVector> vec_muons;
    std::vector<TLorentzVector> vec_ele;
    std::vector<TLorentzVector> vec_leptons;
    
    //generator level
    //before accepance cuts 
    std::vector<TLorentzVector> vec_genmu_bef;
    std::vector<TLorentzVector> vec_genele_bef;
    std::vector<TLorentzVector> vec_genlep_bef;
    std::vector<TLorentzVector> vec_genqg_bef; //quarks or gluons
    std::vector<TLorentzVector> vec_genb_bef; // only b's
    
    //after cuts
    std::vector<TLorentzVector> vec_genmu;
    std::vector<TLorentzVector> vec_genele;
    std::vector<TLorentzVector> vec_genlep;
    std::vector<TLorentzVector> vec_gena;
    std::vector<TLorentzVector> vec_genqg; //quarks or gluons
    std::vector<TLorentzVector> vec_genb; // only b's
    std::vector<TLorentzVector> vec_genbb1_truth;
    std::vector<TLorentzVector> vec_genbb2_truth;

    std::vector<TLorentzVector> vec_zll;
    std::vector<TLorentzVector> vec_zbb;
    std::vector<TLorentzVector> vec_zqq_light;
    std::vector<TLorentzVector> vec_zvv;
   
 
    
   
    for (int imc=0; imc<nmcparticles;imc++)
     {
       if( jentry <3) 
       {
	      cout << " imcparticle " << imc << " : is a " << mc_id[imc] << " , and has a mother at: " << mc_momidx[imc] <<  "("<<mc_mom[imc] << " )" <<"has status   "<<mc_status[imc] << "  and has a 4-vector p = (" << mc_en[imc] << ", " << mc_px[imc] << ", " << mc_py[imc] << ", " << mc_pz[imc] << " ) " << endl;
       }
       
        if(mc_id[imc]==25) 
            { // found the Higgs boson
	            TLorentzVector pH;
	            pH.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);

              double ptH  = pH.Pt();
              double etaH = pH.Eta();
              double phiH = pH.Phi();
              double mH   = pH.M();

              // Make histograms of the pt, eta, phi, mass of the Higgs boson:
              h_ptH->Fill(pH.Pt()); 
              h_etaH->Fill(etaH); h_phiH->Fill(phiH); h_mH->Fill(mH);
              
            }
        if(mc_id[imc]==23) 
            { // found the z boson

              TLorentzVector pz;
	            pz.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
              // Make histograms of the pt, eta, phi, mass of the z boson:
              h_pt_z->Fill(pz.Pt()); 
              h_eta_z->Fill(pz.Eta()); h_phi_z->Fill(pz.Phi()); h_m_z->Fill(pz.M());
              
            }
        if(mc_id[imc]==36) 
            { // found the A

              TLorentzVector pa;
	            pa.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
              // Make histograms of the pt, eta, phi, mass of the A:
              //h_pt_a->Fill(pa.Pt()); 
              h_eta_a->Fill(pa.Eta()); h_phi_a->Fill(pa.Phi()); h_m_a->Fill(pa.M());
              vec_gena.push_back(pa);
            }

           
        if(abs(mc_id[imc])==11) 
            { // electrons   
              TLorentzVector lep1;
              lep1.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
              vec_genele_bef.push_back(lep1);
              vec_genlep_bef.push_back(lep1);
              if (lep1.Pt()>20. && fabs(lep1.Eta())<2.4) 
                { // acceptance cuts
                  vec_genele.push_back(lep1);
                  vec_genlep.push_back(lep1);
                }
            }
        // electrons END
        if(abs(mc_id[imc])==13) 
            { //muons
              TLorentzVector lep2;
              lep2.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
              vec_genmu_bef.push_back(lep2);
              vec_genlep_bef.push_back(lep2);
              if (lep2.Pt()>20. && fabs(lep2.Eta())<2.4) 
                {// acceptance cuts
                  vec_genmu.push_back(lep2);
                  vec_genlep.push_back(lep2);
                }
              
            } // muons END
        
        //quarks or gluons
       
        if( ((abs(mc_id[imc])>=1 && abs(mc_id[imc])<=5) || abs(mc_id[imc])==21) &&  (mc_status[imc]==23 || mc_status[imc]==1))
        {
            TLorentzVector pb;
            pb.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
            vec_genqg_bef.push_back(pb);
            if (abs(mc_id[imc])==5) vec_genb_bef.push_back(pb);
            if (pb.Pt()>20. && fabs(pb.Eta())<2.4) 
            {
              // acceptance cuts :
              vec_genqg.push_back(pb);
              if (abs(mc_id[imc])==5)  vec_genb.push_back(pb);
	    
              
              if (abs(mc_id[imc])==5 && mc_mom[imc]==36 && (mc_momidx[imc]==4))  // mom at place 6 --> corresponds at bb_1
                { 
                  vec_genbb1_truth.push_back(pb);
                }
              if (abs(mc_id[imc])==5 && mc_mom[imc]==36 && (mc_momidx[imc]==5))  // mom at place 7 --> corresponds at bb_1
                { 
                  vec_genbb2_truth.push_back(pb);
                }
		
	          } // end acceptance cuts
         }// if q/g END    

         // z->qqlight/z->bb/z->vv/z->ll
         if ((abs(mc_id[imc])>=1 && abs(mc_id[imc])<=4)  && mc_mom[imc]==23)
         {
            TLorentzVector q;
            q.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
            vec_zqq_light.push_back(q);
         }
         if (abs(mc_id[imc])==5  && mc_mom[imc]==23)
         {
            TLorentzVector b;
            b.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
            vec_zbb.push_back(b);
         }
         if( (abs(mc_id[imc])==12 || abs(mc_id[imc])==14 || abs(mc_id[imc])==16) && mc_mom[imc]==23)
         {
            TLorentzVector vv;
            vv.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
            vec_zvv.push_back(vv);
            
         }
         if( (abs(mc_id[imc])==11 || abs(mc_id[imc])==13 || abs(mc_id[imc])==15) && mc_mom[imc]==23)
         {
            TLorentzVector ll;
            ll.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
            vec_zll.push_back(ll);
         } 
      } // MC particle loop END
   
    
    // Fill lepton histos:
    //sort lepton vector
    if (vec_genlep.size()>0)
    {
      for (int i=0;i<vec_genlep.size()-1;i++)
      {
          int maxIndex = i;  
          for (int j = i + 1; j < vec_genlep.size(); j++) 
            {  
              if (vec_genlep[j].Pt() > vec_genlep[maxIndex].Pt())
          {  
            maxIndex = j;  
          }  
            }   
          TLorentzVector temp;
          temp = vec_genlep[i];  
          vec_genlep[i] = vec_genlep[maxIndex];  
          vec_genlep[maxIndex] = temp; 
      }
    }
    if (vec_genlep.size()>0)
    {
      h_pt_l1->Fill(vec_genlep[0].Pt());
      h_eta_l1->Fill(vec_genlep[0].Eta());
    }
    if (vec_genlep.size()>1) 
    {
      h_pt_l2->Fill(vec_genlep[1].Pt());
      h_eta_l2->Fill(vec_genlep[1].Eta());
    }  

   // if (vec_genlep.size()>=2) 
    //  {
    //    TLorentzVector dilept;
    //    dilept=vec_genlep[0]+vec_genlep[1];
    //    h_pt_dlep->Fill(dilept.Pt());
    //    h_eta_dlep->Fill(dilept.Eta());
    //    h_phi_dlep->Fill(dilept.Phi());
    //    h_m_dlep->Fill(dilept.M());
    //  }
    
    
    
      // find dr between aas
    if (vec_gena.size()>1) h_dr_aa->Fill(getDeltaR(vec_gena[0],vec_gena[1]));

    //sort qluon-quark vector
    if (vec_genqg.size()>0)
    {
      for (int i=0;i<vec_genqg.size()-1;i++)
      {
          int maxIndex = i;  
          for (int j = i + 1; j < vec_genqg.size(); j++) 
            {  
              if (vec_genqg[j].Pt() > vec_genqg[maxIndex].Pt())
          {  
            maxIndex = j;  
          }  
            }   
          TLorentzVector temp;
          temp = vec_genqg[i];  
          vec_genqg[i] = vec_genqg[maxIndex];  
          vec_genqg[maxIndex] = temp; 
      }
    }
    //sort bquark vector
    
    if (vec_genb.size()>0){
      for (int i=0;i<vec_genb.size()-1;i++)
      {
        int maxIndex = i;  
        for (int j = i + 1; j < vec_genb.size(); j++) 
          {  
            if (vec_genb[j].Pt() > vec_genb[maxIndex].Pt())
        {  
          maxIndex = j;  
        }  
          }  
        TLorentzVector temp;
        temp = vec_genb[i];  
        vec_genb[i] = vec_genb[maxIndex];  
        vec_genb[maxIndex] = temp; 
      } 
    }
    //check z->;
    bool is_vv(false);
    bool is_qq_light(false);
    bool is_bb(false);
    bool is_ll(false);
    if (vec_zll.size()>1)
    { zll++;
      is_ll=true;}
    if (vec_zqq_light.size()>1) 
    { zqq_light++;
      is_qq_light=true;}
    if (vec_zbb.size()>1) 
    { zbb++;  
      is_bb=true;}  
    if (vec_zvv.size()>1) 
    { zvv++;  
      is_vv=true;}  
    
    //plot dr and  vector sum pt of bbs of same mom and z->vv event
    if (vec_genbb1_truth.size()>1 && is_vv ){
     h_dr_bb1->Fill(getDeltaR(vec_genbb1_truth[0],vec_genbb1_truth[1]));
     h_pt_bb1_vec->Fill((vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt());
     //h_pt_bb1_scal->Fill(vec_genbb1_truth[0].Pt()+vec_genbb1_truth[1].Pt());
     h_pt_a1->Fill(vec_gena[0].Pt()); 
     h_bb_a_vec->Fill(vec_gena[0].Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt());
    // h_bb_a_scal->Fill(vec_gena[0].Pt(),vec_genbb1_truth[0].Pt()+vec_genbb1_truth[1].Pt());
    }
    if (vec_genbb2_truth.size()>1 && is_vv){
      h_dr_bb2->Fill(getDeltaR(vec_genbb2_truth[0],vec_genbb2_truth[1]));
      h_pt_bb2->Fill((vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt());
      h_pt_a2->Fill(vec_gena[1].Pt()); 
    }
    //plot scalr sum pt of vvs
    if(vec_zvv.size()>1){
      h_pt_vv_scal->Fill(vec_zvv[0].Pt()+vec_zvv[1].Pt());
      h_pt_vv_vec->Fill((vec_zvv[0]+vec_zvv[1]).Pt());

    } 
     // plot eta and pt of 4 b quarks 
    if (vec_genb.size()>0)
    {
      h_pt_b1->Fill(vec_genb[0].Pt());
      h_eta_b1->Fill(vec_genb[0].Eta());
    }
    if (vec_genb.size()>1) 
    {
      h_pt_b2->Fill(vec_genb[1].Pt());
      h_eta_b2->Fill(vec_genb[1].Eta());
    }  
    if (vec_genb.size()>2)
    {
      h_pt_b3->Fill(vec_genb[2].Pt());
      h_eta_b3->Fill(vec_genb[2].Eta());
    }
    if (vec_genb.size()>3)
    {
      h_pt_b4->Fill(vec_genb[3].Pt());
      h_eta_b4->Fill(vec_genb[3].Eta());
    }

    //plot eta and pt of 4 gluons-quarks
    if (vec_genqg.size()>0)
    {
      h_pt_qg1->Fill(vec_genqg[0].Pt());
      h_eta_qg1->Fill(vec_genqg[0].Eta());
    }
    if (vec_genqg.size()>1) 
    {
      h_pt_qg2->Fill(vec_genqg[1].Pt());
      h_eta_qg2->Fill(vec_genqg[1].Eta());
    }  
    if (vec_genqg.size()>2)
    {
      h_pt_qg3->Fill(vec_genqg[2].Pt());
      h_eta_qg3->Fill(vec_genqg[2].Eta());
    }
    if (vec_genqg.size()>3)
    {
      h_pt_qg4->Fill(vec_genqg[3].Pt());
      h_eta_qg4->Fill(vec_genqg[3].Eta());
    }
    
    // Make Histograms: lepton mulciplicity, q/g multiplicity , b-quark multiplicity --> Check effect of the Acceptance cuts (pt/eta)
    h_engen_mult_bef->Fill(vec_genele_bef.size());
    h_mngen_mult_bef->Fill(vec_genmu_bef.size());
    h_leptgen_mult_bef->Fill(vec_genlep_bef.size());
    h_qg_mult_bef->Fill(vec_genqg_bef.size());
    h_bq_mult_bef->Fill(vec_genb_bef.size());
    //after acceptance cuts
    h_engen_mult->Fill(vec_genele.size());
    h_mngen_mult->Fill(vec_genmu.size());
    h_leptgen_mult->Fill(vec_genlep.size());
    h_qg_mult->Fill(vec_genqg.size());
    h_bq_mult->Fill(vec_genb.size());
    
    
    // END GENERATOR LEVEL ANALYSIS




    // RECONSTRUCTION LEVEL ANALYSIS:


    for (int i = 0; i < mn; i++)
    {
      TLorentzVector p_muon;
      p_muon.SetPxPyPzE(mn_px[i], mn_py[i], mn_pz[i], mn_en[i]);

      if (p_muon.Pt() < 20. || std::fabs(p_muon.Eta()) > 2.4)
        continue;

      // id + Isolation
      if (mn_passId[i] && mn_passIso[i])
      {
        vec_muons.push_back(p_muon);
	      vec_leptons.push_back(p_muon);
      }
    }

    for (int i = 0; i < en; i++)
    {
      TLorentzVector p_electron;
      p_electron.SetPxPyPzE(en_px[i], en_py[i], en_pz[i], en_en[i]);

      if (p_electron.Pt() < 20. || std::fabs(p_electron.Eta()) > 2.4)
        continue;

      if (en_passId[i] && en_passIso[i])
      {
        vec_ele.push_back(p_electron);
	      vec_leptons.push_back(p_electron);
      }
    }
    // jets & cross cleaning
    for (int i = 0; i < jet; i++)
    {
      bool overlap = false;
      double dR1 = 0;
      

      TLorentzVector p_jet;
      p_jet.SetPxPyPzE(jet_px[i], jet_py[i], jet_pz[i], jet_en[i]);

      if (p_jet.Pt() < 20. || std::fabs(p_jet.Eta()) > 2.5)
        continue;
      for (int mn_count = 0; mn_count < vec_muons.size(); mn_count++)
      {
        dR1 = getDeltaR(p_jet, vec_muons[mn_count]);
        h_dR_jet_mn_before->Fill(dR1);
        if (dR1 < 0.4){
          overlap = true;
         }
         else {
            h_dR_jet_mn_after->Fill(dR1);
         }
      }

      for (int en_count = 0; en_count < vec_ele.size(); en_count++)
      {
        dR1 = getDeltaR(p_jet, vec_ele[en_count]);
        h_dR_jet_en_before->Fill(dR1);
        if (dR1 < 0.4){
         overlap = true;
        }
        else {
         h_dR_jet_en_after->Fill(dR1);
        }
      }

      if (!overlap) vec_jet.push_back(p_jet);
        
       //bjets
      if (! overlap && jet_btag1[i] > 0.4)
      {
        vec_bjets.push_back(p_jet);
      }
    }

    //fjets
    int index1=0;
    int index2=0;
    int count=0;
    for (int i = 0; i < fjet; i++)
      {  
         bool overlap = false;
         double dR1 = 0;
         TLorentzVector p_fjet;
         p_fjet.SetPxPyPzE(fjet_px[i], fjet_py[i], fjet_pz[i], fjet_en[i]);
         if(fjet_subjet_count[i]<2) continue;
         if (p_fjet.Pt() <30 || std::fabs(p_fjet.Eta()) > 2.5) continue;
         
         for (int mn_count = 0; mn_count < vec_muons.size(); mn_count++)
         {
         dR1 = getDeltaR(p_fjet, vec_muons[mn_count]);
         h_dR_fjet_mn_before->Fill(dR1);
         if (dR1 < 0.8){
            overlap = true;
            }
            else {
               h_dR_fjet_mn_after->Fill(dR1);
            }
         }

         for (int en_count = 0; en_count < vec_ele.size(); en_count++)
         {
         dR1 = getDeltaR(p_fjet, vec_ele[en_count]);
         h_dR_fjet_en_before->Fill(dR1);
         if (dR1 < 0.8){
            overlap = true;
         }
         else {
            h_dR_fjet_en_after->Fill(dR1);
         }
         }

         if (!overlap)
         { 
          count++;
          if (count==1)index1=i;
          if (count==2)index2=i;
          vec_fjet.push_back(p_fjet);
         }
      }
 
    
      //jet f jet cross-cleaning
    for (int j = 0; j < vec_jet.size(); j++)  {
      if(dR_j_fj_min(vec_jet[j],vec_fjet)>0.8) vec_jet_cc.push_back(vec_jet[j]);
    }  
    
     


    //-----End config ------------

    // ----- cut and count -------
   

    // lepton cuts---------------------------------------------------

   // if (vec_ele.size() + vec_muons.size() < 2) continue;
    if ( vec_leptons.size() !=0) continue;
    
    n_event_lepton_test++;

    h_met_pt->Fill(met_pt);
    h_met_phi->Fill(met_phi);
    if(is_qq_light)h_met_qq_pt->Fill(met_pt);
    if(is_bb)h_met_bb_pt->Fill(met_pt);
    if(is_vv)h_met_vv_pt->Fill(met_pt);

    
    // control plots mult
    h_en_mult->Fill(vec_ele.size());
    h_mn_mult->Fill(vec_muons.size());
    //h_lepton_mult->Fill(vec_ele.size() + vec_muons.size());
    h_jet_mult->Fill(vec_jet.size());
    h_bjet_mult->Fill(vec_bjets.size());
    h_fjet_mult->Fill(vec_fjet.size());
    h_jet_cc_mult1->Fill(vec_jet_cc.size());
    // ele kinematics pt, eta, phi  e[1]
    if (vec_ele.size() > 0)
      {
        h_en1_pt->Fill(vec_ele[0].Pt());
        h_en1_eta->Fill(vec_ele[0].Eta());
        h_en1_phi->Fill(vec_ele[0].Phi());
      }
    if (vec_ele.size() > 1)
      { 
        n_2el++;
        h_en2_pt->Fill(vec_ele[1].Pt());
        h_en2_eta->Fill(vec_ele[1].Eta());
        h_en2_phi->Fill(vec_ele[1].Phi());
        TLorentzVector diele=vec_ele[1]+vec_ele[0];
      }
    
    if (vec_muons.size() > 0)
      {
        h_mn1_pt->Fill(vec_muons[0].Pt());
        h_mn1_eta->Fill(vec_muons[0].Eta());
        h_mn1_phi->Fill(vec_muons[0].Phi());
      }
    if (vec_muons.size() > 1)
      {  
        n_2mn++;
        h_mn2_pt->Fill(vec_muons[1].Pt());
        h_mn2_eta->Fill(vec_muons[1].Eta());
        h_mn2_phi->Fill(vec_muons[1].Phi());
        TLorentzVector dimuon=vec_muons[1]+vec_muons[0];
        h_2mn_inv_m1->Fill(dimuon.M());
      }
        
    
     // inv mass cuts
    
   // float dileptonmass=(vec_leptons[0]+vec_leptons[1]).M();
   // if (dileptonmass>100. || dileptonmass < 80.) continue;
    if (met_pt < 80.) continue;
    METCUT++;
    //INVM++;
    h_fjet_mult_after1->Fill(vec_fjet.size());
    h_jet_cc_mult2->Fill(vec_jet_cc.size());
   // if (diele.M()>80 && diele.M()<100)
     // {
     //   n_event_m2el_good++;
     //   h_2en_inv_m->Fill(diele.M());
     //   h_2en_pt->Fill(diele.Pt());
     // }
    //if (dimuon.M()>80 && dimuon.M()<100)
    //  {
    //          n_event_m2mn_good++;
    //          h_2mn_inv_m2->Fill(dimuon.M());
    //          h_2mn_pt->Fill(dimuon.Pt());
    //  }
        
        // Fat-jet analysis:
       
      if (vec_fjet.size()>0)
      { 
         //fj1 kinematics
        h_fjet1_pt->Fill(vec_fjet[0].Pt());
        h_fjet1_eta->Fill(vec_fjet[0].Eta());
        h_fjet1_phi->Fill(vec_fjet[0].Phi());
        //xbb discriminants
        float bbtag_fj1=fjet_btag10[index1]/(fjet_btag10[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);  
        h_fjet1_btagXbb->Fill(bbtag_fj1);
        float bbccqqtagfj1=(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1])/(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);
        h_fjet1_btagXbbXccXqq->Fill(bbccqqtagfj1);
        // fj1 sdmass 
        fj_sd_mass1->Fill(fjet_softdropM[index1]);
        fj_pt_sd_mass1->Fill(vec_fjet[0].Pt(), fjet_softdropM[index1]);
        subcount1->Fill(fjet_subjet_count[index1]);
        bool matched(false);
        bool matched1(false);
        bool matched2(false);
        if(isSignal) 
        {
         if(fjet_bkg_matched_qg_pt(vec_fjet[0],vec_genbb1_truth)>0)
         {
               h_fj1_pt_bb_pt->Fill(vec_fjet[0].Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt());
               matched1=true;
            } 
            else if( fjet_bkg_matched_qg_pt(vec_fjet[0],vec_genbb2_truth)>0)
            {
               h_fj1_pt_bb_pt->Fill(vec_fjet[0].Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt());
               matched2=true;  
            }
         if(matched1||matched2)matched=true;
        }
        else {
         bool matched(false);
         if( fjet_bkg_matched_qg_pt(vec_fjet[1],vec_genqg)>0){
          h_fj1_pt_bb_pt->Fill(vec_fjet[0].Pt(),fjet_bkg_matched_qg_pt(vec_fjet[0],vec_genqg));
          matched=true;
         }
        }
        if(matched){
         h_fj1_pt_mat->Fill(vec_fjet[0].Pt());
         fj_sd_mass1_mat->Fill(fjet_softdropM[index1]);
         subcount1_mat->Fill(fjet_subjet_count[index1]);
        }
	   }
      if (vec_fjet.size()>1)
      {  
          if(verbose) { // debugging
	       cout << "Fatjet1 pT = " << vec_fjet[0].Pt() << " and fatjet2 pT= " << vec_fjet[1].Pt() << endl;
	        }
          //fj2 kinematics
         h_fjet2_pt->Fill(vec_fjet[1].Pt());
         h_fjet2_eta->Fill(vec_fjet[1].Eta());
         h_fjet2_phi->Fill(vec_fjet[1].Phi());
          //fj2 xbb discriminants
         float bbtag_fj2=fjet_btag10[index2]/(fjet_btag10[index2]+fjet_btag13[index2]+fjet_btag14[index2]+fjet_btag15[index2]+fjet_btag16[index2]+fjet_btag17[index2]);  
         h_fjet2_btagXbb->Fill(bbtag_fj2);
         float bbccqqtagfj2=(fjet_btag10[index2]+fjet_btag11[index2]+fjet_btag12[index2])/(fjet_btag10[index2]+fjet_btag11[index2]+fjet_btag12[index2]+fjet_btag13[index2]+fjet_btag14[index2]+fjet_btag15[index2]+fjet_btag16[index2]+fjet_btag17[index2]);
         h_fjet2_btagXbbXccXqq->Fill(bbccqqtagfj2);
         float bbtag_fj1=fjet_btag10[index1]/(fjet_btag10[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);  
         float bbccqqtagfj1=(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1])/(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);
         h_fj1_fj2_btagXbb->Fill(bbtag_fj1,bbtag_fj2);
         h_fj1_fj2_btagXbbXccXqq->Fill(bbccqqtagfj1,bbccqqtagfj2);
          //fj2 sdmass
         fj_sd_mass2->Fill(fjet_softdropM[index2]);
         fj_sd_mass1_2->Fill(fjet_softdropM[index1],fjet_softdropM[index2]);
         fj_pt_sd_mass2->Fill(vec_fjet[1].Pt(), fjet_softdropM[index2]);
         subcount2->Fill(fjet_subjet_count[index2]);
         bool matched(false);
         bool matched1(false);
         bool matched2(false);
         if(isSignal) 
         {
            if(  fjet_bkg_matched_qg_pt(vec_fjet[1],vec_genbb1_truth)>0)
            {
                  h_fj2_pt_bb_pt->Fill(vec_fjet[1].Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt());
                  matched1=true;
               } 
               else if(  fjet_bkg_matched_qg_pt(vec_fjet[1],vec_genbb2_truth)>0)
               {
                  h_fj2_pt_bb_pt->Fill(vec_fjet[1].Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt());
                  matched2=true;  
               }
            if(matched1||matched2)matched=true;
         }
         else{
            if( fjet_bkg_matched_qg_pt(vec_fjet[1],vec_genqg)>0){
            h_fj2_pt_bb_pt->Fill(vec_fjet[1].Pt(),fjet_bkg_matched_qg_pt(vec_fjet[1],vec_genqg));
            matched=true;
            }
         }
         if(matched){
            h_fj2_pt_mat->Fill(vec_fjet[1].Pt());
            fj_sd_mass2_mat->Fill(fjet_softdropM[index2]);
            subcount2_mat->Fill(fjet_subjet_count[index2]);
        }
      }
      
      
      
   
       
        //---------------JET CUTS ------------------------------------------------
        // if (vec_jet.size() >= 1)
        // {
        //   n_event_jet_test++;
        //   //particle multiplicity after cuts
          
        //   if(vec_fjet.size()>=1){
        //     n_fjets++;
        //     h_fjet_mult_after->Fill(vec_fjet.size());
        //     h_jet_mult_after->Fill(vec_jet.size());
        //     h_bjet_mult_after->Fill(vec_bjets.size());
        //     h_mn_mult_after->Fill(vec_muons.size());
        //     h_en_mult_after->Fill(vec_ele.size());
        //     h_lepton_mult_after->Fill(vec_ele.size()+vec_muons.size());
        //     if (dimuon.M()>80 && dimuon.M()<100) h_2mn_inv_m3->Fill(dimuon.M());
        //   }
        // }
	  
	  
	
    
     
	    // AT LEAST 2 FAT JETS
	 if (vec_fjet.size() < 2) continue;
    n_fjets++;
	 h_fjet_mult_after2->Fill(vec_fjet.size());
	 h_jet_mult_after->Fill(vec_jet.size());
    h_jet_cc_mult3->Fill(vec_jet_cc.size());
    h_Ht->Fill(Ht(vec_jet_cc,vec_fjet));
	 // h_met_pt->Fill(met_pt);
    //h_met_phi->Fill(met_phi);
         
    //z-?
    if(is_vv) st4_vv++;
    if(is_qq_light) st4_qq++;
    if(is_bb) st4_bb++;
    if(is_ll) st4_ll++;
	  //fill jet histos
	  if (vec_jet.size() > 0)
	    {
        h_jet1_pt->Fill(vec_jet[0].Pt());
        h_jet1_eta->Fill(vec_jet[0].Eta());
        h_jet1_phi->Fill(vec_jet[0].Phi());
      }
	  if (vec_jet.size() > 1)
            {
	      h_jet2_pt->Fill(vec_jet[1].Pt());
	      h_jet2_eta->Fill(vec_jet[1].Eta());
	      h_jet2_phi->Fill(vec_jet[1].Phi());
            }
	  if (vec_jet.size() > 2)
            {
	      h_jet3_pt->Fill(vec_jet[2].Pt());
	      h_jet3_eta->Fill(vec_jet[2].Eta());
	      h_jet3_phi->Fill(vec_jet[2].Phi());
	      n_3j++;
            }
	  if (vec_jet.size() > 3)
            {
	      h_jet4_pt->Fill(vec_jet[3].Pt());
	      h_jet4_eta->Fill(vec_jet[3].Eta());
	      h_jet4_phi->Fill(vec_jet[3].Phi());
	      n_4j++;   
            }
    
      
            
  } // END EVENT LOOP

  std::cout << "reconstruction level: "  << std::endl;
  std::cout << "number of events: " << nentries << std::endl;
  std::cout << "number of 0 lepton events: " << n_event_lepton_test << std::endl;
  std::cout << "number of events after met_pt cut: " << METCUT << std::endl;
  std::cout << "number of events after at least 2 fjet cuts: " << n_fjets << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and z->vv: " << st4_vv << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and Z->qq_light: " << st4_qq << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and Z->bb: " << st4_bb << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and Z->ll: " << st4_ll << std::endl;
  std::cout << "generator level: "  << std::endl;
  std::cout << "number of events with z->qq_light: " << zqq_light << std::endl;
  std::cout << "number of events with z->bb: " << zbb << std::endl;
  std::cout << "number of events with z->vv: " << zvv << std::endl;
  std::cout << "number of events with z->ll: " << zll << std::endl;
 // tabulate::Table cut_flow_table;
 // cut_flow_table.add_row({"Requirements", "No. of Events", "Passed/Total 100%"});
 //  cut_flow_table.add_row({"Raw events", std::to_string(nentries), "100%"});
 // cut_flow_table.add_row({"N Leptons = 0 or >=2, ", std::to_string(n_event_lepton_test), std::to_string(100 * ((float)n_event_lepton_test) / (float)nentries) + "%"});
 // cut_flow_table.add_row({"N jet >= 3", std::to_string(n_event_jet_test), std::to_string(100 * ((float)n_event_jet_test) / (float)nentries) + "%"});
  // cut_flow_table.add_row({"N fjet >=1 ", std::to_string(n_event_fjet_test), std::to_string(100 * ((float)n_event_fjet_test) / (float)nentries) + "%"});

  //  std::cout << cut_flow_table << std::endl;

  std::ofstream out_file(file_name+".csv");
   out_file<< "Cuts," << "No. of events," << "Efficiency" <<std::endl;
   out_file << "(raw)," << nentries << ","<< nentries/nentries * 100 <<std::endl;
   out_file << "0 lept cut," << n_event_lepton_test << ","<< (float)n_event_lepton_test/(float)nentries*100<<std::endl;
   out_file << "MET PT cut," <<METCUT<< ","<< (float)INVM/(float)nentries*100<<std::endl;
   out_file << "2fjet cut," << n_fjets<< ","<< (float)n_fjets/(float)nentries*100<<std::endl;
   out_file.close();

  TFile ff("histos_" + file_name, "RECREATE");
  // control plots
  h_en_mult->Write();
  h_mn_mult->Write();
  h_jet_mult->Write();
  h_fjet_mult->Write();
  h_bjet_mult->Write();
  h_jet_cc_mult1->Write();
  h_jet_cc_mult2->Write();
  h_jet_cc_mult3->Write();

 
  //----lepton cuts
  h_en1_pt->Write();
  h_en1_eta->Write();
  h_en1_phi->Write();
  h_en2_pt->Write();
  h_en2_eta->Write();
  h_en2_phi->Write();
  //lept_inv_m
  h_2en_inv_m->Write();
  h_2mn_inv_m1->Write();
  h_2mn_inv_m2->Write();
  h_2mn_inv_m3->Write();

  h_mn1_pt->Write();
  h_mn1_eta->Write();
  h_mn1_phi->Write();
  h_mn2_pt->Write();
  h_mn2_eta->Write();
  h_mn2_phi->Write();
  // jet cuts--------

  h_jet1_pt->Write();
  h_jet1_eta->Write();
  h_jet1_phi->Write();
  h_jet2_pt->Write();
  h_jet2_eta->Write();
  h_jet2_phi->Write();
  h_jet3_pt->Write();
  h_jet3_eta->Write();
  h_jet3_phi->Write();
  h_jet4_pt->Write();
  h_jet4_eta->Write();
  h_jet4_phi->Write();
  //bjets
  h_bjet1_pt->Write();
  h_bjet1_eta->Write();
  h_bjet1_phi->Write();
  h_bjet2_pt->Write();
  h_bjet2_eta->Write();
  h_bjet2_phi->Write();
  h_bjet3_pt->Write();
  h_bjet3_eta->Write();
  h_bjet3_phi->Write();
  h_bjet4_pt->Write();
  h_bjet4_eta->Write();
  h_bjet4_phi->Write();
  //fjets
  h_fjet1_pt->Write();
  h_fjet1_eta->Write();
  h_fjet1_phi->Write();
  h_fjet2_pt->Write();
  h_fjet2_eta->Write();
  h_fjet2_phi->Write();


  h_fjet1_btagXbbXccXqq->Write();
  h_fjet2_btagXbbXccXqq->Write();
  h_fjet1_btagXbb->Write();
  h_fjet2_btagXbb->Write();
  h_fj1_fj2_btagXbbXccXqq->Write();
  h_fj1_fj2_btagXbb->Write();
 


  h_mn_mult_after->Write();
  h_en_mult_after->Write();
  h_jet_mult_after->Write();
  h_bjet_mult_after->Write();
  h_fjet_mult_after1->Write();
  h_fjet_mult_after2->Write();


  h_2en_pt->Write();
  h_2mn_pt->Write();
  
  h_higgs_inv_m->Write();
  h_higgs_pt->Write();
  h_higgs_eta->Write();
  h_higgs_phi->Write();
  //deltar
  h_dR_jet_en_before ->Write();
  h_dR_jet_mn_before ->Write(); 
  h_dR_jet_en_after ->Write();
  h_dR_jet_mn_after->Write();

  h_dR_fjet_en_before ->Write();
  h_dR_fjet_mn_before ->Write(); 
  h_dR_fjet_en_after ->Write();
  h_dR_fjet_mn_after->Write();


  h_dR_j_fj->Write();
  h_dR_j_fj_min->Write();
  h_fj1_pt_mat->Write();
  h_fj2_pt_mat->Write();
  h_pt_dR_min->Write();
  h_fj1_pt_bb_pt->Write();
  h_fj2_pt_bb_pt->Write();

  pt_j->Write();
  pt_fj->Write();
  pt_j_matched->Write();
  pt_fj_matched->Write();

  
  fj_sd_mass1_mat->Write();
  fj_sd_mass2_mat->Write();
  fj_sd_mass1->Write();
  fj_sd_mass2->Write();
  fj_pt_sd_mass1->Write();
  fj_pt_sd_mass2->Write();
  fj_sd_mass1_2_mat->Write();
  


  //write gen level
  

  h_ptH->Write(); 
  h_etaH->Write();
  h_phiH->Write();
  h_mH->Write();

  h_pt_dlep->Write();
  h_eta_dlep->Write();
  h_phi_dlep->Write();
  h_m_dlep->Write();

  h_pt_z->Write();
  h_eta_z->Write();
  h_phi_z->Write();
  h_m_z->Write();

  h_pt_a->Write();
  h_pt_a1->Write();
  h_pt_a2->Write();

  h_eta_a->Write();
  h_phi_a->Write();
  h_m_a->Write();
    
  h_engen_mult_bef->Write();
  h_mngen_mult_bef->Write();
  h_leptgen_mult_bef->Write();
  h_qg_mult_bef->Write();
  h_bq_mult_bef->Write();
    
  h_engen_mult->Write();
  h_mngen_mult->Write();
  h_leptgen_mult->Write();
  h_qg_mult->Write();
  h_bq_mult->Write();
  

 
  h_dr_aa->Write();
  h_dr_bb1->Write();
  h_pt_bb1_scal->Write();
  h_pt_bb1_vec->Write();
  h_dr_bb2->Write();
  h_pt_bb2->Write();
 
 
  h_pt_vv_scal->Write();
  h_pt_vv_vec->Write();


  h_pt_b1->Write();
  h_eta_b1->Write();
  h_pt_b2->Write();
  h_eta_b2->Write();
  h_pt_b3->Write();
  h_eta_b3->Write();
  h_pt_b4->Write();
  h_eta_b4->Write();

  h_pt_qg1->Write();
  h_eta_qg1->Write();
  h_pt_qg2->Write();
  h_eta_qg2->Write();
  h_pt_qg3->Write();
  h_eta_qg3->Write();
  h_pt_qg4->Write();
  h_eta_qg4->Write();

  h_pt_l1->Write();
  h_eta_l1->Write();
  h_pt_l2->Write();
  h_eta_l2->Write();

  h_Ht->Write();
  h_met_pt->Write();
  h_met_qq_pt->Write();
  h_met_bb_pt->Write();
  h_met_vv_pt->Write();

  h_met_phi->Write();
  h_bb_a_vec->Write();
  h_bb_a_scal->Write();

  subcount1->Write();
  subcount1_mat->Write();
  subcount2->Write();
  subcount2_mat->Write();


}
//deltar
double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2)
{
  double delta_phi;
  double delta_eta;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_eta = vec_1.Eta() - vec_2.Eta();

  return std::sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
}

bool jet_matched(TLorentzVector jet,std::vector<TLorentzVector> vecb) 
{
  bool matched=false;
  float dR_min=999;
  float dR=0;
   for (int i=0; i<vecb.size(); i++) {
      dR=getDeltaR(jet, vecb[i]);
      if (dR<dR_min) dR_min=dR;
   }
   if (dR_min<0.4) matched=true;
   return matched; 
}
//bool fjet_matched(TLorentzVector fjet) 
//{
  //bool matched1=false;
  //bool matched2=false;
  //bool matched=false;
  //if(vec_genbb1_truth.size()>1)
  //{
   //float dR[vec_genbb1_truth.size().size()];
   //for (int i=0; i<vec_genbb1_truth.size(); i++) 
   //{
     // dR[i]=getDeltaR(fjet, vec_genbb1_truth[i]);
   //}
   //if (dR[0]<0.8 && dR[1]<0.8) matched1=true;
   //}
   //if(vec_genbb2_truth.size()>1)
  //{
   //float dR[vec_genbb2_truth.size()];
   //for (int i=0; i<vec_genbb2_truth.size(); i++)
   // {
     // dR[i]=getDeltaR(fjet, vec_genbb2_truth[i]);
   //}
   //if (dR[0]<0.8 && dR[1]<0.8) matched2=true;
   //}
   //if(matched1||matched2)matched=true;

   //return matched; 
//}

bool fjet_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb) 
{
  bool matched=false;
  float dR[vecb.size()];
   for (int i=0; i<vecb.size(); i++) {
      dR[i]=getDeltaR(fjet, vecb[i]);
      
   }
   if (dR[0]<0.8 && dR[1]<0.8) matched=true;
   return matched; 
}
bool fjet_bkg_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb) 
{
  bool matched=false;
  int n=-1;
 // if(vecb.size()>0){
  float dR;
   for (int i=0; i<vecb.size(); i++) 
   {
      dR=getDeltaR(fjet, vecb[i]);
      if (dR<0.8) n++;
   }
   if (n>0) matched=true;
 // }
   return matched; 
}
float fjet_bkg_matched_qg_pt(TLorentzVector fjet,std::vector<TLorentzVector> vecb)
{
  float qg_pt=-999.;

  std::vector<TLorentzVector> vecmin;
  TLorentzVector vecqgpt;
bool ismatched(false);
  for (int i=0; i<vecb.size(); i++)
   {
ismatched=false;
      float dR=getDeltaR(fjet, vecb[i]);
      if (dR<0.8) {
ismatched=true;
     vecmin.push_back(vecb[i]);
   }
   }
   for (int i=0; i<vecmin.size(); i++)
   {
      vecqgpt+=vecmin[i];
   }
  
if(ismatched) qg_pt=vecqgpt.Pt();
   return qg_pt;
}
double dR_j_fj_min(TLorentzVector jet,std::vector<TLorentzVector> vecfjet)
{
      double dR[vecfjet.size()];
      double dR_min=0;
      int  min_ind=0;
      
      for (int fj = 0; fj < vecfjet.size(); fj++)
      {
         dR[fj] = getDeltaR(jet, vecfjet[fj]);
         if(dR[fj] <dR[min_ind]) min_ind = fj; 
      }
      dR_min=dR[min_ind];
      return dR_min;
}
double Ht(std::vector<TLorentzVector> vecjet,std::vector<TLorentzVector> vecfjet)
{
   double Ht=0;
   if (vecjet.size()>0){
   for(int i=0;i< vecjet.size();i++) Ht+=vecjet[i].Pt();}
   if(vecfjet.size()>0){
   for(int j=0;j< vecfjet.size();j++) Ht+=vecfjet[j].Pt();}
   return Ht;
}

     
