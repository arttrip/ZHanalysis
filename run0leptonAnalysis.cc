#define YEAR_2017

#include <iostream>
#include <map>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/TMVAReader.h"
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/bsmhiggs_fwk/interface/PDFInfo.h"
//#include "UserCode/bsmhiggs_fwk/interface/rochcor2016.h" 
#include "UserCode/bsmhiggs_fwk/interface/RoccoR.h" 
//#include "UserCode/bsmhiggs_fwk/interface/muresolution_run2.h" 
#include "UserCode/bsmhiggs_fwk/interface/LeptonEfficiencySF.h"
//#include "UserCode/bsmhiggs_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/bsmhiggs_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/bsmhiggs_fwk/interface/METUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/XYMETCorrection.h"
//#include "UserCode/bsmhiggs_fwk/interface/BTagUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/EventCategory.h"
#include "UserCode/bsmhiggs_fwk/interface/statWgt.h"

#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"

#include <unistd.h>
using namespace std;

double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2);
bool jet_matched(TLorentzVector jet,std::vector<TLorentzVector> vecb);
double dR_j_fj_min(TLorentzVector jet,std::vector<TLorentzVector> vecfjet);
float fjet_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb) ;
double Ht(std::vector<TLorentzVector> vecjet,std::vector<TLorentzVector> vecfjet);
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec);

int main(int argc, char* argv[])
{
  //##################################################################################
  //##########################    GLOBAL INITIALIZATION     ##########################
  //##################################################################################

  bool verbose(false);




  SmartSelectionMonitor mon;
  //reconstruction level
  //dR between leptons-jets/leptons-fjets/jets-fjets: before and after cross-cleaning
   // -- -- -- Delta R -- -- --
  //dR jets-leptons
  TH1F *h_dR_jet_en_before =(TH1F*) mon.addHistogram( new TH1F("h_dR_jet_en_before", "dR_fjet_en_before", 100, 0, 8));
  TH1F *h_dR_jet_mn_before =(TH1F*) mon.addHistogram( new TH1F("h_dR_jet_mn_before", "dR_fjet_mn_before", 100, 0, 8));
  TH1F *h_dR_jet_en_after =(TH1F*) mon.addHistogram( new TH1F("h_dR_jet_en_after", "dR_fjet_en_after", 100, 0, 8));
  TH1F *h_dR_jet_mn_after =(TH1F*) mon.addHistogram (new TH1F("h_dR_jet_mn_after", "dR_fjet_mn_after", 100, 0, 8));
  //dR fjets-leptons 
  TH1F *h_dR_fjet_en_before =(TH1F*) mon.addHistogram (new TH1F("h_dR_fjet_en_before", "dR_fjet_en_before", 100, 0, 8));
  TH1F *h_dR_fjet_mn_before =(TH1F*) mon.addHistogram( new TH1F("h_dR_fjet_mn_before", "dR_fjet_mn_before", 100, 0, 8));
  TH1F *h_dR_fjet_en_after =(TH1F*) mon.addHistogram( new TH1F("h_dR_fjet_en_after", "dR_fjet_en_after", 100, 0, 8));
  TH1F *h_dR_fjet_mn_after =(TH1F*) mon.addHistogram (new TH1F("h_dR_fjet_mn_after", "dR_fjet_mn_after", 100, 0, 8));
  //dR jets-fjets
  TH1F *h_dR_j_fj =(TH1F*) mon.addHistogram( new TH1F("h_dR_j-fj", "dR_jet-fjet", 100, 0, 8));
  TH1F *h_dR_j_fj_min =(TH1F*) mon.addHistogram( new TH1F("h_dR_j-fj_min", "dR_jet-fjet_min", 100, 0,8));
  
   // fjets
  //fjet mult
  TH1F *h_fjet_mult=(TH1F*) mon.addHistogram( new TH1F("fjet_mult", "fjets", 10, 0, 10));
  TH1F *h_fjet_mult_after1 =(TH1F*) mon.addHistogram( new TH1F("fjet_mult_after1", "fJets", 10, 0, 10));
  TH1F *h_fjet_mult_after2 =(TH1F*) mon.addHistogram( new TH1F("fjet_mult_after2", "fJets", 10, 0, 10));
  //fjet kinematics
  TH1F *h_fjet1_pt =(TH1F*) mon.addHistogram( new TH1F("h_fjet1_pt", "fjet1_pt", 100, 0, 500));
  TH1F *h_fjet1_eta =(TH1F*) mon.addHistogram( new TH1F("h_fjet1_eta", "fjet1_eta", 100, -5, 5));
  TH1F *h_fjet1_phi =(TH1F*) mon.addHistogram( new TH1F("h_fjet1_phi", "fjet1_phi", 100, -TMath::Pi(), TMath::Pi()));
  TH1F *h_fjet2_pt =(TH1F*) mon.addHistogram( new TH1F("h_fjet2_pt", "fjet2_pt", 100, 0, 500));
  TH1F *h_fjet2_eta =(TH1F*) mon.addHistogram( new TH1F("h_fjet2_eta", "fjet2_eta", 100, -5, 5));
  TH1F *h_fjet2_phi =(TH1F*) mon.addHistogram( new TH1F("h_fjet2_phi", "fjet2_phi", 100, -TMath::Pi(), TMath::Pi()));
  //fjet matched pt
  TH1F *h_fj1_pt_mat =(TH1F*) mon.addHistogram( new TH1F("h_fj1_pt_mat", "fj1_pt_mat", 100, 0,500));
  TH1F *h_fj2_pt_mat =(TH1F*) mon.addHistogram( new TH1F("h_fj2_pt_mat", "fj2_pt_mat", 100, 0,500));
  //fjet pt vs pt(bb)
  
  TH2F *h_fj1_pt_bb_pt =(TH2F*) mon.addHistogram( new TH2F("h_fj1_pt_bb_pt", "fj1_pt_bb_pt",100,0,500, 100, 0,500));
  TH2F *h_fj2_pt_bb_pt =(TH2F*) mon.addHistogram( new TH2F("h_fj2_pt_bb_pt", "fj2_pt_bb_pt",100,0,500, 100, 0,500));
  //fjet and fjet matched soft drop mass 
  TH1F *fj_sd_mass1=(TH1F*) mon.addHistogram(new TH1F("fj_sd_mass1","fj_sd_mass1",40,0,80));
  TH1F *fj_sd_mass2=(TH1F*) mon.addHistogram(new TH1F("fj_sd_mass2","fj_sd_mass2",40,0,80));
  TH1F *fj_sd_mass1_mat=(TH1F*) mon.addHistogram(new TH1F("fj_sd_mass1_matched","fj_sd_mass1_matched",40,0,80));
  TH1F *fj_sd_mass2_mat=(TH1F*) mon.addHistogram(new TH1F("fj_sd_mass2_matched","fj_sd_mass2_matched",40,0,80));
  TH2F *fj_pt_sd_mass1 =(TH2F*) mon.addHistogram( new TH2F("fj_pt_sd_mass1", "fj_pt_sd_mass1",100,0,500, 40, 0,80));
  TH2F *fj_pt_sd_mass2 =(TH2F*) mon.addHistogram( new TH2F("fj_pt_sd_mass2", "fj_pt_sd_mass2",100,0,500, 40, 0,80));
  
  //fjet discriminants
  TH1F *h_fjet1_btagXbb =(TH1F*) mon.addHistogram( new TH1F("h_fjet1_btagXbb", "fjet1_btagXbb", 50, 0, 1));
  TH1F *h_fjet1_btagXbbXccXqq =(TH1F*) mon.addHistogram( new TH1F("h_fjet1_btagXbbXccXqq", "fjet1_btagXbbXccXqq", 50, 0, 1));
 
  TH1F *h_fjet2_btagXbb =(TH1F*) mon.addHistogram( new TH1F("h_fjet2_btagXbb", "fjet2_btagXbb", 50, 0, 1);
  TH1F *h_fjet2_btagXbbXccXqq =(TH1F*) mon.addHistogram( new TH1F("h_fjet2_btagXbbXccXqq", "fjet2_btagXbbXccXqq", 50, 0,1));

  TH2F *h_fj1_fj2_btagXbb =(TH2F*) mon.addHistogram( new TH2F("h_fj1_vs_fj2_btagXbb", "fjet_btagXbb",50,0,1, 50, 0, 1));
  TH2F *h_fj1_fj2_btagXbbXccXqq =(TH2F*) mon.addHistogram( new TH2F("h_fj1_fj2_btagXbbXccXqq", "fjet_btagXbbXccXqq",50,0,1,50, 0, 1));

  TH2F *h_fjet1_Xbb_pt = (TH2F*) mon.addHistogram(new TH2F("h_fjet1_Xbb_pt", "fjet1_btagXbb", 100,0,500,50,0, 1));
  TH2F *h_fjet1_XbbXccXqq_pt =(TH2F*) mon.addHistogram( new TH2F("h_fjet1_XbbXccXqq_pt", "fjet1_btagXbb", 100,0,500,50,0, 1));
  TH2F *h_fjet2_Xbb_pt =(TH2F*) mon.addHistogram( new TH2F("h_fjet2_Xbb_pt", "fjet1_btagXbb", 100,0,500,50,0, 1));
  TH2F *h_fjet2_XbbXccXqq_pt =(TH2F*) mon.addHistogram( new TH2F("h_fjet2_XbbXccXqq_pt", "fjet1_btagXbb", 100,0,500,50,0, 1));
 
  //fjet subjet count
  TH1F *subcount1 =(TH1F*) mon.addHistogram (new TH1F("sub_mult1", "sub_mult1", 5, 0, 5));
  TH1F *subcount1_mat =(TH1F*) mon.addHistogram( new TH1F("sub_mult1_mat", "sub_mult1", 5, 0, 5));
  TH1F *subcount2 =(TH1F*) mon.addHistogram( new TH1F("sub_mult2", "sub_mult1", 5, 0, 5));
  TH1F *subcount2_mat =(TH1F*) mon.addHistogram( new TH1F("sub_mult2_mat", "sub_mult1", 5, 0, 5));
  // cc-jets
  //cc-jets mult
  TH1F *h_jet_cc_mult1 =(TH1F*) mon.addHistogram( new TH1F("jet_cc_mult1", "jets", 10, 0, 10));
  TH1F *h_jet_cc_mult2 =(TH1F*) mon.addHistogram( new TH1F("jet_cc_mult2", "jets", 10, 0, 10));
  TH1F *h_jet_cc_mult3 =(TH1F*) mon.addHistogram( new TH1F("jet_cc_mult3", "jets", 10, 0, 10));						   
  //cc jets kinematics
  TH1F *h_jet1_pt =(TH1F*) mon.addHistogram( new TH1F("h_jet1_pt", "jet1_pt", 100, 0, 300));
  TH1F *h_jet1_eta =(TH1F*) mon.addHistogram( new TH1F("h_jet1_eta", "jet1_eta", 100, -5, 5));
  TH1F *h_jet1_phi =(TH1F*) mon.addHistogram( new TH1F("h_jet1_phi", "jet1_phi", 100, -TMath::Pi(), TMath::Pi()));
  TH1F *h_jet2_pt =(TH1F*) mon.addHistogram( new TH1F("h_jet2_pt", "jet2_pt", 100, 0, 300));
  TH1F *h_jet2_eta =(TH1F*) mon.addHistogram( new TH1F("h_jet2_eta", "jet2_eta", 100, -5, 5));
  TH1F *h_jet2_phi =(TH1F*) mon.addHistogram( new TH1F("h_jet2_phi", "jet2_phi", 100, -TMath::Pi(), TMath::Pi()));
   //MET
  //total met
  TH1F *h_met_pt=(TH1F*) mon.addHistogram(new TH1F("met_pt","met_PT",100,0,500));
  TH1F *h_met_phi=(TH1F*) mon.addHistogram(new TH1F("met_phi","met_phi",100,-TMath::Pi(), TMath::Pi()));
  //met from vv/bb/qq_light 				   
  TH1F *h_met_qq_pt=(TH1F*) mon.addHistogram(new TH1F("met_qq_l_pt","met_PT",100,0,500));
  TH1F *h_met_bb_pt=(TH1F*) mon.addHistogram(new TH1F("met_bb_pt","met_PT",100,0,500));
  TH1F *h_met_vv_pt=(TH1F*) mon.addHistogram(new TH1F("met_vv_pt","met_PT",100,0,500));
  //met after trigger requirement
  TH1F *h_met_trig_pt=(TH1F*) mon.addHistogram(new TH1F("met_trig_pt","met_PT",100,0,500));
  TH1F *h_met_vv_trig_phi=(TH1F*) mon.addHistogram(new TH1F("met_vv_trig_phi","met_phi",100,-TMath::Pi(), TMath::Pi())); 
  TH1F *h_met_vv_trig_pt=(TH1F*) mon.addHistogram(new TH1F("met_vv_trig_pt","met_PT",100,0,500));
  						   
 //Ht
  TH1F *h_Ht=(TH1F*) mon.addHistogram(new TH1F("Ht","Ht",300,0,1500));
						   
   //generator level histos
 
  //HIGGS
  TH1F *h_ptH=(TH1F*) mon.addHistogram(new TH1F("H_HIGGS_PT","HIGGS_PT",500,0,500)); 
  TH1F *h_etaH=(TH1F*) mon.addHistogram(new TH1F("H_HIGGS_ETA","HIGGS_ETA",100,-5,5)); 
  TH1F *h_phiH=(TH1F*) mon.addHistogram(new TH1F("H_HIGGS_PHI","HIGGS_PHI",100,-TMath::Pi(), TMath::Pi()));
  TH1F *h_mH=(TH1F*) mon.addHistogram(new TH1F("H_HIGGS_M","HIGGS_M",50,100,150)); 
  //Z boson
  TH1F *h_pt_z=(TH1F*) mon.addHistogram(new TH1F("h_z_PT","z_PT",500,0,500)); 
  TH1F *h_eta_z=(TH1F*) mon.addHistogram(new TH1F("H_z_eta","z_eta",100,-5,5)); 
  TH1F *h_phi_z=(TH1F*) mon.addHistogram(new TH1F("H_z_phi","z_phi",100,-TMath::Pi(), TMath::Pi())); 
  TH1F *h_m_z=(TH1F*) mon.addHistogram(new TH1F("H_4b_m","4b_m",100,50,150)); 
  //A
  TH1F *h_pt_a=(TH1F*) mon.addHistogram(new TH1F("h_A_PT","A_PT",500,0,500));  
  TH1F *h_eta_a=(TH1F*) mon.addHistogram(new TH1F("H_A_eta","A_eta",100,-5,5)); 
  TH1F *h_phi_a=(TH1F*) mon.addHistogram(new TH1F("H_A_phi","A_phi",100,-TMath::Pi(), TMath::Pi())); 
  TH1F *h_m_a=(TH1F*) mon.addHistogram(new TH1F("H_A_m","A_m",500,0,500));
  //1st and 2nd A
  TH1F *h_pt_a1=(TH1F*) mon.addHistogram(new TH1F("h_A_PT1","A_PT",500,0,500)); 
  TH1F *h_pt_a2=(TH1F*) mon.addHistogram(new TH1F("h_A_PT2","A_PT",500,0,500));
   // multiplicities before acceptance cuts plots
  TH1F *h_engen_mult_bef =(TH1F*) mon.addHistogram( new TH1F("en_gen_mult_bef", "electrons_mult", 5, 0, 5));
  TH1F *h_mngen_mult_bef =(TH1F*) mon.addHistogram (new TH1F("mn_gen_mult_bef", "muons_mult", 5, 0, 5));
  TH1F *h_leptgen_mult_bef =(TH1F*) mon.addHistogram( new TH1F("lepton_gen_mult_bef", "leptons_mult", 5, 0, 5));
  TH1F *h_qg_mult_bef =(TH1F*) mon.addHistogram( new TH1F("qg_gen_mult_bef", "quark-gluon_mult", 12, 0, 12));
  TH1F *h_bq_mult_bef=(TH1F*) mon.addHistogram( new TH1F("bquark_gen_mult_bef", "bquark_mult", 10, 0, 10));

  // Multiplicities after cuts
  TH1F *h_engen_mult =(TH1F*) mon.addHistogram( new TH1F("en_gen_mult", "electrons_mult", 5, 0, 5));
  TH1F *h_mngen_mult =(TH1F*) mon.addHistogram( new TH1F("mn_gen_mult", "muons_mult", 5, 0, 5));
  TH1F *h_leptgen_mult =(TH1F*) mon.addHistogram( new TH1F("lepton_gen_mult", "leptons_mult", 5, 0, 5));
  TH1F *h_qg_mult =(TH1F*) mon.addHistogram( new TH1F("qg_gen_mult", "quark-gluon_mult", 10, 0, 10));
  TH1F *h_bq_mult=(TH1F*) mon.addHistogram( new TH1F("bquark_gen_mult", "bquark_mult", 10, 0, 10));
 
  //deltar between aas
  TH1F *h_dr_aa =(TH1F*) mon.addHistogram( new TH1F("h_dR_aa", "dR_aa", 100, 0, 8));
  // deltar between bbs of same mom
  TH1F *h_dr_bb1 =(TH1F*) mon.addHistogram( new TH1F("h_dR_bb1", "dR_bb", 100, 0, 8));
  TH1F *h_dr_bb2 =(TH1F*) mon.addHistogram( new TH1F("h_dR_bb2", "dR_bb", 100, 0, 8));
 
  //vector sum of pt between bbs of same mom
  TH1F *h_pt_bb1_vec=(TH1F*) mon.addHistogram(new TH1F("h_bb1_PT_vec","bb_PT_vec",100,0,500));
  //vec sum pt of vvs
  TH1F *h_pt_vv_vec=(TH1F*) mon.addHistogram(new TH1F("h_vv_PT_vec","vv_PT_vec",100,0,500));						   
  
  

  
  //##################################################################################
  //#############         GET READY FOR THE EVENT LOOP           #####################
  //##################################################################################

 //------Event Counters-------------
  Int_t n_event_lepton_test = 0;
  Int_t n_fjets=0;
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

  //####################################################################################################################
  //###########################################           EVENT LOOP         ###########################################
  //####################################################################################################################

  for( int iev=evStart; iev<evEnd; iev++) 
    {

      if ( verbose ) printf("\n\n Event info %3d: \n",iev);
      GetEntry(iev);

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

      //start generator level loop
      for (int imc=0; imc<nmcparticles;imc++)
     {
       if(verbose && iev <3) 
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
              h_ptH->Fill(pH.Pt(),weight); 
              h_etaH->Fill(etaH,weight); h_phiH->Fill(phiH,weight); h_mH->Fill(mH,weight);
              
            }
        if(mc_id[imc]==23) 
            { // found the z boson

              TLorentzVector pz;
	            pz.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
              // Make histograms of the pt, eta, phi, mass of the z boson:
		    h_pt_z->Fill(pz.Pt(),weight); 
		    h_eta_z->Fill(pz.Eta(),weight); h_phi_z->Fill(pz.Phi(),weight); h_m_z->Fill(pz.M(),weight);
              
            }
        if(mc_id[imc]==36) 
            { // found the A

              TLorentzVector pa;
	            pa.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
              // Make histograms of the pt, eta, phi, mass of the A:
              //h_pt_a->Fill(pa.Pt()); 
              h_eta_a->Fill(pa.Eta(),weight); h_phi_a->Fill(pa.Phi(),weight); h_m_a->Fill(pa.M(),weight);
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
	    // acceptance cuts :
            if (pb.Pt()>20. && fabs(pb.Eta())<2.4) 
            {
              vec_genqg.push_back(pb);
              if (abs(mc_id[imc])==5)  vec_genb.push_back(pb);
	      //bb pair from same a boson 
              
              if (abs(mc_id[imc])==5 && mc_mom[imc]==36 && (mc_momidx[imc]==4))  
                { 
                  vec_genbb1_truth.push_back(pb);
                }
              if (abs(mc_id[imc])==5 && mc_mom[imc]==36 && (mc_momidx[imc]==5))  
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
   
    
   
    
    
    
      // find dr between aas
    if (vec_gena.size()>1) h_dr_aa->Fill(getDeltaR(vec_gena[0],vec_gena[1]),weight);

    //sort qluon-quark vector
    vec_genqg= sort_vec_pt( vec_genqg);
    //sort bquark vector
    vec_genb=sort_vec_pt(vec_genb);
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
     h_dr_bb1->Fill(getDeltaR(vec_genbb1_truth[0],vec_genbb1_truth[1]),weight);
     h_pt_bb1_vec->Fill((vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
    }
    if (vec_genbb2_truth.size()>1 && is_vv){
      h_dr_bb2->Fill(getDeltaR(vec_genbb2_truth[0],vec_genbb2_truth[1]),weight);
      h_pt_bb2->Fill((vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
      h_pt_a2->Fill(vec_gena[1].Pt(),weight); 
    }
    //plot vec sum pt of vvs
    if(is_vv){
      h_pt_vv_vec->Fill((vec_zvv[0]+vec_zvv[1]).Pt(),weight);

    } 
     // plot eta and pt of 4 b quarks 
    if (vec_genb.size()>0)
    {
      h_pt_b1->Fill(vec_genb[0].Pt());
      h_eta_b1->Fill(vec_genb[0].Eta());
 
    
    // Make Histograms: lepton mulciplicity, q/g multiplicity , b-quark multiplicity --> Check effect of the Acceptance cuts (pt/eta)
    h_engen_mult_bef->Fill(vec_genele_bef.size(),weight);
    h_mngen_mult_bef->Fill(vec_genmu_bef.size(),weight);
    h_leptgen_mult_bef->Fill(vec_genlep_bef.size(),weight);
    h_qg_mult_bef->Fill(vec_genqg_bef.size(),weight);
    h_bq_mult_bef->Fill(vec_genb_bef.size(),weight);
    //after acceptance cuts
    h_engen_mult->Fill(vec_genele.size(),weight);
    h_mngen_mult->Fill(vec_genmu.size(),weight);
    h_leptgen_mult->Fill(vec_genlep.size(),weight);
    h_qg_mult->Fill(vec_genqg.size(),weight);
    h_bq_mult->Fill(vec_genb.size(),weight);
    
    
    // END GENERATOR LEVEL ANALYSIS
    
    //start reco level analysis
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
        h_dR_jet_mn_before->Fill(dR1,weight);
        if (dR1 < 0.4){
          overlap = true;
         }
         else {
            h_dR_jet_mn_after->Fill(dR1,weight);
         }
      }

      for (int en_count = 0; en_count < vec_ele.size(); en_count++)
      {
        dR1 = getDeltaR(p_jet, vec_ele[en_count]);
        h_dR_jet_en_before->Fill(dR1,weight);
        if (dR1 < 0.4){
         overlap = true;
        }
        else {
         h_dR_jet_en_after->Fill(dR1,weight);
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
         if(fjet_subjet_count[i]<2)continue;
         if (p_fjet.Pt() <30 || std::fabs(p_fjet.Eta()) > 2.5 ) continue;
         
         for (int mn_count = 0; mn_count < vec_muons.size(); mn_count++)
         {
         dR1 = getDeltaR(p_fjet, vec_muons[mn_count]);
         h_dR_fjet_mn_before->Fill(dR1,weight);
         if (dR1 < 0.8){
            overlap = true;
            }
            else {
               h_dR_fjet_mn_after->Fill(dR1,weight);
            }
         }

         for (int en_count = 0; en_count < vec_ele.size(); en_count++)
         {
         dR1 = getDeltaR(p_fjet, vec_ele[en_count]);
         h_dR_fjet_en_before->Fill(dR1,weight);
         if (dR1 < 0.8){
            overlap = true;
         }
         else {
            h_dR_fjet_en_after->Fill(dR1,weight);
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
       // lepton cuts---------------------------------------------------

   
    if (vec_leptons.size() !=0) continue;
    
    n_event_lepton_test++;

    h_met_pt->Fill(met_pt,weight);
    h_met_phi->Fill(met_phi,weight);
    if(is_qq_light)h_met_qq_pt->Fill(met_pt,weight);
    if(is_bb)h_met_bb_pt->Fill(met_pt,weight);
    if(is_vv)h_met_vv_pt->Fill(met_pt,weight);
    bool hasMETtrigger1=(triggerType>>11)&0x1;
    bool hasMETtrigger2=(triggerType>>12)&0x1;
    

    if(hasMETtrigger1||hasMETtrigger2){
    h_met_trig_pt->Fill(met_pt,weight);
    h_met_trig_phi->Fill(met_phi,weight);
    if(is_vv){
    h_met_vv_trig_pt->Fill(met_pt,weight);
    h_met_vv_trig_phi->Fill(met_phi,weight);}
    }
    
    // control plots mult
    h_en_mult->Fill(vec_ele.size(),weight);
    h_mn_mult->Fill(vec_muons.size(),weight);
    //h_lepton_mult->Fill(vec_ele.size() + vec_muons.size());
    h_jet_mult->Fill(vec_jet.size(),weight);
    h_bjet_mult->Fill(vec_bjets.size(),weight);
    h_fjet_mult->Fill(vec_fjet.size(),weight);
    h_jet_cc_mult1->Fill(vec_jet_cc.size(),weight);
   
        
    
    //metcut 
    if (met_pt < 170.) continue;
    METCUT++;
    //INVM++;
    h_fjet_mult_after1->Fill(vec_fjet.size(),weight);
    h_jet_cc_mult2->Fill(vec_jet_cc.size(),weight);
  
        // Fat-jet analysis:
       
      if (vec_fjet.size()>0)
      { 
         //fj1 kinematics
        h_fjet1_pt->Fill(vec_fjet[0].Pt(),weight);
        h_fjet1_eta->Fill(vec_fjet[0].Eta(),weight);
        h_fjet1_phi->Fill(vec_fjet[0].Phi(),weight);
       
        // fj1 sdmass 
        fj_sd_mass1->Fill(fjet_softdropM[index1],weight);
        fj_pt_sd_mass1->Fill(vec_fjet[0].Pt(), fjet_softdropM[index1],weight);
        subcount1->Fill(fjet_subjet_count[index1],weight);
	
        bool matched(false);
        bool matched1(false);bool matched2(false);
	
        if(isSignal) 
	  {
	    if(fjet_matched(vec_fjet[0],vec_genbb1_truth)>0)
	      {
		h_fj1_pt_bb_pt->Fill(vec_fjet[0].Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
		matched1=true;
	      } 
            else if( fjet_matched(vec_fjet[0],vec_genbb2_truth)>0)
	      {
		h_fj1_pt_bb_pt->Fill(vec_fjet[0].Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
		matched2=true;  
	      }
	    if(matched1||matched2)matched=true;
	  }
        else { // backgrounds:
	  
	  if( fjet_matched(vec_fjet[1],vec_genqg)>0){
	    h_fj1_pt_bb_pt->Fill(vec_fjet[0].Pt(),fjet_matched(vec_fjet[0],vec_genqg),weight);
	    matched=true;
	  }
        }
        if(matched){
	  h_fj1_pt_mat->Fill(vec_fjet[0].Pt(),weight);
	  fj_sd_mass1_mat->Fill(fjet_softdropM[index1],weight);
	  subcount1_mat->Fill(fjet_subjet_count[index1],weight);
        }
	
	
	//xbb discriminants
        float bbtag_fj1=fjet_btag10[index1]/(fjet_btag10[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);  
        h_fjet1_btagXbb->Fill(bbtag_fj1,weight);
        float bbccqqtagfj1=(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1])/(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);
        h_fjet1_btagXbbXccXqq->Fill(bbccqqtagfj1,weight);
        h_fjet1_XbbXccXqq_pt->Fill(vec_fjet[0].Pt(),bbccqqtagfj1);
        h_fjet1_Xbb_pt->Fill(vec_fjet[0].Pt(),bbtag_fj1,weight);
	
	
      }
      
      // fjet2:
      if (vec_fjet.size()>1)
	{  
          if(verbose) { // debugging
	    cout << "Fatjet1 pT = " << vec_fjet[0].Pt() << " and fatjet2 pT= " << vec_fjet[1].Pt() << endl;
	  }
          //fj2 kinematics
	  h_fjet2_pt->Fill(vec_fjet[1].Pt(),weight);
	  h_fjet2_eta->Fill(vec_fjet[1].Eta(),weight);
	  h_fjet2_phi->Fill(vec_fjet[1].Phi(),weight);
	  
          //fj2 sdmass
	  fj_sd_mass2->Fill(fjet_softdropM[index2],weight);
	  fj_pt_sd_mass2->Fill(vec_fjet[1].Pt(), fjet_softdropM[index2],weight);
	  subcount2->Fill(fjet_subjet_count[index2],weight);
	  
	  bool matched(false);
	  bool matched1(false); bool matched2(false);
	  
	  if(isSignal) 
	    {
	      if(  fjet_matched(vec_fjet[1],vec_genbb1_truth)>0)
		{
                  h_fj2_pt_bb_pt->Fill(vec_fjet[1].Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
                  matched1=true;
		} 
	      else if(  fjet_matched(vec_fjet[1],vec_genbb2_truth)>0)
		{
                  h_fj2_pt_bb_pt->Fill(vec_fjet[1].Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
                  matched2=true;  
		}
	      if(matched1||matched2)matched=true;
	    }
	  else
	    {
	      if( fjet_matched(vec_fjet[1],vec_genqg)>0){
		h_fj2_pt_bb_pt->Fill(vec_fjet[1].Pt(),fjet_matched(vec_fjet[1],vec_genqg),weight);
		matched=true;
	      }
	    }
	  if(matched){
	    h_fj2_pt_mat->Fill(vec_fjet[1].Pt(),weight);
	    fj_sd_mass2_mat->Fill(fjet_softdropM[index2],weight);
	    subcount2_mat->Fill(fjet_subjet_count[index2],weight);
	  }
	  
	  
	  //fj2 xbb discriminants
	  float bbtag_fj2=fjet_btag10[index2]/(fjet_btag10[index2]+fjet_btag13[index2]+fjet_btag14[index2]+fjet_btag15[index2]+fjet_btag16[index2]+fjet_btag17[index2]);  
	  h_fjet2_btagXbb->Fill(bbtag_fj2,weight);
	  float bbccqqtagfj2=(fjet_btag10[index2]+fjet_btag11[index2]+fjet_btag12[index2])/(fjet_btag10[index2]+fjet_btag11[index2]+fjet_btag12[index2]+fjet_btag13[index2]+fjet_btag14[index2]+fjet_btag15[index2]+fjet_btag16[index2]+fjet_btag17[index2]);
	  h_fjet2_btagXbbXccXqq->Fill(bbccqqtagfj2,weight);
	  float bbtag_fj1=fjet_btag10[index1]/(fjet_btag10[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);  
	  float bbccqqtagfj1=(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1])/(fjet_btag10[index1]+fjet_btag11[index1]+fjet_btag12[index1]+fjet_btag13[index1]+fjet_btag14[index1]+fjet_btag15[index1]+fjet_btag16[index1]+fjet_btag17[index1]);
	  h_fj1_fj2_btagXbb->Fill(bbtag_fj1,bbtag_fj2,weight);
	  h_fj1_fj2_btagXbbXccXqq->Fill(bbccqqtagfj1,bbccqqtagfj2,weight);
	  h_fjet2_XbbXccXqq_pt->Fill(vec_fjet[1].Pt(),bbccqqtagfj2,weight);
	  h_fjet2_Xbb_pt->Fill(vec_fjet[1].Pt(),bbtag_fj2,weight);
	}
      
	    // AT LEAST 2 FAT JETS
      if (vec_fjet.size() < 2) continue;
      n_fjets++;
      h_fjet_mult_after2->Fill(vec_fjet.size(),weight);
      h_jet_mult_after->Fill(vec_jet.size(),weight);
      h_jet_cc_mult3->Fill(vec_jet_cc.size(),weight);
      h_Ht->Fill(Ht(vec_jet_cc,vec_fjet),weight);
      // h_met_pt->Fill(met_pt);
      //h_met_phi->Fill(met_phi);
      
      //z-?
      if(is_vv) st4_vv++;
      if(is_qq_light) st4_qq++;
      if(is_bb) st4_bb++;
      if(is_ll) st4_ll++;
      //fill jet histos
      if (vec_jet_cc.size() > 0)
	{
	  h_jet1_pt->Fill(vec_jet_cc[0].Pt(),weight);
	  h_jet1_eta->Fill(vec_jet_cc[0].Eta(),weight);
	  h_jet1_phi->Fill(vec_jet_cc[0].Phi(),weight);
	}
      if (vec_jet_cc.size() > 1)
	{
	  h_jet2_pt->Fill(vec_jet_cc[1].Pt(),weight);
	  h_jet2_eta->Fill(vec_jet_cc[1].Eta(),weight);
	  h_jet2_phi->Fill(vec_jet_cc[1].Phi(),weight);
	}
    } // end event loop    

      // in the end: save all to the file
  int nTrial = 0;
  TFile *ofile=TFile::Open(outUrl, "recreate");
  while( !ofile->IsOpen() || ofile->IsZombie() ){
    if(nTrial > 3){
      printf("Output file open failed!");
      if ( outTxtFile_final ) fclose(outTxtFile_final);
      return -1;
    }
    nTrial++;
    usleep(1000000*nTrial);
    ofile=TFile::Open(outUrl, "update");
  }
  mon.Write();
  ofile->Close();
  
  
} // end main



//functions
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

float fjet_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb)
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

    
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec)
{
if (vec.size()>0)
    {
      for (int i=0;i<vec.size()-1;i++)
      {
          int maxIndex = i;  
          for (int j = i + 1; j < vec.size(); j++) 
            {  
              if (vec[j].Pt() > vec[maxIndex].Pt())
          {  
            maxIndex = j;  
          }  
            }   
          TLorentzVector temp;
          temp = vec[i];  
          vec[i] = vec[maxIndex];  
          vec[maxIndex] = temp; 
      }
    }
  return vec;
}
