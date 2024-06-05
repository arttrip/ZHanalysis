#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"
bool run0lep(false);
MVAHandler::MVAHandler()
{
}

//
MVAEvtContainer &MVAHandler::getEvent()
{
  std::cout << " event check " << std::endl;
  return evSummary_;
}

//
void MVAHandler::resetStruct()
{
  evSummary_.isSR1=false;
  evSummary_.isSR2=false;
  evSummary_.is2b=false;
  evSummary_.is3b=false;
  evSummary_.m4b=-1;
  evSummary_.pt4b=-1;
  evSummary_.ptf1=-1;
  evSummary_.sd_mass1=-1;
  evSummary_.ht=-1;
  evSummary_.met=-1;
  evSummary_.xbb1=-2;
  evSummary_.xbbccqq1=-2;
  evSummary_.ptf2=-1;
  evSummary_.sd_mass2=-1;
  evSummary_.xbb2=-2;                                \
  evSummary_.xbbccqq2=-2;
  evSummary_.ptb1=-1;
  evSummary_.n_ad_j=-1;
  evSummary_.btag1=-1;
  evSummary_.btag3=-1;
  evSummary_.drjj=-1;
  evSummary_.dphi_met_j=-1;
  evSummary_.dilep_pt=-1;
  evSummary_.drll=-1;
  evSummary_.dphiHZ=-1;
  evSummary_.dphi_met_l=-1;
  evSummary_.weight=0;
 
  std::cout << " check 2!" << std::endl;

}

//
void MVAHandler::getEntry(
			  bool isSR1,
			  bool isSR2,
			  bool is2b,
			  bool is3b,
			  float m4b,
			  float pt4b,
			  float ptf1,
			  float sd_mass1,
			  float ht,
			  float met,
			  float xbb1,
			  float xbbccqq1,
			  float ptf2,
			  float sd_mass2,
			  float xbb2,		
			  float xbbccqq2,
			  float ptb1,
			  float ptb2,
			  float n_ad_j,
			  float btag1,
			  float btag3,
			  float drjj,
			  float dphi_met_j,
			  float dilep_pt,
			  float drll,
			  float dphiHZ,
			  float dphi_met_l,
			  float weight) 
 
{
  
  resetStruct();
  evSummary_.isSR1=isSR1;
  evSummary_.isSR2=isSR2;
  evSummary_.is2b=is2b;
  evSummary_.is3b=is3b;
  evSummary_.m4b=m4b;
  evSummary_.pt4b=pt4b;
  evSummary_.ptf1=ptf1;
  evSummary_.ptf2=ptf2;
  evSummary_.ptb1=ptb1;
  evSummary_.ptb2=ptb2;
  evSummary_.ht=ht;
  evSummary_.n_ad_j=n_ad_j;
  evSummary_.met=met;
  evSummary_.btag1=btag1;
  evSummary_.btag3=btag3;
  evSummary_.sd_mass1=sd_mass1;
  evSummary_.sd_mass2=sd_mass2;
  evSummary_.xbb1=xbb1;
  evSummary_.xbb2=xbb2;
  evSummary_.xbbccqq1=xbbccqq1;
  evSummary_.xbbccqq2=xbbccqq2;
  evSummary_.drjj=drjj;
  evSummary_.dphi_met_j=dphi_met_j;
  evSummary_.dilep_pt=dilep_pt;
  evSummary_.drll=drll;
  evSummary_.dphiHZ=dphiHZ;
  evSummary_.dphi_met_l=dphi_met_l;

  //weight
  evSummary_.weight = weight;
  std::cout << " check 3 !" << std::endl;
  return ;
}

//
bool MVAHandler::initTree()
{
  
  t1 = new TTree("t1","trMVA");
  t1->SetDirectory(0);
  std::cout << "check init tree"<<std::endl;
  
  t1->Branch("weight",&evSummary_.weight,"weight/F");
  t1->Branch("m4b",&evSummary_.m4b,"m4b/F");
  t1->Branch("pt4b",&evSummary_.pt4b,"pt4b/F");
  t1->Branch("ht",&evSummary_.ht,"ht/F");
  t1->Branch("met",&evSummary_.met,"met/F");
  t1->Branch("ptf1",&evSummary_.ptf1,"ptf1/F");
  t1->Branch("n_ad_j",&evSummary_.n_ad_j,"n_ad_j/F");  //1f+1b jet
  t1->Branch("ptb1",&evSummary_.ptb1,"ptb1/F");
  t1->Branch("btag3",&evSummary_.btag3,"btag3/F");
   
  t1->Branch("sd_mass1",&evSummary_.sd_mass1,"sd_mass1/F");
  
  t1->Branch("xbb1",&evSummary_.xbb1,"xbb1/F");
  
  t1->Branch("xbbccqq1",&evSummary_.xbbccqq1,"xbbccqq1/F");
  t1->Branch("drjj",&evSummary_.drjj,"drjj/F");
  t1->Branch("dphi_met_j",&evSummary_.dphi_met_j,"dphi_met_j/F");
  //2lepton
  if(!run0lep){
    t1->Branch("dilep_pt",&evSummary_.dilep_pt,"dilep_pt/F");
    t1->Branch("drll",&evSummary_.drll,"drll/F");
    t1->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
    t1->Branch("dphi_met_l",&evSummary_.dphi_met_l,"dphi_met_l/F");
 }
 
  
 
  t2 = new TTree("t2","trMVA2");
  t2->SetDirectory(0);
  std::cout << "check init tree2"<<std::endl;

  t2->Branch("weight",&evSummary_.weight,"weight/F");
  t2->Branch("m4b",&evSummary_.m4b,"m4b/F");
  t2->Branch("pt4b",&evSummary_.pt4b,"pt4b/F");
  t2->Branch("ht",&evSummary_.ht,"ht/F");
  t2->Branch("met",&evSummary_.met,"met/F");
  t2->Branch("ptf1",&evSummary_.ptf1,"ptf1/F");
  
  t2->Branch("drjj",&evSummary_.drjj,"drjj/F");
  t2->Branch("dphi_met_j",&evSummary_.dphi_met_j,"dphi_met_j/F");
  //2f jets                                                                                                                                                                                                    
  t2->Branch("ptf2",&evSummary_.ptf2,"ptf2/F");
  t2->Branch("sd_mass1",&evSummary_.sd_mass1,"sd_mass1/F");
  t2->Branch("sd_mass2",&evSummary_.sd_mass2,"sd_mass2/F");
  t2->Branch("xbb1",&evSummary_.xbb1,"xbb1/F");
  t2->Branch("xbb2",&evSummary_.xbb2,"xbb2/F");
  t2->Branch("xbbccqq1",&evSummary_.xbbccqq1,"xbbccqq1/F");
  t2->Branch("xbbccqq2",&evSummary_.xbbccqq2,"xbbccqq2/F");
  //2lepton                                                                                                                                                                                                     
  if(!run0lep){
    t2->Branch("dilep_pt",&evSummary_.dilep_pt,"dilep_pt/F");
    t2->Branch("drll",&evSummary_.drll,"drll/F");
    t2->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
    t2->Branch("dphi_met_l",&evSummary_.dphi_met_l,"dphi_met_l/F");
  }

  t3 = new TTree("t3","trMVA3");
  t3->SetDirectory(0);
  std::cout << "check init tree3"<<std::endl;

  t3->Branch("weight",&evSummary_.weight,"weight/F");
  t3->Branch("m4b",&evSummary_.m4b,"m4b/F");
  t3->Branch("pt4b",&evSummary_.pt4b,"pt4b/F");
  t3->Branch("ht",&evSummary_.ht,"ht/F");
  t3->Branch("met",&evSummary_.met,"met/F");
  t3->Branch("n_ad_j",&evSummary_.n_ad_j,"n_ad_j/F");                          
  t3->Branch("ptb1",&evSummary_.ptb1,"ptb1/F");
  t3->Branch("ptb2",&evSummary_.ptb2,"ptb2/F");
  t3->Branch("btag1",&evSummary_.btag1,"btag1/F");
  t3->Branch("btag3",&evSummary_.btag3,"btag3/F");
  t3->Branch("drjj",&evSummary_.drjj,"drjj/F");
  t3->Branch("dphi_met_j",&evSummary_.dphi_met_j,"dphi_met_j/F");
  if(!run0lep){
    t3->Branch("dilep_pt",&evSummary_.dilep_pt,"dilep_pt/F");
    t3->Branch("drll",&evSummary_.drll,"drll/F");
    t3->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
    t3->Branch("dphi_met_l",&evSummary_.dphi_met_l,"dphi_met_l/F");
  }
  t4 = new TTree("t4","trMVA4");
  t4->SetDirectory(0);
  std::cout << "check init tree4"<<std::endl;

  t4->Branch("weight",&evSummary_.weight,"weight/F");
  t4->Branch("m4b",&evSummary_.m4b,"m4b/F");
  t4->Branch("pt4b",&evSummary_.pt4b,"pt4b/F");
  t4->Branch("ht",&evSummary_.ht,"ht/F");
  t4->Branch("met",&evSummary_.met,"met/F");
  t4->Branch("n_ad_j",&evSummary_.n_ad_j,"n_ad_j/F");
  t4->Branch("ptb1",&evSummary_.ptb1,"ptb1/F");
  t4->Branch("ptb2",&evSummary_.ptb2,"ptb2/F");
  t4->Branch("btag1",&evSummary_.btag1,"btag1/F");
  t4->Branch("btag3",&evSummary_.btag3,"btag3/F");
  t4->Branch("drjj",&evSummary_.drjj,"drjj/F");
  t4->Branch("dphi_met_j",&evSummary_.dphi_met_j,"dphi_met_j/F");
  if(!run0lep){
    t4->Branch("dilep_pt",&evSummary_.dilep_pt,"dilep_pt/F");
    t4->Branch("drll",&evSummary_.drll,"drll/F");
    t4->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
    t4->Branch("dphi_met_l",&evSummary_.dphi_met_l,"dphi_met_l/F");
  }

  //2f jets                                                                                                                               
  std::cout << "branches check!" << std::endl;
  return true;
}

//
void MVAHandler::fillTree()
{

  if(evSummary_.isSR1 && !evSummary_.isSR2)
    {
      // t1->SetDirectory(0);
      t1->Fill();}
  if(evSummary_.isSR2 && !evSummary_.isSR1)
    {
      //t2->SetDirectory(0);
      t2->Fill();}
  if ( evSummary_.is2b  )
    {
      t3->Fill();
    }
  if ( evSummary_.is3b  )
    {
      t4->Fill();
    }
  if ( evSummary_.is3b &&  evSummary_.is2b )
    {
      std::cout<<"ouuuups"<<std::endl;
    }

      std::cout << "tree fill check!" << std::endl;
      
      return;
       
}

void MVAHandler::writeTree(TString mvaout)
{
  TFile *MVAofile=TFile::Open( mvaout, "recreate");  
  t1->SetDirectory(0);
  t1->Write();
  t2->SetDirectory(0);
  t2->Write();
  t3->SetDirectory(0);
  t3->Write();
  t4->SetDirectory(0);
  t4->Write();
  MVAofile->Close();
  std::cout << "write tree check!" << std::endl;
  //fprintf("write tree\n");
  return ;
}
//
MVAHandler::~MVAHandler()
{
}
