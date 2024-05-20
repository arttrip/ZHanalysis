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
  evSummary_.btag3=-1;
  evSummary_.dilep_pt=-1;
  evSummary_.drll=-1;
  evSummary_.dphiHZ=-1;
  evSummary_.weight=0;
 
  std::cout << " check 2!" << std::endl;

}

//
void MVAHandler::getEntry(
			  bool isSR1,
			  bool isSR2,
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
			  float n_ad_j,
			  float btag3,
			  float dilep_pt,
			  float drll,
			  float dphiHZ,
			  float weight) 
 
{
  
  resetStruct();
  evSummary_.isSR1=isSR1;
  evSummary_.isSR2=isSR2;
  evSummary_.m4b=m4b;
  evSummary_.pt4b=pt4b;
  evSummary_.ptf1=ptf1;
  evSummary_.ptf2=ptf2;
  evSummary_.ptb1=ptb1;
  evSummary_.ht=ht;
  evSummary_.n_ad_j=n_ad_j;
  evSummary_.met=met;
  evSummary_.btag3=btag3;
  evSummary_.sd_mass1=sd_mass1;
  evSummary_.sd_mass2=sd_mass2;
  evSummary_.xbb1=xbb1;
  evSummary_.xbb2=xbb2;
  evSummary_.xbbccqq1=xbbccqq1;
  evSummary_.xbbccqq2=xbbccqq2;

  evSummary_.dilep_pt=dilep_pt;
  evSummary_.drll=drll;
  evSummary_.dphiHZ=dphiHZ;

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
  
  //2lepton
  if(run0lep){
    t1->Branch("dilep_pt",&evSummary_.dilep_pt,"dilep_pt/F");
    t1->Branch("drll",&evSummary_.drll,"drll/F");
    t1->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
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
  
   
  //2f jets                                                                                                                                                                                                    
  t2->Branch("ptf2",&evSummary_.ptf2,"ptf2/F");
  t2->Branch("sd_mass1",&evSummary_.sd_mass1,"sd_mass1/F");
  t2->Branch("sd_mass2",&evSummary_.sd_mass2,"sd_mass2/F");
  t2->Branch("xbb1",&evSummary_.xbb1,"xbb1/F");
  t2->Branch("xbb2",&evSummary_.xbb2,"xbb2/F");
  t2->Branch("xbbccqq1",&evSummary_.xbbccqq1,"xbbccqq1/F");
  t2->Branch("xbbccqq2",&evSummary_.xbbccqq2,"xbbccqq2/F");
  //2lepton                                                                                                                                                                                                     
  if(run0lep){
    t2->Branch("dilep_pt",&evSummary_.dilep_pt,"dilep_pt/F");
    t2->Branch("drll",&evSummary_.drll,"drll/F");
    t2->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
  }

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
  if ( evSummary_.isSR1 && evSummary_.isSR2 )
    {
      std::cout << "One event can not be both in 3 and 4 b cat! Please check!" << std::endl;
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
  MVAofile->Close();
  std::cout << "write tree check!" << std::endl;
  //fprintf("write tree\n");
  return ;
}
//
MVAHandler::~MVAHandler()
{
}
