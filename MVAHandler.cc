#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"

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
  // evSummary_.run0lep=false;
  //  evSummary_.m4b=-1;
  // evSummary_.pt4b=-1;
  evSummary_.ptf1=-1;
  //  evSummary_.ptf2=-1;
  // evSummary_.ptb1=-1;
  // evSummary_.ptb2=-1;
  evSummary_.ht=-1;
  // evSummary_.n_ad_j=-1;
  evSummary_.met=-1;
  // evSummary_.btag3=-1;
  evSummary_.sd_mass1=-1;
  // evSummary_.sd_mass2=-1;
  evSummary_.xbb1=-1;
  // evSummary_.xbb2=-1;
  evSummary_.xbbccqq1=-1;
  // evSummary_.xbbccqq2=-1; 
  evSummary_.weight=0;
  // evSummary_.ptl1=-1  // evSummary_.ptl2=-1;
  // evSummary_.drll=-1;
  // evSummary_.dphiHZ=-1;
  std::cout << " check 2!" << std::endl;

}

//
void MVAHandler::getEntry(float ptf1,float ht,float met,float sd_mass1, float xbb1,float xbbccqq1,  float weight) //,float ptl1,float ptl2,float drll,float dphiHZ
 
{
  
  resetStruct();
  //  evSummary_.m4b=m4b;
  // evSummary_.pt4b=pt4b;
  evSummary_.ptf1=ptf1;
  //  evSummary_.ptf2=ptf2;
  //evSummary_.ptb1=ptb1;
  // evSummary_.ptb2=ptb2;
  evSummary_.ht=ht;
  // evSummary_.n_ad_j=n_ad_j;
  evSummary_.met=met;
  // evSummary_.btag3=btag3;
  evSummary_.sd_mass1=sd_mass1;
  // evSummary_.sd_mass2=sd_mass2;
  evSummary_.xbb1=xbb1;
  // evSummary_.xbb2=xbb2;
  evSummary_.xbbccqq1=xbbccqq1;
  // evSummary_.xbbccqq2=xbbccqq2;

  // evSummary_.ptl1=ptl1;
  // evSummary_.ptl2=ptl2;
  // evSummary_.drll=drll;
  // evSummary_.dphiHZ=dphiHZ;

  //weight
  evSummary_.weight = weight;
  std::cout << " check 3 !" << std::endl;
  return ;
}

//
bool MVAHandler::initTree(TString mvaout)
{
  std::cout << " check intree" << std::endl;
  //write mode, to mva tree
  TFile * MVAofile =new TFile( mvaout,"RECREATE");
  MVAofile->cd();
  // MVAofile->mkdir("dir");
  std::cout << " check open file" << std::endl;
  TTree* t1 = new TTree("t1","trMVA");
  std::cout << "check init tree"<<std::endl;
  
  t1->Branch("weight",&evSummary_.weight,"weight/F");
  //  t1->Branch("m4b",&evSummary_.m4b,"m4b/F");
  // t1->Branch("pt4b",&evSummary_.pt4b,"pt4b/F");
  t1->Branch("ht",&evSummary_.ht,"ht/F");
  t1->Branch("met",&evSummary_.met,"met/F");
  t1->Branch("ptf1",&evSummary_.ptf1,"ptf1/F");
  //  t1->Branch("n_ad_j",&evSummary_.n_ad_j,"n_ad_j/F");
  //1f+1b jet
  // t1->Branch("ptb1",&evSummary_.ptb1,"ptb1/F");
  //  t1->Branch("ptb2",&evSummary_.ptb2,"ptb2/F");
  // t1->Branch("btag3",&evSummary_.btag3,"btag3/F");
  //2f jets
  // t1->Branch("ptf2",&evSummary_.ptf2,"ptf2/F");
  t1->Branch("sd_mass1",&evSummary_.sd_mass1,"sd_mass1/F");
  // t1->Branch("sd_mass2",&evSummary_.sd_mass2,"sd_mass2/F");
  t1->Branch("xbb1",&evSummary_.xbb1,"xbb1/F");
  // t1->Branch("xbb2",&evSummary_.xbb2,"xbb2/F");
  t1->Branch("xbbccqq1",&evSummary_.xbbccqq1,"xbbccqq1/F");
  // t1->Branch("xbbccqq2",&evSummary_.xbbccqq2,"xbbccqq2/F");
  //2lepton
  //  if(run0lep){
  // t1->Branch("ptl1",&evSummary_.ptl1,"ptl1/F");
  //  t1->Branch("ptl2",&evSummary_.ptl2,"ptl2/F");
  // t1->Branch("drll",&evSummary_.drll,"drll/F");
  // t1->Branch("dphiHZ",&evSummary_.dphiHZ,"dphiHZ/F");
  // }
  std::cout << "branches check!" << std::endl;
  return true;
}

//
void MVAHandler::fillTree()
{

     
      t1->Fill();

      std::cout << "tree fill check!" << std::endl;
      
      return;
       
    }

void MVAHandler::writeTree()
{
  //TFile *MVAofile=TFile::Open( outURL, "recreate");  
  // toSignal_e->Write(); toSignal_o->Write();
  // to3b_e->Write(); to3b_o->Write();
  // to4b_e->Write(); to4b_o->Write();
  t1->Write();
  MVAofile->Close();
  std::cout << "write tree check!" << std::endl;
  //fprintf("write tree\n");
  return ;
}
//
MVAHandler::~MVAHandler()
{
}
