#include "UserCode/bsmhiggs_fwk/interface/TMVAReader.h"

void TMVAReader::InitTMVAReader()
{
  myreader = new TMVA::Reader( "!Color:!Silent" );
  return ;
}


void TMVAReader::SetupMVAReader( std::string methodName, std::string modelPath )
{
  if ( methodName.find("SR1") != std::string::npos )
  {
     if ( modelPath.find("2lepton") != std::string::npos)
      {
	myreader->AddVariable( "dilep_pt", &dilep_pt );
	myreader->AddVariable( "drll", &drll );
	myreader->AddVariable( "dphiHZ", &dphiHZ );
	myreader->AddVariable( "dphi_met_l", &dphi_met_l );
	myreader->AddVariable( "dphi_met_j", &dphi_met_j );
      }
    myreader->AddVariable( "m4b"    , &m4b );
    myreader->AddVariable( "pt4b"  , &pt4b );
    myreader->AddVariable( "ht", &ht );
    myreader->AddVariable( "met"    ,&met  );
    myreader->AddVariable( "ptf1", &ptf1 );  
    myreader->AddVariable( "sd_mass1"   , &sd_mass1 );
    myreader->AddVariable( "xbb1"   , &xbb1);              
    myreader->AddVariable( "xbbccqq1"   , &xbbccqq1 );
    myreader->AddVariable( "n_ad_j", &n_ad_j );
    myreader->AddVariable( "ptb1"    , &ptb1 );
    myreader->AddVariable( "btag3", &btag3 );
    myreader->AddVariable( "drjj", &drjj );
   
    
  }
  if (( methodName.find("SR2") != std::string::npos )||( methodName.find("SR3") != std::string::npos) )
    {
      if ( methodName.find("SR2") != std::string::npos )
	{
	  if ( modelPath.find("2lepton") != std::string::npos)
	    {
	      myreader->AddVariable( "dilep_pt", &dilep_pt );
	      myreader->AddVariable( "drll", &drll );
	      myreader->AddVariable( "dphiHZ", &dphiHZ );
	      myreader->AddVariable( "dphi_met_l", &dphi_met_l );
	      myreader->AddVariable( "dphi_met_j", &dphi_met_j );
	    }
	  myreader->AddVariable( "m4b"    , &m4b );
	  myreader->AddVariable( "pt4b"  , &pt4b );
	  myreader->AddVariable( "ht", &ht );
	  myreader->AddVariable( "met"    ,&met );
	  myreader->AddVariable( "drjj", &drjj );
	  myreader->AddVariable( "n_ad_j", &n_ad_j );
	  myreader->AddVariable( "ptb1"    , &ptb1 );
	  myreader->AddVariable( "ptb2"    , &ptb2 );
	  myreader->AddVariable( "btag1", &btag1 );
	  myreader->AddVariable( "btag3", &btag3 );
     
	}
      else  if ( methodName.find("SR3") != std::string::npos )
	{
	  if ( modelPath.find("2lepton") != std::string::npos)
	    {
	      myreader->AddVariable( "dilep_pt", &dilep_pt );
	      myreader->AddVariable( "drll", &drll );
	      myreader->AddVariable( "dphiHZ", &dphiHZ );
	      myreader->AddVariable( "dphi_met_l", &dphi_met_l );
	      myreader->AddVariable( "dphi_met_j", &dphi_met_j );
	    }
	  myreader->AddVariable( "m4b"    , &m4b );
	  myreader->AddVariable( "pt4b"  , &pt4b );
	  myreader->AddVariable( "ht", &ht );
	  myreader->AddVariable( "met"    ,&met );
	  myreader->AddVariable( "n_ad_j", &n_ad_j );
	  myreader->AddVariable( "ptb1"    , &ptb1 );
	  myreader->AddVariable( "ptb2"    , &ptb2 );
	  myreader->AddVariable( "btag1", &btag1 );
	  myreader->AddVariable( "btag3", &btag3 );
	  
	}
    }
      else
	{
	  return ;
	}
      
      myreader->BookMVA( methodName.c_str(), modelPath.c_str() );
      return ;
}

float TMVAReader::GenReMVAReader(
				  float thisdilep_pt, float thisdrll,float thisdhiHZ,float thisdphi_met_l,  float thisdphi_met_j,
				  float thism4b,float thispt4b,float thismet,  float thisht, 
				  float thisptf1,float thissd_mass1,float thisxbb1, float thisxbbccqq1 ,
				  float thisdrjj,float thisn_ad_j,float thisptb1,float thisptb2,float thisbtag1,float thisbtag3,
				  
				  std::string methodName
                                )
{
  
  m4b=thism4b; pt4b=thispt4b;met=thismet; ht=thisht;
 ptb1=thisptb1;ptb2=thisptb2;  n_ad_j=thisn_ad_j;drjj=thisdrjj;
ptf1=thisptf1; sd_mass1=thissd_mass1;xbb1=thisxbb1; xbbccqq1=thisxbbccqq1;
 btag1=thisbtag1; btag3=thisbtag3;
dilep_pt=thisdilep_pt; drll=thisdrll;dphiHZ= thisdhiHZ;dphi_met_j=thisdphi_met_j; dphi_met_l=thisdphi_met_l;
  
  float mvaOut = myreader->EvaluateMVA( methodName.c_str() );
  return mvaOut;
  
}

void TMVAReader::CloseMVAReader()
{
  delete myreader;
  return ;
}
