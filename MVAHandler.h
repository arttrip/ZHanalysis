#ifndef mvahandler_h
#define mvahandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <vector>

#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"

#endif

struct MVAEvtContainer
{
  //catagory type
  bool isSR1=false;
  bool isSR2=false;
  float m4b=-1;
  float pt4b=-1;
  float ptf1=-1;
  float sd_mass1=-1;
  float ht=-1;
  float met=-1;
  float xbb1=-2;
  float xbbccqq1=-2;
  float ptf2=-1;
  float sd_mass2=-1;
  float xbb2=-2;				\
  float xbbccqq2=-2;
  float ptb1=-1;
  float n_ad_j=-1;
  float btag3=-1;
  float dilep_pt=-1;
  float drll=-1;
  float dphiHZ=-1;
  float weight=0;
};

class MVAHandler 
{
 public:
  
  MVAHandler();
  ~MVAHandler();

  //current event
  MVAEvtContainer evSummary_;
  MVAEvtContainer &getEvent();

  //read mode, from calculated var
  void resetStruct();
  void getEntry(bool isSR1,
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
		float xbb2,			\
		float xbbccqq2,
		float ptb1,
		float n_ad_j,
		float btag3,
		float dilep_pt,
		float drll,
		float dphiHZ,
	        float weight);
		// float ptf1,float ht,float met,float sd_mass1,
		//float xbb1,float xbbccqq1,
  //float weight); //,float ptl1,float ptl2,float drll,float dphiHZ
              

  //write mode, to mva tree
  TFile* MVAofile;
  //the tree, 2 for 3b 4b separately
  TTree *t1;
  TTree *t2;
  bool initTree();
  void fillTree();
  void writeTree(TString mvaout);
 private:
};
#endif
