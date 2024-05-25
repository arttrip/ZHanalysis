#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
using namespace std;
using namespace TMVA;

void mva2(){
  
  TString tag("1");
  TString lep("2");
	      
  TString outf( "TMVA_sr" +tag+"_"+lep+ "lep.root" );
  TFile *outputFile =  TFile::Open(outf, "RECREATE" );
  TString classif("TMVA"+tag+"_"+lep+"lep");
  TMVA::Factory *factory = new TMVA::Factory(classif , outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification");
 
  TFile *sig = TFile::Open("zh"+lep+"lep_all.root");
  TFile * bkg1 = TFile::Open("ZNN100.root");
  TFile * bkg2 = TFile::Open("ZNN200.root");
  TFile * bkg3 = TFile::Open("ZNN400.root");
  TFile * bkg4 = TFile::Open("ZNN600.root");
  TFile * bkg5 = TFile::Open("ZNN800.root");
  TFile * bkg6 = TFile::Open("ZNN1200.root");
  TFile * bkg7 = TFile::Open("ZNN2500.root");

   TFile * b1 = TFile::Open("DY1.root");
  TFile * b2 = TFile::Open("DY2.root");
  TFile * b3 = TFile::Open("DY3.root");
  TFile * b4 = TFile::Open("DY4.root");
  
  TString tr("t"+tag);
  TTree* sig_t = (TTree*) sig->Get( tr );
  TTree* bkg1_t = (TTree*) bkg1->Get( tr );
  TTree* bkg2_t = (TTree*) bkg2->Get( tr );
  TTree* bkg3_t = (TTree*) bkg3->Get( tr );
  TTree* bkg4_t = (TTree*) bkg4->Get( tr );      
  TTree* bkg5_t = (TTree*) bkg5->Get( tr );
  TTree* bkg6_t = (TTree*) bkg6->Get( tr );
  TTree* bkg7_t = (TTree*) bkg7->Get( tr );
 
  TTree* b1_t = (TTree*) b1->Get( tr );
  TTree* b2_t = (TTree*) b2->Get( tr );
  TTree* b3_t = (TTree*) b3->Get( tr );
  TTree* b4_t = (TTree*) b4->Get( tr );
    
 TMVA::DataLoader *dataloader = new TMVA::DataLoader("MVAnalysis");

 float sigweight(1.0);
 float bgkweight1(1.0);
 float bgkweight2(1.0);
 float bgkweight3(1.0);
 float bgkweight4(1.0);
 float bgkweight5(1.0);
 float bgkweight6(1.0);
 float bgkweight7(1.0);

 float bweight1(1.0);
 float bweight2(1.0);
 float bweight3(1.0);
 float bweight4(1.0);

 TLeaf *bg1w = bkg1_t->GetLeaf("weight");
 bg1w->GetBranch()->GetEntry(1);
 bgkweight1 = bg1w->GetValue();
 TLeaf *bg2w = bkg2_t->GetLeaf("weight");
 bg2w->GetBranch()->GetEntry(1);
 bgkweight2 = bg2w->GetValue();
 TLeaf *bg3w = bkg3_t->GetLeaf("weight");
 bg3w->GetBranch()->GetEntry(1);
 bgkweight3 = bg3w->GetValue();
 TLeaf *bg4w = bkg4_t->GetLeaf("weight");
 bg4w->GetBranch()->GetEntry(1);
 bgkweight4 = bg4w->GetValue();
 TLeaf *bg5w = bkg5_t->GetLeaf("weight");
 bg5w->GetBranch()->GetEntry(1);
 bgkweight5 = bg5w->GetValue();
 TLeaf *bg6w = bkg6_t->GetLeaf("weight");
 bg6w->GetBranch()->GetEntry(1);
 bgkweight6 = bg6w->GetValue();
 TLeaf *bg7w = bkg7_t->GetLeaf("weight");
 bg7w->GetBranch()->GetEntry(1);
 bgkweight7 = bg7w->GetValue();

TLeaf *b1w = b1_t->GetLeaf("weight");
 b1w->GetBranch()->GetEntry(1);
 bweight1 = b1w->GetValue();
 TLeaf *b2w = b2_t->GetLeaf("weight");
 b2w->GetBranch()->GetEntry(1);
 bweight2 = b2w->GetValue();
 TLeaf *b3w = b3_t->GetLeaf("weight");
 b3w->GetBranch()->GetEntry(1);
 bweight3 = b3w->GetValue();
 TLeaf *b4w = b4_t->GetLeaf("weight");
 b4w->GetBranch()->GetEntry(1);
 bweight4 = b4w->GetValue();

 TCut preselectionCut1;
 TCut preselectionCut2;
dataloader->AddSignalTree(sig_t, 1.0 );
 if(lep=="0")
   {
     dataloader->AddBackgroundTree(bkg1_t,bgkweight1);
     dataloader->AddBackgroundTree(bkg2_t,bgkweight2);
     dataloader->AddBackgroundTree(bkg3_t,bgkweight3);
     dataloader->AddBackgroundTree(bkg4_t,bgkweight4);
     dataloader->AddBackgroundTree(bkg5_t,bgkweight5);
     dataloader->AddBackgroundTree(bkg6_t,bgkweight6);
     dataloader->AddBackgroundTree(bkg7_t,bgkweight7);
   }
 else if(lep=="2")
   {
     dataloader->AddBackgroundTree(b1_t,bweight1);
     dataloader->AddBackgroundTree(b2_t,bweight2);
     dataloader->AddBackgroundTree(b3_t,bweight3);
     dataloader->AddBackgroundTree(b4_t,bweight4);

      dataloader->AddVariable("dilep_pt", 'F');
      dataloader->AddVariable("drll", 'F');
      dataloader->AddVariable("dphiHZ", 'F');
   }
  dataloader->AddVariable("m4b", 'F');
    dataloader->AddVariable("pt4b", 'F');
    dataloader->AddVariable("ht", 'F');
    dataloader->AddVariable("met", 'F');
    dataloader->AddVariable("ptf1", 'F');
    dataloader->AddVariable("sd_mass1", 'F');
    dataloader->AddVariable("xbb1", 'F');
    dataloader->AddVariable("xbbccqq1", 'F');
    if(tag=="1")
      {
	dataloader->AddVariable("n_ad_j", 'F');
	dataloader->AddVariable("ptb1", 'F');
	dataloader->AddVariable("btag3", 'F');
	 preselectionCut1 = "!(sd_mass1<0)";
	 preselectionCut2 ="";
      }
    else if(tag=="2")
      {
	 dataloader->AddVariable("ptf2", 'F');
	 dataloader->AddVariable("sd_mass2", 'F');
	 dataloader->AddVariable("xbb2", 'F');
	 dataloader->AddVariable("xbbccqq2", 'F');
	  preselectionCut1 = "!(sd_mass1<0)";
	  preselectionCut2 = "!(sd_mass2<0)";
      }
    	
 

   
    dataloader->PrepareTrainingAndTestTree(preselectionCut1 &&  preselectionCut2,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT","!H:!V:NTrees=100:MinNodeSize=5.0%:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   
    //factory->BookMethod(dataloader,TMVA::Types::kMLP, "MLP", "!V:NCycles=200:HiddenLayers=N+1,N:TestRate=5" );
    //factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
    //  factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD","!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
    //factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
    
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods(); 

    outputFile->Close();
    delete factory;

    TMVA::TMVAGui(outf);
    
  }
