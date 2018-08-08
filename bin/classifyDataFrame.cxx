/////////////////////////////////////////////////////////////////////////////
///  Simplified Higgs vs. background classification                       ///
///                                                                       ///
///  This runs on flat ntuples with the selections already made.          ///
//   classify.cxx runs on our large ttrees making selections then adding  ///
//   the events one by one.                                               ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/TMVAMultiClassGui.h"

void trainAndTest()
{
    // load TMVA 
    TMVA::Tools::Instance();
    auto outputFile = TFile::Open("TMVAOutputCV.root", "RECREATE");

    ///////////////////////////////////////////////////////////////////////
    // The Training/Testing samples
    ///////////////////////////////////////////////////////////////////////

    TString filedir("/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt/");
   
    // signal ntuples 
    std::vector<TString> sigfilenames  = {
                                       TString("H2Mu_VBF_bdt_training_dyAMC-J.root"),
                                       TString("H2Mu_WH_neg_bdt_training_dyAMC-J.root"),
                                       TString("H2Mu_WH_pos_bdt_training_dyAMC-J.root"),
                                       TString("H2Mu_ZH_bdt_training_dyAMC-J.root"),
                                       TString("H2Mu_gg_bdt_training_dyAMC-J.root"),
                                       };
    // bkg ntuples 
    std::vector<TString> bkgfilenames = {
//                                       TString("WW_bdt_training_dyAMC-J.root"),
//                                       TString("WZ_3l_bdt_training_dyAMC-J.root"),
                                       TString("ZJets_AMC_0j_bdt_training_dyAMC-J.root"),
                                       TString("ZJets_AMC_1j_bdt_training_dyAMC-J.root"),
//                                       TString("ZJets_AMC_2j_bdt_training_dyAMC-J.root"),
//                                       TString("ZZ_4l_bdt_training_dyAMC-J.root"),
//                                       TString("ZZ_2l_2q_bdt_training_dyAMC-J.root"),
//                                       TString("ZZ_2l_2v_bdt_training_dyAMC-J.root"),
                                       TString("tW_neg_bdt_training_dyAMC-J.root"),
                                       TString("tW_pos_bdt_training_dyAMC-J.root"),
//                                       TString("tZq_bdt_training_dyAMC-J.root"),
//                                       TString("ttW_bdt_training_dyAMC-J.root"),
//                                       TString("ttZ_bdt_training_dyAMC-J.root"),
                                       TString("tt_ll_AMC_bdt_training_dyAMC-J.root")
                                       };

    ////////////////////////////////////////////////////////////////////////////////////
    // Load signal and bkg ntuples into TMVA Dataloader
    ///////////////////////////////////////////////////////////////////////////////////
    
    TMVA::Factory factory("TMVAClassification", outputFile,
                          "!V:ROC:Transformations=I:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 
    TMVA::DataLoader loader("dataset");

    // define the training variables
    loader.AddVariable("dimu_pt");
    loader.AddVariable("dimu_eta");
    loader.AddVariable("dimu_abs_dEta");
    loader.AddVariable("dimu_abs_dPhi");
    loader.AddVariable("jet1_pt");
    loader.AddVariable("jet1_eta");
    loader.AddVariable("jet2_pt");
    loader.AddVariable("jet2_eta");
    loader.AddVariable("dijet1_mass");
    loader.AddVariable("dijet1_abs_dEta");
    loader.AddVariable("dijet2_mass");
    loader.AddVariable("dijet2_abs_dEta");
    loader.AddVariable("nJets");
    loader.AddVariable("nBMed");
    loader.AddVariable("MET");
    
    // load signal
    for(auto& filename: sigfilenames) 
    {
        auto inputFile = TFile::Open(filedir+filename);
        TTree *tsignal;
        inputFile->GetObject("theNtuple", tsignal);
        loader.AddSignalTree(tsignal,     1.0);   //signal weight  = 1
    }

    // load background
    for(auto& filename: bkgfilenames) 
    {
        auto inputFile = TFile::Open(filedir+filename);
        TTree *tbackground;
        inputFile->GetObject("theNtuple", tbackground);
        loader.AddBackgroundTree(tbackground, 1.0);   //background weight = 1 
    }
    
    // Apply some additional cuts here
    // Use a mass window with the same number of events on the left, in the signal region, and on the right
    // The mass intervals 113.8 to 120 GeV, 120 GeV to 130 GeV, and 130 GeV to 147.8 GeV fix this
    TCut mycuts, mycutb;
    mycuts = "dimu_mass_KaMu >= 113.8 && dimu_mass_KaMu <= 147.8";
    mycutb = mycuts;

    loader.SetSignalWeightExpression("weight");    // use the weight branch from the ntuple as the weight for training
    loader.SetBackgroundWeightExpression("weight");

    // total signal = ~300k
    // total backgr = ~350k
    loader.PrepareTrainingAndTestTree(mycuts, mycutb,
                                       "nTrain_Signal=200000:nTrain_Background=220000:SplitMode=Random:NormMode=NumEvents:!V" ); 
    
    factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG_AWB",
                        "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001" );


    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();    
    
    auto c1 = factory.GetROCCurve(&loader);
    c1->SaveAs("ROC.png");
}

int main(int argc, char* argv[])
{
    trainAndTest();
}
