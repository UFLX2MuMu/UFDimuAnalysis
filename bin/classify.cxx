
//////////////////////////////////////////////////////////////////////////////
///  Simplified Higgs vs. background classification                         //
//                                                                          //     
//  Runs on our h->mumu signal and background ttrees skimmed from           //
//   miniAOD. Further selections are made here and the events are added     //
//   one at a time to the training/testing sets. Testing each event         //
//   for the different selections takes a long time for these large         // 
//   MC samples. See classifyDataFrame.cxx to cut down the time. This runs  //
//   on flat ntuples that have already further skimmed events.              //
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cassert>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
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

// Specific includes for UFDimuAnalysis
#include "Sample.h"
#include "SampleDatabase.cxx"
#include "MuonSelection.h"
#include "EventSelection.h"
#include "CategorySelection.h"
#include "MuonCollectionCleaner.h"
#include "JetCollectionCleaner.h"
#include "EleCollectionCleaner.h"

// Extra tools - AWB 13.03.17
#include "EventTools.h"
#include "TMVA_helper.h"

using namespace TMVA;
const double LUMI = 36814; // pb-1 

// Create a new root output file
TString outdir = "/home/puno/h2mumu/UFDimuAnalysis_v2/bin/tmva_out";
TString out_file_name = outdir+"/baseline_plus_resolution_weight_mu_eta.root";
TString factoryname = "factory";
TString dlname = "dataset";

// use 1/4 of events for testing
int modtest = 4;

// training and spectator variables for TMVA
std::vector<TString> trainVarNames = {
                                       TString("dimu_pt"),
                                       TString("dimu_eta"),
                                       TString("dimu_abs_dEta"),
                                       TString("dimu_abs_dPhi"),
                                       TString("jet1_pt"),
                                       TString("jet1_eta"),
                                       TString("jet2_pt"),
                                       TString("jet2_eta"),
                                       TString("dijet1_mass"),
                                       TString("dijet1_abs_dEta"),
                                       TString("dijet2_mass"),
                                       TString("dijet2_abs_dEta"),
                                       TString("nJets"),
                                       TString("nBMed"),
                                       TString("MET")
                                     };

std::vector<TString> specVarNames = {
                                       TString("samp_ID"),
                                       TString("evt_wgt"),
                                       TString("res_wgt"),
                                       TString("dimu_mass_KaMu")
                                     };

void setUpDatasetVars(DataLoader* dataloader, std::vector<TString>& trainVarNames, std::vector<TString>& specVarNames,
                      std::vector<TString>& varNames, std::vector<Double_t>& varVals)
{
// Look in lib/VarSet.cxx for the full list of available features
// Add any you want to this database

   // Defined in interface/MVA_helper.h
   // TMVA_var(TString name, TString descr, TString unit, TString type, Double_t def_val)
   std::map<TString,TMVA_var> vars;  // All spectator variables
   
   /////////////////////////////////////////////////////////
   ///  Input variables                                  ///
   /////////////////////////////////////////////////////////

   // Muon variables
   vars["mu1_pt"] = TMVA_var( "mu1_pt",        "p_{T}(#mu1)",       "GeV", 'F', -88 ); 
   vars["mu2_pt"] = TMVA_var( "mu2_pt",        "p_{T}(#mu2)",       "GeV", 'F', -88 ); 
   vars["mu1_abs_eta"] = TMVA_var( "mu1_abs_eta",   "|#eta(#mu1)|",         "", 'F', -88 ); 
   vars["mu2_abs_eta"] = TMVA_var( "mu2_abs_eta",   "|#eta(#mu2)|",         "", 'F', -88 ); 

   vars["dimu_pt"]         = TMVA_var( "dimu_pt",       "p_{T}(#mu#mu)",     "GeV", 'F', -88 ); 
   vars["dimu_eta"]   = TMVA_var( "dimu_eta",      "#eta(#mu#mu)",         "", 'F', -88 );   
   vars["dimu_rapid"] = TMVA_var( "dimu_rapid",    "rapid(#mu#mu)",        "", 'F', -88 );  

   vars["dimu_dR"]       = TMVA_var( "dimu_dR",       "dR(#mu#mu)",           "", 'F', -88 ); 
   vars["dimu_abs_dEta"] = TMVA_var( "dimu_abs_dEta", "|d#eta(#mu#mu)|",      "", 'F', -88 ); 
   vars["dimu_abs_dPhi"] = TMVA_var( "dimu_abs_dPhi", "|d#phi(#mu#mu)|",      "", 'F', -88 );   
   vars["dimu_dPhiStar"] = TMVA_var( "dimu_dPhiStar", "d#phi*(#mu#mu)",       "", 'F', -88 );  

   // Jet variables
   vars["jet1_pt"]  = TMVA_var( "jet1_pt",         "p_{T}(jet1)",       "GeV", 'F', -88 ); 
   vars["jet2_pt"]  = TMVA_var( "jet2_pt",         "p_{T}(jet2)",       "GeV", 'F', -88 ); 
   vars["jet1_eta"] = TMVA_var( "jet1_eta",        "#eta(jet1)",           "", 'F', -88 ); 
   vars["jet2_eta"] = TMVA_var( "jet2_eta",        "#eta(jet2)",           "", 'F', -88 ); 

   vars["dijet1_mass"] = TMVA_var( "dijet1_mass",     "1^{st} M(jj)",      "GeV", 'F', -88 ); 
   vars["dijet2_mass"] = TMVA_var( "dijet2_mass",     "2^{nd} M(jj)",      "GeV", 'F', -88 ); 
   vars["dijet3_mass"] = TMVA_var( "dijet3_mass",     "3^{rd} M(jj)",      "GeV", 'F', -88 ); 
   vars["dijet4_mass"] = TMVA_var( "dijet4_mass",     "4^{th} M(jj)",      "GeV", 'F', -88 ); 

   vars["dijet1_abs_dEta"] = TMVA_var( "dijet1_abs_dEta", "1^{st} |d#eta(jj)|",   "", 'F', -88 ); 
   vars["dijet2_abs_dEta"] = TMVA_var( "dijet2_abs_dEta", "2^{nd} |d#eta(jj)|",   "", 'F', -88 ); 
   vars["dijet3_abs_dEta"] = TMVA_var( "dijet3_abs_dEta", "3^{rd} |d#eta(jj)|",   "", 'F', -88 ); 
   vars["dijet4_abs_dEta"] = TMVA_var( "dijet4_abs_dEta", "4^{th} |d#eta(jj)|",   "", 'F', -88 ); 

   // Global event variables
   vars["nJets"]     = TMVA_var( "nJets",       "# of jets",            "", 'I', -88 ); 
   vars["nJetsCent"] = TMVA_var( "nJetsCent",   "# of central jets",    "", 'I', -88 ); 
   vars["nJetsFwd"]  = TMVA_var( "nJetsFwd",    "# of forward jets",    "", 'I', -88 ); 
   vars["nBMed"]     = TMVA_var( "nBMed",       "# of medium b-tags",   "", 'I', -88 ); 

   vars["MET"]      = TMVA_var( "MET",         "MET",               "GeV", 'F', -88 ); 
   vars["MHT"]      = TMVA_var( "MHT",         "MHT",               "GeV", 'F', -88 ); 
   vars["MT_had"]   = TMVA_var( "MT_had",      "M_{T} of jets",     "GeV", 'F', -88 ); 
   vars["mass_had"] = TMVA_var( "mass_had",    "Mass of jets",      "GeV", 'F', -88 ); 

   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   vars["samp_ID"]         = TMVA_var( "samp_ID",        "Sample ID",            "", 'I', -77 );
   vars["evt_wgt"]        = TMVA_var( "evt_wgt",       "Sample weight",        "", 'F', -77 );
   vars["res_wgt"]         = TMVA_var( "res_wgt",        "Resolution weight",    "", 'F', -77 );
   vars["dimu_mass_Roch"]  = TMVA_var( "dimu_mass_Roch", "mass(#mu#mu)",      "GeV", 'F', -77 );
   vars["dimu_mass_KaMu"]  = TMVA_var( "dimu_mass_KaMu", "mass(#mu#mu)",      "GeV", 'F', -77 );

   // ### FOR EVERY VARIABLE ADD THE TRAINING VARS AND SPECTATOR VARS
   for(auto& varName: trainVarNames) // training vars
   {
       if(vars.count(varName) == 0) std::cout << "Training variable " << varName.Data() << " not in var set." << std::endl;
       TMVA_var v = vars[varName];
       dataloader->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader 
       varNames.push_back(varName);
       varVals.push_back(v.def_val);
   }

   for(auto& varName: specVarNames) // spectator vars
   {
       if(vars.count(varName) == 0) std::cout << "Spectator variable " << varName.Data() << " not in var set." << std::endl;
       TMVA_var v = vars[varName];
       dataloader->AddSpectator( v.name, v.descr, v.unit, v.type );
       varNames.push_back(varName);
       varVals.push_back(v.def_val);
   }

   TString trainVarString = "";
   TString specVarString = "";
   TString totalVarString = "";

   for(auto& varName: trainVarNames)
       trainVarString+=varName+"\n";

   for(auto& varName: specVarNames)
       specVarString+=varName+"\n";
       
   for(auto& varName: varNames)
       totalVarString+=varName+"\n";

    std::cout << std::endl;
    std::cout << "//////////// train vars /////////" << std::endl << trainVarString.Data() << std::endl << std::endl;
    std::cout << std::endl;
    std::cout << "//////////// spectator vars /////" << std::endl << specVarString.Data() << std::endl << std::endl;
    std::cout << std::endl;
    std::cout << "//////////// all vars ///////////" << std::endl << totalVarString.Data() << std::endl << std::endl;
    std::cout << std::endl;
}

void classify(TString factoryname, TString dlname, std::vector<TString>& trainVarNames, std::vector<TString>& specVarNames)
{
   // initialize TMVA
   TMVA::Tools::Instance();

   // create the output file
   TFile* out_file = TFile::Open( out_file_name, "RECREATE" );

   /////////////////////////////////////////////////////////////
   ///  Load input samples: MC signal, MC background, data   ///
   /////////////////////////////////////////////////////////////

   std::map<TString, Sample*> samples;
   std::cout << "\nAbout to get the samples" << std::endl;

   // Load samples into our map. "UF" if you are on the UF servers
   // or "CERN" if you at CERN. "ALL" specifies that we want to load the Data
   // and all of the MC samples. Can loop through and remove the ones you don't want 
   // to use if you desire or just grab the ones you care about from the map.
   GetSamples(samples, "UF", "SIGNAL" );       // all of the signal samples
   GetSamples(samples, "UF", "ZJets_AMC-J");   // Drell Yan with MC samples for 0-jets, 1-jets, 2 or more jets
   GetSamples(samples, "UF", "tt_ll_MG");      // ttbar
   GetSamples(samples, "UF", "singleTop");     // t + X
   
   std::map<TString, int> idmap;

   idmap["H2Mu_gg"]     = -1;
   idmap["H2Mu_VBF"]    = -2;
   idmap["H2Mu_ZH"]     = -3;
   idmap["H2Mu_WH_pos"] = -4;
   idmap["H2Mu_WH_neg"] = -5;

   idmap["ZJets_AMC_0j"] = +1;
   idmap["ZJets_AMC_1j"] = +2;
   idmap["ZJets_AMC_2j"] = +3;

   idmap["tt_ll_MG"] =   +10;
   idmap["tW_pos"]   =   +11;
   idmap["tW_neg"]   =   +12;
   idmap["tZq"]      =   +13;

   std::cout << std::endl << "\n========== Got the samples ============" << std::endl;

   //////////////////////////////////////////////////////////////////
   ///  Initialize some objects to help with cuts and cleaning
   //////////////////////////////////////////////////////////////////

   // Objects to help clean valid objects from the net collections in samp->vars
   JetCollectionCleaner   jetCollectionCleaner;
   MuonCollectionCleaner  muonCollectionCleaner;
   EleCollectionCleaner   eleCollectionCleaner;

   // Objects to cut events from the analysis
   // See selection/MuonSelection.h or selection/EventSelection.h
   // Choose from an available option or make your own implementing the interface
   // Use selection.evaluate(samp->vars) to see whether an event passes or fails the cuts
   Run2MuonSelectionCuts   run2MuonSelection;
   Run2EventSelectionCuts  run2EventSelection;

   //////////////////////////////////////////////////////////////////
   ///  Set Up TMVA Factory [Master object for all ML Algos]      ///
   //////////////////////////////////////////////////////////////////
   
   TString fact_settings         = "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification";
   // if you wanted to do multiclass
   // fact_settings = "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=multiclass";


   TMVA::Factory* factory = new TMVA::Factory(factoryname, out_file, fact_settings); 
   TMVA::DataLoader* dataloader = new TMVA::DataLoader(dlname);                     

   std::vector<TString> varNames;     // Holds values of variables for a given factory
   std::vector<Double_t> varVals;      // Holds values of variables for a given factory
   setUpDatasetVars(dataloader, trainVarNames, specVarNames, varNames, varVals);

   //////////////////////////////////////////////////i//////////////////////////
   // Set up function describing inclusive signal to get resolution information
   //////////////////////////////////////////////////////////////////////////////
   
   // Can weight the events by resolution if we want.
   // The inclusive signal distribution will help figure out where in the width the event is.
   // For inclusive ggH signal, we have the following triple gaussian
   // 0.73*Gaus(124.8, 1.52) + 0.23*Gaus(122.8, 4.24) + 0.04*Gaus(126, 2.1)
   TString tripGausExpr = 
        "(0.73*TMath::Gaus(x, 124.8, 1.52, 1) + 0.23*TMath::Gaus(x, 122.8, 4.24, 1) + 0.04*TMath::Gaus(x, 126, 2.1, 1))";
   TF1* tripGaus   = new TF1("tripGaus", tripGausExpr, 113.8, 147.8);
   TF1* tripGausSq = new TF1("tripGaus", tripGausExpr+" * "+tripGausExpr, 113.8, 147.8);
   Double_t tripGausNorm = tripGausSq->Integral(-1000, 1000);
   std::cout << "Triple gaussian has a normalization of " << tripGaus->Integral(-1000, 1000) << std::endl;

   //////////////////////////////////////////////////////////////////
   // ADD EVENTS FROM EACH SAMPLE THAT PASS CUTS TO TRAINING DATASET
   //////////////////////////////////////////////////////////////////

   std::cout << "\n******* About to loop over samples *******" << std::endl;
   std::cout << std::endl;

   for(auto &s: samples)
   {
     Sample* samp = s.second;
     int num_passed_in_sample = 0;

     std::cout << "///////////////////////////////////////////////////////////////////////" << std::endl;
     std::cout << "Looping over " << samp->name << " [" << samp->N << " events]" << std::endl;
     std::cout << "///////////////////////////////////////////////////////////////////////" << std::endl;
     std::cout << std::endl;

     for (UInt_t iEvt = 0; iEvt < samp->N; iEvt++) 
     {
       //////////////////////////////////////////////////////////////////
       // CUT EVENTS WE DONT WANT
       //////////////////////////////////////////////////////////////////
       samp->branches.eventInfo->GetEntry(iEvt); 

       // Prevent double counting in RunF
       if (samp->name == "RunF_1" && samp->vars.eventInfo->run > 278801)
	 continue;
       if (samp->name == "RunF_2" && samp->vars.eventInfo->run < 278802)
	 continue;

       // Only use half of the signal events for training
       // Save the other half of the signal events for use in limit-setting
       if (samp->sampleType == "signal" && (samp->vars.eventInfo->event % 2 == 0)) 
	 continue;

       // Load muon info
       samp->branches.muPairs->GetEntry(iEvt);

       // no muon pairs in this event
       if (samp->vars.muPairs->size() == 0)
	 continue;

       // Only train on events in window from 113.8 to 147.8
       // In data, we have the same number events in [113.8, 120], [120, 130], and [130, 147.8] GeV
       if ( samp->vars.muPairs->at(0).mass_Roch < 113.8 || samp->vars.muPairs->at(0).mass_Roch > 147.8 )
	 continue;

       // Now load all the info needed
       samp->branches.getEntry(iEvt);

       // Set the Higgs dimuon candidate pair
       samp->vars.dimuCand = &(samp->vars.muPairs->at(0));

       // Throw away events that fail the event or muon selection
       if ( !run2EventSelection.evaluate(samp->vars) )
	 continue;

       if ( !run2MuonSelection.evaluate(samp->vars) )
	 continue;

       if ( samp->vars.muons->at(samp->vars.dimuCand->iMu1).isMediumID != 1 ||
	    samp->vars.muons->at(samp->vars.dimuCand->iMu2).isMediumID != 1 )
	 continue;

       // Throw away events that fail MC specific cuts
       if ( samp->name == "ZJets_MG" && samp->vars.lhe_ht > 70 )
	 continue;

       num_passed_in_sample++;

       // Clear vectors for the valid collections
       samp->vars.validMuons.clear();
       samp->vars.validExtraMuons.clear();
       samp->vars.validElectrons.clear();
       samp->vars.validJets.clear();
       samp->vars.validBJets.clear();

       // Load valid collections from samp->vars raw collections
       // These objects pass certain cuts and are considered in the analysis
       muonCollectionCleaner.getValidMuons(samp->vars, samp->vars.validMuons, samp->vars.validExtraMuons);
       jetCollectionCleaner.getValidJets(samp->vars, samp->vars.validJets, samp->vars.validBJets);
       eleCollectionCleaner.getValidElectrons(samp->vars, samp->vars.validElectrons);

       // Clean jets from muons
       CollectionCleaner::cleanByDR(samp->vars.validJets, samp->vars.validMuons, 0.4);

       //////////////////////////////////////////////////////
       ///  End block mostly lifted from bin/example.cxx  ///
       //////////////////////////////////////////////////////

       MuPairInfo& dimu = samp->vars.muPairs->at(0); 
       MuonInfo& mu1    = samp->vars.muons->at(dimu.iMu1);
       MuonInfo& mu2    = samp->vars.muons->at(dimu.iMu2);
       Int_t nJets      = samp->vars.getValue("nJets");
       Int_t nValJets   = samp->vars.getValue("nValJets");
       Int_t nBMed      = samp->vars.getValue("nBMed");
       Int_t nBLoose    = samp->vars.getValue("nBLoose");
       Float_t MET      = samp->vars.getValue("MET");

       if (nJets != nValJets) {
	 std::cout << "\n  * Bizzare event where nJets = " << nJets << ", nValJets = " << nValJets << std::endl;
	 continue;
       }

       //////////////////////////////////////////////////////////////////////////////////
       ///  Weight signal events by cross section and inclusive H2Mu mass resolution  ///
       //////////////////////////////////////////////////////////////////////////////////

       Double_t res_wgt = tripGaus->Eval(dimu.mass_Roch) / tripGausNorm;
       Double_t evt_wgt = 1.0;

       if (samp->sampleType != "data") // Half of signal / background MC events go into training, half into testing classifier
	 evt_wgt = 2.0 * samp->getWeight() * samp->getLumiScaleFactor(LUMI);
       if (samp->sampleType == "signal" < 0)
	 evt_wgt *= 2.0; // Even-numbered signal events were reserved for limit-setting

       for(unsigned int iVar=0; iVar<varNames.size(); iVar++)
       {
           TString vName = varNames[iVar];
	   /////////////////////////////
	   ///  Spectator variables  ///
	   /////////////////////////////
	   
	   if      ( vName == "samp_ID" )
	     varVals.at(iVar) = idmap[samp->name];
	   else if ( vName == "evt_wgt" )
	     varVals.at(iVar) = evt_wgt;
	   else if ( vName == "res_wgt" )
	     varVals.at(iVar) = res_wgt;

           // Other variables automatically provided through lib/VarSet.h via getValue("varname")
	   else 
	     varVals.at(iVar) = samp->vars.getValue(vName.Data());
         }
	  
         if(num_passed_in_sample <= 2) 
         {
             std::cout << std::endl;
             for(unsigned int i=0; i<varNames.size(); i++)
             {
                 std::cout << varNames[i].Data() << ": " << varVals[i] << std::endl;
             }
             std::cout << std::endl;
         }
	 
	 // Weight by expected sample normalization x signal resolution
	 Double_t sig_evt_weight = evt_wgt*res_wgt;
	 Double_t bkg_evt_weight = evt_wgt;
	     
	 // Load values into event
	 if(num_passed_in_sample % modtest == 0)
         {
             if(samp->sampleType == "signal") dataloader->AddSignalTestEvent( varVals, sig_evt_weight );
             if(samp->sampleType == "background") dataloader->AddBackgroundTestEvent( varVals, sig_evt_weight );
         }
         else
         {
             if(samp->sampleType == "signal") dataloader->AddSignalTrainingEvent( varVals, sig_evt_weight );
             if(samp->sampleType == "background") dataloader->AddBackgroundTrainingEvent( varVals, sig_evt_weight );
         }

     } // end loop over events 
   } // End loop over samples 
   std::cout << "******* Made it out of the event loop *******" << std::endl;
   
     
   dataloader->SetWeightExpression( 1.0 );
   dataloader->PrepareTrainingAndTestTree( "", "", "" ); // Default 
   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG_AWB", (std::string)
                        "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001" );
   
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
    
   // Save the output
   out_file->Close();

   std::cout << "==> Wrote root file: " << out_file->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   // delete factory;
   // delete dataloader;
}


int main( int argc, char** argv )
{
   classify(factoryname, dlname, trainVarNames, specVarNames);
   return 0;
}

