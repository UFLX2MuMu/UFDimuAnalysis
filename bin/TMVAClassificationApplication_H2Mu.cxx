/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication_H2Mu                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

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

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TXMLEngine.h"

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
#include "TMVATools.h"
#include "TMVA_helper.h"

//////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////

void TMVAClassificationApplication_H2Mu( TString myMethodList = "" ) 
{

   /////////////////////////////////////////////////////
   // Get a sample to test the classification on

   std::map<TString, Sample*> sampleMap;
   GetSamples(sampleMap, "CERN", "ZJets_MG_incl");
   GetSamples(sampleMap, "CERN", "H2Mu_gg");
   GetSamples(sampleMap, "CERN", "H2Mu_VBF");
   for(auto& i: sampleMap)
       printf("/// %s in Sample Map\n", i.second->name.Data());

   Sample* s = sampleMap["ZJets_MG"];
   //Sample* s = sampleMap["H2Mu_gg"];
   //Sample* s = sampleMap["H2Mu_VBF"];
   printf("/// Using %s for training\n", s->name.Data());

   TMVA::Tools::Instance();

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication_H2Mu" << std::endl << std::endl;

   TString dir    = "classification/";
   TString methodName = "BDTG_UF_v1";

   TString weightfile       = dir+"f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";
   TString weightfile_multi = dir+"f_Opt_v1_multi_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";

   /////////////////////////////////////////////////////
   // Book training and spectator vars into reader

   std::map<TString, Float_t> tmap;
   std::map<TString, Float_t> smap;
   TMVA::Reader* reader = TMVATools::bookVars(methodName, weightfile, tmap, smap);
   
   std::map<TString, Float_t> tmap_multi;
   std::map<TString, Float_t> smap_multi;
   TMVA::Reader* reader_multi = TMVATools::bookVars(methodName, weightfile_multi, tmap_multi, smap_multi);
   
   std::vector<TString> classes; 
   TMVATools::getClassNames(weightfile_multi, classes);
   
   /////////////////////////////////////////////////////
   // Loop over the events in the sample

   TStopwatch sw;
   sw.Start();

   // Objects to help with the cuts and selections
   JetCollectionCleaner      jetCollectionCleaner;
   MuonCollectionCleaner     muonCollectionCleaner;
   EleCollectionCleaner      eleCollectionCleaner;

   Run2MuonSelectionCuts  run2MuonSelection;
   Run2EventSelectionCuts run2EventSelection;

   // set some flags
   bool isData = s->sampleType.EqualTo("data");

   std::cout << Form("/// Processing %s \n", s->name.Data());
   s->setBranchAddresses(2);
   int ngood = 0;
   for(unsigned int i=0; i<s->N; i++)
   {
      //std::cout << Form("  /// GetEntry muon info %s \n", s->name.Data());
      // only load essential information for the first set of cuts 
      s->branches.muPairs->GetEntry(i);
      s->branches.muons->GetEntry(i);
      s->branches.eventInfo->GetEntry(i);

      // loop and find a good dimuon candidate
      if(s->vars.muPairs->size() != 1) continue;
     
      // the dimuon candidate and the muons that make up the pair
      //std::cout << Form("  /// Set dimuon info %s \n", s->name.Data());
      MuPairInfo& dimu = s->vars.muPairs->at(0); 
      s->vars.dimuCand = &dimu;
      MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
      MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

      ////////////////////////////////////////////////////////////////////
      // CUTS  ----------------------------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      //std::cout << Form("  /// Cuts %s \n", s->name.Data());
      if(dimu.mass < 120 || dimu.mass > 130)
      {
          continue;
      }
      if(!run2EventSelection.evaluate(s->vars))
      {
          continue;
      }
      if(!mu1.isMediumID || !mu2.isMediumID)
      {
          continue;
      }
      if(!run2MuonSelection.evaluate(s->vars))
      {
          continue;
      }

      //std::cout << Form("  /// Get remaining branches %s \n", s->name.Data());
      s->branches.getEntry(i);
      ngood++;

      // clear vectors for the valid collections
      s->vars.validMuons.clear();
      s->vars.validExtraMuons.clear();
      s->vars.validElectrons.clear();
      s->vars.validJets.clear();
      s->vars.validBJets.clear();

      //std::cout << Form("  /// Load valid collections %s \n", s->name.Data());
      // load valid collections from s->vars raw collections
      jetCollectionCleaner.getValidJets(s->vars, s->vars.validJets, s->vars.validBJets);
      muonCollectionCleaner.getValidMuons(s->vars, s->vars.validMuons, s->vars.validExtraMuons);
      eleCollectionCleaner.getValidElectrons(s->vars, s->vars.validElectrons);

      // Clean jets and electrons from muons, then clean remaining jets from remaining electrons
      CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.4);
      CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
      //CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);

      std::cout << std::endl;
      Float_t valb = TMVATools::getClassifierScore(reader, methodName, tmap, s->vars);
      printf("!!! %d) BDT_Prediction: %f\n", s->vars.eventInfo->event, valb);

      std::vector<float> vals = TMVATools::getMulticlassScores(reader_multi, methodName, tmap_multi, s->vars);
      for(unsigned int j=0; j<vals.size(); j++)
      {
          printf("!!! %d) %s: %f\n", s->vars.eventInfo->event, classes[j].Data(), vals[j]);
      }
      std::cout << std::endl;
    
      if(ngood > 10) break;
   }
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   delete reader;
    
   std::cout << "==> TMVAClassificationApplication_H2Mu is done!" << std::endl << std::endl;
}

//////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList; 
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(","); 
      methodList += regMethod;
   }
   TMVAClassificationApplication_H2Mu(methodList);
   return 0;
}
