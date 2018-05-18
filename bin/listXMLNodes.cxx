/////////////////////////////////////////////////////////////////////////////
//                           listXMLNodes.cxx                              //
//=========================================================================//
//                                                                         //
//       List the nodes created by our decision tree autocategorizer.      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"

#include "EventTools.h"
#include "TMVATools.h"
#include "PUTools.h"
#include "SignificanceMetrics.hxx"

#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TROOT.h"

#include "SampleDatabase.cxx"
#include "MuonInfo.h"
#include "EleInfo.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include <unordered_map>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
   TString xmlfile; 

   for(int i=1; i<argc; i++) 
   {    
       std::stringstream ss;  
       ss << argv[i];
       if(i==1) xmlfile = TString(ss.str().c_str());
   }    

   float reductionFactor = 1;

   
   XMLCategorizer categorySelection(xmlfile);

   std::map<TString, Sample*> samples;
   GetSamples(samples, "CERN", "H2Mu_VH");

   for(auto& sample: samples)
   {
      Sample* s = sample.second;
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      /////////////////////////////////////////////////////
      // Book training and spectator vars into reader

      TString dir    = "classification/";
      TString methodName = "BDTG_UF_v1";

      // sig vs bkg and multiclass (ggf, vbf, ... drell yan, ttbar) weight files
      TString weightfile = dir+"f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";
      TString weightfile_multi = dir+"f_Opt_v1_multi_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";

      TMVA::Reader* reader = 0; 
      std::map<TString, Float_t> tmap;
      std::map<TString, Float_t> smap;

      TMVA::Reader* reader_multi = 0; 
      std::map<TString, Float_t> tmap_multi;
      std::map<TString, Float_t> smap_multi;

      reader       = TMVATools::bookVars(methodName, weightfile, tmap, smap);
      reader_multi = TMVATools::bookVars(methodName, weightfile_multi, tmap_multi, smap_multi);

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelection;

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {
         // only load essential information for the first set of cuts 
         s->branches.muPairs->GetEntry(i);
         s->branches.muons->GetEntry(i);
         s->branches.eventInfo->GetEntry(i);

         // loop and find a good dimuon candidate
         if(s->vars.muPairs->size() < 1) continue;
         bool found_good_dimuon = false;

         // find the first good dimuon candidate and fill info
         for(auto& dimu: (*s->vars.muPairs))
         {
            // Reset the categorizer in preparation for the next event
            categorySelection.reset();

            // the dimuon candidate and the muons that make up the pair
            s->vars.dimuCand = &dimu;
            MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
            MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

            // normal selections
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

            found_good_dimuon = true;

            // Load the rest of the information needed for run2 categories
            s->branches.getEntry(i);

            // clear vectors for the valid collections
            s->vars.validMuons.clear();
            s->vars.validExtraMuons.clear();
            s->vars.validElectrons.clear();
            s->vars.validJets.clear();
            s->vars.validBJets.clear();

            // load valid collections from s->vars raw collections
            jetCollectionCleaner.getValidJets(s->vars, s->vars.validJets, s->vars.validBJets);
            muonCollectionCleaner.getValidMuons(s->vars, s->vars.validMuons, s->vars.validExtraMuons);
            eleCollectionCleaner.getValidElectrons(s->vars, s->vars.validElectrons);

            // Clean jets and electrons from muons, then clean remaining jets from remaining electrons
            CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.4);
            CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
            CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);

            // Figure out which category the event belongs to
            categorySelection.evaluate(s->vars);
            EventTools::outputEvent(s->vars, categorySelection);
            return 0;

            // ouput pt, mass info etc for the event
            if(s->vars.validExtraMuons.size() + s->vars.validElectrons.size() >= 1)
            {
               EventTools::outputEvent(s->vars, categorySelection);
               std::cout << "  @@@@ M_mu_mu " << std::endl;
               for(unsigned int i=0; i<s->vars.validExtraMuons.size(); i++) 
               {
                   for(unsigned int j=i+1; j<s->vars.validExtraMuons.size(); j++) 
                   {
                       TLorentzVector dimu = s->vars.validExtraMuons[i] + s->vars.validExtraMuons[j];
                       std::cout << Form("  m%d,%d = %f \n", i,j, dimu.M());
                   }
               }
               std::cout << std::endl;
               std::cout << "  @@@@ M_e_e " << std::endl;
               for(unsigned int i=0; i<s->vars.validElectrons.size(); i++) 
               {
                   for(unsigned int j=i+1; j<s->vars.validElectrons.size(); j++) 
                   {
                       TLorentzVector diele = s->vars.validElectrons[i] + s->vars.validElectrons[j];
                       std::cout << Form("  m%d,%d = %f \n", i,j, diele.M());
                   }
               }
               std::cout << std::endl;
            }
            //------------------------------------------------------------------
            ////////////////////////////////////////////////////////////////////

            if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

         }
      }
      std::cout << Form("  /// Done processing %s \n", s->name.Data());
   }
}
