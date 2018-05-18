/////////////////////////////////////////////////////////////////////////////
//                           synchronize.cxx                               //
//=========================================================================//
//                                                                         //
//   Synchronize with the other groups in the analysis.                    //
//   Output information about selected events, output lists of events      //
//   in each category to csv, read in lists of events from csv if there    //
//   are many events to investigate.                                       //
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

#include "TMVATools.h"
#include "EventTools.h"
#include "PUTools.h"

#include "SignificanceMetrics.hxx"
#include "SampleDatabase.cxx"
#include "ThreadPool.hxx"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TBranchElement.h"
#include "TROOT.h"
//#include "TStreamerInfo.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

UInt_t getNumCPUs()
{
  SysInfo_t s;
  gSystem->GetSysInfo(&s);
  UInt_t ncpu  = s.fCpus;
  return ncpu;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    gROOT->SetBatch();
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    int nthreads = 1;        // number of threads to use in parallelization

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> nthreads;
    }   
    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    float luminosity = 36814;      // pb-1
    float reductionFactor = 1;     // reduce the number of events you run over in case you want to debug or some such thing

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    GetSamples(samples, "CERN");

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING: SetBranchAddresses-------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    // Add to the vector we loop over to categorize
    
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto &i : samples)
    {

        //if(i.second->sampleType != "data" && i.second->name != "H2Mu_gg" && i.second->name != "ZZ_2l_2v") continue;
        //if(i.second->sampleType == "data" && i.second->name != "RunD" && i.second->name != "RunH") continue;
        if(i.second->sampleType != "data" && i.second->name != "H2Mu_gg") continue;
        if(i.second->sampleType == "data") continue;

        // Output some info about the current file
        std::cout << "  /// Using sample " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filenames[0] << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        i.second->setBranchAddresses("");
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    ///////////////////////////////////////////////////////////////////
    // Get Histo Information from Input -------------------------------
    ///////////////////////////////////////////////////////////////////

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    ///////////////////////////////////////////////////////////////////
    // Define Task for Parallelization -------------------------------
    ///////////////////////////////////////////////////////////////////

    auto synchSample = [reductionFactor](Sample* s)
    {
      // set the muon calibration type for synchronization
      TString pf_roch_or_kamu = "PF";
      float dimuMassMin = 60;

      // TMVA Classifier
      TString dir    = "classification/";
      TString methodName = "BDTG_UF_v1";
      TString weightfile = dir+"f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";

      TMVA::Reader* reader = 0; 
      std::map<TString, Float_t> tmap;
      std::map<TString, Float_t> smap;
      reader = TMVATools::bookVars(methodName, weightfile, tmap, smap);

      std::vector<std::pair<int,long long int>> eventsToCheck;
      //loadEventsFromFile("synchcsv/xCheck.txt", eventsToCheck);

      //eventsToCheck.push_back(std::pair<int, long long int>(1, 212289));
      //eventsToCheck.push_back(std::pair<int, long long int>(1, 215743));
      //eventsToCheck.push_back(std::pair<int, long long int>(1, 232626));
      //eventsToCheck.push_back(std::pair<int, long long int>(1, 241209));

      ///////////////////////////////////////////////////////////////////
      // INIT Cuts and Categories ---------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts    synchMuonSelection;
      Run2EventSelectionCuts   synchEventSelection;
      synchEventSelection.cDimuMassMin = dimuMassMin;
      //SynchEventSelectionCuts  synchEventSelection;

      Categorizer* categorySelection = new XMLCategorizer("tree_categorization_final.xml");
      //Categorizer* categorySelection = new CategorySelectionSynch();

      ///////////////////////////////////////////////////////////////////
      // INIT HISTOGRAMS TO FILL ----------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram to fill in the category
      TString hkey = s->name;

      // Different categories for the analysis
      for(auto &c : categorySelection->categoryMap)
      {
          // histogram settings
          float min = 0;
          float max = 2;
          int hbins = 2;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle("count");
          c.second.histoList->Add(c.second.histoMap[hkey]);

          c.second.eventsMap[hkey] = std::vector< std::pair<int, long long int> >();
      }


      ///////////////////////////////////////////////////////////////////
      // LOOP OVER EVENTS -----------------------------------------------
      ///////////////////////////////////////////////////////////////////

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
          // Reset the flags in preparation for the next event
          categorySelection->reset();

          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu; 
          MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);
          s->vars.setCalibrationType(pf_roch_or_kamu);

          // Load the rest of the information needed for run2 categories
          s->branches.getEntry(i);
          s->vars.setCalibrationType(pf_roch_or_kamu); // reloaded the branches, need to set mass,pt to correct calibrations again

          std::pair<int, long long int> e(s->vars.eventInfo->run, s->vars.eventInfo->event); // create a pair that identifies the event uniquely

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

          // Clean jets from muons
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.4);
          CollectionCleaner::cleanByDR(s->vars.validBJets, s->vars.validMuons, 0.4);
          //CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
          //CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);
          
          s->vars.bdt_out = TMVATools::getClassifierScore(reader, methodName, tmap, s->vars); // set tmva's bdt score
          
          if(EventTools::eventInVector(e, eventsToCheck) || true) // Adrian gave a list of events to look at for synch purposes
             EventTools::outputEvent(s->vars);

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          if(!synchMuonSelection.evaluate(s->vars)) 
          {
              if(EventTools::eventInVector(e, eventsToCheck) || true)
                  std::cout << "   !!! Fail Muon Selection !!!" << std::endl;

              continue; 
          }
          if(!synchEventSelection.evaluate(s->vars))
          { 
              if(EventTools::eventInVector(e, eventsToCheck) || true)
                  std::cout << "   !!! Fail Event Selection !!!" << std::endl;

              continue; 
          }
          if(!mu1.isMediumID || !mu2.isMediumID)
          { 
              if(EventTools::eventInVector(e, eventsToCheck) || true)
                  std::cout << "   !!! Fail Medium ID !!!" << std::endl;

              continue; 
          }

          // avoid double counting in RunF
          if(s->name == "RunF_1" && s->vars.eventInfo->run > 278801)
          {
              continue;
          }
          if(s->name == "RunF_2" && s->vars.eventInfo->run < 278802)
          {
              continue;
          }


          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true; 

          // Figure out which category the event belongs to
          categorySelection->evaluate(s->vars);

          // Look at each category
          for(auto &c : categorySelection->categoryMap)
          {
              // skip categories
              if(c.second.hide) continue;
              if(!c.second.inCategory) continue;

              c.second.histoMap[hkey]->Fill(1, 1);
              c.second.eventsMap[hkey].push_back(e);
          } // end category loop

          if(EventTools::eventInVector(e, eventsToCheck) || true)
            categorySelection->outputResults();

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop
      } // end event loop

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return categorySelection;

    }; // done defining synchSample

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(nthreads);
    std::vector< std::future<Categorizer*> > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(synchSample, s));

   ///////////////////////////////////////////////////////////////////
   // Gather all the Histos into one Categorizer----------------------
   ///////////////////////////////////////////////////////////////////

    //Categorizer* cAll = new CategorySelectionRun1();
    Categorizer* cAll = new XMLCategorizer("tree_categorization_final.xml");

    // get histos from all categorizers and put them into one
    for(auto && categorizer: results)  // loop through each Categorizer object, one per sample
    {
        for(auto& category: categorizer.get()->categoryMap) // loop through each category
        {
            // category.first is the category name, category.second is the Category object
            if(category.second.hide) continue;
            for(auto& h: category.second.histoMap) // loop through each histogram in the category
            {
                // std::cout << Form("%s: %f", h.first.Data(), h.second->Integral()) << std::endl;
                // h.first is the sample name, category.histoMap<samplename, TH1D*>
                Sample* s = samples[h.first];

                cAll->categoryMap[category.first].histoMap[h.first] = h.second;
                cAll->categoryMap[category.first].histoList->Add(h.second);
                std::cout << Form("%s: %d\n", h.second->GetName(), (int)h.second->Integral());
            }
            for(auto& v: category.second.eventsMap)
            {
                TString filename = "csv/synchcsv/"+v.first+"_"+category.second.name+".csv";
                EventTools::outputEventsToFile(v.second, filename);
                v.second.clear();
            }
        }
        std::cout << Form("\n");
    }

   ///////////////////////////////////////////////////////////////////
   // Save All of the Histos-----------------------------------------
   ///////////////////////////////////////////////////////////////////
   TList* list = new TList();

    std::cout << "About to loop through cAll" << std::endl;
    for(auto &c : cAll->categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;
        for(auto& item: c.second.histoMap)
        {
            list->Add(item.second);
        }
    }
    std::cout << std::endl;
    TString savename = Form("rootfiles/synchronize.root");

    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;

    TFile* savefile = new TFile(savename, "RECREATE");
    savefile->cd();
    list->Write();
    savefile->Close();
 
    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
