// output the events in a small window around the higgs mass = 125 GeV and their variables
// so that we can train a bdt to make our categories for us

// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "MuonSelection.h"
#include "EventSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"
#include "SampleDatabase.cxx"

#include "EventTools.h"
#include "PUTools.h"
#include "ThreadPool.hxx"

#include "TLorentzVector.h"
#include "TSystem.h"
#include <sstream>
#include <map>
#include <vector>
#include <utility>

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
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    int nthreads = 10;

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        ss >> nthreads;
    }   

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float reductionFactor = 1;
    float luminosity = 36814;      // pb-1
    GetSamples(samples, "CERN");

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto &i : samples)
    {
        if(i.second->sampleType == "data") continue;

        // Output some info about the current file
        std::cout << "  /// Preprocessing " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        i.second->setBranchAddresses(2);
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    auto outputSampleInfo = [luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      int bin = -1;
      float interval = 0.5;
      int nbins = 20;
      float massmin = 120;
      float massmax = 130;

      std::map<TString, double> vars;
      vars["bin"] = -999;
      vars["is_signal"] = -999;
      vars["weight"] = -999;

      // loop through all features in the varset (Sample->vars)
      // and initialize values
      for(auto& item: s->vars.varMap)
          vars[item.first.c_str()] = -999;

      // !!!! output first line of csv to file
      std::ofstream file("csv/"+s->name+"_bdt_training.csv", std::ofstream::out);
      file << EventTools::outputMapKeysCSV(vars).Data() << std::endl;

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts     run2MuonSelection;
      Run2EventSelectionCuts80X run2EventSelectionMC;

      Categorizer* categorySelection = new CategorySelectionRun1();

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {
        // only load essential information for the first set of cuts 
        s->branches.muPairs->GetEntry(i);
        s->branches.muons->GetEntry(i);

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
          MuonInfo& mu1 = s->vars.muons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.muons->at(dimu.iMu2);

          // only want to train on events that are in the higgs mass window
          if(!(dimu.mass_PF > 110 && dimu.mass_PF < 160))
          {
              continue;
          }
          // usual cuts
          if(!mu1.isTightID || !mu2.isTightID)
          { 
              continue; 
          }
          if(!run2EventSelectionMC.evaluate(s->vars))
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          // Load the rest of the information needed for run2 categories
          s->branches.jets->GetEntry(i);
          s->branches.mht->GetEntry(i);
          s->branches.nVertices->GetEntry(i);
          s->branches.electrons->GetEntry(i);

          s->branches.gen_wgt->GetEntry(i);
          s->branches.nPU->GetEntry(i);
          s->branches.pu_wgt->GetEntry(i);
          s->branches.eff_wgt->GetEntry(i);

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
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.3);
          CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.3);
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.3);

          // low stats in these categories. don't include them in the training set
          categorySelection->evaluate(s->vars);
          if(categorySelection->categoryMap["c_2_Jet_VBF_Tight"].inCategory) 
          {
              //EventTools::outputEvent(s->vars, (*categorySelection));
              continue;
          }
          if(categorySelection->categoryMap["c_2_Jet_GGF_Tight"].inCategory)
          {
              //EventTools::outputEvent(s->vars, (*categorySelection));
              continue;
          }

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // bin = -1 for events outside signal region, 
          // [0,nbins) for events inside signal region
          if(dimu.mass_PF < massmin || dimu.mass_PF >= massmax) bin = -1;
          else
          {
              float diff = dimu.mass_PF - massmin;
              bin =  diff/interval;
          }

          // nonfeature variables
          vars["bin"] = bin;
          vars["is_signal"] = s->sampleType.Contains("signal")?1:0;
          vars["weight"] = s->getScaleFactor(luminosity)*s->getWeight();

          // put the actual value for each feature into our map
          // then output the map to CSV later
          for(auto& item: s->vars.varMap)
              vars[item.first.c_str()] = s->vars.getValue(item.first);
/*
          // muon features
          vars["dimu_pt"] = s->vars.getValue("dimu_pt");
          vars["mu0_pt"]  = s->vars.getValue("mu0_pt");
          vars["mu1_pt"]  = s->vars.getValue("mu1_pt");
          vars["mu0_eta"] = s->vars.getValue("mu0_eta");
          vars["mu1_eta"] = s->vars.getValue("mu1_eta");

          // jet variables
          vars["jj_jet0_pt"]      = s->vars.getValue("jj_jet0_pt");
          vars["jj_jet1_pt"]      = s->vars.getValue("jj_jet1_pt");
          vars["jj_jet0_eta"]     = s->vars.getValue("jj_jet0_eta");
          vars["jj_jet1_eta"]     = s->vars.getValue("jj_jet1_eta"); 

          vars["m_jj"]         = s->vars.getValue("m_jj");
          vars["dEta_jj"]      = s->vars.getValue("dEta_jj");
          vars["dEta_jj_mumu"] = s->vars.getValue("dEta_jj_mumu");

          vars["bjet0_pt"]  = s->vars.getValue("bjet0_pt");  
          vars["bjet1_pt"]  = s->vars.getValue("bjet1_pt"); 
          vars["bjet0_eta"] = s->vars.getValue("bjet0_eta"); 
          vars["bjet1_eta"] = s->vars.getValue("bjet1_eta"); 

          vars["m_bb"]    = s->vars.getValue("m_bb"); 
          vars["dEta_bb"] = s->vars.getValue("dEta_bb"); 

          vars["nValJets"]          = s->vars.getValue("nValJets"); 
          vars["nValBTags"]         = s->vars.getValue("nValBTags"); 
          vars["nExtraMu"]   = s->vars.getValue("nExtraMu"); 
          vars["nEle"]     = s->vars.getValue("nEle"); 
          vars["nExtraLep"] = s->vars.getValue("nExtraLep"); 

          vars["MET"] = s->vars.getValue("MET"); 

          vars["extra_muon0_pt"]  = s->vars.getValue("extra_muon0_pt"); 
          vars["extra_muon1_pt"]  = s->vars.getValue("extra_muon1_pt"); 
          vars["extra_muon0_eta"] = s->vars.getValue("extra_muon0_eta"); 
          vars["extra_muon1_eta"] = s->vars.getValue("extra_muon1_eta"); 

          vars["electron0_pt"]  = s->vars.getValue("electron0_pt");  
          vars["electron1_pt"]  = s->vars.getValue("electron1_pt"); 
          vars["electron0_eta"] = s->vars.getValue("electron0_eta"); 
          vars["electron1_eta"] = s->vars.getValue("electron1_eta"); 

          vars["mT_b_MET"] = s->vars.getValue("mT_b_MET");

          vars["zep"]      = s->vars.getValue("zep");
          vars["dimu_dPhiStar"] = s->vars.getValue("dimu_dPhiStar");
          vars["dPhi"]     = s->vars.getValue("dPhi");
*/
          if(false)
            EventTools::outputEvent(s->vars, (*categorySelection));

          // !!!! output event info to file
          file << EventTools::outputMapValuesCSV(vars).Data() << std::endl;

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimucand loop
      } // end event loop

      file.close();
      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return 0;
    }; // end sample lambda function

    ThreadPool pool(nthreads);
    std::vector< std::future<int> > results;

    for(auto& s: samplevec)
        results.push_back(pool.enqueue(outputSampleInfo, s));
}
