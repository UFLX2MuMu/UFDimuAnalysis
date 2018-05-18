#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"

#include "EventTools.h"
#include "PUTools.h"
#include "SignificanceMetrics.hxx"

#include "TSystem.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TROOT.h"

#include "SampleDatabase.cxx"
#include "ThreadPool.hxx"

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
  std::map<TString, Sample*> sampleMap;
  std::vector<Sample*> sampleVec;
  //GetSamples(sampleMap, "CERN", "ZJets_MG");

  // Get DY AMC sample so that we can compare that to MG as well
  GetSamples(sampleMap, "CERN", "ZJets_AMC");

  int nthreads = 4;        // number of threads to use in parallelization
  float luminosity = 36814;
  float reductionFactor=1;
  std::map<TString, TH1D*> histMap;

  TString xname = "nJetsCent";
  int xbins = 11;
  float xmin = 0;
  float xmax = 11;

  // use pf, roch, or kamu values for selections, categories, and fill?
  TString pf_roch_or_kamu = "PF";
  if(xname.Contains("PF")) pf_roch_or_kamu = "PF";
  else if(xname.Contains("Roch")) pf_roch_or_kamu = "Roch";
  else if(xname.Contains("KaMu")) pf_roch_or_kamu = "KaMu";

  for(auto &s : sampleMap)
  {    
      //if(s.second->sampleType != "signal" && s.second->name != "ZJets_AMC" && s.second->name != "tt_ll_AMC" && s.second->sampleType != "data") continue;
      // Output some info about the current file
      std::cout << "  /// Using sample " << s.second->name << std::endl;
      std::cout << std::endl;
      std::cout << "    sample name:       " << s.second->name << std::endl;
      std::cout << "    sample file:       " << s.second->filenames[0] << std::endl;
      std::cout << "    pileup file:       " << s.second->pileupfile << std::endl;
      std::cout << "    nOriginal:         " << s.second->nOriginal << std::endl;
      std::cout << "    N:                 " << s.second->N << std::endl;
      std::cout << "    nOriginalWeighted: " << s.second->nOriginalWeighted << std::endl;
      std::cout << std::endl;

      s.second->setBranchAddresses(""); // tell the ttree to load the variable values into sample->vars
      sampleVec.push_back(s.second);   // add the sample to the vector of samples which we will run over

      // maps sample to histogram
      histMap[s.second->name] = new TH1D(s.second->name, s.second->name+Form(";%s;", xname.Data()), xbins, xmin, xmax);
  }    

  std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
  std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
  std::cout << "@@@ nSamples used  : " << sampleVec.size() << std::endl;


  auto makeHistoForSample = [&histMap, xname, pf_roch_or_kamu, reductionFactor, luminosity](Sample* s)
  {
    std::cout << Form("  /// Processing %s \n", s->name.Data());
    Run2MuonSelectionCuts  run2MuonSelection;
    Run2EventSelectionCuts run2EventSelection;

    JetCollectionCleaner      jetCollectionCleaner;
    MuonCollectionCleaner     muonCollectionCleaner;
    EleCollectionCleaner      eleCollectionCleaner;


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
          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu; 
          MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

          s->vars.setCalibrationType(pf_roch_or_kamu); // set mass calibrations

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

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

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // Load the rest of the information needed for run2 categories
          s->branches.getEntry(i);
          s->vars.setCalibrationType(pf_roch_or_kamu); // reloaded the branches, need to set mass,pt to correct calibrations again

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

          double varvalue = -999;
          varvalue = s->vars.getValue(xname.Data()); 
          histMap[s->name]->Fill(varvalue, s->getWeight());

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop //
    } // end event loop //

    histMap[s->name]->Scale(s->getLumiScaleFactor(luminosity));

    std::cout << Form("  /// Done processing %s \n", s->name.Data());
    return 0;
  };

  ///////////////////////////////////////////////////////////////////
  // PARALLELIZE BY SAMPLE -----------------------------------------
  ///////////////////////////////////////////////////////////////////

  ThreadPool pool(nthreads);
  std::vector< std::future<int> > results;

  TStopwatch timerWatch;
  timerWatch.Start();

  for(auto &s : sampleVec)
      results.push_back(pool.enqueue(makeHistoForSample, s));

  for(auto && result: results)
      std::cout << result.get();

  std::cout << std::endl;

  TList* ht = new TList();
  TList* stacklist = new TList();
  for(auto& h: histMap)
  {
      std::cout << Form("%s: %f\n", h.first.Data(), h.second->Integral());
  }

  ht->Add(histMap["ZJets_AMC_0j"]);
  ht->Add(histMap["ZJets_AMC_1j"]);
  ht->Add(histMap["ZJets_AMC_2j"]);

  TH1D* hNetHT   = DiMuPlottingSystem::addHists(ht, "Net_DY_Jets", "Net_DY_Jets");

  ht->Add(histMap["ZJets_AMC"]);
  std::cout << std::endl;

  TCanvas* stack = DiMuPlottingSystem::stackedHistogramsAndRatio(ht, "DY_Compare", "DY_Compare", xname, "Num Entries", true, false, "Inclusive/Stitched");
  stacklist->Add(stack);

  TString savename = Form("rootfiles/compare_drell_yan_%s.root", xname.Data());
  std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;

  TFile* f = new TFile(savename, "RECREATE");
  f->cd();

  stacklist->Write();
  hNetHT->Write();
  ht->Write();

  f->Close();

  timerWatch.Stop();
  std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;
  return 0;
}
