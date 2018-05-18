/////////////////////////////////////////////////////////////////////////////
//                           fakes.cxx                                     //
//=========================================================================//
//                                                                         //
// Study non prompt muons by looking at same sign pairs.                   //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "SampleDatabase.cxx"
#include "ThreadPool.hxx"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"

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
    gROOT->SetBatch();
    int varNumber = 0;
    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> varNumber;
    }   

    float nthreads = 7;
    float reductionFactor = 1;
    std::map<TString, Sample*> samples;
    std::vector<Sample*> samplevec;
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    ///////////////////////////////////////////////////////////////////
    // Plot Settings Depending on Input--------------------------------
    ///////////////////////////////////////////////////////////////////
    
    TString massname = "dimu_mass";
    float massmin = 60;
    float massmax = 200;
    int massbins = 140;

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    GetSamples(samples, "CERN");

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING: SetBranchAddresses-------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    // Add only data samples to the vector we loop over to categorize

    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);

    for(auto &i : samples)
    {
        if(i.second->sampleType != "data") continue;
        //if(i.second->name == "RunH") continue;

        // Output some info about the current file
        std::cout << "  /// Using sample " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        i.second->setBranchAddresses(1);
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    //std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    
    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "massmin      : " << massmin << std::endl;
    std::cout << "massmax      : " << massmax << std::endl;
    std::cout << "massbins     : " << massbins << std::endl;
    std::cout << std::endl;
  
    auto makePlotsForSample = [massname, massbins, massmin, massmax, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      // cuts to apply to the events in the sample
      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelection;

      // turn off charge matching
      run2EventSelection.cutset.cuts[0].on = false;

      // turn off isolation
      run2MuonSelection.cutset.cuts[2].on = false;
      run2MuonSelection.cutset.cuts[5].on = false;

      TH1D* tight_histo  = new TH1D(s->name+"_tight_muon_id", "tight muon id", massbins, massmin, massmax);
      TH1D* medium_histo = new TH1D(s->name+"_medium_muon_id", "loose muon id", massbins, massmin, massmax);

      std::map<TString, TH1D*> returnMap;
      returnMap[s->name+"_tight_muon_id"] = tight_histo;
      returnMap[s->name+"_medium_muon_id"] = medium_histo;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {
        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////
        
        // only load essential information for the first set of cuts 
        s->branches.muPairs->GetEntry(i);
        s->branches.muons->GetEntry(i);
        s->branches.eventInfo->GetEntry(i);

        // require exactly two muons
        if(s->vars.muons->size() != 2) continue;
        if(s->vars.muPairs->size() < 1) continue;
        bool found_good_dimuon = false;

        // find the first good dimuon candidate and fill info
        for(auto& dimu: (*s->vars.muPairs))
        {
          // set the appropriate dimuon for the cuts
          s->vars.dimuCand = &dimu;
          MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////
          
          // require same sign muons to check out fakes
          if(mu1.charge != mu2.charge)
          { 
              continue; 
          }
          // require that the muons pass the basic run2 selections except for the ones we turned off
          // (opposite sign and isolation are off)
          if(!run2EventSelection.evaluate(s->vars))
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }
          // at minimum require medium isolation
          if(mu1.iso() > 0.12 || mu2.iso() > 0.12)
          { 
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


          found_good_dimuon = true;

          // if it is a medium muon and outside the blinded region then add to the medium histogram
          if(mu1.isMediumID && mu2.isMediumID)
          { 
              if(dimu.mass_PF > 140 || dimu.mass_PF < 110) medium_histo->Fill(dimu.mass_PF);
          }

          // if it is a tight muon with appropriate isolation and outside the blinded region then add to the tight histogram
          if(mu1.isTightID && mu2.isTightID)
          { 
              if(mu1.iso() > 0.12 || mu2.iso() > 0.12) ;
              else if(dimu.mass_PF > 140 || dimu.mass_PF < 110) tight_histo->Fill(dimu.mass_PF);
              
          }
          if(found_good_dimuon) break;

        } // end dimucand loop
      } // end events loop
      
      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return returnMap;

    };

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(nthreads);
    std::vector< std::future< std::map<TString, TH1D*> > > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makePlotsForSample, s));

    // get the tight and medium histograms from the different runs so we can add them together
    TList* tight_list = new TList();
    TList* medium_list = new TList();
    for(auto && result: results)
    {
        for(auto& item: result.get())
        {
            if(item.first.Contains("tight")) tight_list->Add(item.second);
            else if(item.first.Contains("medium")) medium_list->Add(item.second);
        }
    }

    // ad them together
    TH1D* net_tight_histo = DiMuPlottingSystem::addHists(tight_list, "Tight_Muon_ID_Data", "Tight ID");
    TH1D* net_medium_histo = DiMuPlottingSystem::addHists(medium_list, "Medium_Muon_ID_Data", "Medium ID");

    // Create the stack and ratio plot    
    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/fakes_"+massname+"_data_8_0_X.root", "RECREATE");
    savefile->cd();

    net_tight_histo->Write();
    net_medium_histo->Write();

    savefile->Write();
    savefile->Close();

    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
