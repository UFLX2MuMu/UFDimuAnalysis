/////////////////////////////////////////////////////////////////////////////
//                           outputForCategorizer.cxx                      //
//=========================================================================//
//                                                                         //
// output the events in a small window around the higgs mass = 125 GeV     //
// and their variables so that we can train a categorizer to make our      //
// categories for us.                                                      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


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
#include "TMVATools.h"
#include "PUTools.h"
#include "ThreadPool.hxx"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TRandom3.h"
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

    TString whichDY = "dyAMC-J";
    //TString whichDY = "dyMG";
    GetSamples(samples, "UF", "ALL_"+whichDY);

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
        //if(i.second->sampleType == "data") continue;
        //if(i.second->name != "H2Mu_gg") continue;

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

        i.second->setBranchAddresses("");
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    auto outputSampleInfo = [whichDY, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      TString dir    = "classification/";
      //TString methodName = "BDTG_default";
      TString methodName = "BDTG_UF_v1";

      // classification and multiclassification
      //TString weightfile = dir+"f_Opt1_all_sig_all_bkg_ge0j_BDTG_default.weights.xml";                       // >= 0j, use as inclusive
      //TString weightfile = dir+"f_Opt3_oneClass_all_sig_all_bkg_eq2j_eq0b_met80_BDTG_default.weights.xml";   // == 2j, 0b, met<80
      TString weightfile = dir+"f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";           // TMVA binary classification 1/2 sig for training
      //TString weightfile = dir+"binaryclass_amc.weights.xml";           // TMVA binary classification 1/2 sig for training, dy_amc@nlo

      //TString weightfile_multi = dir+"f_Opt2_all_sig_all_bkg_ge0j_eq0b_BDTG_default.weights.xml";              // assumes 0b jets, but use as inclusive
      //TString weightfile_multi = dir+"f_Opt3_all_sig_all_bkg_eq2j_eq0b_met80_BDTG_default.weights.xml";        // == 2j, 0b, met < 80
      //TString weightfile_multi = dir+"f_Opt3_half_all_sig_all_bkg_ge0j_eq0b_met80_BDTG_default.weights.xml";   // >= 0 jets, 0b, met < 80
      //TString weightfile_multi = dir+"f_Opt_v1_multi_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";   // TMVA multiclass 1/2 sig for training


      /////////////////////////////////////////////////////
      // Book training and spectator vars into reader

      std::map<TString, Float_t> tmap;
      std::map<TString, Float_t> smap;
      TMVA::Reader* reader = TMVATools::bookVars(methodName, weightfile, tmap, smap);

      //std::map<TString, Float_t> tmap_multi;
      //std::map<TString, Float_t> smap_multi;
      //TMVA::Reader* reader_multi = TMVATools::bookVars(methodName, weightfile_multi, tmap_multi, smap_multi);

      bool isData = s->sampleType == "data";
      bool isSignal = s->sampleType == "signal";

      TRandom3 r;

      int bin = -1;
      float interval = 0.5;
      int nbins = 20;
      float massmin = 120;
      float massmax = 130;

      //TFile* tfile = new TFile("rootfiles/bdt/"+s->name+"_bdt_training.root", "RECREATE");
      //tfile->cd();
      //TNtuple* ntuple = new TNtuple("theNtuple", "theNtuple");

      std::map<TString, double> vars;
      vars["bin"] = -999;
      vars["is_signal"] = -999;
      vars["weight"] = -999;

      // loop through all features in the varset (Sample->vars)
      // and initialize values
      for(auto& item: s->vars.varMap)
          vars[item.first.c_str()] = -999;

      for(auto& item: s->vars.varMapI)
          vars[item.first.c_str()] = -999;

      // !!!! output first line of csv to file
      std::ofstream file(Form("csv/bdtcsv/%s_bdt_training_%s.csv", s->name.Data(), whichDY.Data()), std::ofstream::out);
      file << EventTools::outputMapKeysCSV(vars).Data() << std::endl;

      // ntuple requires list of variables to be separated by ":" rather than ","
      TString ntuplevars = EventTools::outputMapKeysCSV(vars).ReplaceAll(",", ":");
      //std::cout << Form("\n   /// Vars to be saved in ntuple \n\n");
      //std::cout << Form("%s\n", ntuplevars.Data());
      TNtuple* ntuple = new TNtuple(s->name.Data(), s->name.Data(), ntuplevars.Data());

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelection;

      Categorizer* categorySelection = new CategorySelectionRun1();

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        // We are stitching together zjets_ht from 70-inf. We use the inclusive for
        // ht from 0-70.
        if(!isData)
        {
            s->branches.lhe_ht->GetEntry(i);
            if(s->name == "ZJets_MG" && s->vars.lhe_ht >= 70) continue;
        }

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

          // set pt calibration type: PF, Roch, or KaMu
          s->vars.setCalibrationType("Roch");

          MuonInfo& mu1 = s->vars.muons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.muons->at(dimu.iMu2);

          // only use the odd signal events for training
          // as to separate training and evaluation events
          if(isSignal && (s->vars.eventInfo->event % 2 == 0))
          {
              continue;
          }

          // only train on events that are in the mass window
          if(!(dimu.mass > 110 && dimu.mass < 160))
          {
              continue;
          }
          // only keep data in the sidebands
          if(isData && dimu.mass > massmin && dimu.mass < massmax)
          {
              continue;
          }

          // usual cuts
          if(!mu1.isMediumID || !mu2.isMediumID)
          { 
              continue; 
          }
          if(!run2EventSelection.evaluate(s->vars))
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
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

          // Load the rest of the information needed
          s->branches.getEntry(i);

          // set pt calibration type: PF, Roch, or KaMu again
          s->vars.setCalibrationType("Roch");

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
          CollectionCleaner::cleanByDR(s->vars.validBJets, s->vars.validMuons, 0.4);
          CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
          //CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // bin = -1 for events outside signal region, 
          // [0,nbins) for events inside signal region
          if(dimu.mass < massmin || dimu.mass >= massmax) bin = -1;
          else
          {
              float diff = dimu.mass - massmin;
              bin =  diff/interval;
          }

          // nonfeature variables
          vars["bin"] = bin;
          if(isData) vars["is_signal"] = -1;
          else vars["is_signal"] = s->sampleType.Contains("signal")?1:0;
          if(isData) vars["weight"] = 1;
          else vars["weight"] = s->getLumiScaleFactor(luminosity)*s->getWeight();

          // only using half of the signal events. need to boost the weight to get the correct estimate
          // of the expected signal
          if(isSignal) vars["weight"] = vars["weight"]*2;

          // set tmva's bdt_score
          s->vars.bdt_out = TMVATools::getClassifierScore(reader, methodName, tmap, s->vars);

          // load multi results into varset
          //std::vector<float> bdt_multi_scores = TMVATools::getMulticlassScores(reader_multi, methodName, tmap_multi, s->vars);
          //s->vars.bdt_ggh_out = bdt_multi_scores[0];
          //s->vars.bdt_vbf_out = bdt_multi_scores[1];
          //s->vars.bdt_vh_out  = bdt_multi_scores[2];
          //s->vars.bdt_ewk_out = bdt_multi_scores[3];
          //s->vars.bdt_top_out = bdt_multi_scores[4];

          s->vars.bdt_ggh_out = -999;
          s->vars.bdt_vbf_out = -999;
          s->vars.bdt_vh_out  = -999;
          s->vars.bdt_ewk_out = -999;
          s->vars.bdt_top_out = -999;

          //std::cout << i << " !!! SETTING JETS " << std::endl;
          //s->vars.setJets();    // jets sorted and paired by mjj, turn this off to simply take the leading two jets
          s->vars.setVBFjets();   // jets sorted and paired by vbf criteria

          // put the actual value for each feature into our map
          // then output the map to CSV later
          for(auto& item: s->vars.varMap)
          {
              // item.first is the name of the feature
              //std::cout << item.first.c_str() << std::endl;
              vars[item.first.c_str()] = s->vars.getValue(item.first);
          }
          for(auto& item: s->vars.varMapI)
          {
              // item.first is the name of the feature
              //std::cout << item.first.c_str() << std::endl;
              vars[item.first.c_str()] = s->vars.getValue(item.first);
          }

          if(false)
            EventTools::outputEvent(s->vars, (*categorySelection));

          // !!!! output event info to file
          file << EventTools::outputMapValuesCSV(vars).Data() << std::endl;

          std::vector<Float_t> varvalues;
          //std::cout << Form("\n   /// Vars to be saved \n\n");
          for(auto& item: vars)
          {
              //std::cout << Form("%s: %f\n", item.first.Data(), item.second);
              varvalues.push_back((Float_t)item.second);
          }

          // fill the ntuple
          ntuple->Fill(&varvalues[0]);

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimucand loop
      } // end event loop
      delete reader;
      //delete reader_multi;
      file.close();

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return ntuple;
    }; // end sample lambda function

    ThreadPool pool(nthreads);
    std::vector< std::future<TNtuple*> > results;

    for(auto& s: samplevec)
        results.push_back(pool.enqueue(outputSampleInfo, s));

    std::vector<TNtuple*> ntuples;
    for(auto && result: results)  // get each ntuple and save it to a tfile
    {
        TNtuple* ntuple = result.get();
        ntuples.push_back(ntuple);
    }

    for(auto& ntuple: ntuples)
    {
        TString sname = ntuple->GetName();
        ntuple->SetName("theNtuple");
        ntuple->SetTitle("theNtuple");
        TString filename = Form("rootfiles/bdt/%s_bdt_training_%s.root", sname.Data(), whichDY.Data());
        std::cout << Form("  /// Writing ntuple to %s \n", filename.Data());
        TFile* tfile = new TFile(filename,"RECREATE");
        tfile->cd();
        ntuple->Write();
        tfile->Write();
        tfile->Close();
    }
}
