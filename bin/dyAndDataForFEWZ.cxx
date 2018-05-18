/////////////////////////////////////////////////////////////////////////////
//                            dyAndDataForFEWZ.cxx                         //
//=========================================================================//
//                                                                         //
// Make drell yan and data histograms with the FEWZ selection so that      //
// we can overlay and compare to FEWZ later. FEWZ histos are parsed        //
// from txt to TH1D in ../python/fewz.                                     //
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

// set bins, min, max, and the var to plot based upon terminal varNumber: varNumber and binning
void initPlotSettings(int varNumber, int binning, int& bins, float& min, float& max, TString& varname)
{
    // dimu_mass
    if(varNumber <= 0)
    {
        if(binning == 0)
        {
            bins = 150;
            min = 50;
            max = 200;
        }
        else if(binning == -1)
        {
            bins = 50;
            min = 110;
            max = 160;
        }
        else if(binning == 1)
        {
            bins = 50;
            min = 110;
            max = 160;
        }
        else if(binning == -2)
        {
            bins = 200;
            min = 110;
            max = 310;
        }
        else if(binning == 2)
        {
            bins = 100;
            min = 110;
            max = 310;
        }
        else
        {
            bins = 150;
            min = 50;
            max = 200;
        }

        if(varNumber == 0)  varname = "dimu_mass_PF";
        if(varNumber == -1) varname = "dimu_mass_Roch";
        if(varNumber == -2) varname = "dimu_mass_KaMu";
        if(varNumber == -3) varname = "dimu_mass_gen";
    }

    // dimu_pt 
    if(varNumber == 1)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt";
    }

    // dimu_pt_PF 
    if(varNumber == 100)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt";
    }

    // dimu_pt_Roch 
    if(varNumber == 101)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt_Roch";
    }

    // dimu_pt_KaMu 
    if(varNumber == 102)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt_KaMu";
    }

    // dimu_pt_gen 
    if(varNumber == 103)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt_gen";
    }

    // mu_pt
    if(varNumber == 2)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt";
    }

    // mu_pt_PF
    if(varNumber == 200)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt_PF";
    }

    // mu_pt_Roch
    if(varNumber == 201)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt_Roch";
    }
 
    // mu_pt_KaMu
    if(varNumber == 202)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt_KaMu";
    }
 
    // mu_pt_gen
    if(varNumber == 203)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt_gen";
    }
 
    // mu_eta
    if(varNumber == 3)
    {
        bins = 100;
        min = -2.5;
        max = 2.5;
        varname = "mu_eta";
    }

    // mu_eta_gen
    if(varNumber == 300)
    {
        bins = 100;
        min = -2.5;
        max = 2.5;
        varname = "mu_eta_gen";
    }

    // NPV
    if(varNumber == 4)
    {
        bins = 50;
        min = 0;
        max = 50;
        varname = "NPV";
    }

    // jet_pt
    if(varNumber == 5)
    {   
        bins = 200;
        min = 0;
        max = 200;
        varname = "jet_pt";
    }   

    // jet_eta 
    if(varNumber == 6)
    {   
        bins = 100;
        min = -5; 
        max = 5;
        varname = "jet_eta";
    }   

    // nValJets
    if(varNumber == 7)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "nValJets";
    }   

    // m_jj
    if(varNumber == 8)
    {   
        bins = 200;
        min = 0; 
        max = 2000;
        varname = "m_jj";
    }   

    // dEta_jj
    if(varNumber == 9)
    {   
        bins = 100;
        min = -10; 
        max = 10;
        varname = "dEta_jj";
    }   
    // N_valid_muons
    if(varNumber == 10) 
    {   
        bins = 11; 
        min = 0;  
        max = 11; 
        varname = "N_valid_muons";
    }   

    // nExtraMu
    if(varNumber == 11) 
    {   
        bins = 11; 
        min = 0;  
        max = 11; 
        varname = "nExtraMu";
    }   

    // extra_muon_pt
    if(varNumber == 12) 
    {   
        bins = 200;
        min = 0;
        max = 150;
        varname = "extra_muon_pt";
    }   
 
    // extra_muon_eta
    if(varNumber == 13)
    {
        bins = 100;
        min = -3;
        max = 3;
        varname = "extra_muon_eta";
    }

    // nEle
    if(varNumber == 14)
    {
        bins = 11;
        min = 0;
        max = 11;
        varname = "nEle";
    }
    // electron_pt
    if(varNumber == 15)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "electron_pt";
    }

    // electron_eta
    if(varNumber == 16)
    {
        bins = 100;
        min = -3;
        max = 3;
        varname = "electron_eta";
    }

    // nExtraLep
    if(varNumber == 17)
    {
        bins = 11;
        min = 0;
        max = 11;
        varname = "nExtraLep";
    }

    // nValBTags
    if(varNumber == 18)
    {
        bins = 11;
        min = 0;
        max = 11;
        varname = "nValBTags";
    }

    // bjet_pt
    if(varNumber == 19)
    {
        bins = 200;
        min = 0;
        max = 200;
        varname = "bjet_pt";
    }
    // bjet_eta 
    if(varNumber == 20)
    {
        bins = 100;
        min = -5;
        max = 5;
        varname = "bjet_eta";
    }

    // m_bb
    if(varNumber == 21)
    {
        bins = 200;
        min = 0;
        max = 2000;
        varname = "m_bb";
    }

    // mT_b_MET
    if(varNumber == 22)
    {
        bins = 200;
        min = 0;
        max = 2000;
        varname = "mT_b_MET";
    }

    // MHT
    if(varNumber == 23)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "MHT";
    }

    // dEta_jj_mumu
    if(varNumber == 24)
    {
        bins = 100;
        min = -10;
        max = 10;
        varname = "dEta_jj_mumu";
    }
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

    int varNumber = 0;        // the variable to plot, 0 is dimu_mass for instance
    int binning = 0;          // binning = 1 -> plot dimu_mass from 110 to 160 for limit setting
    int genSet = 0;           // value determines which values to use for cuts, and filling: gen or reco
    bool rebin = true;        // rebin the histograms so that the ratio plots have small errors
    int nthreads = 10;        // number of threads to use in parallelization

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> varNumber;
        if(i==2) ss >> binning;
        if(i==3) ss >> genSet;
        if(i==3) ss >> rebin;
        if(i==4) ss >> nthreads;
    }   
    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    float luminosity = 36814;       // pb-1
    float reductionFactor = 1;      // reduce the number of events you run over in case you want to debug or some such thing

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
    
    for(auto &i : samples)
    {
        // only want to plot Drell Yan and data here
        if(i.second->name == "ZJets_AMC" || i.second->sampleType == "data") ;
        else continue;

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

        i.second->setBranchAddresses(2);    // tell the TTree to load the variable info into s->vars
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    ///////////////////////////////////////////////////////////////////
    // Get Histo Information from Input -------------------------------
    ///////////////////////////////////////////////////////////////////

    // histogram settings for the variable we are going to plot in each sample x category
    TString varname;
    int bins;
    float min;
    float max;

    // set nbins, min, max, and var to plot based upon terminal varNumbers: varNumber and binning
    initPlotSettings(varNumber, binning, bins, min, max, varname);

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "min         : " << min << std::endl;
    std::cout << "max         : " << max << std::endl;
    std::cout << "bins        : " << bins << std::endl;
    std::cout << "binning     : " << binning << std::endl;
    std::cout << "genSet      : " << genSet << std::endl;
    std::cout << "rebin       : " << rebin << std::endl;
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////
    // Define Task for Parallelization -------------------------------
    ///////////////////////////////////////////////////////////////////

    auto makeHistoForSample = [varNumber, varname, binning, bins, min, max, rebin, luminosity, reductionFactor](Sample* s)
    {
      bool isblinded = true;
      if(binning < 0) isblinded = false;

      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      ///////////////////////////////////////////////////////////////////
      // INIT Cuts and Categories ---------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      FEWZCompareCuts        eventSelection;
      Run2MuonSelectionCuts  muonSelection;

      Categorizer* categorySelection = 0;
      categorySelection = new CategorySelectionFEWZ();

      // set some flags
      bool isData = s->sampleType.EqualTo("data");

      // use pf, roch, or kamu values for selections, categories, and fill?
      TString pf_roch_or_kamu = "";
      if(varname.Contains("PF")) pf_roch_or_kamu = "PF";
      else if(varname.Contains("Roch")) pf_roch_or_kamu = "Roch";
      else if(varname.Contains("KaMu")) pf_roch_or_kamu = "KaMu";

      ///////////////////////////////////////////////////////////////////
      // INIT HISTOGRAMS TO FILL ----------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram to fill in the category
      TString hkey = s->name;

      // Different categories for the analysis
      // categoryMap maps a TString to a Category object, one key,value pair for each category
      // in the categorySelection Categorizer object
      for(auto &c : categorySelection->categoryMap)
      {
          //number of bins for the histogram
          int hbins = bins;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;

          // Set up the histogram for the category and variable to plot
          // Each Category has a histoMap that maps a TString to a TH1D*
          // it also has TLists to keep data, sig, and bkg organized separately
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(varname);
          c.second.histoList->Add(c.second.histoMap[hkey]);                                      // need them ordered by xsec for the stack and ratio plot
          if(s->sampleType.EqualTo("data")) c.second.dataList->Add(c.second.histoMap[hkey]);     // data histo
          if(s->sampleType.EqualTo("signal")) c.second.signalList->Add(c.second.histoMap[hkey]); // signal histos
          if(s->sampleType.EqualTo("background")) c.second.bkgList->Add(c.second.histoMap[hkey]);// bkg histos
      }


      ///////////////////////////////////////////////////////////////////
      // LOOP OVER EVENTS -----------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // filter the events into the different categories to make our plots
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
          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu; 
          MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

          // Selection cuts and categories use standard values e.g. mu.pt
          // to use PF, Roch, or KaMu for cuts and categories set the pt calibration type 
          s->vars.setCalibrationType(pf_roch_or_kamu);

          ///////////////////////////////////////////////////////////////////
          // APPLY CUTS  ----------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          if(!eventSelection.evaluate(s->vars))
          { 
              continue; 
          }
          if(!mu1.isMediumID || !mu2.isMediumID)
          { 
              continue; 
          }
          if(!muonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          // avoid double counting for RunF
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

          // Load the rest of the information
          s->branches.jets->GetEntry(i);
          s->branches.mht->GetEntry(i);
          s->branches.nVertices->GetEntry(i);
          s->branches.electrons->GetEntry(i);
          //eventInfoBranch->GetEntry(i);

          if(!isData)
          {
              s->branches.nPU->GetEntry(i);
              s->branches.getEntryWeightsMC(i); // gen_wgt, pu_wgt and hlt, iso, mu_id scale factors
              s->branches.genParents->GetEntry(i);
              s->branches.genMuons->GetEntry(i);
              s->branches.genDimuons->GetEntry(i);
          }

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

          //std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

          ///////////////////////////////////////////////////////////////////
          // CATEGORIZE THE EVENT  ------------------------------------------
          ///////////////////////////////////////////////////////////////////

          // Figure out which category the event belongs to
          categorySelection->evaluate(s->vars);

          // Fill the histogram for the variable we are plotting for each sample x category
          for(auto &c : categorySelection->categoryMap)
          {
              // skip categories that we decided to hide (usually some intermediate categories)
              if(c.second.hide) continue;
              if(!c.second.inCategory) continue;

              ///////////////////////////////////////////////////////////////////
              // FILL HISTOGRAM APPROPRIATELY  ----------------------------------
              ///////////////////////////////////////////////////////////////////

              if(varname.Contains("dimu_mass")) 
              {
                  float varvalue = dimu.mass;
                  if(isData && varvalue > 120 && varvalue < 130 && isblinded) continue; // blind signal region

                  // if the event is in the current category then fill the category's histogram for the given sample and variable
                  c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
                  //std::cout << "    " << c.first << ": " << varvalue;
                  continue;
              }

              if(varname.Contains("dimu_pt"))
              {
                  // if the event is in the current category then fill the category's histogram for the given sample and variable
                  c.second.histoMap[hkey]->Fill(dimu.pt, s->getWeight());
                  continue;
              }
              // mu_pt is a substring of dimu_pt so we need the else if
              else if(varname.Contains("mu_pt"))
              {
                  c.second.histoMap[hkey]->Fill(mu1.pt, s->getWeight());
                  c.second.histoMap[hkey]->Fill(mu2.pt, s->getWeight());
                  continue;
              }

              // recoMu_Eta
              if(varname.EqualTo("mu_eta"))
              {
                  c.second.histoMap[hkey]->Fill(mu1.eta, s->getWeight());
                  c.second.histoMap[hkey]->Fill(mu2.eta, s->getWeight());
                  continue;
              }

              // NPV
              if(varname.EqualTo("NPV"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.nVertices, s->getWeight());
                   continue;
              }

              // jet_pt
              if(varname.EqualTo("jet_pt"))
              {
                   for(auto& jet: s->vars.validJets)
                       c.second.histoMap[hkey]->Fill(jet.Pt(), s->getWeight());
                   continue;
              }

              // jet_eta
              if(varname.EqualTo("jet_eta"))
              {
                   for(auto& jet: s->vars.validJets)
                       c.second.histoMap[hkey]->Fill(jet.Eta(), s->getWeight());
                   continue;
              }

              // nValJets
              if(varname.EqualTo("nValJets"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.validJets.size(), s->getWeight());
                   continue;
              }

              // m_jj
              if(varname.EqualTo("m_jj"))
              {
                   if(s->vars.validJets.size() >= 2)
                   {
                       TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                       c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                   }
                   continue;
              }

              // dEta_jj
              if(varname.EqualTo("dEta_jj"))
              {
                   if(s->vars.validJets.size() >= 2)
                   {
                       float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
                       c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                   }
                   continue;
              }

              // N_valid_muons
              if(varname.EqualTo("N_valid_muons"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.validMuons.size(), s->getWeight());
                   continue;
              }

              // nExtraMu
              if(varname.EqualTo("nExtraMu"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons.size(), s->getWeight());
                   continue;
              }

              // extra_muon_pt
              if(varname.EqualTo("extra_muon_pt"))
              {
                  for(auto& mu: s->vars.validExtraMuons)
                      c.second.histoMap[hkey]->Fill(mu.Pt(), s->getWeight());
                  continue;
              }

              // extra_muon_eta
              if(varname.EqualTo("extra_muon_eta"))
              {
                  for(auto& mu: s->vars.validExtraMuons)
                      c.second.histoMap[hkey]->Fill(mu.Eta(), s->getWeight());
                  continue;
              }

              // nEle
              if(varname.EqualTo("nEle"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.validElectrons.size(), s->getWeight());
                   continue;
              }

              // electron_pt
              if(varname.EqualTo("electron_pt"))
              {
                  for(auto& e: s->vars.validElectrons)
                      c.second.histoMap[hkey]->Fill(e.Pt(), s->getWeight());
                  continue;
              }
              // electron_eta
              if(varname.EqualTo("electron_eta"))
              {
                  for(auto& e: s->vars.validElectrons)
                      c.second.histoMap[hkey]->Fill(e.Eta(), s->getWeight());
                  continue;
              }

              // nExtraLep
              if(varname.EqualTo("nExtraLep"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.validElectrons.size() + s->vars.validExtraMuons.size(), s->getWeight());
                   continue;
              }

              // nValBTags
              if(varname.EqualTo("nValBTags"))
              {
                   c.second.histoMap[hkey]->Fill(s->vars.validBJets.size(), s->getWeight());
                   continue;
              }

              // bjet_pt
              if(varname.EqualTo("bjet_pt"))
              {
                  for(auto& bjet: s->vars.validBJets)
                      c.second.histoMap[hkey]->Fill(bjet.Pt(), s->getWeight());
                  continue;
              }

              // bjet_eta 
              if(varname.EqualTo("bjet_eta"))
              {
                  for(auto& bjet: s->vars.validBJets)
                      c.second.histoMap[hkey]->Fill(bjet.Eta(), s->getWeight());
                  continue;
              }
              // m_bb
              if(varname.EqualTo("m_bb"))
              {
                   if(s->vars.validBJets.size() >= 2)
                   {
                       TLorentzVector dijet = s->vars.validBJets[0] + s->vars.validBJets[1];
                       c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                   }
                   continue;
              }

              // mT_b_MET
              //if(varname.EqualTo("mT_b_MET"))
              //{
              //     if(s->vars.validBJets.size() > 0)
              //     {
              //         TLorentzVector mht(s->vars.mht.px, s->vars.mht.py, 0, s->vars.met.sumEt);
              //         TLorentzVector bjet = s->vars.validBJets[0];
              //         TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
              //         TLorentzVector bmet_t = met + bjet_t;

              //         c.second.histoMap[hkey]->Fill(bmet_t.M(), s->getWeight());
              //     }
              //     continue;
              //}

              // MHT
              if(varname.EqualTo("MHT"))
              {
                  c.second.histoMap[hkey]->Fill(s->vars.mht->pt, s->getWeight());
              }

              // dEta_jj_mumu
              if(varname.EqualTo("dEta_jj_mumu"))
              {
                   if(s->vars.validJets.size() >= 2)
                   {
                       TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                       float dEta = dijet.Eta() - s->vars.dimuCand->eta;
                       c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                   }
                   continue;
              }

          } // end category loop

          if(false)
            // ouput pt, mass info etc for the event
            EventTools::outputEvent(s->vars, *categorySelection);

          // Reset the flags in preparation for the next event
          categorySelection->reset();

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop
      } // end event loop

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection->categoryMap)
      {
          c.second.histoMap[hkey]->Scale(s->getLumiScaleFactor(luminosity));
      }

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return categorySelection;

    }; // done defining makeHistoForSample

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(nthreads);
    std::vector< std::future<Categorizer*> > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makeHistoForSample, s));

   ///////////////////////////////////////////////////////////////////
   // Gather all the Histos into one Categorizer----------------------
   ///////////////////////////////////////////////////////////////////

    Categorizer* cAll = 0;
    cAll = new CategorySelectionFEWZ();

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

                if(s->sampleType.EqualTo("signal"))          cAll->categoryMap[category.first].signalList->Add(h.second);
                else if(s->sampleType.EqualTo("background")) cAll->categoryMap[category.first].bkgList->Add(h.second);
                else                                         cAll->categoryMap[category.first].dataList->Add(h.second);
            }
        }
    }

   ///////////////////////////////////////////////////////////////////
   // Save All of the Histos-----------------------------------------
   ///////////////////////////////////////////////////////////////////

    TList* varstacklist = new TList();   // list to save all of the stacks
    TList* signallist = new TList();     // list to save all of the signal histos
    TList* bglist = new TList();         // list to save all of the background histos
    TList* datalist = new TList();       // list to save all of the data histos
    TList* netlist = new TList();        // list to save all of the net histos

    std::cout << "About to loop through cAll" << std::endl;
    for(auto &c : cAll->categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;

        // we need the data histo, the net signal, and the net bkg dimu mass histos for the datacards
        TH1D* hNetBkg    = dps->addHists(c.second.bkgList,    c.first+"_Net_Bkg",    "Net Background");
        TH1D* hNetData   = dps->addHists(c.second.dataList,   c.first+"_Net_Data",   "Data");

        TList* stacklist = (TList*)c.second.bkgList->Clone();
        stacklist->Add(hNetData);

        // Create the stack and ratio plot    
        TString cname = c.first+"_stack";
        //stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin = false, bool fit = true,
                                  //TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
        // stack signal, bkg, and data
        TCanvas* stack = dps->stackedHistogramsAndRatio(stacklist, cname, cname, varname, "Num Entries", rebin);
        varstacklist->Add(stack);

        // lists will contain signal, bg, and data histos for every category
        bglist->Add(c.second.bkgList);
        datalist->Add(c.second.dataList);

        netlist->Add(hNetBkg);
        netlist->Add(hNetData);
       
        stack->SaveAs("imgs/"+varname+"_"+cname+".png");
    }
    std::cout << std::endl;

    bool isblinded = binning >= 0; // binning == -1 means that we want dimu_mass from 110-160 unblinded
    TString blinded = "blinded";
    if(!isblinded) blinded = "UNBLINDED";

    if(varname.Contains("dimu_mass")) varname=blinded+"_"+varname;
    TString savename = Form("rootfiles/dy_data_fewz_%s__%s_%d_%d_%d.root", blinded.Data(), varname.Data(), (int)min, 
                            (int)max, (int)luminosity);

    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;

    TFile* savefile = new TFile(savename, "RECREATE");

    TDirectory* stacks        = savefile->mkdir("stacks");
    TDirectory* signal_histos = savefile->mkdir("signal_histos");
    TDirectory* bg_histos     = savefile->mkdir("bg_histos");
    TDirectory* data_histos   = savefile->mkdir("data_histos");
    TDirectory* net_histos    = savefile->mkdir("net_histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    stacks->cd();
    varstacklist->Write();

    signal_histos->cd();
    signallist->Write();

    bg_histos->cd();
    bglist->Write();

    data_histos->cd();
    datalist->Write();

    net_histos->cd();
    netlist->Write();

    savefile->Close();
 
    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
