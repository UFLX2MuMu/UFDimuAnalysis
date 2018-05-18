/////////////////////////////////////////////////////////////////////////////
//                           outputCounts.cxx                              //
//=========================================================================//
//                                                                         //
// output # signal, # background, S/B  and significance for                //
// the different categories to CSV. Useful to evaluate the goodness        //
// of categorization. Using it to compare Decision Tree created categories //
// to the Run1 categories at the moment. Works with different significance //
// metrics in SignificanceMetrics.hxx.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#include "Sample.h"
#include "SignificanceMetrics.hxx"
#include "DiMuPlottingSystem.h"
#include "MuonSelection.h"
#include "EventSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"
#include "SampleDatabase.cxx"

#include "TMVATools.h"
#include "EventTools.h"
#include "PUTools.h"
#include "ThreadPool.hxx"

#include "TLorentzVector.h"
#include "TSystem.h"
#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include <iomanip>

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

std::pair<int,int> fwhm(TH1D* h, int start=2)
{
    int maxbin = -999;
    int leftbin = -999;
    int rightbin = -999;
    double max = -999;
    double half_max = -999; 

    // find max
    for(int i=start; i<=h->GetNbinsX(); i++)
    {
        double val = h->GetBinContent(i);
        if(val > max)
        {
            max = val;
            maxbin = i;
        } 
    }

    half_max = max/2;

    // find half-max bin on left of max
    for(int i=start; i<maxbin; i++)
    {
        double val = h->GetBinContent(i);
        if(val >= half_max )
        {
            leftbin = i;
            break;
        } 
    }

    // find half-max bin on right of max
    for(int i=maxbin; i<=h->GetNbinsX(); i++)
    {
        double val = h->GetBinContent(i);
        if(val <= half_max)
        {
            rightbin = i;
            break;
        } 
    }

    double left = h->GetBinCenter(leftbin);
    double right = h->GetBinCenter(rightbin);
    double full_width_half_max = right - left;

    return std::pair<int, int>(leftbin, rightbin);
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

Categorizer* getHistos(TString xmlfile, int nthreads, float reductionFactor, float luminosity, TString whichDY)
{
    std::map<TString, Sample*> samples;
    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    GetSamples(samples, "CERN", "ALL_"+whichDY);

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

    auto outputSampleInfo = [xmlfile, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      bool isData = s->sampleType == "data";

      float interval = 0.5;
      int nbins = 20;
      float massmin = 120;
      float massmax = 130;


      /////////////////////////////////////////////////////
      // Load TMVA classifiers
    
      TString dir    = "classification/";
      //TString methodName = "BDTG_default";
      TString methodName = "BDTG_UF_v1";


      // sig vs bkg and multiclass (ggf, vbf, ... drell yan, ttbar) weight files
      TString weightfile = dir+"f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";
      //TString weightfile = dir+"binaryclass_amc.weights.xml";

      /////////////////////////////////////////////////////
      // Book training and spectator vars into reader

      TMVA::Reader* reader = 0;
      std::map<TString, Float_t> tmap;
      std::map<TString, Float_t> smap;

      //TMVA::Reader* reader_multi = 0;
      //std::map<TString, Float_t> tmap_multi;
      //std::map<TString, Float_t> smap_multi;

      // load tmva binary classification and multiclass classifiers
      reader       = TMVATools::bookVars(methodName, weightfile, tmap, smap);
      //reader_multi = TMVATools::bookVars(methodName, weightfile_multi, tmap_multi, smap_multi);

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelection;

      Categorizer* categorySelection = 0;
      if(xmlfile.Contains(".xml")) categorySelection = new XMLCategorizer(xmlfile); 
      else categorySelection = new CategorySelectionRun1();

      TString hkeyn;
      TString hkeyw;
      for(auto &c : categorySelection->categoryMap)
      {   
            //number of bins for the histogram
            int hbins = 21;

            // c.second is the category object, c.first is the category name
            TString hnamen = c.first+"_"+s->name+"_num";
            TString hnamew = c.first+"_"+s->name+"_weights";
            hkeyn = s->name+"_num";
            hkeyw = s->name+"_weights";

            // The number of events histo
            c.second.histoMap[hkeyn] = new TH1D(hnamen, hnamen, hbins, -1, 20);
            c.second.histoMap[hkeyn]->GetXaxis()->SetTitle("bin");
            c.second.histoList->Add(c.second.histoMap[hkeyn]);                                        
            if(s->sampleType.Contains("data")) c.second.dataList->Add(c.second.histoMap[hkeyn]);      
            if(s->sampleType.Contains("signal")) c.second.signalList->Add(c.second.histoMap[hkeyn]);  
            if(s->sampleType.Contains("background")) c.second.bkgList->Add(c.second.histoMap[hkeyn]); 

            // The sum of weights histo
            c.second.histoMap[hkeyw] = new TH1D(hnamew, hnamew, hbins, -1, 20);
            c.second.histoMap[hkeyw]->GetXaxis()->SetTitle("bin");
            c.second.histoList->Add(c.second.histoMap[hkeyw]);                                        
            if(s->sampleType.Contains("data")) c.second.dataList->Add(c.second.histoMap[hkeyw]);      
            if(s->sampleType.Contains("signal")) c.second.signalList->Add(c.second.histoMap[hkeyw]);  
            if(s->sampleType.Contains("background")) c.second.bkgList->Add(c.second.histoMap[hkeyw]); 
      }   

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
          MuonInfo& mu1 = s->vars.muons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.muons->at(dimu.iMu2);

          // only want to train on events that are in the higgs mass window
          if(!(dimu.mass_PF > 110 && dimu.mass_PF < 160))
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

          // avoid double counting for RunF
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

          //std::cout << i << " !!! SETTING JETS " << std::endl;
          //s->vars.setJets();    // jets sorted and paired by mjj, turn this off to simply take the leading two jets
          s->vars.setVBFjets();   // jets sorted and paired by vbf criteria
          s->vars.bdt_out = TMVATools::getClassifierScore(reader, methodName, tmap, s->vars); // set tmva's bdt score

          categorySelection->evaluate(s->vars);

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // bin = -1 for events outside signal region, 
          // [0,nbins) for events inside signal region
          int bin = -1;
          if(dimu.mass_PF < massmin || dimu.mass_PF >= massmax) bin = -1;
          else
          {
              float diff = dimu.mass_PF - massmin;
              bin =  diff/interval;
          }

          // Look at each category
          for(auto &c : categorySelection->categoryMap)
          {
              // if the event is in the current category then fill the category's histogram for the given sample and variable
              if(c.second.inCategory) c.second.histoMap[hkeyn]->Fill(bin);
              if(c.second.inCategory) c.second.histoMap[hkeyw]->Fill(bin, s->getWeight());
          } // end category loop

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimucand loop
      } // end event loop

      delete reader;

      for(auto& c : categorySelection->categoryMap)
      {
          for(auto& h: c.second.histoMap)
          {
              if(h.first.Contains("weights") && !isData)
              {
                  h.second->Scale(s->getLumiScaleFactor(luminosity));
                  //if(s->sampleType == "signal") std::cout << Form("    %s, %s, %f\n", c.first.Data(), h.first.Data(), h.second->Integral());
              }
          }
      }

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      delete s;
      return categorySelection;
    }; // end sample lambda function

    ThreadPool pool(nthreads);
    std::vector< std::future<Categorizer*> > results;

    for(auto& s: samplevec)
        results.push_back(pool.enqueue(outputSampleInfo, s));

   ///////////////////////////////////////////////////////////////////
   // Gather all the Histos into one Categorizer----------------------
   ///////////////////////////////////////////////////////////////////

    Categorizer* cAll = 0;
    if(xmlfile.Contains(".xml")) cAll = new XMLCategorizer(xmlfile);
    else cAll = new CategorySelectionRun1();

    // get histos from all categorizers and put them into one
    for(auto && categorizer: results)  // loop through each Categorizer object, one per sample
    {
        for(auto& category: categorizer.get()->categoryMap) // loop through each category
        {
            // category.first is the category name, category.second is the Category object
            if(category.second.hide) continue;
            for(auto& h: category.second.histoMap) // loop through each histogram in the category
            {
                cAll->categoryMap[category.first].histoMap[h.first] = h.second;
                cAll->categoryMap[category.first].histoList->Add(h.second);

                if(h.first.Contains("H2Mu"))      cAll->categoryMap[category.first].signalList->Add(h.second);
                else if(h.first.Contains("Run"))  cAll->categoryMap[category.first].dataList->Add(h.second);
                else                              cAll->categoryMap[category.first].bkgList->Add(h.second);
            }
        }
    }

    return cAll;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    TString xmlfile = "";
    int nthreads = 10;
    float reductionFactor = 1;
    float luminosity = 36814;      // pb-1
    //TString whichDY = "dyAMC";
    TString whichDY = "dyAMC-J";

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) xmlfile = ss.str().c_str();
        if(i==2) ss >> nthreads;
    }   

    TStopwatch timerWatch;
    timerWatch.Start();

    ///////////////////////////////////////////////////////////////////
    // Categorize Events, Make Histograms -----------------------------
    ///////////////////////////////////////////////////////////////////

    Categorizer* cAll = getHistos(xmlfile, nthreads, reductionFactor, luminosity, whichDY);

    ///////////////////////////////////////////////////////////////////
    // Make Tables, Save Histos ---------------------------------------
    ///////////////////////////////////////////////////////////////////

    TList* varstacklist = new TList();   // list to save all of the stacks
    TList* signallist = new TList();     // list to save all of the signal histos
    TList* bglist = new TList();         // list to save all of the background histos
    TList* datalist = new TList();       // list to save all of the data histos
    TList* netlist = new TList();        // list to save all of the net histos

    TString xcategoryString = "run1";
    if(xmlfile.Contains(".xml")) 
    {    
        Ssiz_t i = xmlfile.Last('/');
        xcategoryString = xmlfile(i+1, xmlfile.Length()); // get the name of the xmlfile without all the /path/to/dir/ business 
        xcategoryString = xcategoryString.ReplaceAll(".xml", ""); 
    }    
    xcategoryString = "_"+xcategoryString;

    TString csvfilename = Form("csv/sigcsv/significance%s_categories_%s.csv", xcategoryString.Data(), whichDY.Data());
    std::cout << std::endl << "  /// Exporting counts and significnace to " << csvfilename << " ..." << std::endl;
    std::cout << std::endl;

    std::ofstream file(csvfilename, std::ofstream::out);
    TString outtitles = "category,limit,s0,S_div_sqrtB,S_div_B,signal,signal_fwhm,background_fwhm";
    TString stdouttitles = outtitles;
    std::cout << stdouttitles << std::endl;

    file << outtitles.Data() << std::endl;

    // map of run1 limits
    std::map<TString, double> mlim;
    mlim["c_01_Jet_Loose_BB"] =  8.59;
    mlim["c_01_Jet_Loose_BE"] =  15.9;
    mlim["c_01_Jet_Loose_BO"] =  8.90;
    mlim["c_01_Jet_Loose_EE"] =  34.4;
    mlim["c_01_Jet_Loose_OE"] =  18.9;
    mlim["c_01_Jet_Loose_OO"] =  16.1;
    mlim["c_01_Jet_Tight_BB"] =  5.58;
    mlim["c_01_Jet_Tight_BE"] =  10.0;
    mlim["c_01_Jet_Tight_BO"] =  5.55;
    mlim["c_01_Jet_Tight_EE"] =  22.1;
    mlim["c_01_Jet_Tight_OE"] =  12.7;
    mlim["c_01_Jet_Tight_OO"] =  9.75;
    mlim["c_2_Jet_VBF_Tight"] =  6.7;
    mlim["c_2_Jet_GGF_Tight"] =  10.2;
    mlim["c_2_Jet_VBF_Loose"] =  8.72;
    mlim["c_ALL"]             =  2.13;

    // map of bdt category limits
    mlim["c0"] =  19.3;
    mlim["c1"] =  19.8;
    mlim["c2"] =  15.3;
    mlim["c3"] =  14.8;
    mlim["c4"] =  8.47;
    mlim["c5"] =  7.66;
    mlim["c6"] =  7.28;
    mlim["c7"] =  6.72;
    mlim["c8"] =  6.25;
    mlim["c9"] =  7.59;
    mlim["c10"] =  6.53;
    mlim["c11"] =  6.34;
    mlim["c12"] =  4.67;
    mlim["root"] = 1.85;
    SignificanceMetric* s0 = new PoissonSignificance(0, 10);

    double net_significance = 0;
    double net_s_sqrt_b = 0;

    // loop through categories to calculate significance, other info for the table of category info
    for(auto &c : cAll->categoryMap)
    {
        double lim = mlim[c.second.name];
        TString cname = c.second.name;
        double sig0 = -999;
        double signal = -999;
        double signal_fwhm = -999;
        double bkg_fwhm = -999;
        double s_over_b = -999;
        double s2_over_b = -999;

        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;

        TIter isig(c.second.signalList);
        TIter ibkg(c.second.bkgList);
        TIter idata(c.second.dataList);
        TObject* obj = 0;

        TList* signalListNum = new TList();
        TList* signalListW = new TList();
        TList* bkgListNum = new TList();
        TList* bkgListW = new TList();
        TList* dataListNum = new TList();
        TList* dataListW = new TList();

        // filter num/weight histos to appropriate lists
        TH1D* hvbf = 0;
        TH1D* hggf = 0;
        TH1D* hdy = 0;
        TH1D* httbar = 0;
        while(obj = isig())
        {
            TH1D* h = (TH1D*) obj;
            if(TString(h->GetName()).Contains("num")) signalListNum->Add(h);
            if(TString(h->GetName()).Contains("weights")) signalListW->Add(h);
        }
        while(obj = ibkg())
        {
            TH1D* h = (TH1D*) obj;
            if(TString(h->GetName()).Contains("num")) bkgListNum->Add(h);
            if(TString(h->GetName()).Contains("weights")) bkgListW->Add(h);
        }

        while(obj = idata())
        {
            TH1D* h = (TH1D*) obj;
            if(TString(h->GetName()).Contains("num")) dataListNum->Add(h);
            if(TString(h->GetName()).Contains("weights")) dataListW->Add(h);
        }

        // group the MC together
        TList* mcListNum = (TList*)signalListNum->Clone(); 
        TList* mcListW = (TList*)signalListW->Clone();

        TIter in(bkgListNum);
        while(obj = in())
        {
            TH1D* h = (TH1D*) obj;
            mcListNum->Add(h);
        }
        TIter iw(bkgListW);
        while(obj = iw())
        {
            TH1D* h = (TH1D*) obj;
            mcListW->Add(h);
        }
        TList* groupedListNum = DiMuPlottingSystem::groupMC(mcListNum, c.second.name);
        TList* groupedListW = DiMuPlottingSystem::groupMC(mcListW, c.second.name);

        // num histos
        TH1D* hNetSignalNum = DiMuPlottingSystem::addHists(signalListNum, c.second.name+"_Num_Signal", "Num Signal");
        TH1D* hNetBkgNum    = DiMuPlottingSystem::addHists(bkgListNum,    c.second.name+"_Num_Bkg",    "Num Background");
        TH1D* hNetDataNum   = DiMuPlottingSystem::addHists(dataListNum,    c.second.name+"_Num_Data",  "Num Data");

        // weighted histos
        TH1D* hNetSignalW = DiMuPlottingSystem::addHists(signalListW, c.second.name+"_Net_Signal", "Net Signal");
        TH1D* hNetBkgW    = DiMuPlottingSystem::addHists(bkgListW,    c.second.name+"_Net_Bkg",    "Net Background");
        TH1D* hNetDataW   = DiMuPlottingSystem::addHists(dataListW,   c.second.name+"_Net_Data",   "Net Data");

        netlist->Add(hNetSignalW);
        netlist->Add(hNetBkgW);
        netlist->Add(hNetDataW);
        netlist->Add(hNetSignalNum);
        netlist->Add(hNetBkgNum);
        netlist->Add(hNetDataNum);

        // we want the counts inside the signal region
        
        // Full Width Half Max interval for net signal
        std::pair<int, int> bounds = fwhm(hNetSignalW);
        int leftbin = bounds.first;
        int rightbin = bounds.second;

        // we want the counts inside the signal region
        signal      = hNetSignalW->Integral(1, hNetSignalW->GetNbinsX());
        signal_fwhm = hNetSignalW->Integral(leftbin, rightbin);
        bkg_fwhm = hNetBkgW->Integral(leftbin, rightbin);

        s_over_b = signal_fwhm/bkg_fwhm;
        s2_over_b = signal_fwhm*signal_fwhm/bkg_fwhm;
        sig0 = s0->significance2(hNetSignalW, hNetBkgW, hNetDataW, hNetSignalNum, hNetBkgNum, hNetDataNum);

        if(c.second.isTerminal) net_significance += sig0;
        if(c.second.isTerminal) net_s_sqrt_b += s2_over_b;

        TString outtitles = "category,limit,s0,S_div_sqrtB,S_div_B,signal,signal_fwhm,background_fwhm";
        TString outstring = Form("%s,%f,%f,%f,%f,%f,%f,%f", cname.Data(), lim, TMath::Sqrt(sig0), 
                                 TMath::Sqrt(s2_over_b), s_over_b, signal, signal_fwhm, bkg_fwhm);
        file << outstring.Data() << std::endl;

        

        TString stdoutstring = outstring;
        std::cout << stdoutstring << std::endl;
    }
    file.close();  
    net_significance = TMath::Sqrt(net_significance);
    net_s_sqrt_b = TMath::Sqrt(net_s_sqrt_b);

    std::cout << std::endl;
    std::cout << "  ////////////////////////////////////////// " << std::endl;
    std::cout << "  ### NET SIGNIFICANCE: " << net_significance << std::endl;
    std::cout << "  ### NET S/SQRT(B)   : " << net_s_sqrt_b << std::endl;
    std::cout << "  ////////////////////////////////////////// " << std::endl;

    TString savename = Form("rootfiles/significance%s_categories_%s.root", xcategoryString.Data(), whichDY.Data());

    std::cout << std::endl;
    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;

    TFile* savefile = new TFile(savename, "RECREATE");
    savefile->cd();
    netlist->Write();
    savefile->Close();

    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}

