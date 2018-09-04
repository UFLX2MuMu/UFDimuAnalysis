/////////////////////////////////////////////////////////////////////////////
//                           masscalibration.cxx                           //
//=========================================================================//
//                                                                         //
// plot dimu-mass mean and resolution in data and drell yan MC             //
// vs phi+, phi-, eta+, eta-. Used to investigate different pt corrections //
// Compare Particle Flow to KaMu and Rochester.                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "MassCalibration.h"
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

struct Settings
{
    TString peak = "Z";          // Peak to fit: Z, JPsi, Upsilon
    TString xname = "dimu_pt";  // xvar to plot against 
    TString muID = "medium";     // muon id to consider
    TString whichDY = "dyAMC";   // use drell yan AMC or Madgraph for DY MC
    TString whichHiggs = "SIGNALX";
 
    int nthreads = 10;           // # of threads to use when plotting
    float reductionFactor = 1;  // reduce the # of events to run over for debugging

    float ymean_min = 90.5;     // yrange for mean vs x plot
    float ymean_max = 91.5;
    float yres_min = 0.6;       // yrange for resolution vs x plot
    float yres_max = 1.8;

    float massmin = 82.2;       // the window around the peak to fit the voigtians
    float massmax = 100.2;
    int   massbins = 62;

    float xmin = -3.15;         // xmin,max,bins defaulted to phi
    float xmax = 3.15;
    int   xbins = 12;
    std::vector<Float_t> binning; // variable binning for xvar if we want it

    float fitsig = 1;           // width to fit around the mean of the distribution
    float voigt_gamma = 2.5;    // theoretical width of the peak

    float leadPtMin = 26;       // min lead muon pt
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void combineDataRunInfo(std::map<TString, MassCalibration*>& masscalmap)
{
    // get an item from the map, so that we have the masscal histo and plot settings
    MassCalibration* zfirst = 0;
    for(auto& zc: masscalmap)
        zfirst = zc.second;

    bool variableBinning = zfirst->variableBinning;

    std::cout << Form("  /// Running combineDataRunInfo \n");

    TString massname = "Net_Data_mass_";
    TString massname_pf = massname+"PF";
    TString massname_roch = massname+"Roch";
    TString massname_kamu = massname+"KaMu";

    MassCalibration* masscal_net_pf = 0;
    MassCalibration* masscal_net_roch = 0;
    MassCalibration* masscal_net_kamu = 0;

    if(!variableBinning)
    {
        masscal_net_pf = new MassCalibration(zfirst->xname, massname_pf, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->xmin, zfirst->xmax, zfirst->xbins, zfirst->voigt_gamma);

        masscal_net_roch = new MassCalibration(zfirst->xname, massname_roch, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->xmin, zfirst->xmax, zfirst->xbins, zfirst->voigt_gamma);

        masscal_net_kamu = new MassCalibration(zfirst->xname, massname_kamu, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->xmin, zfirst->xmax, zfirst->xbins, zfirst->voigt_gamma);
    }
    else
    {
        masscal_net_pf = new MassCalibration(zfirst->xname, massname_pf, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->binning, zfirst->voigt_gamma);

        masscal_net_roch = new MassCalibration(zfirst->xname, massname_roch, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->binning, zfirst->voigt_gamma);

        masscal_net_kamu = new MassCalibration(zfirst->xname, massname_kamu, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->binning, zfirst->voigt_gamma);
    }

    for(unsigned int i=0; i<zfirst->histos.size(); i++)
    {
        TList* add_hist_list_pf = new TList();
        TList* add_hist_list_roch = new TList();
        TList* add_hist_list_kamu = new TList();

	std::cout << Form("  i=%d \n",i);

        for(auto& zc: masscalmap)
        {
            if(!zc.second->massname.Contains("data")) continue;

	    if(zc.second->massname.Contains("PF")) add_hist_list_pf->Add(zc.second->histos[i]);
            else if(zc.second->massname.Contains("Roch")) add_hist_list_roch->Add(zc.second->histos[i]);
            else if(zc.second->massname.Contains("KaMu")) add_hist_list_kamu->Add(zc.second->histos[i]);
        }
 
        if(add_hist_list_pf->GetSize() != 0)
        {
	    std::cout << Form(" %d \n",add_hist_list_pf->GetSize());
            TH1D* net_hist_i_pf = DiMuPlottingSystem::addHists(add_hist_list_pf, masscal_net_pf->histos[i]->GetName(), masscal_net_pf->histos[i]->GetTitle()); 
            TH1D* temp_pf = masscal_net_pf->histos[i];
            masscal_net_pf->histos[i] = net_hist_i_pf;
            delete temp_pf;
        }

        if(add_hist_list_roch->GetSize() != 0)
        {
            TH1D* net_hist_i_roch = DiMuPlottingSystem::addHists(add_hist_list_roch, masscal_net_roch->histos[i]->GetName(), masscal_net_roch->histos[i]->GetTitle()); 
            TH1D* temp_roch = masscal_net_roch->histos[i];
            masscal_net_roch->histos[i] = net_hist_i_roch;
            delete temp_roch;
        }

        if(add_hist_list_kamu->GetSize() != 0)
        {
            TH1D* net_hist_i_kamu = DiMuPlottingSystem::addHists(add_hist_list_kamu, masscal_net_kamu->histos[i]->GetName(), masscal_net_kamu->histos[i]->GetTitle()); 
            TH1D* temp_kamu = masscal_net_kamu->histos[i];
            masscal_net_kamu->histos[i] = net_hist_i_kamu;
            delete temp_kamu;
        }
    }

    std::cout << Form("%d \n %d",masscal_net_pf->histos.size(),masscal_net_pf->histos[0]->Integral());

    if(masscal_net_pf->histos.size() > 0 ) masscalmap[massname_pf] = masscal_net_pf;
    if(masscal_net_roch->histos.size() > 0 ) masscalmap[massname_roch] = masscal_net_roch;
    if(masscal_net_kamu->histos.size() > 0 ) masscalmap[massname_kamu] = masscal_net_kamu;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void initPlotSettings(Settings& settings) 
{

    settings.ymean_min = 90.5;
    settings.ymean_max = 91.5;
    settings.yres_min = 0.6;
    settings.yres_max = 1.8;

    settings.massmin = 82.2;
    settings.massmax = 100.2;

    settings.massbins = 62;
    settings.xbins = 12;
    settings.fitsig = 1;
    settings.voigt_gamma = 2.5;

    if(settings.peak == "JPsi")
    {
        settings.massmin = 3.0;
        settings.massmax = 3.2;
        settings.ymean_min = 3.08;
        settings.ymean_max = 3.12;
        settings.yres_min = 0.033;
        settings.yres_max = 0.043;
        //settings.voigt_gamma = 0.0929;
        settings.voigt_gamma = 0;
        settings.fitsig = 1.5;
    }
    else if(settings.peak == "Upsilon")
    {
        settings.massmin = 9.3;
        settings.massmax = 9.6;
        settings.ymean_min = 9.40;
        settings.ymean_max = 9.55;
        settings.yres_min = 0.05;
        settings.yres_max = 0.20;
        //settings.voigt_gamma = 0.0540;
        settings.voigt_gamma = 0;
        settings.fitsig = 1.5;
    }

    if(settings.peak == "Higgs")
    {
	settings.whichDY = settings.whichHiggs;
        settings.massmin = 113;
	settings.massmax = 137;
	settings.ymean_min = 124;
	settings.ymean_max = 126;
	settings.yres_min = 0;
	settings.yres_max = 4;
        
	settings.fitsig =1;
	settings.voigt_gamma = 0;
    }

    if(settings.xname.Contains("Dphi"))
    {
        settings.xmin = -6.3;
        settings.xmax = 6.3;
        settings.xbins = 24;
    }

    else if(settings.xname.Contains("phi"))
    {
        settings.xmin = -3.15;
        settings.xmax = 3.15;
    }
    
    if(settings.xname.Contains("Deta"))
    {   
        settings.yres_max = 1.55*settings.yres_max;
        settings.xmin = -4.82;
        settings.xmax = 4.82;
	settings.xbins =24;
    }

    else if(settings.xname.Contains("eta"))
    {
        settings.yres_max = 1.55*settings.yres_max;
        if(settings.peak == "JPsi") settings.yres_max = 4.0*settings.yres_max;
        if(settings.peak == "Upsilon") settings.yres_max = 2.25*settings.yres_max;
        settings.xmin = -2.41;
        settings.xmax = 2.41;
    }

    if(settings.xname.Contains("dimu_pt"))
    {
        settings.ymean_min = 0.997*settings.ymean_min;
        settings.ymean_max = 1.007*settings.ymean_max;
//	settings.ymean_min = 90.6;
//	settings.ymean_max = 92.5;

        settings.yres_min = 0.33*settings.yres_min;
        settings.yres_max = 2.6*settings.yres_max;
        settings.xmin = 0;
        settings.xmax = 500;
//	settings.xmax = 120;
        settings.xbins=50;
        settings.binning = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120,130,140,150,200,500};
    }

    else if(settings.xname.Contains("pt"))
    {
        settings.ymean_min = 0.997*settings.ymean_min;
        settings.ymean_max = 1.007*settings.ymean_max;
        settings.yres_min = 0.33*settings.yres_min;
        settings.yres_max = 2.6*settings.yres_max;
        settings.xmin = 10;
        settings.xmax = 210;
        settings.xbins = 50;
        settings.binning = {10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 510};
    }

}

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
    Settings settings;

    gROOT->SetBatch();

    ///////////////////////////////////////////////////////////////
    // Parse Arguments -------------------------------------------
    ///////////////////////////////////////////////////////////////

    for(int i=1; i<argc; i++) 
    {    
        std::stringstream ss;  
        TString in = argv[i];
        TString option = in(0, in.First("="));
        option = option.ReplaceAll("--", ""); 
        TString value  = in(in.First("=")+1, in.Length());
        value = value.ReplaceAll("\"", ""); 
        ss << value.Data();

        if(option=="x")                    settings.xname = value;
        else if(option=="peak")            settings.peak = value;
        else if(option=="leadPtMin")       ss >> settings.leadPtMin;
        else if(option=="nthreads")        ss >> settings.nthreads;
        else if(option=="reductionFactor") ss >> settings.reductionFactor;
        else if(option=="whichDY")         settings.whichDY = value;
        else
        {
            std::cout << Form("!!! %s is not a recognized option.", option) << std::endl;
        }
    }

    std::map<TString, Sample*> samples;
    std::vector<Sample*> samplevec;
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    initPlotSettings(settings);

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    GetSamples(samples, "CERN", "CERN_2017_collection");  //+settings.whichDY

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
        if(!(i.second->sampleType == "data" || i.second->name.Contains("ZJets_AMC_2017") || i.second->name == "ZJets_MG" )) continue;

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
    std::cout << "@@@ nCPUs used     : " << settings.nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "peak           : " << settings.peak << std::endl;
    std::cout << "xvar           : " << settings.xname << std::endl;
    std::cout << "xmin           : " << settings.xmin << std::endl;
    std::cout << "xmax           : " << settings.xmax << std::endl;
    std::cout << "xbins          : " << settings.xbins << std::endl;
    std::cout << std::endl;
    std::cout << "massmin        : " << settings.massmin << std::endl;
    std::cout << "massmax        : " << settings.massmax << std::endl;
    std::cout << "massbins       : " << settings.massbins << std::endl;
    std::cout << std::endl;
    std::cout << "muID           : " << settings.muID << std::endl;
    std::cout << "reductionFactor: " << settings.reductionFactor << std::endl;
    std::cout << std::endl;

    for(auto& i: settings.binning)
        std::cout << i << ",";
    std::cout << std::endl << std::endl;


  
    auto makePlotsForSample = [settings](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      bool isData = s->sampleType == "data";
      int muID = 0;
      int varNumber = 0;
      
      // integers are faster to compare than strings
      if     (settings.muID == "tight")       muID = 0;
      else if(settings.muID == "medium")      muID = 1;
      else if(settings.muID == "mediumx")     muID = -1;
      else if(settings.muID == "medium2016x") muID = -2;
      else muID = 0;

      if     (settings.xname == "phi_plus")  varNumber = 0;
      else if(settings.xname == "phi_minus") varNumber = 1;
      else if(settings.xname == "eta_plus")  varNumber = 2;
      else if(settings.xname == "eta_minus") varNumber = 3;
      else if(settings.xname == "pt_plus")   varNumber = 4;
      else if(settings.xname == "pt_minus")  varNumber = 5;
      else if(settings.xname == "dimu_pt")   varNumber = 6;
      else if(settings.xname == "Dphi")	     varNumber = 7;
      else if(settings.xname == "Dphi_given")  varNumber = 8;
      else if(settings.xname == "Deta")	     varNumber = 9;
      else if(settings.xname == "Dphi_barrel") varNumber = 10;
      else if(settings.xname == "Dphi_overlap") varNumber = 11;
      else if(settings.xname == "Dphi_endcap") varNumber = 12;

      else varNumber = 0;

      // cuts to apply to the events in the sample
      MuonSelectionCuts_2017  run2MuonSelection;
      EventSelectionCuts_2017 run2EventSelection(settings.leadPtMin, settings.massmin*0.5);

//      Run2MuonSelectionCuts  run2MuonSelection(5 , 2.4 , 0.25);
//      Run2EventSelectionCuts run2EventSelection(settings.leadPtMin, settings.massmin*0.5);



      // use tight isolation
      // run2MuonSelection.cMaxRelIso = 0.25;

      // will probably want one for each type of mass
      MassCalibration* masscal_pf   = 0; 
      MassCalibration* masscal_roch =  0; 
//      MassCalibration* masscal_kamu =  0; 

      if(settings.binning.size()==0) 
      {
          masscal_pf = new MassCalibration(settings.xname, s->name+"_mass_PF",   settings.fitsig, settings.massmin, settings.massmax, settings.massbins, 
                                             settings.xmin, settings.xmax, settings.xbins, settings.voigt_gamma);

          masscal_roch = new MassCalibration(settings.xname, s->name+"_mass_Roch", settings.fitsig, settings.massmin, settings.massmax, settings.massbins, 
                                             settings.xmin, settings.xmax, settings.xbins, settings.voigt_gamma);

//          masscal_kamu = new MassCalibration(settings.xname, s->name+"_mass_KaMu", settings.fitsig, settings.massmin, settings.massmax, settings.massbins, 
//                                             settings.xmin, settings.xmax, settings.xbins, settings.voigt_gamma);
      }
      else
      {
          masscal_pf = new MassCalibration(settings.xname, s->name+"_mass_PF",   settings.fitsig, settings.massmin, settings.massmax, settings.massbins, 
                                           settings.binning, settings.voigt_gamma);

          masscal_roch = new MassCalibration(settings.xname, s->name+"_mass_Roch", settings.fitsig, settings.massmin, settings.massmax, settings.massbins, 
                                             settings.binning, settings.voigt_gamma);

//          masscal_kamu = new MassCalibration(settings.xname, s->name+"_mass_KaMu", settings.fitsig, settings.massmin, settings.massmax, settings.massbins, 
//                                             settings.binning, settings.voigt_gamma);
      }

      //std::cout << Form("  VGAMMA: pf = %f,  roch = %f, kamu = %f \n", masscal_pf->voigt_gamma, masscal_roch->voigt_gamma, masscal_kamu->voigt_gamma);

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      for(unsigned int i=0; i<s->N/settings.reductionFactor; i++)
      {
        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////
        
        // only load essential information for the first set of cuts 
        s->branches.muPairs->GetEntry(i);
        s->branches.muons->GetEntry(i);
        s->branches.eventInfo->GetEntry(i);
	
        // need a dimuon candidate to fill the zmass
        if(s->vars.muPairs->size() < 1) continue;
        if(s->vars.muons->size() < 2) continue;
        bool found_good_dimuon = false;

        // find the first good dimuon candidate and fill info
        for(auto& dimu: (*s->vars.muPairs))
//	for(int pairnumber=0 ; pairnumber < s->vars.muPairs->size() ; pairnumber++)
        {
//	  auto& dimu = s->vars.muPairs->at(pairnumber);
          s->vars.dimuCand = &dimu;
          MuonInfo& mu1 = s->vars.muons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.muons->at(dimu.iMu2);

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          if(dimu.mass < settings.massmin || dimu.mass > settings.massmax)
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

                // complicated muon ID selection
          bool isTight = false;
          bool isMedium = false;

          if(mu1.isTightID && mu2.isTightID)
              isTight = true; 

          if(mu1.isMediumID && mu2.isMediumID)
              isMedium = true;

          // tight
          if(muID == 0)
          {
              if(isTight) ;
              else continue;
          }
          // medium
          if(muID == 1)
          {
              if(isMedium) ;
              else continue;
          }
          // exclusive medium
          if(muID == -1)
          {
              if(!isTight && isMedium) ;
              else continue;
          }

          found_good_dimuon = true;

          // mc branches
          if(s->sampleType != "data")
          {
              s->branches.nPU->GetEntry(i);
              s->branches.getEntryWeightsMC(i); // gen_wgt, pu_wgt and hlt, iso, mu_id scale factors
          }

          if(varNumber == 0) 
          {
              float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
          
              masscal_pf->fill(phi_plus, dimu.mass, s->getWeight());
              masscal_roch->fill(phi_plus, dimu.mass_Roch, s->getWeight());
//              masscal_kamu->fill(phi_plus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 1) 
          {
              float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;

              masscal_pf->fill(phi_minus, dimu.mass, s->getWeight());
              masscal_roch->fill(phi_minus, dimu.mass_Roch, s->getWeight());
//              masscal_kamu->fill(phi_minus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 2) 
          {
              float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;
	      float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;
	      if(eta_minus > 1.2) continue;
	      if(eta_minus < -1.2) continue;


              masscal_pf->fill(eta_plus, dimu.mass, s->getWeight());
              masscal_roch->fill(eta_plus, dimu.mass_Roch, s->getWeight());
//              masscal_kamu->fill(eta_plus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 3) 
          {
	      float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;
              float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;
	      if(eta_plus < 1.2) continue;

              masscal_pf->fill(eta_minus, dimu.mass, s->getWeight());
              masscal_roch->fill(eta_minus, dimu.mass_Roch, s->getWeight());
//              masscal_kamu->fill(eta_minus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 4) 
          {
              float pt_plus_pf = (mu1.charge==1)?mu1.pt:mu2.pt;
              float pt_plus_roch = (mu1.charge==1)?mu1.pt_Roch:mu2.pt_Roch;
              float pt_plus_kamu = (mu1.charge==1)?mu1.pt_KaMu:mu2.pt_KaMu;

              masscal_pf->fill(pt_plus_pf, dimu.mass, s->getWeight());
              masscal_roch->fill(pt_plus_roch, dimu.mass_Roch, s->getWeight());
//              masscal_kamu->fill(pt_plus_kamu, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 5) 
          {
              float pt_minus_pf = (mu1.charge==-1)?mu1.pt:mu2.pt;
              float pt_minus_roch = (mu1.charge==-1)?mu1.pt_Roch:mu2.pt_Roch;
              float pt_minus_kamu = (mu1.charge==-1)?mu1.pt_KaMu:mu2.pt_KaMu;

              masscal_pf->fill(pt_minus_pf, dimu.mass, s->getWeight());
              masscal_roch->fill(pt_minus_roch, dimu.mass_Roch, s->getWeight());
//              masscal_kamu->fill(pt_minus_kamu, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 6) 
          {
              float dimu_pt_pf = dimu.pt;
              float dimu_pt_roch = dimu.pt_Roch;
//              float dimu_pt_kamu = dimu.pt_KaMu;
	      float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;
              float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;
//              if(std::abs(mu1.eta-mu2.eta) < 2.5) continue;

//	      if( i%10000 == 1)
//	      {
//		std::cout<< "sample " << s->name << std::endl;
//		std::cout<< "i = " << i << std::endl;
//		std::cout<< "dimu.mass = " << dimu.mass << ", " << dimu.mass_Roch << std::endl;
//		std::cout<< "dimu.pt = " << dimu_pt_pf << ", " << dimu_pt_roch << std::endl; 
//		std::cout<< "weight = " << s->getWeight() << std::endl;
//	      }


              masscal_pf->fill(dimu_pt_pf, dimu.mass, 1);   // s->getWeight()
              masscal_roch->fill(dimu_pt_roch, dimu.mass_Roch, 1);
//              masscal_kamu->fill(dimu_pt_kamu, dimu.mass_KaMu, s->getWeight());
          }
	
 	  if(varNumber == 7)
	  { 
	     float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
	     float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;
	     float dPhi = phi_plus - phi_minus;
	     if(dPhi > M_PI)  dPhi = dPhi - 2*M_PI;
	     if(dPhi <= -M_PI) dPhi = dPhi + 2*M_PI;
	     float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;
             float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;  
//	     if(std::abs(mu1.eta-mu2.eta) >2.5) continue;


	     masscal_pf->fill(dPhi, dimu.mass_PF, s->getWeight());
	     masscal_roch->fill(dPhi, dimu.mass_Roch, s->getWeight());
//	     masscal_kamu->fill(dPhi, dimu.mass_KaMu, s->getWeight());
	  }

	  if(varNumber == 8)
	  {
	     masscal_pf->fill(dimu.dPhi, dimu.mass, s->getWeight());
             masscal_roch->fill(dimu.dPhi, dimu.mass_Roch, s->getWeight());
//             masscal_kamu->fill(dimu.dPhi, dimu.mass_KaMu, s->getWeight());
	  }

	  if(varNumber == 9)
          {
             float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;
             float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;

             masscal_pf->fill(eta_plus - eta_minus, dimu.mass, s->getWeight());
             masscal_roch->fill(eta_plus - eta_minus, dimu.mass_Roch, s->getWeight());
//             masscal_kamu->fill(eta_plus - eta_minus, dimu.mass_KaMu, s->getWeight());
          }


	  if(varNumber == 10)
          {
             float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
             float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;

	     if(std::abs(mu1.eta)>0.9 || std::abs(mu2.eta)>0.9) continue;

             masscal_pf->fill(phi_plus - phi_minus, dimu.mass_PF, s->getWeight());
             masscal_roch->fill(phi_plus - phi_minus, dimu.mass_Roch, s->getWeight());
//             masscal_kamu->fill(phi_plus - phi_minus, dimu.mass_KaMu, s->getWeight());
          }

	  if(varNumber == 11)
          {
             float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
             float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;

             if(std::abs(mu1.eta)<0.9 || std::abs(mu1.eta)>1.2 || std::abs(mu2.eta)<0.9 || std::abs(mu2.eta)>1.2) continue;

             masscal_pf->fill(phi_plus - phi_minus, dimu.mass_PF, s->getWeight());
             masscal_roch->fill(phi_plus - phi_minus, dimu.mass_Roch, s->getWeight());
//             masscal_kamu->fill(phi_plus - phi_minus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 12)
          {
             float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
             float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;

             if(std::abs(mu1.eta)<1.2 || std::abs(mu2.eta)<1.2) continue;

             masscal_pf->fill(phi_plus - phi_minus, dimu.mass_PF, s->getWeight());
             masscal_roch->fill(phi_plus - phi_minus, dimu.mass_Roch, s->getWeight());
//             masscal_kamu->fill(phi_plus - phi_minus, dimu.mass_KaMu, s->getWeight());
          }



          if(found_good_dimuon) break;

        } // end dimucand loop
      } // end events loop
      
      std::vector<MassCalibration*>* returnVector = new std::vector<MassCalibration*>();
      returnVector->push_back(masscal_pf);
      returnVector->push_back(masscal_roch);
//      returnVector->push_back(masscal_kamu);

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return returnVector;

    };

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(settings.nthreads);
    std::vector< std::future< std::vector<MassCalibration*>* > > results;

    TStopwatch timerWatch;
    timerWatch.Start();


    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makePlotsForSample, s));

    // put all of the masscalibration pointers into this map
    std::map<TString, MassCalibration*> masscalMap;  

    // Gather results for each sample
    for(auto && result: results)
    {
        // loop over each sample's PF, Roch, KaMu masscal
        for(auto masscal: (*result.get()))
        {
           masscalMap[masscal->massname] = masscal;
        }
    }


//now we have a map containing every masscal result (data samples, correction types, maybe tuple numbers)
//if different tuple files from one sample are taken as different samples, combine them here
    std::map<TString, bool> have_sum;
    std::array<TString, 6>  names_of_samples = {"ZJets_AMC","data_B","data_C","data_D","data_E","data_F"};
    for(int iname = 0 ; iname < names_of_samples.size() ; iname ++)
    {
	have_sum[names_of_samples[iname] + "_PF"] = false;
	have_sum[names_of_samples[iname] + "_Roch"] = false;
	have_sum[names_of_samples[iname] + "_KaMu"] = false;
    }
    std::map<TString, MassCalibration*> map_for_cal_sum;
 //   std::map<TString, TList*> histos_for_each_sum;
    for(auto& item: masscalMap)
    {
	TString masscal_key = item.first;
	MassCalibration* masscalib = item.second;
	MassCalibration* temp_masscal = 0;
	
	TString what_sample_it_from = "no";
	TString what_calib = "nay";

	if( !masscal_key.Contains("tuple") ) continue;
	
	if( masscal_key.Contains("ZJets") ) what_sample_it_from = "ZJets_AMC";
	else if ( masscal_key.Contains("data_B") ) what_sample_it_from = "data_B";
	else if ( masscal_key.Contains("data_C") ) what_sample_it_from = "data_C";
	else if ( masscal_key.Contains("data_D") ) what_sample_it_from = "data_D";
	else if ( masscal_key.Contains("data_E") ) what_sample_it_from = "data_E";
	else if ( masscal_key.Contains("data_F") ) what_sample_it_from = "data_F";

	if( masscal_key.Contains("PF") ) what_calib = "_PF";
	else if( masscal_key.Contains("Roch") ) what_calib = "_Roch";
	else if( masscal_key.Contains("KaMu") ) what_calib = "_KaMu";

	if( have_sum[what_sample_it_from + what_calib] == false )
	{
		if(masscalib->variableBinning)
			temp_masscal = new MassCalibration(masscalib->xname, "sum_" + what_sample_it_from + "_2017_mass" + what_calib, masscalib->fitsig,
                                         	   	   masscalib->massmin, masscalib->massmax, masscalib->massbins,
                                         	   	   masscalib->binning, masscalib->voigt_gamma);
		else
			temp_masscal = new MassCalibration(masscalib->xname, "sum_" + what_sample_it_from + "_2017_mass" + what_calib, masscalib->fitsig, 
                                         		   masscalib->massmin, masscalib->massmax, masscalib->massbins, 
                                         		   masscalib->xmin, masscalib->xmax, masscalib->xbins, masscalib->voigt_gamma);

		for(int ihis = 0; ihis < masscalib->histos.size() ; ihis++)
		{
			temp_masscal->histos[ihis] = (TH1D*) masscalib->histos[ihis]->Clone();
		}
		masscalMap[ "sum_" + what_sample_it_from + "_2017_mass" + what_calib ] = temp_masscal;
		have_sum[what_sample_it_from + what_calib] = true;
	}

	else
	{
		temp_masscal = masscalMap["sum_" + what_sample_it_from + "_2017_mass" + what_calib ];
		for(int ihis = 0; ihis < masscalib->histos.size() ; ihis++)
		{
			temp_masscal->histos[ihis]->Add(masscalib->histos[ihis],1.0); 
		}
	}
	delete masscalib;
	masscalMap.erase(masscal_key);
    }
//summation of different files in each sample complete






    // add up the individual data runs into net data and add to the map
    combineDataRunInfo(masscalMap);


    // Now get the plots from the map and organize them into lists.
    // Then, we can save them into the different root file directories easily.
    TList* overlays_list = new TList();
    TList* histos_list = new TList();
    TList* plots_list = new TList();
    std::map<TString, TList*> overlayMap;

    for(auto& item: masscalMap)
    {
        TString masscal_key = item.first;
        MassCalibration* masscal = item.second;
        std::cout << masscal_key << ": " << masscal->histos[0]->GetName() << ": " << masscal->histos[1]->Integral() << std::endl;
        masscal->plot();

        TGraph* mean_graph = masscal->mean_vs_x;
        TGraph* res_graph = masscal->resolution_vs_x;
        
        for(auto& h: masscal->histos)
            histos_list->Add(h);

        plots_list->Add(mean_graph);
        plots_list->Add(res_graph);

        // strip extraneous naming information to initialize the overlay map
        TString overlay_name = masscal->massname;
        overlay_name = overlay_name.ReplaceAll("_mass_PF", "");
        overlay_name = overlay_name.ReplaceAll("_mass_KaMu", "");
        overlay_name = overlay_name.ReplaceAll("_mass_Roch", "");

        // Name by mass type only
        //if(masscal->massname.Contains("PF")) 
        //{        
        //    mean_graph->SetTitle("PF");
        //    res_graph->SetTitle("PF");
        //}
        //else if(masscal->massname.Contains("Roch")) 
        //{        
        //    mean_graph->SetTitle("Roch");
        //    res_graph->SetTitle("Roch");
        //}
        //else if(masscal->massname.Contains("KaMu")) 
        //{        
        //    mean_graph->SetTitle("KaMu");
        //    res_graph->SetTitle("KaMu");
        //}

        // name by sample_MASSTYPE 
        TString title = masscal_key;
        title = title.ReplaceAll("_mass", "")+"_"+settings.muID;
        title = title.ReplaceAll("-", "n");

        mean_graph->SetTitle(title);
        res_graph->SetTitle(title);
       
        // overlay pf, roch, kamu mean for sample
        if(overlayMap.count(overlay_name+"_Mean") == 0) 
        {
            overlayMap[overlay_name+"_Mean"] = new TList();
            overlayMap[overlay_name+"_Mean"]->Add(mean_graph);
        }
        else 
        {
            overlayMap[overlay_name+"_Mean"]->Add(mean_graph);
        }

        // overlay pf, roch, kamu res for the sample
        if(overlayMap.count(overlay_name+"_Resolution") == 0) 
        {
            overlayMap[overlay_name+"_Resolution"] = new TList();
            overlayMap[overlay_name+"_Resolution"]->Add(res_graph);
        }
        else 
        {
            overlayMap[overlay_name+"_Resolution"]->Add(res_graph);
        }

        // overlay data and mc for pf
        if(masscal_key.Contains("PF") && (masscal_key.Contains("Net_Data") || masscal_key.Contains("ZJets")))
        {
            if(overlayMap.count("Data_MC_Mean_PF") == 0) 
            {
                overlayMap["Data_MC_Mean_PF"] = new TList();
                overlayMap["Data_MC_Mean_PF"]->Add((TGraph*)mean_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Mean_PF"]->Add((TGraph*)mean_graph->Clone());
            }

            if(overlayMap.count("Data_MC_Resolution_PF") == 0) 
            {
                overlayMap["Data_MC_Resolution_PF"] = new TList();
                overlayMap["Data_MC_Resolution_PF"]->Add((TGraph*)res_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Resolution_PF"]->Add((TGraph*)res_graph->Clone());
            }
        }
        // overlay data and mc for roch
        if(masscal_key.Contains("Roch") && (masscal_key.Contains("Net_Data") || masscal_key.Contains("ZJets")))
        {
            if(overlayMap.count("Data_MC_Mean_Roch") == 0) 
            {
                overlayMap["Data_MC_Mean_Roch"] = new TList();
                overlayMap["Data_MC_Mean_Roch"]->Add((TGraph*)mean_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Mean_Roch"]->Add((TGraph*)mean_graph->Clone());
            }

            if(overlayMap.count("Data_MC_Resolution_Roch") == 0) 
            {
                overlayMap["Data_MC_Resolution_Roch"] = new TList();
                overlayMap["Data_MC_Resolution_Roch"]->Add((TGraph*)res_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Resolution_Roch"]->Add((TGraph*)res_graph->Clone());
            }
        }
        // overlay data and mc for kamu
        if(masscal_key.Contains("KaMu") && (masscal_key.Contains("Net_Data") || masscal_key.Contains("ZJets")))
        {
            if(overlayMap.count("Data_MC_Mean_KaMu") == 0) 
            {
                overlayMap["Data_MC_Mean_KaMu"] = new TList();
                overlayMap["Data_MC_Mean_KaMu"]->Add((TGraph*)mean_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Mean_KaMu"]->Add((TGraph*)mean_graph->Clone());
            }

            if(overlayMap.count("Data_MC_Resolution_KaMu") == 0) 
            {
                overlayMap["Data_MC_Resolution_KaMu"] = new TList();
                overlayMap["Data_MC_Resolution_KaMu"]->Add((TGraph*)res_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Resolution_KaMu"]->Add((TGraph*)res_graph->Clone());
            }
        }


    }


    // overlay the mean and resolution vs x plots
    for(auto& item: overlayMap)
    {   
        TString yname;
        float ymin = -999;
        float ymax = -999;

        if(item.first.Contains("Mean"))
        {
            yname = "Fit Mean (GeV)";
            ymin = settings.ymean_min;
            ymax = settings.ymean_max;
        }
        else
        {
            yname = "Fit Resolution (GeV)";
            ymin = settings.yres_min;
            ymax = settings.yres_max;
        }

        TCanvas* c = DiMuPlottingSystem::overlay(item.second, ymin, ymax, item.first+"_"+settings.muID, item.first+"_"+settings.muID, settings.xname, yname, false);
        overlays_list->Add(c);
        c->SaveAs("rootplots/dimu_pt/"+settings.peak+"_"+item.first+"_"+settings.muID+"_"+settings.xname+"_"+settings.whichDY+".png");
    }   

    TString savename = Form("rootplots/dimu_pt/masscalibration_%s_%s_%s_data_MC_%s_leadPt%d.root", settings.peak.Data(), settings.xname.Data(), 
                                                                                           settings.muID.Data(), settings.whichDY.Data(), (int) settings.leadPtMin);
    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile(savename, "RECREATE");
    savefile->cd();
    TDirectory* overlaydir = savefile->mkdir("overlays");
    TDirectory* plotdir = savefile->mkdir("plots");
    TDirectory* histodir = savefile->mkdir("histos");

    // Create the stack and ratio plot    
    overlaydir->cd();
    overlays_list->Write();

    plotdir->cd();
    plots_list->Write();
 
    histodir->cd();
    histos_list->Write();

    savefile->Write();
    savefile->Close();

    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
