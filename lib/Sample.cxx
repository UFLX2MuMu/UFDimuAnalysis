///////////////////////////////////////////////////////////////////////////
//                             Sample.cxx                                //
//=======================================================================//
//                                                                       //
// Load info from analysis ttrees and variables into this datastructure. //
// Keep track of xsec, name, and files.                                  //
// TTree vars are loaded into objects in Sample.vars.                    //
// This object also has functionality to set the appropriate             //
// weights for MC based on gen weight, efficiency, and PU.               //
// can set the normalization based upon xsec and lumi.                   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "TCanvas.h"

///////////////////////////////////////////////////////////////////////////
// _______________________Constructor/Destructor_________________________//
///////////////////////////////////////////////////////////////////////////

Sample::Sample() 
{
    lumiWeights = 0;
    xsec = -999; 
    lumi = -999;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Sample::Sample(std::vector<TString> ifilenames, TString iname, TString isampleType)
{
    filenames = ifilenames;
    name = iname;
    if (isampleType != "NONE")
      sampleType = isampleType;

    treename = TString("dimuons/tree");
    chain = new TChain(treename);
    for (int i = 0; i < filenames.size(); i++) {
      chain->Add(filenames.at(i));
      //std::cout << i+1 << ", ";
    }
    N = chain->GetEntries();

    lumiWeights = 0;
    xsec = -999; 
    lumi = -999;

    setBranchAddresses();
    calculateNoriginal();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Sample::Sample(TString iname, TString isampleType)
{
    name = iname;
    if (isampleType != "NONE")
      sampleType = isampleType;

    treename = TString("dimuons/tree");

    lumiWeights = 0;
    xsec = -999; 
    lumi = -999;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Sample::~Sample() {
  // free pointed to memory!
  if (chain != 0) {
    delete chain;
  }
  if (files.size() !=0) {
    files.clear();
  }
  if (lumiWeights !=0) {
    delete lumiWeights;
  }
}

///////////////////////////////////////////////////////////////////////////
// _______________________Other Functions________________________________//
///////////////////////////////////////////////////////////////////////////

void Sample::setBranchAddresses(TString options)
{
      // Link only to the branches we need to save a lot of time
      // run 1 category info 
      branches.muons      = chain->GetBranch("muons");
      branches.muPairs    = chain->GetBranch("muPairs");
      branches.jets       = chain->GetBranch("jets");
      branches.jetPairs   = chain->GetBranch("jetPairs");
      branches.mht        = chain->GetBranch("mht");
      branches.met        = chain->GetBranch("met");

      branches.eventInfo  = chain->GetBranch("event");
//      branches.nVertices  = chain->GetBranch("nVertices");
//      branches.nJets      = chain->GetBranch("nJets");
//      branches.nJetsCent  = chain->GetBranch("nJetsCent");
//      branches.nJetsFwd   = chain->GetBranch("nJetsFwd");
//      branches.nBLoose    = chain->GetBranch("nBLoose");
//      branches.nBMed      = chain->GetBranch("nBMed");
//      branches.nBTight    = chain->GetBranch("nBTight");

      branches.electrons = chain->GetBranch("eles");

      // branches for MC samples
      if(!sampleType.EqualTo("data"))
      {
          branches.nPU     = chain->GetBranch("nPU");
          branches.gen_wgt = chain->GetBranch("GEN_wgt");
          branches.pu_wgt  = chain->GetBranch("PU_wgt");

          if(options == "PU_up")
              branches.pu_wgt  = chain->GetBranch("PU_wgt_up");

          if(options == "PU_down")
              branches.pu_wgt  = chain->GetBranch("PU_wgt_down");

          if(options == "JES_up")
          {
              branches.jets      = chain->GetBranch("jets_JES_up");
              branches.jetPairs  = chain->GetBranch("jetPairs_JES_up");
              branches.mht       = chain->GetBranch("mht_JES_up");
              branches.met       = chain->GetBranch("met_JES_up");
//              branches.nJets     = chain->GetBranch("nJets_JES_up");
//              branches.nJetsCent = chain->GetBranch("nJetsCent_JES_up");
//              branches.nJetsFwd  = chain->GetBranch("nJetsFwd_JES_up");
//              branches.nBLoose   = chain->GetBranch("nBLoose_JES_up");
//              branches.nBMed     = chain->GetBranch("nBMed_JES_up");
//              branches.nBTight   = chain->GetBranch("nBTight_JES_up");
          }

          if(options == "JES_down")
          {
              branches.jets      = chain->GetBranch("jets_JES_down");
              branches.jetPairs  = chain->GetBranch("jetPairs_JES_down");
              branches.mht       = chain->GetBranch("mht_JES_down");
              branches.met       = chain->GetBranch("met_JES_down");
//              branches.nJets     = chain->GetBranch("nJets_JES_down");
//              branches.nJetsCent = chain->GetBranch("nJetsCent_JES_down");
//              branches.nJetsFwd  = chain->GetBranch("nJetsFws_JES_down");   // typo in analyzer collection
//              branches.nBLoose   = chain->GetBranch("nBLoose_JES_down");
//              branches.nBMed     = chain->GetBranch("nBMed_JES_down");
//              branches.nBTight   = chain->GetBranch("nBTight_JES_down");
          }

//          branches.eff_wgt = chain->GetBranch("IsoMu_eff_3");
//          branches.lhe_ht  = chain->GetBranch("LHE_HT");

          branches.isoMu_SF_3 = chain->GetBranch("IsoMu_eff_3");   //these three should SF_3, but for now (2017 ntuple v3), they are eff_3
          branches.muID_SF_3  = chain->GetBranch("MuID_eff_3");
          branches.muIso_SF_3 = chain->GetBranch("MuIso_eff_3");

          branches.gen_wgt->SetAddress(&vars.gen_wgt);
          branches.nPU->SetAddress(&vars.nPU);
          branches.pu_wgt->SetAddress(&vars.pu_wgt);
//          branches.eff_wgt->SetAddress(&vars.eff_wgt);
//          branches.lhe_ht->SetAddress(&vars.lhe_ht);

          branches.isoMu_SF_3->SetAddress(&vars.isoMu_SF_3); 
          branches.muID_SF_3->SetAddress(&vars.muID_SF_3); 
          branches.muIso_SF_3->SetAddress(&vars.muIso_SF_3); 

          branches.genParents = chain->GetBranch("genParents");
          branches.genMuons   = chain->GetBranch("genMuons");
          branches.genDimuons = chain->GetBranch("genMuPairs");

          branches.genParents->SetAddress(&vars.genParents);
          branches.genMuons->SetAddress(&vars.genMuons);
          branches.genDimuons->SetAddress(&vars.genDimuons);
      }

      branches.muons->SetAddress(&vars.muons);
      branches.muPairs->SetAddress(&vars.muPairs);
      branches.jets->SetAddress(&vars.jets);
      branches.jetPairs->SetAddress(&vars.jetPairs);
      branches.mht->SetAddress(&vars.mht);
      branches.met->SetAddress(&vars.met);

      branches.eventInfo->SetAddress(&vars.eventInfo);
//      branches.nVertices->SetAddress(&vars.nVertices);
//      branches.nJets->SetAddress(&vars.nJets);
//      branches.nJetsCent->SetAddress(&vars.nJetsCent);
//      branches.nJetsFwd->SetAddress(&vars.nJetsFwd);
//      branches.nBLoose->SetAddress(&vars.nBLoose);
//      branches.nBMed->SetAddress(&vars.nBMed);
//      branches.nBTight->SetAddress(&vars.nBTight);

      branches.electrons->SetAddress(&vars.electrons);
}
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Sample::calculateNoriginal()
{
  gROOT->SetBatch(1); // Don't draw TCanvas
  TCanvas c;
  // Calculate the number of original events using the meta data tree
  // Calculate the weighted number of original events as well
  TChain* metadata = new TChain("dimuons/metadata");
  for (auto f_name : filenames) {
    metadata->Add(f_name);
  }
  metadata->Draw("sumEventWeights>>eweights_"+name);
  TH1F* sumEventWeightsHist = (TH1F*) gDirectory->Get("eweights_"+name); 
  
  // There are many ttrees combined so the histogram has a numEvents entry for each
  // file. The total number is equal to the average number times the total number of entries.
  nOriginalWeighted = sumEventWeightsHist->GetEntries()*sumEventWeightsHist->GetMean();
  
  metadata->Draw("originalNumEvents>>nevents_"+name);
  TH1F* nEventsHist = (TH1F*) gDirectory->Get("nevents_"+name); 
  nOriginal = nEventsHist->GetEntries()*nEventsHist->GetMean();
  if(metadata !=0) delete metadata;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int Sample::getEntry(int i)
{
    branches.getEntry(i);
    return i;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int Sample::getEntry(int i, TEntryList* list)
{
    int chainnum = list->GetEntry(i);
    branches.getEntry(chainnum);
    return chainnum;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

double Sample::getWeight()
{
// Assumes getEntry has already been called to load the appropriate values.
// Gets the weight for the histogram depending on the sample type 
    if(sampleType == "data") return 1.0;
    else if(lumiWeights == 0) return 1.0*vars.gen_wgt*vars.pu_wgt*vars.sf();
//    else return 1.0*vars.gen_wgt*lumiWeights->weight(vars.nPU)*vars.sf();
    //else return 1.0*vars.gen_wgt;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

float Sample::getLumiScaleFactor(float luminosity)
{
// Scale the MC histograms based upon the data luminosity, the number of events
// that the CMSSW analyzer looked at, and the xsec for the process
    if(sampleType.EqualTo("data")) return 1.0;
    else return luminosity*xsec/nOriginalWeighted;
}
