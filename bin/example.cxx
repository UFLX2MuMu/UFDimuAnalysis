/////////////////////////////////////////////////////////////////////////////
//                           example.cxx                                   //
//=========================================================================//
//                                                                         //
// Serves as a basic example showing some functionality of the             //
// UFDimuAnalysis framework. I also use this to prototype new              //
// functionality before integrating it into the other executables.         //
//                                                                         //
// Need to change "UF" to "CERN" in getSamples(...) in the code below      // 
// if you are running this at CERN.                                        //
//                                                                         //
// Set MAIN=example in the makefile then run `make` to compile this        //
// into an executable. Then run via ./example.                             //
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

double function()
{
    return 1.0;
}

int main(int argc, char* argv[])
{

  std::cout << "Inside example.cxx" << std::endl;
  
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    std::map<TString, Sample*> samples;
    std::cout << "\nAbout to get the samples" << std::endl;
    
    // load samples into our map. "UF" if you are on the UF servers
    // or "CERN" if you at CERN. "ALL" specifies that we want to load the Data
    // and all of the MC samples. Can loop through and remove the ones you don't want 
    // to use if you desire or just grab the ones you care about from the map.
    GetSamples(samples, "CERN_2017", "ALL");

    std::cout << std::endl << "\nGot the samples" << std::endl;

    // Grab the gluon gluon fusion h to mu mu sample
    Sample* s = samples["H2Mu_gg"];

    std::cout << "Using sample " << s->name << std::endl;

    // Get branches, set addresses
    // tells the TTree that it should load the event information into s->vars
    s->setBranchAddresses("");

    std::cout << "Setting up Selection, Collection Cleaning, and Categorizer objects" << std::endl;

    // Objects to help clean valid objects from the net collections
    // in Sample->vars
    JetCollectionCleaner      jetCollectionCleaner;
    MuonCollectionCleaner     muonCollectionCleaner;
    EleCollectionCleaner      eleCollectionCleaner;

    // Objects to cut events from the analysis
    // see ../selection/MuonSelection.h or ../selection/EventSelection.h
    // choose from an available option or make your own implementing the interface
    // use selection.evaluate(s->vars) to see whether an event passes or fails the cuts
    Run2MuonSelectionCuts  run2MuonSelection;
    Run2EventSelectionCuts run2EventSelectionMC;

    // object to categorize the event appropriately. can make your own categorizer or use one of those in
    // ../selection/CategorySelection.h. use one defined there or make your own implementing the interface.

    // categorizer object has a map of the different categories in Categorizer.categoryMap<TString, Category> 
    // which maps the category name to the Category object
    // Then each category object can store histograms in category.histoMap<TString, TH1D*>
    // use categorySelection.evaluate(s->vars) to see which category the event falls into
    CategorySelectionRun1 categorySelection;

    // Set up the histograms for each cateory, so we can fill them later
    // c.first is the category name, c.second is the category object
    std::cout << "Setting up histograms for the different categories" << std::endl;
    for(auto &c : categorySelection.categoryMap)
    {
        // Each category has a map to keep track of histos you want to fill
        // let's make dimuon_mass histos for each category
        // we will fill them after applying cuts and evaluating the categorization
        c.second.histoMap["dimu_mass"] = new TH1D("GGF_H2Mu_Mass_"+c.first, "GGF_H2Mu_Mass_"+c.first, 120, 65, 185);
        c.second.histoMap["dimu_mass"]->GetXaxis()->SetTitle("dimu_mass");
    }

    std::cout << "Looping over the events in the sample" << std::endl;
    for(unsigned int i=0; i<25; i++)
    {
       // load info from the ttree into s->vars
       // s->branches.object (load info) <-> s->vars.object (access info)
       s->branches.muPairs->GetEntry(i);
       s->branches.muons->GetEntry(i);

       if(s->vars.muPairs->size() > 0)
       {
           // reset the categories so we get the correct categorization for this event
           categorySelection.reset();

           // let's only fill the histograms for events with one dimuon candidate
           if(s->vars.muPairs->size() != 1) continue;

           // Set aliases for the dimuon candidate and its muons so we don't have to type as much
           // access objects and their info in s->vars
           MuPairInfo& dimu = s->vars.muPairs->at(0); 
           s->vars.dimuCand = &dimu;
           MuonInfo& mu1 = s->vars.muons->at(dimu.iMu1);
           MuonInfo& mu2 = s->vars.muons->at(dimu.iMu2);

           // Load the rest of the information we might need
           s->branches.jets->GetEntry(i);
           s->branches.mht->GetEntry(i);
           s->branches.eventInfo->GetEntry(i);
           s->branches.nVertices->GetEntry(i);
           s->branches.electrons->GetEntry(i);

           // load more information if we have a mc sample
           if(s->sampleType != "data")
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

           // Can now categorize the event with our categorizer 
           categorySelection.evaluate(s->vars);

           // see which categories the event fell into
           std::cout << std::endl << "// categorizing event " << i << " ..." << std::endl;
           categorySelection.outputResults();

           // Look at each category, if the event belongs to that category fill the histogram for the sample x category
           for(auto &c : categorySelection.categoryMap)
           {
               // c.first is the category name, c.second is the category object    
               // only fill histograms for the categories the event fell into
               if(!c.second.inCategory) continue;
                   c.second.histoMap["dimu_mass"]->Fill(dimu.mass, s->getWeight());  // fill the histogram for this category
                                                                                // weight the event appropriately 
           }

           std::cout << std::endl << "// printing some object info for event " << i << " ..." << std::endl;
           // output information about some objects
           // see ../lib/VarSet.h to see what other objects we have access to
           // see ../lib/analyzer_objects/ to look at the information you can access
           // for the different objects like pt, eta, etc
           for(auto& dimu: (*s->vars.muPairs))
           {
               std::cout << i << " dimu: " << dimu.outputInfo() << std::endl;
               // some fields you can access: dimu.mass, dimu.pt, dimu.eta, etc
           }
           for(auto& m: (*s->vars.muons))
           {
               std::cout << i << " muon: " << m.outputInfo() << std::endl;
               // some fields you can access: mu.pt, mu.eta, etc
           }
       }
    }

    std::cout << std::endl;
    std::cout << "Save histograms and print out histo info ... " << std::endl;
    TFile* f = new TFile("example.root", "RECREATE");
    f->cd();

    // c.first is the category name, c.second is the category object
    for(auto &c : categorySelection.categoryMap)
    {
        double luminosity = 36814;
        TH1D* h = c.second.histoMap["dimu_mass"];
        h->Scale(s->getLumiScaleFactor(luminosity));  // scale histogram based upon data luminosity
                                                  // and the xsec of the process (should loop over all the events
                                                  // in the sample for this to be accurate)
        std::cout << Form("category: %s, histo: %s, integral: %f \n", c.first.Data(), "mass", h->Integral());
        h->Write();
    }
    f->Write();
    f->Close();
    std::cout << "Done." << std::endl << std::endl;
}
