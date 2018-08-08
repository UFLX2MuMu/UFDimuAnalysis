///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// VarSet.h                                                              //
// ======================================================================//
// TTree variables are loaded into these objects. Each sample has        //
// its own VarSet data structure.                                        //
// We also set up a map<TString,FUNCTION> to get feature values based    //
// upon the name of the feature.                                         //
// Has muons, dimuons, jets, electrons, gen muons, etc...                //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_VARSET
#define ADD_VARSET

#include "DataFormats.h"
#include "EventInfo.h"
#include "MuonInfo.h"
#include "MuPairInfo.h"
#include "JetPairInfo.h"
#include "EleInfo.h"
#include "MhtInfo.h"
#include "MetInfo.h"
#include "SlimJetInfo.h"
#include "GenParentInfo.h"
#include "GenMuonInfo.h"
#include "GenMuPairInfo.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <map>

const int N_JETS      = 4;
const int N_JET_PAIRS = 4;

class VarSet
{
    public:
        VarSet();
        ~VarSet(){};

        //////////////////////////////////////////////////////////////////////////////
        // TMVA Binary/Multi classification vars -------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        double bdt_out = -999;
        double bdt_ggh_out = -999;
        double bdt_vbf_out = -999;
        double bdt_vh_out = -999;
        double bdt_ewk_out = -999;
        double bdt_top_out = -999;

        //////////////////////////////////////////////////////////////////////////////
        // Vars loaded from Analyzer TTrees ------------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        // reco weights for mc
        float eff_wgt;            // use this if you don't match MC to trigger
                                  // scales mc to trigger passed data based on pt,eta

        float isoMu_SF_3;  // match MC to trigger, account for discrepancy from there
        float isoMu_SF_4;  // scale factors change for different eras hence the 3 and 4
        float muID_SF_3;   // match MC to muID, adjust to data from there
        float muID_SF_4;   // ...
        float muIso_SF_3;  // cut mc based on iso, adjust to data from there
        float muIso_SF_4;  // ...

        float sf()
        {
            // average the scale factors for the different eras and multiply all of them
            // to get the net scale factor
            return 0.5*(isoMu_SF_3 + isoMu_SF_4)*0.5*(muID_SF_3 + muID_SF_4)*0.5*(muIso_SF_3 + muIso_SF_4);
        }

        float pu_wgt;      // weight mc based upon PU to match data PU distribution
        

        // reco info
        Int_t nVertices;
        Int_t nJets;
        Int_t nJetsCent;
        Int_t nJetsFwd;
        Int_t nBLoose;
        Int_t nBMed;
        Int_t nBTight;
        EventInfo* eventInfo = 0;
        MuPairInfo* dimuCand = 0; // this is a pointer to one of the dimu candidates in the vector
                                  // we don't want to copy the object for ~40 million events
                                  // the other objects/primitives are loaded via TBranch->Get();
        MhtInfo* mht = 0;
        MetInfo* met = 0;

        std::vector<MuonInfo>* muons = 0;
        std::vector<MuPairInfo>* muPairs = 0;
        std::vector<EleInfo>* electrons = 0;
        std::vector<SlimJetInfo>* jets = 0;
        std::vector<JetPairInfo>* jetPairs = 0;

        // gen info
        std::vector<GenParentInfo>* genParents = 0;
        std::vector<GenMuonInfo>* genMuons = 0;
        std::vector<GenMuPairInfo>* genDimuons = 0;

        int nPU;
        float lhe_ht;

        // gen weights
        int gen_wgt;

        //////////////////////////////////////////////////////////////////////////////
        // Cleaned Collections ------------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        std::vector<TLorentzVector> validMuons;
        std::vector<TLorentzVector> validExtraMuons;
        std::vector<TLorentzVector> validElectrons;
        std::vector<TLorentzVector> validJets;
        std::vector<TLorentzVector> validBJets;

        //////////////////////////////////////////////////////////////////////////////
        // Map String to Feature value (double) -------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        // map variable string to variable value via function pointers
        // uses functions below
        std::unordered_map<std::string, double(VarSet::*)()> varMap;
        std::unordered_map<std::string, double(VarSet::*)(int)> varMapI;

        // get the the value for some variable in one of the structs above
        // by name (string). Use the varMap to get the appropriate function
        // of those below.
        // Now we can easily output the values for training TString->Value
        // And evaluate the XML categories via getValue(varTString) > cut
        double getValue(const std::string& name) 
        {
	  // std::cout << "\n  * Getting value for variable with name " << name << std::endl;
          // Must call the function for this particular instance
          if (varMap[name])
            return (this->*varMap[name])();
	  else if (varMapI[name]) {
	    std::string iStr = &(name.substr(0, name.find("_")).back());
	    int iObj = -99;
	    std::stringstream convert(iStr);
	    convert >> iObj;
            return (this->*varMapI[name])(iObj - 1); // Vector indices are 1 lower
	  }
	  else return -999;
        }
	
        // see if the var string is in the map
        bool checkForVar(const std::string& name){return (varMap[name] || varMapI[name]); }

        //////////////////////////////////////////////////////////////////////////////
        // Set systematic to vary ---------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        //void setSystematicVariation(TString which = "")
        //{
        //    if(which == "jes_up") 
        //    {
        //        jets = jets_up;
        //    }
        //    else if(which == "jes_down")
        //    {
        //        jets = jets_down;
        //    }
        //    else if(which == "pu_up")
        //    {
        //        pu_wgt = pu_wgt_up;
        //        // jets = jets;
        //    }
        //    else if(which == "pu_down")
        //    {
        //        pu_wgt = pu_wgt_down;
        //        // jets = jets;
        //    }
        //    else
        //    {
        //        // pu_wgt = pu_wgt;
        //        // jets = jets;
        //    }
        //}


        //////////////////////////////////////////////////////////////////////////////
        // Set dimu mass, mu_pt to PF, Roch, or KaMu ---------------------------------
        //////////////////////////////////////////////////////////////////////////////

        void setCalibrationType(TString ctype)
        {
          if(ctype == "PF") 
          {
              dimuCand->mass = dimuCand->mass_PF;
              muons->at(dimuCand->iMu1).pt = muons->at(dimuCand->iMu1).pt_PF;
              muons->at(dimuCand->iMu2).pt = muons->at(dimuCand->iMu2).pt_PF;
          }
          else if(ctype == "Roch") 
          {
              dimuCand->mass = dimuCand->mass_Roch;
              muons->at(dimuCand->iMu1).pt = muons->at(dimuCand->iMu1).pt_Roch;
              muons->at(dimuCand->iMu2).pt = muons->at(dimuCand->iMu2).pt_Roch;
          }
          else if(ctype == "KaMu") 
          {
              dimuCand->mass = dimuCand->mass_KaMu;
              muons->at(dimuCand->iMu1).pt = muons->at(dimuCand->iMu1).pt_KaMu;
              muons->at(dimuCand->iMu2).pt = muons->at(dimuCand->iMu2).pt_KaMu;
          }
        }

        //////////////////////////////////////////////////////////////////////////////
        // Initialize VBF Jets -------------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        // cuts for vbf jets
        double cLeadPtMin = 40;
        double cDijetMassMinVBFT = 650;
        double cDijetDeltaEtaMinVBFT = 3.5;

        // set the above cuts if they change
        void setVBFcuts(double leadPtMin, double dijetMassMinVBFT, double dijetDeltaEtaMinVBFT)
        {
            cLeadPtMin = leadPtMin;
            cDijetMassMinVBFT = dijetMassMinVBFT;
            cDijetDeltaEtaMinVBFT = dijetDeltaEtaMinVBFT;
        }

        // index for vbf jets
        int vbf_j0 = -999;
        int vbf_j1 = -999;
 
        // get jets that represent vbf jets
        void setVBFjets()
        {
            vbf_j0 = 0;
            vbf_j1 = 1;
            double mjj_max = -999;

            for(unsigned int i=0; i<validJets.size(); i++) 
            {
                if(!(validJets[i].Pt() > cLeadPtMin)) break;
                for(unsigned int j=i+1; j<validJets.size(); j++) 
                {
                    double dEtajj = TMath::Abs(validJets[i].Eta() - validJets[j].Eta());
                    double mjj = (validJets[i] + validJets[j]).M();
                    if(mjj > mjj_max)
                    {
                        vbf_j0 = i;
                        vbf_j1 = j;
                    }
                    if(mjj > cDijetMassMinVBFT && dEtajj > cDijetDeltaEtaMinVBFT)
                    {
                        vbf_j0 = i;
                        vbf_j1 = j;
                        return; 
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////
        // Initialize Standard Jets --------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////

        // index for our "standard" jets
        int j0 = 0;
        int j1 = 1;

        void setJets()
        {
        // standard jets will be the two corresponding to the max mjj value 
            double mjj_max = -999;
            j0 = 0;
            j1 = 1;

            for(unsigned int i=0; i<validJets.size(); i++) 
            {
                for(unsigned int j=i+1; j<validJets.size(); j++) 
                {
                    double dEtajj = TMath::Abs(validJets[i].Eta() - validJets[j].Eta());
                    double mjj = (validJets[i] + validJets[j]).M();
                    if(mjj > mjj_max)
                    {
                        j0 = i;
                        j1 = j;
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////
        // Functions that return feature values via getValue(TString) ----------------
        //////////////////////////////////////////////////////////////////////////////

        // Map initialized in cxx file -----------------------------------------------
        double _nVertices(){ return nVertices;            };

        // load the bdt score manually beforehand
        double bdt_score()    { return bdt_out;                };
        double bdt_ggh_score(){ return bdt_ggh_out;            };
        double bdt_vbf_score(){ return bdt_vbf_out;            };
        double bdt_vh_score() { return bdt_vh_out;             };
        double bdt_ewk_score(){ return bdt_ewk_out;            };
        double bdt_top_score(){ return bdt_top_out;            };

        // Object counting
	double _nJets()       { return nJets;                  };
	double _nJetsCent()   { return nJetsCent;              };
	double _nJetsFwd()    { return nJetsFwd;               };
	double _nBLoose()     { return nBLoose;                };
	double _nBMed()       { return nBMed;                  };
	double _nBTight()     { return nBTight;                };
        double nValJets()     { return validJets.size();       };
        // double nValJetsCent() { return nJetsCent;              };
        // double nValJetsFwd()  { return nJetsFwd;               };
        double nValBTags()    { return validBJets.size();      };
        double nExtraLep()    { return validExtraMuons.size() + validElectrons.size(); };
        double nExtraMu()     { return validExtraMuons.size(); };
        double nEle()         { return validElectrons.size();  };

	// Dimuon variables
        double dimu_mass()        { return dimuCand->mass;             };
        double dimu_mass_Roch()   { return dimuCand->mass_Roch;        };
        double dimu_mass_KaMu()   { return dimuCand->mass_KaMu;        };

        double massErr_PF()  { return   dimuCand->massErr_PF; }
        double massErr_Roch(){ return dimuCand->massErr_Roch; }
        double massErr_KaMu(){ return dimuCand->massErr_KaMu; }

        double dimu_pt()          { return dimuCand->pt;               };
        double dimu_eta()         { return dimuCand->eta;              };
        double dimu_abs_eta()     { return TMath::Abs(dimuCand->eta);  };
        double dimu_rapid()       { return dimuCand->rapid;            };
        double dimu_dR()          { return dimuCand->dR;               };
        double dimu_dEta()        { return dimuCand->dEta;             };
        double dimu_abs_dEta()    { return TMath::Abs(dimuCand->dEta); };
        double dimu_dPhi()        { return dimuCand->dPhi;             };
        double dimu_abs_dPhi()    { return TMath::Abs(dimuCand->dPhi); };
        double dimu_dPhiStar() {  // Phi separation in the parent's rest frame
	  double phi_star = 0;
	  double mu_dPhi = TMath::Abs(muons->at(dimuCand->iMu1).phi - muons->at(dimuCand->iMu2).phi);
	  if(mu_dPhi > TMath::Pi()) mu_dPhi = 2*TMath::Pi() - mu_dPhi;
	  double phiACOP = TMath::Pi() - mu_dPhi;
	  double thetaStarEta = TMath::ACos(TMath::TanH((muons->at(dimuCand->iMu1).eta - muons->at(dimuCand->iMu2).eta)/2));
	  phi_star = TMath::Tan(phiACOP/2)*TMath::Sin(thetaStarEta);
	  return phi_star;
	}
        double dimu_avg_abs_eta() { return ( TMath::Abs(muons->at(dimuCand->iMu1).eta) +
					     TMath::Abs(muons->at(dimuCand->iMu2).eta) ) / 2.; };
        double dimu_min_abs_eta() { return std::min( TMath::Abs(muons->at(dimuCand->iMu1).eta), 
						     TMath::Abs(muons->at(dimuCand->iMu2).eta) );   };
        double dimu_max_abs_eta() { return std::max( TMath::Abs(muons->at(dimuCand->iMu1).eta), 
						     TMath::Abs(muons->at(dimuCand->iMu2).eta) );   };
	
        // Muon variables from dimuon pair
        double mu1_pt()      { return muons->at(dimuCand->iMu1).pt;              };
        double mu2_pt()      { return muons->at(dimuCand->iMu2).pt;              };
        double mu1_eta()     { return muons->at(dimuCand->iMu1).eta;             };
	double mu2_eta()     { return muons->at(dimuCand->iMu2).eta;             };
	double mu1_abs_eta() { return TMath::Abs(muons->at(dimuCand->iMu1).eta); };
	double mu2_abs_eta() { return TMath::Abs(muons->at(dimuCand->iMu2).eta); };

	// Dijet variables
        double dijet_mass(int i)        { return (jetPairs->size() > i) ? jetPairs->at(i).mass : 0; };
        double dijet_pt(int i)          { return (jetPairs->size() > i) ? jetPairs->at(i).pt   : 0; };
        double dijet_dEta(int i)        { return (jetPairs->size() > i) ? jetPairs->at(i).dEta : -10; };
        double dijet_abs_dEta(int i)    { return (jetPairs->size() > i) ? TMath::Abs(jetPairs->at(i).dEta) : -1; };
        double dijet_min_abs_eta(int i) { 
	  if (jetPairs->size() > i)
	    return std::min( TMath::Abs(jets->at(jetPairs->at(i).iJet1).eta), TMath::Abs(jets->at(jetPairs->at(i).iJet2).eta) );
	  else return 0;
	}
        double dijet_max_abs_eta(int i) { 
	  if (jetPairs->size() > i)
	    return std::max( TMath::Abs(jets->at(jetPairs->at(i).iJet1).eta), TMath::Abs(jets->at(jetPairs->at(i).iJet2).eta) );
	  else return 0;
	}

        // Jet variables
	double jet_pt(int i)      { return (validJets.size() > i) ? validJets[i].Pt() : 0; };
	double jet_eta(int i)     { return (validJets.size() > i) ? validJets[i].Eta() : -5; };
	double jet_abs_eta(int i) { return (validJets.size() > i) ? TMath::Abs(validJets[i].Eta()) : -1; };

        // Global event variables
        double MET()      { return met->pt;       };
        double MHT()      { return mht->pt;       };
        double MT_had()   { return mht->MT_had;   };
        double mass_had() { return mht->mass_had; };


	// AMC variables
	double jj_jet0_pt()  { return (validJets.size()>=1)?validJets[j0].Pt():-999;  };
        double jj_jet1_pt()  { return (validJets.size()>=2)?validJets[j1].Pt():-999;  };
	double jj_jet0_eta() { return (validJets.size()>=1)?validJets[j0].Eta():-999; };
        double jj_jet1_eta() { return (validJets.size()>=2)?validJets[j1].Eta():-999; };

	double m_jj()         { return (validJets.size()>=2)?(validJets[j0]+validJets[j1]).M():-999; };
        double dEta_jj()      { return (validJets.size()>=2)?TMath::Abs(validJets[j0].Eta()-validJets[j1].Eta()):-999; };
        double dEta_jj_mumu() { return (validJets.size()>=2)?TMath::Abs((validJets[j0]+validJets[j1]).Eta()-dimuCand->eta):-999; };
        double dPhi_jj_mumu() { 
	  if (validJets.size() < 2) return -999;
	  double dphi = TMath::Abs((validJets[j0] + validJets[j1]).Phi() - dimuCand->phi);
	  if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        };
        double zep() { 
	  if (validJets.size() < 2) return -999; 
	  double meanEta = (validJets[j0].Eta() +validJets[j1].Eta())/2;
	  return validJets[j0].Eta() - meanEta;
        };

        // vbf jet variables
        double vbf_jet0_pt()  { return (validJets.size()>=1)?validJets[vbf_j0].Pt():-999;  };
        double vbf_jet1_pt()  { return (validJets.size()>=2)?validJets[vbf_j1].Pt():-999;  };
        double vbf_jet0_eta() { return (validJets.size()>=1)?validJets[vbf_j0].Eta():-999; };
        double vbf_jet1_eta() { return (validJets.size()>=2)?validJets[vbf_j1].Eta():-999; };

        double vbf_m_jj()         { return (validJets.size()>=2)?(validJets[vbf_j0]+validJets[vbf_j1]).M():-999; };
        double vbf_dEta_jj()      { return (validJets.size()>=2)?TMath::Abs(validJets[vbf_j0].Eta()-validJets[vbf_j1].Eta()):-999; };
        double vbf_dEta_jj_mumu() { return (validJets.size()>=2)?TMath::Abs((validJets[vbf_j0]+validJets[vbf_j1]).Eta()-dimuCand->eta):-999; };

        // bjet variables
        double bjet0_pt() { return (validBJets.size()>=1)?validBJets[0].Pt():-999;  };
        double bjet1_pt() { return (validBJets.size()>=2)?validBJets[1].Pt():-999;  };
        double bjet0_eta(){ return (validBJets.size()>=1)?validBJets[0].Eta():-999; };
        double bjet1_eta(){ return (validBJets.size()>=2)?validBJets[1].Eta():-999; };

        double m_bb()   { return (validBJets.size()>=2)?(validBJets[0]+validBJets[1]).M():-999; };
        double dEta_bb(){ return (validBJets.size()>=2)?TMath::Abs(validBJets[0].Eta()-validBJets[1].Eta()):-999; };

        double vbf_dPhi_jj_mumu()    
        { 
            if(validJets.size() < 2) return -999;
            double dphi = TMath::Abs((validJets[vbf_j0] + validJets[vbf_j1]).Phi() - dimuCand->phi);
            if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        };

        double vbf_zep() { 
	  if (validJets.size() < 2) return -999; 
	  double meanEta = (validJets[vbf_j0].Eta() +validJets[vbf_j1].Eta())/2;
	  return validJets[vbf_j0].Eta() - meanEta;
        };

        double extra_muon0_pt() { return (validExtraMuons.size()>=1)?validExtraMuons[0].Pt():-999; };
        double extra_muon1_pt() { return (validExtraMuons.size()>=2)?validExtraMuons[1].Pt():-999; };
        double extra_muon0_eta(){ return (validExtraMuons.size()>=1)?validExtraMuons[0].Eta():-999; };
        double extra_muon1_eta(){ return (validExtraMuons.size()>=2)?validExtraMuons[1].Eta():-999; };

        double electron0_pt() { return (validElectrons.size()>=1)?validElectrons[0].Pt():-999; };
        double electron1_pt() { return (validElectrons.size()>=2)?validElectrons[1].Pt():-999; };
        double electron0_eta(){ return (validElectrons.size()>=1)?validElectrons[0].Eta():-999; };
        double electron1_eta(){ return (validElectrons.size()>=2)?validElectrons[1].Eta():-999; };

        double mT_b_MET()
        { 
            if(validBJets.size() < 1) return -999;
            TLorentzVector metv(met->pt*TMath::Cos(met->phi), met->pt*TMath::Cos(met->phi), 0, met->pt);                  
            TLorentzVector bjet = validBJets[0];
            TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
            TLorentzVector bmet_t = metv + bjet_t;
            return bmet_t.M();
        };

};

#endif
