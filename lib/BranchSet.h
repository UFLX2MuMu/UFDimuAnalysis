///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// BranchSet.h                                                           //
// ======================================================================//
// Link TTree info to objects in VarSet.h                                //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_BRANCHSET
#define ADD_BRANCHSET

#include "TBranch.h"

class BranchSet
{
    public:
        BranchSet(){};
        ~BranchSet(){};

//        TBranch* nVertices = 0;
//        TBranch* nJets     = 0;
//        TBranch* nJetsCent = 0;
//        TBranch* nJetsFwd  = 0;
//	TBranch* nBLoose   = 0;
//	TBranch* nBMed     = 0;
//	TBranch* nBTight   = 0;
        TBranch* eventInfo = 0;
        TBranch* mht = 0;
        TBranch* met = 0;

        TBranch* muons     = 0;
        TBranch* muPairs   = 0;
        TBranch* electrons = 0;
        TBranch* jets      = 0;
        TBranch* jetPairs  = 0;

//        TBranch* eff_wgt = 0;
        TBranch* pu_wgt  = 0;
        TBranch* nPU     = 0;
        TBranch* gen_wgt = 0;
//        TBranch* lhe_ht  = 0;

        TBranch* isoMu_SF_3 = 0;
        TBranch* muID_SF_3  = 0; 
        TBranch* muIso_SF_3 = 0; 

        TBranch* genParents = 0;
        TBranch* genMuons   = 0;
        TBranch* genDimuons = 0;

        void getEntry(int i)
        {
//            if(nVertices != 0) nVertices->GetEntry(i);
//            if(nJets     != 0) nJets->GetEntry(i);
//            if(nJetsCent != 0) nJetsCent->GetEntry(i);
//            if(nJetsFwd  != 0) nJetsFwd->GetEntry(i);
//            if(nBLoose   != 0) nBLoose->GetEntry(i);
//            if(nBMed     != 0) nBMed->GetEntry(i);
//            if(nBTight   != 0) nBTight->GetEntry(i);
            if(eventInfo != 0) eventInfo->GetEntry(i);
            if(mht != 0) mht->GetEntry(i);
            if(met != 0) met->GetEntry(i);

            if(muons     != 0) muons->GetEntry(i);
            if(muPairs   != 0) muPairs->GetEntry(i);
            if(electrons != 0) electrons->GetEntry(i);
            if(jets      != 0) jets->GetEntry(i);
            if(jetPairs  != 0) jetPairs->GetEntry(i);

//            if(eff_wgt != 0) eff_wgt->GetEntry(i);
            if(pu_wgt  != 0) pu_wgt->GetEntry(i);
            if(nPU     != 0) nPU->GetEntry(i);
            if(gen_wgt != 0) gen_wgt->GetEntry(i);
//            if(lhe_ht  != 0) lhe_ht->GetEntry(i);

            if(isoMu_SF_3 != 0) isoMu_SF_3->GetEntry(i);
            if(muID_SF_3  != 0)  muID_SF_3->GetEntry(i);
            if(muIso_SF_3 != 0) muIso_SF_3->GetEntry(i);

            if(genParents != 0) genParents->GetEntry(i);
            if(genMuons   != 0) genMuons->GetEntry(i);
            if(genDimuons != 0) genDimuons->GetEntry(i);
        }

        void getEntryReco(int i)
        {
//            if(nVertices != 0) nVertices->GetEntry(i);
//            if(nJets     != 0) nJets->GetEntry(i);
//            if(nJetsCent != 0) nJetsCent->GetEntry(i);
//            if(nJetsFwd  != 0) nJetsFwd->GetEntry(i);
//            if(nBLoose   != 0) nBLoose->GetEntry(i);
//            if(nBMed     != 0) nBMed->GetEntry(i);
//            if(nBTight   != 0) nBTight->GetEntry(i);
            if(eventInfo != 0) eventInfo->GetEntry(i);
            if(mht != 0) mht->GetEntry(i);
            if(met != 0) met->GetEntry(i);

            if(muons     != 0) muons->GetEntry(i);
            if(muPairs   != 0) muPairs->GetEntry(i);
            if(electrons != 0) electrons->GetEntry(i);
            if(jets      != 0) jets->GetEntry(i);
            if(jetPairs  != 0) jetPairs->GetEntry(i);
        }

//        void getEntryGenHT(int i)
//        {
//            if(lhe_ht  != 0) lhe_ht->GetEntry(i);
//        }

        void getEntryGenCollections(int i)
        {
            if(genParents != 0) genParents->GetEntry(i);
            if(genMuons   != 0) genMuons->GetEntry(i);
            if(genDimuons != 0) genDimuons->GetEntry(i);
        }

        void getEntryWeightsMC(int i)
        {
//            if(eff_wgt != 0) eff_wgt->GetEntry(i);
            if(pu_wgt  != 0) pu_wgt->GetEntry(i);
            if(gen_wgt != 0) gen_wgt->GetEntry(i);

            if(isoMu_SF_3 != 0) isoMu_SF_3->GetEntry(i);
            if(muID_SF_3  != 0)  muID_SF_3->GetEntry(i);
            if(muIso_SF_3 != 0) muIso_SF_3->GetEntry(i);
        }
};

#endif
