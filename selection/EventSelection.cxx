///////////////////////////////////////////////////////////////////////////
//                             EventSelection.cxx                        //
//=======================================================================//
//                                                                       //
//        Some selections to cut events. Implements Cut.h.               //
//        Cut based on trigger matching, or other event info.            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "EventSelection.h"
#include "ParticleTools.h"
#include "TMath.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run2EventSelectionCuts______________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Event selection based on the run1 h2mu analysis, still using these
// for the run2 analysis modified for our hlt trigger matching
// criteria

// The 80X samples didn't have trigger matching, so the data requires
// hlt trigger matching but the mc doesn't, then the mc is scaled

Run2EventSelectionCuts::Run2EventSelectionCuts()
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = 26; 
    // sublead pt condition is met in muon selection
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(3, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run2EventSelectionCuts::Run2EventSelectionCuts(float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = trigMuPtMin;          // >
    // sublead pt condition is met in muon selection

    cDimuMassMin = dimuMassMin;          // >
    cutset.cuts = std::vector<CutInfo>(3, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run2EventSelectionCuts::evaluate(VarSet& vars)
{
    // check the mass cut first since it's the most likely to cut the event
    // saves computation time so that we don't have to check the others
    if(!(vars.dimuCand->mass > cDimuMassMin) && cutset.cuts[2].on) return false;

    // Pt cuts on leading muon
    MuonInfo& mu1 = vars.muons->at(vars.dimuCand->iMu1);
    MuonInfo& mu2 = vars.muons->at(vars.dimuCand->iMu2);

    if(!cutset.cuts[1].on) ;
    else if((mu1.isHltMatched[2] || mu1.isHltMatched[3]) && mu1.pt > cTrigMuPtMin) ;
    else if((mu2.isHltMatched[2] || mu2.isHltMatched[3]) && mu2.pt > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    // Require oppositely charged muons in the event
    if(!(mu1.charge != mu2.charge) && cutset.cuts[0].on) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run2EventSelectionCuts::string()
{
    return TString("Run2_Event_Selection");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run2EventSelectionCuts::makeCutSet()
{
    cutset.cuts[0].name = "recoMu1.charge != recoMu2.charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("hltMatchedMu.pt > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("mass > %6.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.concatCuts();
}













///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________EventSelectionCuts_2017______________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Event selection based on the run1 h2mu analysis, still using these
// for the run2 analysis modified for our hlt trigger matching
// criteria

// The 80X samples didn't have trigger matching, so the data requires
// hlt trigger matching but the mc doesn't, then the mc is scaled

EventSelectionCuts_2017::EventSelectionCuts_2017()
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = 30; 
    // sublead pt condition is met in muon selection
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(3, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

EventSelectionCuts_2017::EventSelectionCuts_2017(float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = trigMuPtMin;          // >
    // sublead pt condition is met in muon selection

    cDimuMassMin = dimuMassMin;          // >
    cutset.cuts = std::vector<CutInfo>(3, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool EventSelectionCuts_2017::evaluate(VarSet& vars)
{
    // check the mass cut first since it's the most likely to cut the event
    // saves computation time so that we don't have to check the others
    if(!(vars.dimuCand->mass > cDimuMassMin) && cutset.cuts[2].on) return false;

    // Pt cuts on leading muon
    MuonInfo& mu1 = vars.muons->at(vars.dimuCand->iMu1);
    MuonInfo& mu2 = vars.muons->at(vars.dimuCand->iMu2);

    if(!cutset.cuts[1].on) ;
    else if((mu1.isHltMatched[2] || mu1.isHltMatched[3]) && mu1.pt > cTrigMuPtMin) ;
    else if((mu2.isHltMatched[2] || mu2.isHltMatched[3]) && mu2.pt > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    // Require oppositely charged muons in the event
    if(!(mu1.charge != mu2.charge) && cutset.cuts[0].on) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString EventSelectionCuts_2017::string()
{
    return TString("Event_Selection_2017");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void EventSelectionCuts_2017::makeCutSet()
{
    cutset.cuts[0].name = "recoMu1.charge != recoMu2.charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("hltMatchedMu.pt > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("mass > %6.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.concatCuts();
}










///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________SynchEventSelectionCuts______________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Event selection based on the run1 h2mu analysis, still using these
// for the run2 analysis modified for our hlt trigger matching
// criteria

// The 80X samples didn't have trigger matching, so the data requires
// hlt trigger matching but the mc doesn't, then the mc is scaled

SynchEventSelectionCuts::SynchEventSelectionCuts()
{
    cTrigMuPtMin = 26; 
    cDimuMassMin = 0;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

SynchEventSelectionCuts::SynchEventSelectionCuts(float trigMuPtMin, float dimuMassMin)
{
    cTrigMuPtMin = trigMuPtMin;          // >
    cDimuMassMin = dimuMassMin;          // >
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchEventSelectionCuts::evaluate(VarSet& vars)
{
    // check the mass cut first since it's the most likely to cut the event
    // saves computation time so that we don't have to check the others
    if(!(vars.dimuCand->mass > cDimuMassMin)) return false;

    // Pt cuts on leading muon
    MuonInfo& mu1 = vars.muons->at(vars.dimuCand->iMu1);
    MuonInfo& mu2 = vars.muons->at(vars.dimuCand->iMu2);

    // 2 and 3 should be isoMu24 and isoTkMu24
    if(!( (mu1.isHltMatched[2] && mu1.pt > cTrigMuPtMin) || (mu1.isHltMatched[3] && mu1.pt > cTrigMuPtMin) 
          || (mu2.isHltMatched[2] && mu2.pt > cTrigMuPtMin) || (mu2.isHltMatched[3] && mu2.pt > cTrigMuPtMin) ))
        return false;

    // Require oppositely charged muons and zero extra muons in the event
    if(!( (mu1.charge != mu2.charge) && vars.validExtraMuons.size() == 0 )) return false;

    // require that there are zero electrons
    if(!(vars.validElectrons.size() == 0)) return false;
    if(!(vars.validExtraMuons.size() == 0)) return false;
    if(!(vars.validBJets.size() == 0)) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString SynchEventSelectionCuts::string()
{
    return TString("Synch_Event_Selection");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void SynchEventSelectionCuts::makeCutSet()
{
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________FEWZCompareCuts________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Yet another set of event selection criteria. These are used
// in FEWZ (matrix element based "generator" used to model Drell Yan) 
// comparison studies

FEWZCompareCuts::FEWZCompareCuts()
{
// Default values for the fewz comparison selection cuts

    useReco = false;
    cLeadPtMin = 20;   
    cSubleadPtMin = 10;   
    cMaxEta = 2.4;     
    cDimuMassMin = 110; 
    cDimuMassMax = 310; 
    cMaxRelIso = 0.12;
    cutset.cuts = std::vector<CutInfo>(9, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

FEWZCompareCuts::FEWZCompareCuts(bool useReco)
{
// Default values for the fewz comparison selection cuts
    this->useReco = useReco;
    cLeadPtMin = 20;   
    cSubleadPtMin = 10;   
    cMaxEta = 2.4;     
    cDimuMassMin = 110; 
    cDimuMassMax = 310; 
    cMaxRelIso = 0.12;      
    cutset.cuts = std::vector<CutInfo>(9, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

FEWZCompareCuts::FEWZCompareCuts(bool useReco, float leadPtMin, float subleadPtMin, float maxEta, float dimuMassMin, float dimuMassMax, float maxRelIso)
{
// Custom selection cuts for the fewz comparison
    this->useReco = useReco;
    cLeadPtMin = leadPtMin;
    cSubleadPtMin = subleadPtMin;
    cMaxEta = maxEta;
    cDimuMassMin = dimuMassMin;
    cDimuMassMax = dimuMassMax;
    cMaxRelIso = maxRelIso; 
    cutset.cuts = std::vector<CutInfo>(9, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool FEWZCompareCuts::evaluate(VarSet& vars)
{
    float leadPt = -999;
    float subleadPt = -999; 
    float eta0 = -999;
    float eta1 = -999;
    float dimu_mass = -999;
    int charge0 = -999;
    int charge1 = -999;

    // the indices for the muons of the dimu pair
    int mu1 = vars.dimuCand->iMu1;
    int mu2 = vars.dimuCand->iMu2;

    leadPt = TMath::Max(vars.muons->at(mu1).pt, vars.muons->at(mu2).pt);
    subleadPt = TMath::Min(vars.muons->at(mu1).pt, vars.muons->at(mu2).pt);
    eta0 = vars.muons->at(mu1).eta;
    eta1 = vars.muons->at(mu2).eta;
    dimu_mass = vars.dimuCand->mass;
    charge0 = vars.muons->at(mu1).charge;
    charge1 = vars.muons->at(mu2).charge;

    // if a muon fails any of the criterea that are turned on then reutrn false
    if(cutset.cuts[0].on)
    {
        if(leadPt <= cLeadPtMin) return false;
    }
    if(cutset.cuts[1].on)
    {
        if(subleadPt <= cSubleadPtMin) return false;
    }
    if(cutset.cuts[2].on)
    {
        if(TMath::Abs(eta0) >= cMaxEta) return false;
    }
    if(cutset.cuts[3].on)
    {
        if(TMath::Abs(eta1) >= cMaxEta) return false;
    }
    if(cutset.cuts[4].on)
    {
        if(dimu_mass <= cDimuMassMin) return false;
    }
    if(cutset.cuts[5].on)
    {
        if(dimu_mass >= cDimuMassMax) return false;
    }
    if(cutset.cuts[6].on)
    {
        if(charge0 == charge1) return false;
    }
    if(cutset.cuts[7].on)
    {
        if(!(vars.muons->at(mu1).iso() <= cMaxRelIso) ) return false;
    }
    if(cutset.cuts[8].on)
    {
        if(!(vars.muons->at(mu2).iso() <= cMaxRelIso) ) return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString FEWZCompareCuts::string()
{
    return TString("FEWZ_Comparison_Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void FEWZCompareCuts::makeCutSet()
{
    //cLeadPtMin 
    //cSubleadPtMin
    //cMaxEta1
    //cMaxEta2
    //cDimuMassMin
    //cDimuMassMax
    //charge
    //iso0
    //iso1

    cutset.cuts[0].name = Form("recoLead.pt > %6.2f", cLeadPtMin);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cLeadPtMin;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("recoSublead.pt > %6.2f", cSubleadPtMin);
    cutset.cuts[1].bins = 200;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cSubleadPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("recoMu1.eta < %6.3f", cMaxEta);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 3;
    cutset.cuts[2].cutvalue = &cMaxEta;
    cutset.cuts[2].ismin = false;

    cutset.cuts[3].name = Form("recoMu2.eta < %6.3f", cMaxEta);
    cutset.cuts[3].bins = 50;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 3;
    cutset.cuts[3].cutvalue = &cMaxEta;
    cutset.cuts[3].ismin = false;

    cutset.cuts[4].name = Form("dimuCand.mass > %6.2f", cDimuMassMin);
    cutset.cuts[4].bins = 100;
    cutset.cuts[4].min = 110;
    cutset.cuts[4].max = 310;
    cutset.cuts[4].cutvalue = &cDimuMassMin;
    cutset.cuts[4].ismin = true;

    cutset.cuts[5].name = Form("dimuCand.mass < %6.2f", cDimuMassMax);
    cutset.cuts[5].bins = 100;
    cutset.cuts[5].min = 110;
    cutset.cuts[5].max = 310;
    cutset.cuts[5].cutvalue = &cDimuMassMax;
    cutset.cuts[5].ismin = false;

    cutset.cuts[6].name = "recoMu1.charge != recoMu2.charge";
    cutset.cuts[6].bins = 2;
    cutset.cuts[6].min = 0;
    cutset.cuts[6].max = 2;
    cutset.cuts[6].ismin = true;

    cutset.cuts[7].name = Form("recoMu1.iso < %5.3f");
    cutset.cuts[7].bins = 50;
    cutset.cuts[7].min = 0;
    cutset.cuts[7].max = 1;
    cutset.cuts[7].cutvalue = &cMaxRelIso;
    cutset.cuts[7].ismin = false;

    cutset.cuts[8].name = Form("recoMu2.iso < %5.3f");
    cutset.cuts[8].bins = 50;
    cutset.cuts[8].min = 0;
    cutset.cuts[8].max = 1;
    cutset.cuts[8].cutvalue = &cMaxRelIso;
    cutset.cuts[8].ismin = false;
}

