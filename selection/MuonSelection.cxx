///////////////////////////////////////////////////////////////////////////
//                             MuonSelection.cxx                         //
//=======================================================================//
//                                                                       //
//        Eliminate events by cutting on individual muon info.           //
//        The dimuon candidate muons must pass these criteria.           //
//        Usually cuts based on muon eta, pt, and iso.                   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "MuonSelection.h"
#include "ParticleTools.h"
#include "TMath.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run2MuonSelection______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// muon selection based on run 1, updated for run 2

Run2MuonSelectionCuts::Run2MuonSelectionCuts()
{
// Default muon selection values for the synchronization exercise

    cMinPt = 10;       // >, originally 15
    cMaxEta = 2.4;     // <, originally 2.1
    cMaxRelIso = 0.25; // <
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run2MuonSelectionCuts::Run2MuonSelectionCuts(float minPt, float maxEta, float maxRelIso)
{
// Custom muon selection values for the synchronization exercise

    cMinPt = minPt;         // >
    cMaxEta = maxEta;       // <
    cMaxRelIso = maxRelIso; // <
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run2MuonSelectionCuts::evaluate(VarSet& vars)
{
    // cuts for first muon
    MuonInfo& mu1 = vars.muons->at(vars.dimuCand->iMu1); 
    if(cutset.cuts[0].on)
        if(!(mu1.pt > cMinPt))               return false;

    if(cutset.cuts[1].on)
        if(!(TMath::Abs(mu1.eta) < cMaxEta)) return false;

    if(cutset.cuts[2].on)
        if(!(mu1.iso() <= cMaxRelIso))       return false;


    // cuts for the second muon
    MuonInfo& mu2 = vars.muons->at(vars.dimuCand->iMu2); 
    if(cutset.cuts[3].on)
        if(!(mu2.pt > cMinPt))               return false;

    if(cutset.cuts[4].on)
        if(!(TMath::Abs(mu2.eta) < cMaxEta)) return false;

    if(cutset.cuts[5].on)
        if(!(mu2.iso() <= cMaxRelIso))       return false;
    

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run2MuonSelectionCuts::string()
{
    return TString("Run2_Muon_Selection_Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run2MuonSelectionCuts::makeCutSet()
{
    // recoMuon1 cuts
    cutset.cuts[0].name = Form("recoMu1.pt > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cMinPt;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("|recoMu1.eta| < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 3;
    cutset.cuts[1].cutvalue = &cMaxEta;
    cutset.cuts[1].ismin = false;

    cutset.cuts[2].name = Form("recoMu1.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 1;
    cutset.cuts[2].cutvalue = &cMaxRelIso;
    cutset.cuts[2].ismin = false;

     // recoMuon2 cuts
    cutset.cuts[3].name = Form("recoMu2.pt > %4.1f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;
    cutset.cuts[3].cutvalue = &cMinPt;
    cutset.cuts[3].ismin = true;

    cutset.cuts[4].name = Form("|recoMu2.eta| < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 3;
    cutset.cuts[4].cutvalue = &cMaxEta;
    cutset.cuts[4].ismin = false;

    cutset.cuts[5].name = Form("recoMu2.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[5].bins = 50;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 1;
    cutset.cuts[5].cutvalue = &cMaxRelIso;
    cutset.cuts[5].ismin = false;

    // combine the names and cuts in tstring form, for all cut sin the set  
    // to set strings with all of the cut information in one
    cutset.concatCuts();
}




















///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________MuonSelection_2017______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// muon selection based on run 1, updated for run 2

MuonSelectionCuts_2017::MuonSelectionCuts_2017()
{
// Default muon selection values for the synchronization exercise

    cMinPt = 20;       // >, originally 15
    cMaxEta = 2.4;     // <, originally 2.1
    cMaxRelIso = 0.25; // <
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

MuonSelectionCuts_2017::MuonSelectionCuts_2017(float minPt, float maxEta, float maxRelIso)
{
// Custom muon selection values for the synchronization exercise

    cMinPt = minPt;         // >
    cMaxEta = maxEta;       // <
    cMaxRelIso = maxRelIso; // <
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool MuonSelectionCuts_2017::evaluate(VarSet& vars)
{
    // cuts for first muon
    MuonInfo& mu1 = vars.muons->at(vars.dimuCand->iMu1); 
    if(cutset.cuts[0].on)
        if(!(mu1.pt > cMinPt))               return false;

    if(cutset.cuts[1].on)
        if(!(TMath::Abs(mu1.eta) < cMaxEta)) return false;

    if(cutset.cuts[2].on)
        if(!(mu1.iso() <= cMaxRelIso))       return false;


    // cuts for the second muon
    MuonInfo& mu2 = vars.muons->at(vars.dimuCand->iMu2); 
    if(cutset.cuts[3].on)
        if(!(mu2.pt > cMinPt))               return false;

    if(cutset.cuts[4].on)
        if(!(TMath::Abs(mu2.eta) < cMaxEta)) return false;

    if(cutset.cuts[5].on)
        if(!(mu2.iso() <= cMaxRelIso))       return false;
    

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString MuonSelectionCuts_2017::string()
{
    return TString("Muon_Selection_2017_Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void MuonSelectionCuts_2017::makeCutSet()
{
    // recoMuon1 cuts
    cutset.cuts[0].name = Form("recoMu1.pt > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cMinPt;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("|recoMu1.eta| < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 3;
    cutset.cuts[1].cutvalue = &cMaxEta;
    cutset.cuts[1].ismin = false;

    cutset.cuts[2].name = Form("recoMu1.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 1;
    cutset.cuts[2].cutvalue = &cMaxRelIso;
    cutset.cuts[2].ismin = false;

     // recoMuon2 cuts
    cutset.cuts[3].name = Form("recoMu2.pt > %4.1f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;
    cutset.cuts[3].cutvalue = &cMinPt;
    cutset.cuts[3].ismin = true;

    cutset.cuts[4].name = Form("|recoMu2.eta| < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 3;
    cutset.cuts[4].cutvalue = &cMaxEta;
    cutset.cuts[4].ismin = false;

    cutset.cuts[5].name = Form("recoMu2.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[5].bins = 50;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 1;
    cutset.cuts[5].cutvalue = &cMaxRelIso;
    cutset.cuts[5].ismin = false;

    // combine the names and cuts in tstring form, for all cut sin the set  
    // to set strings with all of the cut information in one
    cutset.concatCuts();
}


