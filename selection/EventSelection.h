///////////////////////////////////////////////////////////////////////////
//                             EventSelection.h                          //
//=======================================================================//
//                                                                       //
//        Some selections to cut events. Implements Cut.h.               //
//        Cut based on trigger matching, or other event info.            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_EVENTSELECTIONCUTS
#define ADD_EVENTSELECTIONCUTS

#include "Cut.h"

class Run2EventSelectionCuts : public Cut
{
// if data apply hlt trigger matching
// if MC do not apply hlt trigger matching
    public:
        Run2EventSelectionCuts();
        Run2EventSelectionCuts(float trigMuPtMin, float dimuMassMin);

        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >

        void makeCutSet();
        bool evaluate(VarSet& vars);
        TString string();
};


class EventSelectionCuts_2017 : public Cut
{
// if data apply hlt trigger matching
// if MC do not apply hlt trigger matching
    public:
        EventSelectionCuts_2017();
        EventSelectionCuts_2017(float trigMuPtMin, float dimuMassMin);

        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >

        void makeCutSet();
        bool evaluate(VarSet& vars);
        TString string();
};



class SynchEventSelectionCuts : public Cut
{
    public:
        SynchEventSelectionCuts();
        SynchEventSelectionCuts(float trigMuPtMin, float dimuMassMin);

        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >

        void makeCutSet();
        bool evaluate(VarSet& vars);
        TString string();
};

class FEWZCompareCuts : public Cut
{
    public:
        FEWZCompareCuts();
        FEWZCompareCuts(bool useReco);
        FEWZCompareCuts(bool useReco, float leadPtMin, float subleadPtMin, float maxEta, float dimuMassMin, float dimuMassMax, float maxRelIso);

        bool useReco;

        float cLeadPtMin;           
        float cSubleadPtMin;        
        float cMaxEta;              
        float cDimuMassMin;         
        float cDimuMassMax;         
        float cMaxRelIso;         

        void makeCutSet();
        bool evaluate(VarSet& vars);
        bool evaluate(_MuonInfo& recoMu, int m);
        TString string();
};
#endif
