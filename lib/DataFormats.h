#ifndef DATAFORMATS
#define DATAFORMATS

#include "TString.h"

double static constexpr MASS_MUON     = 0.105658367;     //GeV/c2
double static constexpr MASS_ELECTRON = 0.0005109989461; //GeV/c2
double static constexpr MASS_TAU      = 1.77682;         //GeV/c2

struct _EventInfo{

  // data fields
  int run;
  int lumi;
  long long int event;
  int bx;
  int orbit;

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      return TString("run/I:lumi/I:event/L:bx/I:orbit/I");
  };

  void init()
  {
      run = -999;
      lumi = -999;
      event = -999;
      bx = -999;
      orbit = -999;
  };

};

// vertex info
struct _VertexInfo{

  // data fields
  const static unsigned int arraySize = 20;
  int nVertices;
  int isValid[arraySize];
  float x[arraySize];	
  float y[arraySize];	
  float z[arraySize];	
  float xErr[arraySize];	
  float yErr[arraySize];	
  float zErr[arraySize];	
  float chi2[arraySize];
  int ndf[arraySize];
  float normChi2[arraySize];

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r =   TString("nVertices/I:isValid[N]/I:")+
                    TString("x[N]/F:y[N]/F:z[N]/F:xErr[N]/F:yErr[N]/F:zErr[N]/F:")+
                    TString("chi2[N]/F:ndf[N]/I:normChi2[N]/F");
      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      return r;
  };

  void init()
  {
      nVertices = -999;
      for (unsigned int i=0;i<arraySize;i++) 
      {
          isValid[i]  = -999;
          x[i]        = -999;
          y[i]        = -999;
          z[i]        = -999;
          xErr[i]     = -999;
          yErr[i]     = -999;
          zErr[i]     = -999;
          chi2[i]     = -999;
          ndf[i]      = -999;
          normChi2[i] = -999;
      }
  };
};


struct _DimuCandInfo{

  float recoCandMass;
  float recoCandPt;
  float recoCandEta;
  float recoCandY;
  float recoCandPhi;

  float recoCandMassPF;
  float recoCandPtPF;
  float recoCandEtaPF;
  float recoCandYPF;
  float recoCandPhiPF;

  float angleDiMuons;

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      return TString("recoCandMass/F:recoCandPt/F:recoCandEta/F:recoCandY/F:recoCandPhi/F:")+
             TString("recoCandMassPF/F:recoCandPtPF/F:recoCandEtaPF/F:recoCandYPF/F:recoCandPhiPF/F:angleDiMuons/F");
  };

  void init()
  {
      recoCandMass = -999;
      recoCandPt   = -999;
      recoCandEta  = -999;
      recoCandY    = -999;
      recoCandPhi  = -999;

      recoCandMassPF = -999;
      recoCandPtPF   = -999;
      recoCandEtaPF  = -999;
      recoCandYPF    = -999;
      recoCandPhiPF  = -999;

      angleDiMuons = -999;
  };

};

// muon info
struct _MuonInfo{

  const static unsigned int arraySize = 10;
  const static unsigned int triggerArraySize = 8;
  int nMuons;
  int nMuonPairs;

  int isTracker[arraySize];
  int isStandAlone[arraySize];
  int isGlobal[arraySize];

  int isTightMuon[arraySize];
  int isMediumMuon[arraySize];
  int isLooseMuon[arraySize];

  int charge[arraySize];
  float pt[arraySize];
  float ptErr[arraySize];
  float eta[arraySize];
  float phi[arraySize];

  float trkPt[arraySize];
  float trkPtErr[arraySize];
  float trketa[arraySize];
  float trkPhi[arraySize];

  float d0_BS[arraySize];
  float dz_BS[arraySize];

  float d0_PV[arraySize];
  float dz_PV[arraySize];

  float trackIsoSumPt[arraySize];
  float trackIsoSumPtCorr[arraySize];
  float hcalIso[arraySize];
  float ecalIso[arraySize];
  float relCombIso[arraySize];

  // PF information
  int isPFMuon[arraySize];

  float pfPt[arraySize];
  float pfEta[arraySize];
  float pfPhi[arraySize];
  
  float sumChargedHadronPtR03[arraySize];   // sum-pt of charged Hadron 
  float sumChargedParticlePtR03[arraySize]; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03[arraySize];   // sum pt of neutral hadrons
  float sumPhotonEtR03[arraySize];          // sum pt of PF photons
  float sumPUPtR03[arraySize];              // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04[arraySize]; 
  float sumChargedParticlePtR04[arraySize];
  float sumNeutralHadronEtR04[arraySize];  
  float sumPhotonEtR04[arraySize];
  float sumPUPtR04[arraySize];

  int isHltMatched[arraySize][triggerArraySize];
  float hltPt[arraySize][triggerArraySize];
  float hltEta[arraySize][triggerArraySize];
  float hltPhi[arraySize][triggerArraySize];

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r =  TString("nMuons/I:nMuonPairs/I:")+
                   TString("isTracker[N]/I:isStandAlone[N]/I:isGlobal[N]/I:")+
                   TString("isTightMuon[N]/I:isMediumMuon[N]/I:isLooseMuon[N]/I:")+
                   TString("charge[N]/I:pt[N]/F:ptErr[N]/F:eta[N]/F:phi[N]/F:")+
                   TString("trkPt[N]/F:trkPtErr[N]/F:trkEta[N]/F:trkPhi[N]/F:")+
                   TString("d0_BS[N]/F:dz_BS[N]/F:")+
                   TString("d0_PV[N]/F:dz_PV[N]/F:")+
                   TString("trackIsoSumPt[N]/F:")+
                   TString("trackIsoSumPtCorr[N]/F:")+
                   TString("hcalIso[N]/F:")+
                   TString("ecalIso[N]/F:")+
                   TString("relCombIso[N]/F:")+
                   TString("isPFMuon[N]/I:")+
                   TString("pfPt[N]/F:")+
                   TString("pfEta[N]/F:")+
                   TString("pfPhi[N]/F:")+
                   TString("sumChargedHadronPtR03[N]/F:")+
                   TString("sumChargedParticlePtR03[N]/F:")+
                   TString("sumNeutralHadronEtR03[N]/F:")+
                   TString("sumPhotonEtR03[N]/F:")+
                   TString("sumPUPtR03[N]/F:")+
                   TString("sumChargedHadronPtR04[N]/F:")+
                   TString("sumChargedParticlePtR04[N]/F:")+
                   TString("sumNeutralHadronEtR04[N]/F:")+
                   TString("sumPhotonEtR04[N]/F:")+
                   TString("sumPUPtR04[N]/F:")+
                   TString("isHltMatched[N][T]/I:")+
                   TString("hltPt[N][T]/F:")+
                   TString("hltEta[N][T]/F:")+
                   TString("hltPhi[N][T]/F");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      r.ReplaceAll("[T]",Form("[%d]", triggerArraySize));
      return r;
  };

  void init()
  {
    nMuons = -999;
    nMuonPairs = -999;

    for(unsigned int i=0; i<arraySize; i++)
    {
      isTracker[i]    = -999;
      isStandAlone[i] = -999;
      isGlobal[i]     = -999;

      isTightMuon[i]    = -999;
      isMediumMuon[i]   = -999;
      isLooseMuon[i]    = -999;

      charge[i] = -999;
      pt[i]     = -999;
      eta[i]    = -999;
      phi[i]    = -999;

      d0_BS[i] = -999;
      dz_BS[i] = -999;

      d0_PV[i] = -999;
      dz_PV[i] = -999;

      trackIsoSumPt[i]     = -999;
      trackIsoSumPtCorr[i] = -999;
      hcalIso[i]           = -999;
      ecalIso[i]           = -999;
      relCombIso[i]        = -999;

      isPFMuon[i] = -999;

      pfPt[i]  = -999;
      pfEta[i] = -999;
      pfPhi[i] = -999;

      sumChargedHadronPtR03[i]   = -999;
      sumChargedParticlePtR03[i] = -999;
      sumNeutralHadronEtR03[i]   = -999;
      sumPhotonEtR03[i]          = -999;
      sumPUPtR03[i]              = -999;

      sumChargedHadronPtR04[i]   = -999;
      sumChargedParticlePtR04[i] = -999;
      sumNeutralHadronEtR04[i]   = -999;
      sumPhotonEtR04[i]          = -999;
      sumPUPtR04[i]              = -999;

      for (unsigned int iTrigger=0;iTrigger<triggerArraySize;iTrigger++) 
      {
        isHltMatched[i][iTrigger] = -999;
        hltPt[i][iTrigger] = -999;
        hltEta[i][iTrigger] = -999;
        hltPhi[i][iTrigger] = -999;
      }
    }

  };

};

// electron info
struct _ElectronInfo{

  const static unsigned int arraySize = 10;
  int nElectrons;

  // electron cut based IDs
  int isTightElectron[arraySize];
  int isMediumElectron[arraySize];
  int isLooseElectron[arraySize];
  int isVetoElectron[arraySize];
  int passConversionVeto[arraySize];

  int charge[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];

  float d0_PV[arraySize];
  float dz_PV[arraySize];
  float missingInnerHits[arraySize];
 
  int isPFElectron[arraySize]; 

  float sumChargedHadronPtR03[arraySize];   // sum-pt of charged Hadron 
  float sumNeutralHadronEtR03[arraySize];   // sum pt of neutral hadrons
  float sumPhotonEtR03[arraySize];          // sum pt of PF photons
  float sumPUPtR03[arraySize];              // sum pt of charged Particles not from PV  (for Pu corrections)

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r = TString("nElectrons/I:")+
                  TString("isTightElectron[N]/I:isMediumElectron[N]/I:isLooseElectron[N]/I:isVetoElectron[N]/I:")+
                  TString("passConversionVeto[N]/I:charge[N]/I:pt[N]/F:eta[N]/F:phi[N]/F:")+
                  TString("d0_PV[N]/F:dz_PV[N]/F:missingInnerHits[N]/F:")+
                  TString("isPFElectron[N]/I:")+
                  TString("sumChargedHadronPtR03[N]/F:")+
                  TString("sumNeutralHadronEtR03[N]/F:")+
                  TString("sumPhotonEtR03[N]/F:")+
                  TString("sumPUPtR03[N]/F");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      return r;
  };

  void init()
  {
    nElectrons = -999;

    for(unsigned int i=0; i<arraySize; i++)
    {
        isTightElectron[i]    = -999;
        isMediumElectron[i]   = -999;
        isLooseElectron[i]    = -999;
        isVetoElectron[i]     = -999;
        passConversionVeto[i] = -999;

        charge[i] = -999;
        pt[i]     = -999;
        eta[i]    = -999;
        phi[i]    = -999;

        d0_PV[i] = -999;
        dz_PV[i] = -999;

        missingInnerHits[i] = -999;

        isPFElectron[i]            = -999;
        sumChargedHadronPtR03[i]   = -999;
        sumNeutralHadronEtR03[i]   = -999;
        sumPhotonEtR03[i]          = -999;
        sumPUPtR03[i]              = -999;
    }

  };

};

// tau info
struct _TauInfo{

  const static unsigned int arraySize = 10;
  const static unsigned int idArraySize = 16;
  int nTaus;

  // tau IDs
  // accessed like so in CMSSW, pat::Tau::tauID("name")
  float tauID[arraySize][idArraySize];

  int charge[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];

  float d0_PV[arraySize];
  float dz_PV[arraySize];
 
  int isPFTau[arraySize]; 

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r = TString("nTaus/I:tauID[N][I]/F:")+
                  TString("charge[N]/I:pt[N]/F:eta[N]/F:phi[N]/F:")+
                  TString("d0_PV[N]/F:dz_PV[N]/F:")+
                  TString("isPFTau[N]/I");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      r.ReplaceAll("[I]",Form("[%d]", idArraySize));
      return r;
  };

  void init()
  {
      nTaus = -999;
      for(unsigned int i=0; i<arraySize; i++)
      {
          for(unsigned int id=0; id<idArraySize; id++)
          { 
              tauID[i][id] = -999;
          }

          charge[i] = -999;
          pt[i] = -999;
          eta[i] = -999;
          phi[i] = -999;

          d0_PV[i] = -999;
          dz_PV[i] = -999;
 
          isPFTau[i] = -999; 
      };
  };
};

//MET
struct _MetInfo{
  float px;
  float py;
  float pt;
  float phi;
  float sumEt;

  static TString getVarString()
  {
      return TString("px/F:py/F:pt/F:phi/F:sumEt/F");
  };

  void init()
  {
      px = -999;
      py = -999;
      pt = -999;
      phi = -999;
      sumEt = -999;
  };

};


// pf Jets
struct _PFJetInfo{

  const static unsigned int arraySize = 10;
  int   nJets;
  float px[arraySize];
  float py[arraySize];
  float pz[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];
  float mass[arraySize];
  int   charge[arraySize];
  float isB[arraySize];
  int   partonFlavour[arraySize];
  /////// Energy Fractions //////
  //Charged Hadron
  float chf[arraySize];
  //NeutralHadron
  float nhf[arraySize];
  //Charged EM
  float cef[arraySize];
  //Neutral EM
  float nef[arraySize];
  //Mu
  float muf[arraySize];
  // HF Hadron Fraction
  float hfhf[arraySize];
  // HF EM Fraction
  float hfef[arraySize];
  /////// Multiplicities //////
  // Total Charged Mult
  int cm[arraySize];
  //Charged Hadron Mult
  int chm[arraySize];
  //Total Neutral Mult
  int nm[arraySize];
  //NeutralHadron Mult
  int nhm[arraySize];
  //Charged EM Mult
  int cem[arraySize];
  //Neutral EM Mult
  int nem[arraySize];
  //Mu Mult
  int mum[arraySize];
  // HF Hadron Mult
  int hfhm[arraySize];
  // HF EM Mult
  int hfem[arraySize];
  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  float jecFactor[arraySize];
  // Jet Energy Correction Uncertainty
  float jecUnc[arraySize];
  // Gen Jet Values
  int genMatched[arraySize];
  float genPx[arraySize];
  float genPy[arraySize];
  float genPz[arraySize];
  float genPt[arraySize];
  float genEta[arraySize];
  float genPhi[arraySize];
  float genMass[arraySize];
  ///// Gen Jet Energy Fractions ///////
  // EM Fraction
  float genEMF[arraySize];
  // Had Fraction
  float genHadF[arraySize];
  // Invisible Fraction
  float genInvF[arraySize];
  // Auxiliary Fraction (Undecayed Sigmas, etc.)
  float genAuxF[arraySize];
  // PUID
  float puid[arraySize];

  static TString getVarString()
  {
       TString r = TString("nJets/I:px[N]/F:py[N]/F:pz[N]/F:pt[N]/F:eta[N]/F:")+
                   TString("phi[N]/F:mass[N]/F:charge[N]/I:isB[N]/F:partonFlavour[N]/I:chf[N]/F:")+
                   TString("nhf[N]/F:cef[N]/F:nef[N]/F:muf[N]/F:hfhf[N]/F:hfef[N]/F:")+
                   TString("cm[N]/I:chm[N]/I:nm[N]/I:nhm[N]/I:cem[N]/I:nem[N]/I:mum[N]/I:hfhm[N]/I:")+
                   TString("hfem[N]/I:jecFactor[N]/F:jecUnc[N]/F:genMatched[N]/I:genPx[N]/F:genPy[N]/F:")+
                   TString("genPz[N]/F:genPt[N]/F:genEta[N]/F:genPhi[N]/F:genMass[N]/F:genEMF[N]/F:")+
                   TString("genHadF[N]/F:genInvF[N]/F:genAux[N]/F:puid[N]/F");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      return r;
  };

  void init()
  {

      nJets = -999;

      for(unsigned int i=0; i<arraySize; i++)
      { 
          px[i] = -999;
          py[i] = -999;
          pz[i] = -999;
          pt[i] = -999;
          eta[i] = -999;
          phi[i] = -999;
          mass[i] = -999;
          charge[i] = -999;
          isB[i] = -999;
          partonFlavour[i] = -999;
          /////// Energy Fractions //////
          //Charged Hadron
          chf[i] = -999;
          //NeutralHadron
          nhf[i] = -999;
          //Charged EM
          cef[i] = -999;
          //Neutral EM
          nef[i] = -999;
          //Mu
          muf[i] = -999;
          // HF Hadron Fraction
          hfhf[i] = -999;
          // HF EM Fraction
          hfef[i] = -999;
          /////// Multiplicities //////
          // Total Charged Mult
          cm[i] = -999;
          //Charged Hadron Mult
          chm[i] = -999;
          //Total Neutral Mult
          nm[i] = -999;
          //NeutralHadron Mult
          nhm[i] = -999;
          //Charged EM Mult
          cem[i] = -999;
          //Neutral EM Mult
          nem[i] = -999;
          //Mu Mult
          mum[i] = -999;
          // HF Hadron Mult
          hfhm[i] = -999;
          // HF EM Mult
          hfem[i] = -999;
          // Jet Correction Factor--Above momentum is already corrected!!
          // This factor will return momentum to uncorrected value!!
          jecFactor[i] = -999;
          // Jet Energy Correction Uncertainty
          jecUnc[i] = -999;
          // Gen Jet Values
          genMatched[i] = -999;
          genPx[i] = -999;
          genPy[i] = -999;
          genPz[i] = -999;
          genPt[i] = -999;
          genEta[i] = -999;
          genPhi[i] = -999;
          genMass[i] = -999;
          ///// Gen Jet Energy Fractions ///////
          // EM Fraction
          genEMF[i] = -999;
          // Had Fraction
          genHadF[i] = -999;
          // Invisible Fraction
          genInvF[i] = -999;
          // Auxiliary Fraction (Undecayed Sigmas, etc.)
          genAuxF[i] = -999;
          // PUID
          puid[i] = -999;
      }

  };
};

// generator level jets
struct _GenJetInfo{
  const static unsigned int arraySize = 10;
  int nJets;
  float px[arraySize];
  float py[arraySize];
  float pz[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];
  float mass[arraySize];
  int   charge[arraySize];

  static TString getVarString()
  {
    TString r = TString("nJets/I:px[N]/F:py[N]/F:pz[N]/F:pt[N]/F:eta[N]/F:phi[N]/F:mass[N]/F:charge[N]/I");
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  };

  void init()
  {
      nJets = -999;
      for(unsigned int i=0; i<arraySize; i++)
      {
          px[i] = -999;
          py[i] = -999;
          pz[i] = -999;
          pt[i] = -999;
          eta[i] = -999;
          phi[i] = -999;
          mass[i] = -999;
          charge[i] = -999;
      }

  };

};

// used for gen particles that aren't composite
// not sure why we have two objects for gen particles ...
struct _TrackInfo{
  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      return TString("charge/I:pt/F:ptErr/F:eta/F:phi/F");
  };

  void init()
  {
      charge = -999;
      pt = -999;
      ptErr = -999;
      eta = -999;
      phi = -999;  
  };

};

// generator level composite Candidate
struct _genPartInfo{
  int charge;
  float mass;
  float pt;
  float eta;  // pseudo rapidity
  float y;    // rapidity
  float phi;  // phi

  static TString getVarString()
  {
      return TString("charge/I:mass/F:pt/F:eta/F:y/F:phi/F");
  };

  void init()
  {
      charge = -999;
      mass = -999;
      pt = -999;
      eta = -999;
      y = -999;
      phi = -999;  
  };

};
#endif
