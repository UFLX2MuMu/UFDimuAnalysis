///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// VarSet.cxx                                                            //
// ======================================================================//
// TTree variables are loaded into these objects. Each sample has        //
// its own VarSet data structure.                                        //
// We also set up a map<TString,FUNCTION> to get feature values based    //
// upon the name of the feature.                                         //
// Has muons, dimuons, jets, electrons, gen muons, etc...                //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

// most of the methods are short and defined in the header

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "VarSet.h"

///////////////////////////////////////////////////////////////////////////
// _______________________Constructor/Destructor_________________________//
///////////////////////////////////////////////////////////////////////////

VarSet::VarSet() 
{

  TString str;

  varMap["nVertices"]     = &VarSet::_nVertices;

  varMap["bdt_score"]     = &VarSet::bdt_score;
  varMap["bdt_vbf_score"] = &VarSet::bdt_vbf_score;
  varMap["bdt_ggh_score"] = &VarSet::bdt_ggh_score;
  varMap["bdt_vh_score"]  = &VarSet::bdt_vh_score;
  varMap["bdt_ewk_score"] = &VarSet::bdt_ewk_score;
  varMap["bdt_top_score"] = &VarSet::bdt_top_score;

  // Object counting
  varMap["nJets"]        = &VarSet::_nJets;
  varMap["nJetsCent"]    = &VarSet::_nJetsCent;
  varMap["nJetsFwd"]     = &VarSet::_nJetsFwd;
  varMap["nBLoose"]      = &VarSet::_nBLoose;
  varMap["nBMed"]        = &VarSet::_nBMed;
  varMap["nBTight"]      = &VarSet::_nBTight;
  varMap["nValJets"]     = &VarSet::nValJets;
  // varMap["nValJetsCent"] = &VarSet::nValJetsCent;
  // varMap["nValJetsFwd"]  = &VarSet::nValJetsFwd;
  varMap["nValBTags"]    = &VarSet::nValBTags;
  varMap["nExtraLep"]    = &VarSet::nExtraLep;
  varMap["nExtraMu"]     = &VarSet::nExtraMu;
  varMap["nEle"]         = &VarSet::nEle;
    
  // Dimuon variables
  varMap["dimu_mass"]        = &VarSet::dimu_mass;
  varMap["dimu_mass_Roch"]   = &VarSet::dimu_mass_Roch;
  varMap["dimu_mass_KaMu"]   = &VarSet::dimu_mass_KaMu;
  varMap["massErr_PF"]       = &VarSet::massErr_PF;
  varMap["massErr_Roch"]     = &VarSet::massErr_Roch;
  varMap["massErr_KaMu"]     = &VarSet::massErr_KaMu;

  varMap["dimu_pt"]          = &VarSet::dimu_pt;
  varMap["dimu_eta"]         = &VarSet::dimu_eta;
  varMap["dimu_abs_eta"]     = &VarSet::dimu_abs_eta;
  varMap["dimu_rapid"]       = &VarSet::dimu_rapid;
  varMap["dimu_dR"]          = &VarSet::dimu_dR;
  varMap["dimu_dEta"]        = &VarSet::dimu_dEta;
  varMap["dimu_abs_dEta"]    = &VarSet::dimu_abs_dEta;
  varMap["dimu_dPhi"]        = &VarSet::dimu_dPhi;
  varMap["dimu_abs_dPhi"]    = &VarSet::dimu_abs_dPhi;
  varMap["dimu_dPhiStar"]    = &VarSet::dimu_dPhiStar;
  varMap["dimu_avg_abs_eta"] = &VarSet::dimu_avg_abs_eta;
  varMap["dimu_min_abs_eta"] = &VarSet::dimu_min_abs_eta;
  varMap["dimu_max_abs_eta"] = &VarSet::dimu_max_abs_eta;
  
  // Muon variables from dimuon pair
  varMap["mu1_pt"]      = &VarSet::mu1_pt;
  varMap["mu2_pt"]      = &VarSet::mu2_pt;
  varMap["mu1_eta"]     = &VarSet::mu1_eta;
  varMap["mu2_eta"]     = &VarSet::mu2_eta;
  varMap["mu1_abs_eta"] = &VarSet::mu1_abs_eta;
  varMap["mu2_abs_eta"] = &VarSet::mu2_abs_eta;

  // Dijet variables
  for (int i = 1; i <= N_JET_PAIRS; i++) {
    str.Form("dijet%d_", i);
    varMapI[(str+"mass").Data()]        = &VarSet::dijet_mass;
    varMapI[(str+"pt").Data()]          = &VarSet::dijet_pt;
    varMapI[(str+"dEta").Data()]        = &VarSet::dijet_dEta;
    varMapI[(str+"abs_dEta").Data()]    = &VarSet::dijet_abs_dEta;
    varMapI[(str+"min_abs_eta").Data()] = &VarSet::dijet_min_abs_eta;
    varMapI[(str+"max_abs_eta").Data()] = &VarSet::dijet_max_abs_eta;
  }

  // Jet variables
  for (int i = 1; i <= N_JETS; i++) {
    str.Form("jet%d_", i);
    varMapI[(str+"pt").Data()]      = &VarSet::jet_pt;
    varMapI[(str+"eta").Data()]     = &VarSet::jet_eta;
    varMapI[(str+"abs_eta").Data()] = &VarSet::jet_abs_eta;
  }

  // Global event variables
  varMap["MET"]      = &VarSet::MET;
  varMap["MHT"]      = &VarSet::MHT;
  varMap["MT_had"]   = &VarSet::MT_had;
  varMap["mass_had"] = &VarSet::mass_had;  


  // AMC variables
  varMap["jj_jet0_pt"]  = &VarSet::jj_jet0_pt;
  varMap["jj_jet1_pt"]  = &VarSet::jj_jet1_pt;
  varMap["jj_jet0_eta"] = &VarSet::jj_jet0_eta;
  varMap["jj_jet1_eta"] = &VarSet::jj_jet1_eta;

  varMap["m_jj"]          = &VarSet::m_jj;
  varMap["dEta_jj"]       = &VarSet::dEta_jj;
  varMap["dEta_jj_mumu"]  = &VarSet::dEta_jj_mumu;
  varMap["dPhi_jj_mumu"]  = &VarSet::dPhi_jj_mumu;
  varMap["zep"]           = &VarSet::zep;
  
  varMap["vbf_jet0_pt"]   = &VarSet::vbf_jet0_pt;
  varMap["vbf_jet1_pt"]   = &VarSet::vbf_jet1_pt;
  varMap["vbf_jet0_eta"]  = &VarSet::vbf_jet0_eta;
  varMap["vbf_jet1_eta"]  = &VarSet::vbf_jet1_eta;
  
  varMap["vbf_m_jj"]          = &VarSet::vbf_m_jj;
  varMap["vbf_dEta_jj"]       = &VarSet::vbf_dEta_jj;
  varMap["vbf_dEta_jj_mumu"]  = &VarSet::vbf_dEta_jj_mumu;
  varMap["vbf_dPhi_jj_mumu"]  = &VarSet::vbf_dPhi_jj_mumu;
  varMap["vbf_zep"]           = &VarSet::vbf_zep;
  
  varMap["bjet0_pt"]   = &VarSet::bjet0_pt;
  varMap["bjet1_pt"]   = &VarSet::bjet1_pt;
  varMap["bjet0_eta"]  = &VarSet::bjet0_eta;
  varMap["bjet1_eta"]  = &VarSet::bjet1_eta;
    
  varMap["m_bb"]     = &VarSet::m_bb;
  varMap["dEta_bb"]  = &VarSet::dEta_bb;
  
  varMap["extra_muon0_pt"]  = &VarSet::extra_muon0_pt; 
  varMap["extra_muon1_pt"]  = &VarSet::extra_muon1_pt; 
  varMap["extra_muon0_eta"] = &VarSet::extra_muon0_eta; 
  varMap["extra_muon1_eta"] = &VarSet::extra_muon1_eta; 
  
  varMap["electron0_pt"]  = &VarSet::electron0_pt;  
  varMap["electron1_pt"]  = &VarSet::electron1_pt; 
  varMap["electron0_eta"] = &VarSet::electron0_eta; 
  varMap["electron1_eta"] = &VarSet::electron1_eta; 
  
  varMap["mT_b_MET"] = &VarSet::mT_b_MET;
  
}
