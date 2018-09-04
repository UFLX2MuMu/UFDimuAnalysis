/////////////////////////////////////////////////////////////////////////////
//                           SampleDatabase.cxx                            //
//=========================================================================//
//                                                                         //
// Load TTrees from lxplus (location = "CERN") or from UF HPC/IHPEA        //
// (location = "UF").                                                      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "DiMuPlottingSystem.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

std::vector<Sample*>& GetSamples(std::map<TString, Sample*>& samples, TString location, TString select="ALL", bool info_only=false) {
  
  std::cout << "\n======== Getting samples: " << select << "\n" << std::endl;
  
  ///////////////////////////////////////////////////////////////////
  // SAMPLES---------------------------------------------------------
  ///////////////////////////////////////////////////////////////////
  
  TString in_dir;
  TString data_dir;
  if (location == "UF")
    in_dir = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/"; 
  else if (location == "UF_DoubleMu")
  {
    in_dir = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/"; 
    data_dir = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/doubleMu_simplified/"; 
  }
  else if (location == "CERN")
    in_dir = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17";
  else if (location == "CERN_hiM")
    in_dir = "root://eoscms.cern.ch//store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM";
  else
    std::cout << "\n\nInput location is " << location << ", not UF, CERN, or CERN_hiM.  NOT AN OPTION!!!\n\n" << std::endl;
  
  TString in_dir_hiM = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/hiM_simplified/";

  // ================================================================
  // Data -----------------------------------------------------------
  // ================================================================

  

  if (select.Contains("CERN_2017_collection"))
  {
    TString in_dir_CERN = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17";
    std::cout << "Adding files for  2017 MC and data from CERN" << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    TString tup;
    TString tuple = "_tuple_";

    //ggH
//    for(int ii=0 ; ii<79 ; ii++)
//    {
//      TString tuple = ("/tuple_%d.root",ii)
//      in_file = Form( in_dir_CERN + "/GluGlu_HToMuMu_M125_13TeV_amcatnloFXFX_pythia8/H2Mu_gg/180720_145916/0000" + "/tuple_%d.root",ii);
//      in_files.push_back(in_file);
//    }
//    samples["H2mu_gg_2017"] = new Sample(in_files, "H2Mu_gg_2017","signal");
//    samples["H2mu_gg_2017"]->xsec = 48.58 * 0.00022;
//    std::cout << "H2Mu_gg_2017 file added" << in_files.size() <<std::endl;
//
//    in_files.clear();
//    in_files.shrink_to_fit();
//
//    //VBF
//    in_file.Form( "%s/VBFH_HToMuMu_M125_13TeV_amcatnloFXFX_pythia8/H2Mu_VBF/NTuple_0.root", in_dir_CERN.Data());
//    in_files.push_back(in_file);
//    samples["H2mu_VBF_2017"] = new Sample(in_files, "H2Mu_VBF_2017","signal");
//    samples["H2mu_VBF_2017"]->xsec = 1;
//    std::cout << "H2Mu_VBF_2017 file added" << in_files.size() <<std::endl;
//
//    in_files.clear();
//    in_files.shrink_to_fit();
//    //ZH
//    in_file.Form( "%s/ZH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_ZH/NTuple_0.root", in_dir_CERN.Data());
//    in_files.push_back(in_file);
//    samples["H2mu_ZH_2017"] = new Sample(in_files, "H2Mu_ZH_2017","signal");
//    samples["H2mu_ZH_2017"]->xsec = 1;
//    std::cout << "H2Mu_ZH_2017 file added" << in_files.size() <<std::endl;
//
//    in_files.clear();
//    in_files.shrink_to_fit();
//    //Wpos
//    in_file.Form( "%s/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/NTuple_0.root", in_dir_CERN.Data());
//    in_files.push_back(in_file);
//    samples["H2mu_WH_pos_2017"] = new Sample(in_files, "H2Mu_WH_pos_2017","signal");
//    samples["H2mu_WH_pos_2017"]->xsec = 1;
//    std::cout << "H2Mu_WH_pos_2017 file added" << in_files.size() <<std::endl;
//
//    in_files.clear();
//    in_files.shrink_to_fit();
//    //Wneg
//    in_file.Form( "%s/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_neg/NTuple_0.root", in_dir_CERN.Data());
//    in_files.push_back(in_file);
//    samples["H2mu_WH_neg_2017"] = new Sample(in_files, "H2Mu_WH_neg_2017","signal");
//    samples["H2mu_WH_neg_2017"]->xsec = 1;
//    std::cout << "H2Mu_WH_neg_2017 file added" << in_files.size() <<std::endl;
//
//    in_files.clear();
//    in_files.shrink_to_fit();
    //ttbar
//    in_file.Form( "%s/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/tt/NTuple_100-200.root", in_dir_CERN.Data());
//    in_files.push_back(in_file);
//    samples["tt_ll_2017"] = new Sample(in_files, "tt_ll_2017","bakcground");
//    samples["tt_ll_2017"]->xsec = 1;
//    std::cout << "tt_ll_2017 file added" << in_files.size() <<std::endl;

//    in_files.clear();
//    in_files.shrink_to_fit();
    //DY
    for(int ii=1 ; ii<552 ; ii++)
    {
	if (ii > 100) continue;      //for testing
	if (ii == 403) continue;
	tup.Form(tuple+"%d",ii);
        in_file.Form( in_dir + "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/180720_151319/0000" + "/tuple_%d.root", ii);
        in_files.push_back(in_file);
	samples["ZJets_AMC_2017"+tup] = new Sample(in_files, "ZJets_AMC_2017" + tup,"background");
	in_files.clear();
	in_files.shrink_to_fit();
    }
//    samples["ZJets_AMC_2017"] = new Sample(in_files, "ZJets_AMC_2017","background");
//    samples["ZJets_AMC_2017"]->xsec = 6025.2;
    std::cout << "ZJets_AMC_2017 file added" << in_files.size() <<std::endl;

    in_files.clear();
    in_files.shrink_to_fit();
    //data
    for(int ii=1 ; ii < 98; ii++)
    {
	if(ii > 20) continue;
	if(ii == 30) continue;
	tup.Form(tuple+"%d",ii);
        in_file.Form( in_dir + "/SingleMuon/SingleMu_2017B/180720_145352/0000/" + "tuple_%d.root", ii);
        in_files.push_back(in_file);
	samples["data_B_2017"+tup] = new Sample(in_files, "data_B_2017" + tup, "data");
	in_files.clear();
    	in_files.shrink_to_fit();

    }
//    samples["data_B_2017"] = new Sample(in_files, "data_B_2017" , "data");
//    samples["data_B_2017"]->xsec = 1;
    std::cout << "data_B_2017 file added" << in_files.size() << std::endl;
    in_files.clear();
    in_files.shrink_to_fit();

    for(int ii=1 ; ii < 209; ii++)
    {
	if(ii > 40) continue;
	tup.Form(tuple+"%d",ii);
        in_file.Form( in_dir + "/SingleMuon/SingleMu_2017C/180720_145447/0000/" + "tuple_%d.root", ii);
        in_files.push_back(in_file);
	samples["data_C_2017"+tup] = new Sample(in_files, "data_C_2017" + tup, "data");
        in_files.clear();
        in_files.shrink_to_fit();
    }
//    samples["data_C_2017"] = new Sample(in_files, "data_C_2017" , "data");
//    samples["data_C_2017"]->xsec = 1;
    std::cout << "data_C_2017 file added" << in_files.size() << std::endl;
    in_files.clear();
    in_files.shrink_to_fit();
 
    for(int ii=1 ; ii < 110; ii++)
    {
	if(ii > 20) continue;
	tup.Form(tuple+"%d",ii);
        in_file.Form( in_dir + "/SingleMuon/SingleMu_2017D/180720_145545/0000/" + "tuple_%d.root", ii);
        in_files.push_back(in_file);
	samples["data_D_2017"+tup] = new Sample(in_files, "data_D_2017" + tup, "data");
        in_files.clear();
        in_files.shrink_to_fit();
    }
//    samples["data_D_2017"] = new Sample(in_files, "data_D_2017" , "data");
//    samples["data_D_2017"]->xsec = 1;
    std::cout << "data_D_2017 file added" << in_files.size() << std::endl;
    in_files.clear();
    in_files.shrink_to_fit();

    for(int ii=1 ; ii < 173; ii++)
    {
	if(ii > 40) continue;
	if( ii == 92) continue;
	tup.Form(tuple+"%d",ii);
        in_file.Form( in_dir + "/SingleMuon/SingleMu_2017E/180720_145706/0000/" + "tuple_%d.root", ii);
        in_files.push_back(in_file);
	samples["data_E_2017"+tup] = new Sample(in_files, "data_E_2017" + tup, "data");
        in_files.clear();
        in_files.shrink_to_fit();
    }
//    samples["data_E_2017"] = new Sample(in_files, "data_E_2017" , "data");
//    samples["data_E_2017"]->xsec = 1;
    std::cout << "data_E_2017 file added" << in_files.size() << std::endl;
    in_files.clear();
    in_files.shrink_to_fit();

    for(int ii=1 ; ii < 235; ii++)
    {
	if(ii > 40) continue;
	if (ii == 20 or ii == 28 or ii == 42 or ii == 81 or ii == 140 or ii == 170 or ii == 176 or ii == 185 or ii == 188) continue;
	tup.Form(tuple+"%d",ii);
        in_file.Form( in_dir + "/SingleMuon/SingleMu_2017F/180720_145810/0000/" + "tuple_%d.root", ii);
        in_files.push_back(in_file);
	samples["data_F_2017"+tup] = new Sample(in_files, "data_F_2017" + tup, "data");
        in_files.clear();
        in_files.shrink_to_fit();
    }
//    samples["data_F_2017"] = new Sample(in_files, "data_F_2017" , "data");
//    samples["data_F_2017"]->xsec = 1;
    std::cout << "data_F_2017 file added" << in_files.size() << std::endl;
    in_files.clear();
    in_files.shrink_to_fit();

  }



}
  
  // ================================================================
  // H2Mu_gg ---------------------------------------------------------
  // ================================================================

