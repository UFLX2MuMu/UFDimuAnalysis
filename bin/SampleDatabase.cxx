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
    in_dir = "root://eoscms.cern.ch//store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13";
  else if (location == "CERN_hiM")
    in_dir = "root://eoscms.cern.ch//store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM";
  else
    std::cout << "\n\nInput location is " << location << ", not UF, CERN, or CERN_hiM.  NOT AN OPTION!!!\n\n" << std::endl;
  
  TString in_dir_hiM = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/hiM_simplified/";

  // ================================================================
  // Data -----------------------------------------------------------
  // ================================================================
  
  /////// INDIVIDUAL ERAS /////////////////////////////////////////////////////
  std::vector< std::tuple< TString, float, int > > eras;
  // Era tuple has name, luminosity, and number of files
  // Very rough lumi splitting between eras; needs to be updated - AWB 01.02.17
  eras.push_back( std::make_tuple("B",   5800, 1) );
  eras.push_back( std::make_tuple("C",   2600, 1) );
  eras.push_back( std::make_tuple("D",   4300, 1) );
  eras.push_back( std::make_tuple("E",   4100, 1) );
  eras.push_back( std::make_tuple("F_1", 1600, 1) );
  eras.push_back( std::make_tuple("F_2", 1600, 1) );
  eras.push_back( std::make_tuple("G",   7800, 1) );
  eras.push_back( std::make_tuple("H",   9014, 2) );
  
  std::vector<TString> in_files_all_data;
  for (auto era: eras) {
    if (select != "DATA" && !select.Contains("ALL") && select != "Run"+std::get<0>(era))
      continue;
    std::cout << "Adding files for Run" << std::get<0>(era) << " ...." << std::endl;
    
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"data/SingleMuon_SingleMu_2016"+std::get<0>(era)+".root") );
    }
    else if (location == "UF_DoubleMu") {
      in_files.push_back( TString(data_dir+"DoubleMu_2016"+std::get<0>(era)+".root") );
    } else {
      for (int i = 0; i < std::get<2>(era); i++) {
	TString jDat = "";
	if ( std::get<2>(era) > 1) {
	  if (i == 0) jDat = "_1";
	  if (i == 1) jDat = "_2";
	}

	// if (std::get<0>(era) != "F_1" || location != "CERN_hiM") {
	in_file.Form( "%s/SingleMuon/SingleMu_2016%s%s/NTuple_0.root", in_dir.Data(), std::get<0>(era).Data(), jDat.Data() );
	in_files.push_back(in_file);
	// } else {  // Load input files for F_1 manually because of buggy event in CERN_hiM sample - AWB 28.03.17
	// for (int j = 1; j <= 37; j++) { 
	//   if (j ==  9) continue;  // Leave out tuple_9.root  with buggy event 1026584286 - AWB 28.03.17
	//   if (j == 34) continue;  // Leave out tuple_34.root with buggy event 1160552942 - AWB 28.03.17
	//   if (j == 35) continue;  // Leave out tuple_35.root with buggy event 1310428464 - AWB 28.03.17
	//   if (j == 36) continue;  // Leave out tuple_36.root with buggy event   90777759 - AWB 28.03.17
	//   if (j == 37) continue;  // Leave out tuple_37.root with buggy event 2179193035 - AWB 28.03.17
	//   in_file.Form( "%s/SingleMuon/SingleMu_2016F_1/170315_104802/0000/tuple_%d.root", in_dir.Data(), j );
	//   in_files.push_back(in_file);
	// }
	// }
	  
      } // End loop: for (int i = 0; i < std::get<2>(era); i++)
    } 
    
    Sample* data_sample = (!info_only)?new Sample(in_files, "Run"+std::get<0>(era), "data"):new Sample("Run"+std::get<0>(era), "data");
    data_sample->lumi = std::get<1>(era);
    data_sample->xsec = 9999;
    samples["Run"+std::get<0>(era)] = data_sample;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
    in_files_all_data.insert( in_files_all_data.end(), in_files.begin(), in_files.end() );
  }

  // if (select == "DATA" || select.Contains("ALL") || select == "RunAll") {
  //   std::cout << "Adding files for RunAll ...." << std::endl;
  //   Sample* data_sample_all = new Sample(in_files_all_data, "RunAll", "data");
  //   data_sample_all->lumi = 36814;
  //   data_sample_all->xsec = 9999;
  //   samples["RunAll"] = data_sample_all;
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files_all_data.size() << " files added." << std::endl;
  // }
  

  
  // ================================================================
  // H2Mu_gg ---------------------------------------------------------
  // ================================================================
  
  if (select.Contains("SIGNALX") || (select.Contains("ALL") && !select.Contains("1")) || select == "MC" || select == "SIGNAL" || select == "H2Mu_gg") {
    std::cout << "Adding files for H2Mu_gg ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_gg.root") );
    } else {
      in_file.Form( "%s/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_gg/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_gg"] = (!info_only)?new Sample(in_files, "H2Mu_gg", "signal"):new Sample("H2Mu_gg", "signal");
    samples["H2Mu_gg"]->xsec = 0.009618; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select.Contains("SIGNALX") || select.Contains("ALL120") || select == "MC120" || select == "SIGNAL120" || select == "H2Mu_gg_120") {
    std::cout << "Adding files for H2Mu_gg_120 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir_hiM+"signal/GluGlu_HToMuMu_M120_13TeV_powheg_pythia8_H2Mu_gg_120.root") );
    } 
    else {
      in_file.Form( "%s/GluGlu_HToMuMu_M120_13TeV_powheg_pythia8/H2Mu_gg_120/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_gg_120"] = (!info_only)?new Sample(in_files, "H2Mu_gg_120", "signal"):new Sample("H2Mu_gg_120", "signal");
    samples["H2Mu_gg_120"]->xsec = 0.009618; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select.Contains("SIGNALX") || select.Contains("ALL130") || select == "MC130" || select == "SIGNAL130" || select == "H2Mu_gg_130") {
    std::cout << "Adding files for H2Mu_gg_130 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir_hiM+"signal/GluGlu_HToMuMu_M130_13TeV_powheg_pythia8_H2Mu_gg_130.root") );
    } 
    else {
      in_file.Form( "%s/GluGlu_HToMuMu_M130_13TeV_powheg_pythia8/H2Mu_gg_130/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_gg_130"] = (!info_only)?new Sample(in_files, "H2Mu_gg_130", "signal"):new Sample("H2Mu_gg_130", "signal");
    samples["H2Mu_gg_130"]->xsec = 0.009618; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  // ================================================================
  // H2Mu_VBF ---------------------------------------------------------
  // ================================================================
  
  if (select.Contains("SIGNALX") || (select.Contains("ALL") && !select.Contains("1")) || select == "MC" || select == "SIGNAL" || select == "H2Mu_VBF") {
    std::cout << "Adding files for H2Mu_VBF ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/VBF_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_VBF.root") );
    } else {
      in_file.Form( "%s/VBF_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_VBF/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_VBF"] = (!info_only)?new Sample(in_files, "H2Mu_VBF", "signal"):new Sample("H2Mu_VBF", "signal");
    samples["H2Mu_VBF"]->xsec = 0.0008208; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 

  // if (select.Contains("SIGNALX") || select.Contains("ALL120") || select == "MC120" || select == "SIGNAL120" || select == "H2Mu_VBF_120") {
  //   std::cout << "Adding files for H2Mu_VBF_120 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/VBF_HToMuMu_M120_13TeV_powheg_pythia8_H2Mu_VBF_120.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/VBF_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_VBF_120/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_VBF_120"] = (!info_only)?new Sample(in_files, "H2Mu_VBF_120", "signal"):new Sample("H2Mu_VBF_120", "signal");
  //   samples["H2Mu_VBF_120"]->xsec = 0.0008208; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("SIGNALX") || select.Contains("ALL130") || select == "MC130" || select == "SIGNAL130" || select == "H2Mu_VBF_130") {
  //   std::cout << "Adding files for H2Mu_VBF_130 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/VBF_HToMuMu_M130_13TeV_powheg_pythia8_H2Mu_VBF_130.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/VBF_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_VBF_130/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_VBF_130"] = (!info_only)?new Sample(in_files, "H2Mu_VBF_130", "signal"):new Sample("H2Mu_VBF_130", "signal");
  //   samples["H2Mu_VBF_130"]->xsec = 0.0008208; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

 
  // ================================================================
  // H2Mu_VH ---------------------------------------------------------
  // ================================================================
  
  if (select.Contains("SIGNALX") || (select.Contains("ALL") && !select.Contains("1")) || select == "MC" 
      || select == "SIGNAL" || select == "H2Mu_VH" || select == "H2Mu_ZH") {

    std::cout << "Adding files for H2Mu_ZH ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/ZH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_ZH.root") );
    } else {
      in_file.Form( "%s/ZH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_ZH/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_ZH"] = (!info_only)?new Sample(in_files, "H2Mu_ZH", "signal"):new Sample("H2Mu_ZH", "signal");
    samples["H2Mu_ZH"]->xsec = 0.0002136;     // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 

  // if (select.Contains("SIGNALX") || select.Contains("ALL120") || select == "MC120" || select == "SIGNAL120" || select == "H2Mu_ZH_120") {
  //   std::cout << "Adding files for H2Mu_ZH_120 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/ZH_HToMuMu_M120_13TeV_powheg_pythia8_H2Mu_ZH_120.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/ZH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_ZH_120/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_ZH_120"] = (!info_only)?new Sample(in_files, "H2Mu_ZH_120", "signal"):new Sample("H2Mu_ZH_120", "signal");
  //   samples["H2Mu_ZH_120"]->xsec = 0.0002136; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("SIGNALX") || select.Contains("ALL130") || select == "MC130" || select == "SIGNAL130" || select == "H2Mu_ZH_130") {
  //   std::cout << "Adding files for H2Mu_ZH_130 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/ZH_HToMuMu_M130_13TeV_powheg_pythia8_H2Mu_ZH_130.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/ZH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_ZH_130/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_ZH_130"] = (!info_only)?new Sample(in_files, "H2Mu_ZH_130", "signal"):new Sample("H2Mu_ZH_130", "signal");
  //   samples["H2Mu_ZH_130"]->xsec = 0.0002136; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("SIGNALX") || (select.Contains("ALL") && !select.Contains("1"))  || select == "MC" || 
      select == "SIGNAL" || select == "H2Mu_VH" || select == "H2Mu_WH"  || select == "H2Mu_WH_pos") {

    std::cout << "Adding files for H2Mu_WH_pos ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_WH_pos.root") );
    } else {
      in_file.Form( "%s/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_WH_pos"] = (!info_only)?new Sample(in_files, "H2Mu_WH_pos", "signal"):new Sample("H2Mu_WH_pos", "signal");
    samples["H2Mu_WH_pos"]->xsec = 0.0001858; // 0.851*0.0002176;     // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 
 
  // if (select.Contains("SIGNALX") || select.Contains("ALL120") || select == "MC120" || select == "SIGNAL120" || select == "H2Mu_WH_pos_120") {
  //   std::cout << "Adding files for H2Mu_WH_pos_120 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/WPlusH_HToMuMu_M120_13TeV_powheg_pythia8_H2Mu_WH_pos_120.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/WH_pos_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos_120/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_WH_pos_120"] = (!info_only)?new Sample(in_files, "H2Mu_WH_pos_120", "signal"):new Sample("H2Mu_WH_pos_120", "signal");
  //   samples["H2Mu_WH_pos_120"]->xsec = 0.0001858; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("SIGNALX") || select.Contains("ALL130") || select == "MC130" || select == "SIGNAL130" || select == "H2Mu_WH_pos_130") {
  //   std::cout << "Adding files for H2Mu_WH_pos_130 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/WPlusH_HToMuMu_M130_13TeV_powheg_pythia8_H2Mu_WH_pos_130.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/WH_pos_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos_130/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_WH_pos_130"] = (!info_only)?new Sample(in_files, "H2Mu_WH_pos_130", "signal"):new Sample("H2Mu_WH_pos_130", "signal");
  //   samples["H2Mu_WH_pos_130"]->xsec = 0.0001858; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("SIGNALX") || (select.Contains("ALL") && !select.Contains("1")) || select == "MC" || 
      select == "SIGNAL" || select == "H2Mu_VH" || select == "H2Mu_WH" || select == "H2Mu_WH_neg") {

    std::cout << "Adding files for H2Mu_WH_neg ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_WH_neg.root") );
    } else {
      in_file.Form( "%s/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_neg/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_WH_neg"] = (!info_only)?new Sample(in_files, "H2Mu_WH_neg", "signal"):new Sample("H2Mu_WH_neg", "signal");
    samples["H2Mu_WH_neg"]->xsec = 0.0001164; //0.5331*0.0002176;     // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 

  // if (select.Contains("SIGNALX") || select.Contains("ALL120") || select == "MC120" || select == "SIGNAL120" || select == "H2Mu_WH_neg_120") {
  //   std::cout << "Adding files for H2Mu_WH_neg_120 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/WMinusH_HToMuMu_M120_13TeV_powheg_pythia8_H2Mu_WH_neg_120.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/WH_neg_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_neg_120/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_WH_neg_120"] = (!info_only)?new Sample(in_files, "H2Mu_WH_neg_120", "signal"):new Sample("H2Mu_WH_neg_120", "signal");
  //   samples["H2Mu_WH_neg_120"]->xsec = 0.0001164; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("SIGNALX") || select.Contains("ALL130") || select == "MC130" || select == "SIGNAL130" || select == "H2Mu_WH_neg_130") {
  //   std::cout << "Adding files for H2Mu_WH_neg_130 ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir_hiM+"signal/WMinusH_HToMuMu_M130_13TeV_powheg_pythia8_H2Mu_WH_neg_130.root") );
  //   } 
  //   else {
  //     in_file.Form( "%s/WH_neg_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_neg_130/NTuple_0.root", in_dir.Data() );
  //     in_files.push_back(in_file);
  //   }
  //   samples["H2Mu_WH_neg_130"] = (!info_only)?new Sample(in_files, "H2Mu_WH_neg_130", "signal"):new Sample("H2Mu_WH_neg_130", "signal");
  //   samples["H2Mu_WH_neg_130"]->xsec = 0.0001164; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // ================================================================
  // H2Mu_ttH ---------------------------------------------------------
  // ================================================================
  
  if (select == "MC" || select == "SIGNAL" || select == "H2Mu_ttH") {
    std::cout << "Adding files for H2Mu_ttH ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_ttH.root") );
    } else {
      in_file.Form( "%s/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/ttH/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["H2Mu_ttH"] = new Sample(in_files, "H2Mu_ttH", "signal");
    samples["H2Mu_ttH"]->xsec = 0.2151;     // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 
 
 
  // ================================================================
  // DYJetsToLL -----------------------------------------------------
  // ================================================================
  
  float DY_xsec_m50 = 5765.4; // pb // old value = 6025.2
  float DY_m100to200_factor = 1.235; 
  
  if ((select.Contains("ALL") && select.Contains("AMC") && !select.Contains("AMC-J")) || select == "MC" || select == "BACKGROUND" 
       || select == "ZJets" || select == "ZJets_AMC") {
    std::cout << "Adding files for ZJets_AMC ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_AMC.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_AMC"] = new Sample(in_files, "ZJets_AMC", "background");
    samples["ZJets_AMC"]->xsec = DY_xsec_m50; 
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // if ((select.Contains("ALL") && select.Contains("AMC-J")) || select == "MC" || select == "BACKGROUND" 
  //      || select == "ZJets" || select == "ZJets_AMC" || select == "ZJets_AMC-J") {
  //   std::cout << "Adding files for ZJets_AMC_0j ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"dy/DYToLL_0J_13TeV-amcatnloFXFX-pythia8_ZJets_AMC_0j_A.root") );
  //   } else {
  //   //  in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/NTuple_0.root", in_dir.Data() );
  //   //  in_files.push_back(in_file);
  //   }
  //   samples["ZJets_AMC_0j"] = new Sample(in_files, "ZJets_AMC_0j", "background");
  //   samples["ZJets_AMC_0j"]->xsec = 4754*0.96; 
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // if ((select.Contains("ALL") && select.Contains("AMC-J")) || select == "MC" || select == "BACKGROUND" || 
  //      select == "ZJets" || select == "ZJets_AMC" || select == "ZJets_AMC-J") {
  //   std::cout << "Adding files for ZJets_AMC_1j ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"dy/DYToLL_1J_13TeV-amcatnloFXFX-pythia8_ZJets_AMC_1j_A.root") );
  //   } else {
  //   //  in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/NTuple_0.root", in_dir.Data() );
  //   //  in_files.push_back(in_file);
  //   }
  //   samples["ZJets_AMC_1j"] = new Sample(in_files, "ZJets_AMC_1j", "background");
  //   samples["ZJets_AMC_1j"]->xsec = 888.9*0.86*0.985*0.995; //888.9*0.95; 
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // if ((select.Contains("ALL") && select.Contains("AMC-J")) || select == "MC" || select == "BACKGROUND" || 
  //      select == "ZJets" || select == "ZJets_AMC" || select == "ZJets_AMC-J") {
  //   std::cout << "Adding files for ZJets_AMC_2j ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"dy/DYToLL_2J_13TeV-amcatnloFXFX-pythia8_ZJets_AMC_2j.root") );
  //   } else {
  //   //  in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/NTuple_0.root", in_dir.Data() );
  //   //  in_files.push_back(in_file);
  //   }
  //   samples["ZJets_AMC_2j"] = new Sample(in_files, "ZJets_AMC_2j", "background");
  //   samples["ZJets_AMC_2j"]->xsec = 348.8*0.88*0.975*0.992; 
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL") && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG" || select == "ZJets_MG_incl") {
    std::cout << "Adding files for ZJets_MG ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG"] = new Sample(in_files, "ZJets_MG", "background");
    samples["ZJets_MG"]->xsec = DY_xsec_m50;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_70_100 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_70_100.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_70_100/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_70_100"] = new Sample(in_files, "ZJets_MG_HT_70_100", "background");
    samples["ZJets_MG_HT_70_100"]->xsec = 0.98*178.952; // old value = 206.184;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_100_200 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_100_200.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_100_200_A/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);

      // if (location != "CERN_hiM") {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_100_200_B/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      // } else {  // Load input files for HT_100_200_B manually because of buggy event in CERN_hiM sample - AWB 28.03.17
      // 	for (int i = 1; i <= 27; i++) { 
      // 	  if (i ==  9) continue;  // Leave out tuple_9.root  with buggy event 81999086 - AWB 28.03.17
      // 	  if (i == 27) continue;  // Leave out tuple_27.root with buggy event 22590445 - AWB 28.03.17
      // 	  if (i == 1 || i == 4 || i == 11 || i == 14) continue;  // Failed jobs, no root files
      // 	  in_file.Form( "%s/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_100_200_B/170315_223309/0000/tuple_%d.root", in_dir.Data(), i );
      // 	  std::cout << "  * Adding file " << in_file << std::endl;
      // 	  in_files.push_back(in_file);
      // 	}
      // }
    }
    samples["ZJets_MG_HT_100_200"] = new Sample(in_files, "ZJets_MG_HT_100_200", "background");
    samples["ZJets_MG_HT_100_200"]->xsec = 0.96*181.302;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_200_400 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_200_400.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_200_400_A/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      in_file.Form( "%s/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_200_400_B/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_200_400"] = new Sample(in_files, "ZJets_MG_HT_200_400", "background");
    samples["ZJets_MG_HT_200_400"]->xsec = 0.96*50.4177;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_400_600 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_400_600.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_400_600_A/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      in_file.Form( "%s/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_400_600_B/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_400_600"] = new Sample(in_files, "ZJets_MG_HT_400_600", "background");
    samples["ZJets_MG_HT_400_600"]->xsec = 0.96*6.98394;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_600_800 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_600_800.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_600_800/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_600_800"] = new Sample(in_files, "ZJets_MG_HT_600_800", "background");
    samples["ZJets_MG_HT_600_800"]->xsec = 0.96*1.68141;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_800_1200 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_800_1200.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_800_1200/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_800_1200"] = new Sample(in_files, "ZJets_MG_HT_800_1200", "background");
    samples["ZJets_MG_HT_800_1200"]->xsec = 0.96*0.775392;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL")  && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_1200_2500 ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_1200_2500.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_1200_2500/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_1200_2500"] = new Sample(in_files, "ZJets_MG_HT_1200_2500", "background");
    samples["ZJets_MG_HT_1200_2500"]->xsec = 0.96*0.186222;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if ((select.Contains("ALL") && select.Contains("MG")) || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "Adding files for ZJets_MG_HT_2500_inf ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG_HT_2500_inf.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG_HT_2500_inf/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_MG_HT_2500_inf"] = new Sample(in_files, "ZJets_MG_HT_2500_inf", "background");
    samples["ZJets_MG_HT_2500_inf"]->xsec = 0.96*0.004385;
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_hiM") {
    std::cout << "Adding files for ZJets_hiM ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_hiM.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_hiM/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_hiM"] = new Sample(in_files, "ZJets_hiM", "background");
    samples["ZJets_hiM"]->xsec = DY_xsec_m50 * DY_m100to200_factor; // 7117
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_hiM_SpringPU") {
    std::cout << "Adding files for ZJets_hiM_SpringPU ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_hiM_SpringPU.root") );
    } else {
      in_file.Form( "%s/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_hiM_SpringPU/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ZJets_hiM_SpringPU"] = new Sample(in_files, "ZJets_hiM_SpringPU", "background");
    samples["ZJets_hiM_SpringPU"]->xsec = DY_xsec_m50 * DY_m100to200_factor; // 7117
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }


  // ================================================================
  // TTJets ---------------------------------------------------------
  // ================================================================
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "ttbar" || select == "tt_ll_AMC") {
    std::cout << "Adding files for tt_ll_AMC ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"ttjets/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_tt_ll_AMC.root") );
    } else {
      in_file.Form( "%s/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/tt_ll_AMC/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["tt_ll_AMC"] = new Sample(in_files, "tt_ll_AMC", "background");
    samples["tt_ll_AMC"]->xsec = 85.656*0.9; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select == "MC" || select == "BACKGROUND" || select == "ttbar" || select == "tt_ll_MG") {
    std::cout << "Adding files for tt_ll_MG ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"ttjets/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_tt_ll_MG_2.root") );
    } else {
      in_file.Form( "%s/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/tt_ll_MG_1/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      in_file.Form( "%s/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/tt_ll_MG_2/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["tt_ll_MG"] = new Sample(in_files, "tt_ll_MG", "background");
    samples["tt_ll_MG"]->xsec = 85.656; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }


  // ================================================================
  // SingleTop ------------------------------------------------------
  // ================================================================
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tW" || select == "tW_pos") {
    std::cout << "Adding files for tW_pos ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"singletop/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_tW_pos.root") );
    } else {
      in_file.Form( "%s/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/tW_pos_1/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      in_file.Form( "%s/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/tW_pos_2/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["tW_pos"] = new Sample(in_files, "tW_pos", "background");
    samples["tW_pos"]->xsec = 35.85; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tW" || select == "tW_neg") {
    std::cout << "Adding files for tW_neg ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"singletop/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_tW_neg.root") );
    } else {

      in_file.Form( "%s/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/tW_neg_1/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);

      // if (location != "CERN_hiM") {
      in_file.Form( "%s/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/tW_neg_2/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      // } else { // Load input files manually because of buggy event in CERN_hiM sample - AWB 28.03.17
      // 	for (int i = 1; i <= 2; i++) { 
      // 	  if (i == 1) continue;  // Leave out tuple_1.root with buggy event 4906108 - AWB 28.03.17
      // 	  if (i == 2) continue;  // Leave out tuple_2.root with buggy event 3586447 - AWB 28.03.17
      // 	  in_file.Form( "%s/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/tW_neg_2/170315_105904/0000/tuple_%d.root", in_dir.Data(), i);
      // 	  in_files.push_back(in_file);
      // 	}
      // }

    }
    samples["tW_neg"] = new Sample(in_files, "tW_neg", "background");
    samples["tW_neg"]->xsec = 35.85; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tZq") {
    std::cout << "Adding files for tZq ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"singletop/tZq_ll_4f_13TeV-amcatnlo-pythia8_tZq.root") );
    } else {
      in_file.Form( "%s/tZq_ll_4f_13TeV-amcatnlo-pythia8/tZq/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["tZq"] = new Sample(in_files, "tZq", "background");
    samples["tZq"]->xsec = 0.0758; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //if (select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tZW") {
  //  std::cout << "Adding files for tZW ..." << std::endl;
  //  std::vector<TString> in_files;
  //  TString in_file;
  //  if (location.Contains("UF")) {
  //    in_files.push_back( TString(in_dir+"singletop/ST_tWll_5f_LO_13TeV-MadGraph-pythia8_tZW.root") );
  //  } else {
  //    in_file.Form( "%s/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/tZW/NTuple_0.root", in_dir.Data() );
  //    in_files.push_back(in_file);
  //  }
  //  samples["tZW"] = new Sample(in_files, "tZW", "background");
  //  samples["tZW"]->xsec = -999; // pb
  //   if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  //}

    // ================================================================
    // TTX ---------------------------------------------------------
    // ================================================================
    
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "ttX" || select == "ttV" || select == "ttW") {
    std::cout << "Adding files for ttW ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"ttv/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_ttW.root") );
    } else {
      in_file.Form( "%s/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/ttW_1/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
      in_file.Form( "%s/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/ttW_2/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ttW"] = new Sample(in_files, "ttW", "background");
    samples["ttW"]->xsec = 0.2043; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "ttX" || select == "ttV" || select == "ttZ") {
    std::cout << "Adding files for ttZ ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"ttv/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_ttZ.root") );
    } else {
      in_file.Form( "%s/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/ttZ/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ttZ"] = new Sample(in_files, "ttZ", "background");
    samples["ttZ"]->xsec = 0.2529; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select == "MC" || select == "BACKGROUND" || select == "ttX" || select == "ttH" ) {
    std::cout << "Adding files for H2Mu_ttH ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"signal/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_ttH.root") );
    } else {
      in_file.Form( "%s/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/ttH/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["ttH"] = new Sample(in_files, "ttH", "background");
    samples["ttH"]->xsec = 0.2151;     // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }
 
  // ================================================================
  // Diboson ------------------------------------------------------
  // ================================================================
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "WW") {
    std::cout << "Adding files for WW ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"diboson/WWTo2L2Nu_13TeV-powheg_WW.root") );
    } else {
      in_file.Form( "%s/WWTo2L2Nu_13TeV-powheg/WW/NTuple_0.root", in_dir.Data() );
      in_files.push_back(in_file);
    }
    samples["WW"] = new Sample(in_files, "WW", "background");
    samples["WW"]->xsec = 12.46; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "WZ" || select == "WZ_2l") {
    std::cout << "Adding files for WZ_2l ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"diboson/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_WZ_2l.root") );
    } else {
     in_file.Form( "%s/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/WZ_2l/NTuple_0.root", in_dir.Data() );
     in_files.push_back(in_file);
    }
    samples["WZ_2l"] = new Sample(in_files, "WZ_2l", "background");
    samples["WZ_2l"]->xsec = 4.409; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
   }

  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "WZ" || select == "WZ_3l") {
    std::cout << "Adding files for WZ_3l ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"diboson/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_WZ_3l_AMC.root") );
    } else {
     in_file.Form( "%s/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WZ_3l_AMC/NTuple_0.root", in_dir.Data() );
     in_files.push_back(in_file);
    }
    samples["WZ_3l"] = new Sample(in_files, "WZ_3l", "background");
    samples["WZ_3l"]->xsec = 2.113; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "ZZ" || select == "ZZ_2l_2v") {
    std::cout << "Adding files for ZZ_2l_2v ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"diboson/ZZTo2L2Nu_13TeV_powheg_pythia8_ZZ_2l_2v.root") );
    } else {
     in_file.Form( "%s/ZZTo2L2Nu_13TeV_powheg_pythia8/ZZ_2l_2v/NTuple_0.root", in_dir.Data() );
     in_files.push_back(in_file);
    }
    samples["ZZ_2l_2v"] = new Sample(in_files, "ZZ_2l_2v", "background");
    samples["ZZ_2l_2v"]->xsec = 0.564; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "ZZ" || select == "ZZ_2l_2q") {
    std::cout << "Adding files for ZZ_2l_2q ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"diboson/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_ZZ_2l_2q.root") );
    } else {
     in_file.Form( "%s/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/ZZ_2l_2q/NTuple_0.root", in_dir.Data() );
     in_files.push_back(in_file);
    }
    samples["ZZ_2l_2q"] = new Sample(in_files, "ZZ_2l_2q", "background");
    samples["ZZ_2l_2q"]->xsec = 3.22; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "ZZ" || select == "ZZ_4l") {
    std::cout << "Adding files for ZZ_4l ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location.Contains("UF")) {
      in_files.push_back( TString(in_dir+"diboson/ZZTo4L_13TeV-amcatnloFXFX-pythia8_ZZ_4l_AMC.root") );
    } else {
     in_file.Form( "%s/ZZTo4L_13TeV-amcatnloFXFX-pythia8/ZZ_4l_AMC/NTuple_0.root", in_dir.Data() );
     in_files.push_back(in_file);
    }
    samples["ZZ_4l"] = new Sample(in_files, "ZZ_4l", "background");
    samples["ZZ_4l"]->xsec = 1.212; // pb
     if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  // ================================================================
  // Triboson ------------------------------------------------------
  // ================================================================
  
  // if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "WWW") {
  //   std::cout << "Adding files for WWW ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"triboson/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WWW.root") );
  //   } else {
  //   in_file.Form( "%s/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/WWW/NTuple_0.root", in_dir.Data() );
  //   in_files.push_back(in_file);
  //   }
  //   samples["WWW"] = new Sample(in_files, "WWW", "background");
  //   samples["WWW"]->xsec = -999; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "WWZ") {
  //   std::cout << "Adding files for WWZ ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"triboson/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WWZ.root") );
  //   } else {
  //   in_file.Form( "%s/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/WWZ/NTuple_0.root", in_dir.Data() );
  //   in_files.push_back(in_file);
  //     // }
  //   }
  //   samples["WWZ"] = new Sample(in_files, "WWZ", "background");
  //   samples["WWZ"]->xsec = -999; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "WZZ") {
  //   std::cout << "Adding files for WZZ ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"triboson/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WZZ.root") );
  //   } else {
  //   in_file.Form( "%s/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/WZZ/NTuple_0.root", in_dir.Data() );
  //   in_files.push_back(in_file);
  //   }
  //   samples["WZZ"] = new Sample(in_files, "WZZ", "background");
  //   samples["WZZ"]->xsec = -999; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select.Contains("ALL") || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "ZZZ") {
  //   std::cout << "Adding files for ZZZ ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location.Contains("UF")) {
  //     in_files.push_back( TString(in_dir+"triboson/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_ZZZ.root") );
  //   } else {
  //   in_file.Form( "%s/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/ZZZ/NTuple_0.root", in_dir.Data() );
  //   in_files.push_back(in_file);
  //   }
  //   samples["ZZZ"] = new Sample(in_files, "ZZZ", "background");
  //   samples["ZZZ"]->xsec = -999; // pb
  //    if(!location.Contains("UF")) std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

}
