///////////////////////////////////////////////////////////////////////////
//                       CategorySelection.cxx                           //
//=======================================================================//
//                                                                       //
// Categorizer objects to categorize each event. Define the cuts for the //
// categories and the evaluate function to determine the category        //
// structure. We keep track of the different categories in the           //
// in the categorizer via categoryMap<TString, Category>. Each category  //
// tracks its historams with histoMap<TString, TH1D*> and some TLists.   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "CategorySelection.h"
#include "ParticleTools.h"
#include "EventTools.h"
#include <sstream>


void CategoryNode::theMiracleOfChildBirth()
{
    left = new CategoryNode(this, 0, 0, this->key+"_left", -999, "", -999, -999) ; 
    right = new CategoryNode(this, 0, 0, this->key+"_right", -999, "", -999, -999) ;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________XMLCategorizer_________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void XMLCategorizer::initCategoryMap()
{
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void XMLCategorizer::evaluate(VarSet& vars)
{
    evaluateRecursive(vars, rootNode); 
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void XMLCategorizer::evaluateRecursive(VarSet& vars, CategoryNode* cnode)
{
    // if it was filtered into this node then it's in this category
    categoryMap[cnode->key].inCategory = true;

    // return if terminal node
    if(cnode->left == 0 || cnode->right == 0)
    { 
        return;
    }

    // if not terminal node, filter to correct daughter node
    const char* splitVarName = cnode->splitVarName.Data();
    double splitValue   = cnode->splitVal;
    double xvalue = vars.getValue(splitVarName); 
    
    if(xvalue <= splitValue) 
        evaluateRecursive(vars, cnode->left);
    else 
        evaluateRecursive(vars, cnode->right);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

XMLCategorizer::XMLCategorizer()
{
    rootNode = new CategoryNode(0, 0, 0, "", -999, "", -999, -999);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

XMLCategorizer::XMLCategorizer(TString xmlfile)
{
    rootNode = new CategoryNode(0, 0, 0, "", -999, "", -999, -999);
    loadFromXML(xmlfile);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void XMLCategorizer::loadFromXML(TString filename)
{
    //std::cout << Form("/// Loading categories from %s... \n\n", filename.Data());
    // First create the engine.
    TXMLEngine* xml = new TXMLEngine;

    // Now try to parse xml file.
    XMLDocPointer_t xmldoc = xml->ParseFile(filename);
    if (xmldoc==0)
    {
        delete xml;
        return;  
    }

    // Get access to main node of the xml file.
    XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
   
    // Recursively connect nodes together.
    loadFromXMLRecursive(xml, mainnode, rootNode);
   
    // Release memory before exit
    xml->FreeDoc(xmldoc);

    delete xml;

    int i=0;
    for(auto& c: categoryMap)
    {
        if(c.second.isTerminal)
        {
            c.second.name = Form("c%d", i);
            i++;
        }
        else if(c.second.key == "root")
        {
            c.second.name = c.second.key;
        }
        else
        {
            c.second.name = c.second.key;
            c.second.hide = true;
        }
    }
      
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void XMLCategorizer::loadFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t xnode, CategoryNode* cnode)
{
   // Get the split information from xml.
    XMLAttrPointer_t attr = xml->GetFirstAttr(xnode);
    std::vector<std::string> splitInfo(4);
    for(unsigned int i=0; i<splitInfo.size(); i++)
    {
        splitInfo[i] = xml->GetAttrValue(attr); 
        attr = xml->GetNextAttr(attr);  
    }

    // Convert strings into numbers.
    std::stringstream converter;
    int splitVar;
    TString splitVarName;
    double splitVal;
    double significanceSquared;  

    converter << splitInfo[0];
    converter >> splitVar;
    converter.str("");
    converter.clear();

    converter << splitInfo[1];
    splitVarName = TString(converter.str().c_str());
    converter.str("");
    converter.clear();

    converter << splitInfo[2];
    converter >> splitVal;
    converter.str("");
    converter.clear();

    converter << splitInfo[3];
    converter >> significanceSquared;
    converter.str("");
    converter.clear();

    //std::cout << "svar: " << splitVar << ", svar_name: " << splitVarName << ", split_val: " << splitVal << ", sig2: " << significanceSquared << std::endl;
    cnode->splitVar = splitVar;
    cnode->splitVarName = splitVarName;
    cnode->splitVal = splitVal;
    cnode->significanceSquared = significanceSquared;

    // Get the xml daughters of the current xml node. 
    XMLNodePointer_t xleft = xml->GetChild(xnode);
    XMLNodePointer_t xright = xml->GetNext(xleft);

    // If there are no daughters we are done.
    if(xleft == 0 || xright == 0)
    {
        TString sig = Form("%5.4f", TMath::Sqrt(significanceSquared));
        sig.ReplaceAll(" ", "");
        sig.ReplaceAll(".", "p");
        sig.ReplaceAll("-", "n");
        cnode->key = "T_"+sig+"_"+cnode->key;

        //cnode->output();
        categoryMap[cnode->key] = Category(cnode->key);
        categoryMap[cnode->key].isTerminal = true;
        return;
    }
    else
    {
        TString scut = Form("%s_%7.3f", splitVarName.Data(), splitVal); 
        scut.ReplaceAll(" ", "");
        scut.ReplaceAll(".", "p");
        scut.ReplaceAll("-", "n");
        if(cnode->mother == 0) cnode->key = "root";
        //cnode->output();

        categoryMap[cnode->key] = Category(cnode->key);

        cnode->theMiracleOfChildBirth();
        CategoryNode* cleft = cnode->left; 
        CategoryNode* cright = cnode->right; 

        if(cnode->mother==0)
        {
            cleft->key = "lt_"+scut;
            cright->key = "gt_"+scut;
        }
        else
        {
            cleft->key = cnode->key+"_lt_"+scut;
            cright->key = cnode->key+"_gt_"+scut;
        }

        loadFromXMLRecursive(xml, xleft, cleft);
        loadFromXMLRecursive(xml, xright, cright);
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelectionBDT___________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelectionBDT::initCategoryMap()
{
// Initialize the categories

    // Category(name), Category(name, hide?), Category(name, hide?, isTerminal?)
    // defaults are hide = false, isTerminal = false

    categoryMap["cAll"] = Category("cAll", false, false);
    categoryMap["c0"] = Category("c0", false, true);
    categoryMap["c1"] = Category("c1", false, true);
    categoryMap["c2"] = Category("c2", false, true);
    categoryMap["c3"] = Category("c3", false, true);
    categoryMap["c4"] = Category("c4", false, true);
    categoryMap["c5"] = Category("c5", false, true);
    categoryMap["c6"] = Category("c6", false, true);
    categoryMap["c7"] = Category("c7", false, true);
    categoryMap["c8"] = Category("c8", false, true);
    categoryMap["c9"] = Category("c9", false, true);
    categoryMap["c10"] = Category("c10", false, true);
    categoryMap["c11"] = Category("c11", false, true);
    categoryMap["c12"] = Category("c12", false, true);
    categoryMap["c13"] = Category("c13", false, true);
    categoryMap["c14"] = Category("c14", false, true);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionBDT::CategorySelectionBDT()
{
    initCategoryMap();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionBDT::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    double bdt_score = vars.getValue("bdt_score");
    double max_eta = vars.getValue("dimu_max_abs_eta");

    // Inclusive set of events
    categoryMap["cAll"].inCategory = true;

    if( bdt_score < -0.400 ) 
        categoryMap["c0"].inCategory = true;

    else if( bdt_score >= -0.400 && bdt_score < 0.050 && max_eta >= 1.900 )
        categoryMap["c1"].inCategory = true;

    else if( bdt_score >= -0.400 && bdt_score < 0.050 && max_eta < 1.900  && max_eta >=0.9) 
        categoryMap["c2"].inCategory = true;

    else if( bdt_score >= -0.400 && bdt_score < 0.050 && max_eta < 0.9 ) 
        categoryMap["c3"].inCategory = true;

    else if( bdt_score >= 0.050 && bdt_score < 0.250 && max_eta >= 1.9) 
        categoryMap["c4"].inCategory = true;

    else if( bdt_score >= 0.050 && bdt_score < 0.250 && max_eta >= 0.900 && max_eta < 1.9) 
        categoryMap["c5"].inCategory = true;

    else if( bdt_score >= 0.050 && bdt_score < 0.250 && max_eta < 0.900 ) 
        categoryMap["c6"].inCategory = true;

    else if( bdt_score >= 0.250 && bdt_score < 0.400 && max_eta >= 1.900 )
        categoryMap["c7"].inCategory = true;

    else if( bdt_score >= 0.250 && bdt_score < 0.400 && max_eta < 1.900 && max_eta >= 0.900 ) 
        categoryMap["c8"].inCategory = true;

    else if( bdt_score >= 0.250 && bdt_score < 0.400 && max_eta < 0.900 ) 
        categoryMap["c9"].inCategory = true;

    else if( bdt_score < 0.650 && bdt_score >= 0.400 && max_eta >= 1.900 ) 
        categoryMap["c10"].inCategory = true;

    else if( bdt_score < 0.650 && bdt_score >= 0.400 && max_eta < 1.900 && max_eta >= 0.900 ) 
        categoryMap["c11"].inCategory = true;

    else if( bdt_score < 0.650 && bdt_score >= 0.400 && max_eta < 0.900 ) 
        categoryMap["c12"].inCategory = true;

    else if( bdt_score < 0.730 && bdt_score >= 0.650 ) 
        categoryMap["c13"].inCategory = true;

    else if( bdt_score >= 0.730 )
        categoryMap["c14"].inCategory = true;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelectionRun1__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelectionRun1::initCategoryMap()
{
// Initialize the categories

    // Category(name), Category(name, hide?), Category(name, hide?, isTerminal?)
    // defaults are hide = false, isTerminal = false

    categoryMap["c_ALL"] = Category("c_ALL");

    // intermediate categories to make things easier
    categoryMap["c_2_Jet"] = Category("c_2_Jet");
    categoryMap["c_01_Jet"] = Category("c_01_Jet");

    categoryMap["c_2_Jet_VBF_Tight"] = Category("c_2_Jet_VBF_Tight", false, true);
    categoryMap["c_2_Jet_VBF_Loose"] = Category("c_2_Jet_VBF_Loose", false, true);
    categoryMap["c_2_Jet_GGF_Tight"] = Category("c_2_Jet_GGF_Tight", false, true);
    categoryMap["c_01_Jet_Tight"] = Category("c_01_Jet_Tight");
    categoryMap["c_01_Jet_Loose"] = Category("c_01_Jet_Loose");

    // intermediate categories to make things easier, don't plot these, hence the hide=true
    categoryMap["c_BB"] = Category("c_BB", true);
    categoryMap["c_BO"] = Category("c_BO", true);
    categoryMap["c_BE"] = Category("c_BE", true);
    categoryMap["c_OO"] = Category("c_OO", true);
    categoryMap["c_OE"] = Category("c_OE", true);
    categoryMap["c_EE"] = Category("c_EE", true);

    categoryMap["c_01_Jet_Tight_BB"] = Category("c_01_Jet_Tight_BB", false, true);
    categoryMap["c_01_Jet_Tight_BO"] = Category("c_01_Jet_Tight_BO", false, true);
    categoryMap["c_01_Jet_Tight_BE"] = Category("c_01_Jet_Tight_BE", false, true);
    categoryMap["c_01_Jet_Tight_OO"] = Category("c_01_Jet_Tight_OO", false, true);
    categoryMap["c_01_Jet_Tight_OE"] = Category("c_01_Jet_Tight_OE", false, true);
    categoryMap["c_01_Jet_Tight_EE"] = Category("c_01_Jet_Tight_EE", false, true);

    categoryMap["c_01_Jet_Loose_BB"] = Category("c_01_Jet_Loose_BB", false, true);
    categoryMap["c_01_Jet_Loose_BO"] = Category("c_01_Jet_Loose_BO", false, true);
    categoryMap["c_01_Jet_Loose_BE"] = Category("c_01_Jet_Loose_BE", false, true);
    categoryMap["c_01_Jet_Loose_OO"] = Category("c_01_Jet_Loose_OO", false, true);
    categoryMap["c_01_Jet_Loose_OE"] = Category("c_01_Jet_Loose_OE", false, true);
    categoryMap["c_01_Jet_Loose_EE"] = Category("c_01_Jet_Loose_EE", false, true);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionRun1::CategorySelectionRun1()
{
// initialize the default cut values in the constructor

    initCategoryMap();

    // Preselection
    cLeadPtMin = 40;
    cSubleadPtMin = 30;
    cMETMax = 40;

    // VBF Tight
    cDijetMassMinVBFT = 650;
    cDijetDeltaEtaMinVBFT = 3.5;

    // VBF Loose

    // GGF Tight
    cDijetMassMinGGFT = 250;
    cDimuPtMinGGFT = 50;

    // 01Tight
    cDimuPtMin01T = 25;

    // 01Loose
    
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionRun1::CategorySelectionRun1(float leadPtMin, float subleadPtMin, float METMax, float dijetMassMinVBFT, float dijetDeltaEtaMinVBFT, 
                                             float dijetMassMinGGFT, float dimuPtMinGGFT, float dimuPtMin01T)
{
// Initialize custom cut values

    initCategoryMap();

    // Preselection
    cLeadPtMin = leadPtMin;
    cSubleadPtMin = subleadPtMin;
    cMETMax = METMax;

    // VBF Tight
    cDijetMassMinVBFT = dijetMassMinVBFT;
    cDijetDeltaEtaMinVBFT = dijetDeltaEtaMinVBFT;

    // VBF Loose

    // GGF Tight
    cDijetMassMinGGFT = dijetMassMinGGFT;
    cDimuPtMinGGFT = dimuPtMinGGFT;

    // 01Tight
    cDimuPtMin01T = dimuPtMin01T;

    // 01Loose
}


///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionRun1::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // Inclusive category, all events that passed the selection cuts
    categoryMap["c_ALL"].inCategory = true;

    // Geometric Categories
    // Barrel Barrel
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < 0.8) 
        categoryMap["c_BB"].inCategory = true;

    // Overlap Overlap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta)>=0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta)<1.6 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta)>=0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta)<1.6) 
        categoryMap["c_OO"].inCategory = true;

    // Endcap Endcap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= 1.6 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= 1.6) 
        categoryMap["c_EE"].inCategory = true;

    // Barrel Overlap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < 1.6) 
        categoryMap["c_BO"].inCategory = true;

    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < 1.6) 
        categoryMap["c_BO"].inCategory = true;

    // Barrel Endcap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= 1.6) 
        categoryMap["c_BE"].inCategory = true;

    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= 1.6) 
        categoryMap["c_BE"].inCategory = true;

    // Overlap Endcap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < 1.6 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= 1.6) 
        categoryMap["c_OE"].inCategory = true;

    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= 0.8 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < 1.6 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= 1.6) 
        categoryMap["c_OE"].inCategory = true;

    // jet category selection
    if(vars.validJets.size() >= 2)
    {
        //if(vars.validJets[0].Pt() > cLeadPtMin && vars.validJets[1].Pt() > cSubleadPtMin) // No MET for now
        if(vars.validJets[0].Pt() > cLeadPtMin && vars.validJets[1].Pt() > cSubleadPtMin && vars.met->pt < cMETMax)
        {
            categoryMap["c_2_Jet"].inCategory = true;
            double mjj_max = -1;

            for(unsigned int i=0; i<vars.validJets.size(); i++)
            {
                if(!(vars.validJets[i].Pt() > cLeadPtMin)) break;
                for(unsigned int j=i+1; j<vars.validJets.size(); j++)
                {
                    double dEtajj = vars.validJets[i].Eta() - vars.validJets[j].Eta();
                    double mjj = (vars.validJets[i] + vars.validJets[j]).M();
                    if(mjj > mjj_max) mjj_max = mjj;
                    if(mjj > cDijetMassMinVBFT && TMath::Abs(dEtajj) > cDijetDeltaEtaMinVBFT)
                    { 
                        categoryMap["c_2_Jet_VBF_Tight"].inCategory = true; 
                        return; 
                    }
                }
            }

            if(mjj_max > cDijetMassMinGGFT && vars.dimuCand->pt > cDimuPtMinGGFT)
            { 
                categoryMap["c_2_Jet_GGF_Tight"].inCategory = true; 
                return; 
            }
            else
            { 
                categoryMap["c_2_Jet_VBF_Loose"].inCategory = true; 
                return; 
            }
        }
    }
    if(!categoryMap["c_2_Jet"].inCategory) // fails 2jet preselection enters 01 categories
    {
        categoryMap["c_01_Jet"].inCategory = true;
        if(vars.dimuCand->pt > cDimuPtMin01T){ categoryMap["c_01_Jet_Tight"].inCategory = true;}
        else{ categoryMap["c_01_Jet_Loose"].inCategory = true; }

        // Geometric categories for 01_Jet categories
        // tight
        if(categoryMap["c_01_Jet_Tight"].inCategory && categoryMap["c_BB"].inCategory) categoryMap["c_01_Jet_Tight_BB"].inCategory = true;
        if(categoryMap["c_01_Jet_Tight"].inCategory && categoryMap["c_BO"].inCategory) categoryMap["c_01_Jet_Tight_BO"].inCategory = true;
        if(categoryMap["c_01_Jet_Tight"].inCategory && categoryMap["c_BE"].inCategory) categoryMap["c_01_Jet_Tight_BE"].inCategory = true;
        if(categoryMap["c_01_Jet_Tight"].inCategory && categoryMap["c_OO"].inCategory) categoryMap["c_01_Jet_Tight_OO"].inCategory = true;
        if(categoryMap["c_01_Jet_Tight"].inCategory && categoryMap["c_OE"].inCategory) categoryMap["c_01_Jet_Tight_OE"].inCategory = true;
        if(categoryMap["c_01_Jet_Tight"].inCategory && categoryMap["c_EE"].inCategory) categoryMap["c_01_Jet_Tight_EE"].inCategory = true;

        // loose
        if(categoryMap["c_01_Jet_Loose"].inCategory && categoryMap["c_BB"].inCategory) categoryMap["c_01_Jet_Loose_BB"].inCategory = true;
        if(categoryMap["c_01_Jet_Loose"].inCategory && categoryMap["c_BO"].inCategory) categoryMap["c_01_Jet_Loose_BO"].inCategory = true;
        if(categoryMap["c_01_Jet_Loose"].inCategory && categoryMap["c_BE"].inCategory) categoryMap["c_01_Jet_Loose_BE"].inCategory = true;
        if(categoryMap["c_01_Jet_Loose"].inCategory && categoryMap["c_OO"].inCategory) categoryMap["c_01_Jet_Loose_OO"].inCategory = true;
        if(categoryMap["c_01_Jet_Loose"].inCategory && categoryMap["c_OE"].inCategory) categoryMap["c_01_Jet_Loose_OE"].inCategory = true;
        if(categoryMap["c_01_Jet_Loose"].inCategory && categoryMap["c_EE"].inCategory) categoryMap["c_01_Jet_Loose_EE"].inCategory = true;
    }

}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelectionSynch_________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelectionSynch::initCategoryMap()
{
// Initialize the categories
    categoryMap["c_ALL"] = Category("c_ALL");

    // intermediate categories to make things easier
    categoryMap["c_2_Jet"] = Category("c_2_Jet", true);
    categoryMap["c_01_Jet"] = Category("c_01_Jet", true);

    categoryMap["c_2_Jet_VBF_Tight"] = Category("c_2_Jet_VBF_Tight");
    categoryMap["c_2_Jet_VBF_Loose"] = Category("c_2_Jet_VBF_Loose");
    categoryMap["c_2_Jet_GGF_Tight"] = Category("c_2_Jet_GGF_Tight");
    categoryMap["c_01_Jet_Tight"] = Category("c_01_Jet_Tight");
    categoryMap["c_01_Jet_Loose"] = Category("c_01_Jet_Loose");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionSynch::CategorySelectionSynch()
{
// initialize the default cut values in the constructor

    initCategoryMap();

    // Preselection
    cLeadPtMin = 40;
    cSubleadPtMin = 30;
    cMETMax = 40;

    // VBF Tight
    cDijetMassMinVBFT = 650;
    cDijetDeltaEtaMinVBFT = 3.5;

    // VBF Loose

    // GGF Tight
    cDijetMassMinGGFT = 250;
    cDimuPtMinGGFT = 50;

    // 01Tight
    cDimuPtMin01T = 25;

    // 01Loose
    
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionSynch::CategorySelectionSynch(float leadPtMin, float subleadPtMin, float METMax, float dijetMassMinVBFT, float dijetDeltaEtaMinVBFT, 
                                             float dijetMassMinGGFT, float dimuPtMinGGFT, float dimuPtMin01T)
{
// Initialize custom cut values

    initCategoryMap();

    // Preselection
    cLeadPtMin = leadPtMin;
    cSubleadPtMin = subleadPtMin;
    cMETMax = METMax;

    // VBF Tight
    cDijetMassMinVBFT = dijetMassMinVBFT;
    cDijetDeltaEtaMinVBFT = dijetDeltaEtaMinVBFT;

    // VBF Loose

    // GGF Tight
    cDijetMassMinGGFT = dijetMassMinGGFT;
    cDimuPtMinGGFT = dimuPtMinGGFT;

    // 01Tight
    cDimuPtMin01T = dimuPtMin01T;

    // 01Loose
}


///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionSynch::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // Inclusive category, all events that passed the selection cuts
    categoryMap["c_ALL"].inCategory = true;

    // jet category selection
    if(vars.validJets.size() >= 2)
    {
        //if(leadJet.Pt() > cLeadPtMin && subleadJet.Pt() > cSubleadPtMin && vars.met->pt < cMETMax && vars.validBJets.size() == 0)
        if(vars.validJets[0].Pt() > cLeadPtMin && vars.validJets[1].Pt() > cSubleadPtMin && vars.validBJets.size() == 0) // No MET for now
        {
            categoryMap["c_2_Jet"].inCategory = true;
            double mjj_max = -1;

            for(unsigned int i=0; i<vars.validJets.size(); i++)
            {
                if(!(vars.validJets[i].Pt() > cLeadPtMin)) break;
                for(unsigned int j=i+1; j<vars.validJets.size(); j++)
                {
                    double dEtajj = vars.validJets[i].Eta() - vars.validJets[j].Eta();
                    double mjj = (vars.validJets[i] + vars.validJets[j]).M();
                    if(mjj > mjj_max) mjj_max = mjj;
                    if(mjj > cDijetMassMinVBFT && TMath::Abs(dEtajj) > cDijetDeltaEtaMinVBFT)
                    { 
                        categoryMap["c_2_Jet_VBF_Tight"].inCategory = true; 
                        return; 
                    }
                }
            }

            if(mjj_max > cDijetMassMinGGFT && vars.dimuCand->pt > cDimuPtMinGGFT)
            { 
                categoryMap["c_2_Jet_GGF_Tight"].inCategory = true; 
                return; 
            }
            else
            { 
                categoryMap["c_2_Jet_VBF_Loose"].inCategory = true; 
                return; 
            }
        }
    }
    if(!categoryMap["c_2_Jet"].inCategory) // fails 2jet preselection enters 01 categories
    {
        categoryMap["c_01_Jet"].inCategory = true;
        if(vars.dimuCand->pt > cDimuPtMin01T){ categoryMap["c_01_Jet_Tight"].inCategory = true;}
        else{ categoryMap["c_01_Jet_Loose"].inCategory = true; }
    }

}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelectionFEWZ__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelectionFEWZ::initCategoryMap()
{
// Initialize the categories
    categoryMap["c_Wide"] = Category("c_Wide");
    categoryMap["c_Narrow"] = Category("c_Narrow");

    // intermediate categories to make things easier, don't plot these, hence the true
    categoryMap["c_Central_Central"] = Category("c_Central_Central", true);
    categoryMap["c_Central_Not_Central"] = Category("c_Central_Not_Central", true);

    categoryMap["c_1Jet"] = Category("c_1Jet");
    categoryMap["c_Central_Central_Wide"] = Category("c_Central_Central_Wide");
    categoryMap["c_Central_Not_Central_Wide"] = Category("c_Central_Not_Central_Wide");
    categoryMap["c_1Jet_Wide"] = Category("c_1Jet_Wide");
    categoryMap["c_Central_Central_Narrow"] = Category("c_Central_Central_Narrow");
    categoryMap["c_Central_Not_Central_Narrow"] = Category("c_Central_Not_Central_Narrow");
    categoryMap["c_1Jet_Narrow"] = Category("c_1Jet_Narrow");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ(bool useRecoMu, bool useRecoJets)
{
// Standard values for the FEWZ categorization

    initCategoryMap();
    this->useRecoMu = useRecoMu;
    this->useRecoJets = useRecoJets;

    // init cut values
    cMassSplit = 160;
    cEtaCentralSplit = 0.8;
    cJetPtMin = 30;
    cJetEtaMax = 4.7;
}


///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ()
{
// Standard values for the FEWZ categorization

    initCategoryMap();
    useRecoMu = false;
    useRecoJets = false;

    // init cut values
    cMassSplit = 160;
    cEtaCentralSplit = 0.8;
    cJetPtMin = 30;
    cJetEtaMax = 4.7;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ(bool useRecoMu, bool useRecoJets, float massSplit, float etaCentralSplit, float jetPtMin, float jetEtaMax)
{
// Custom values for the cuts

    initCategoryMap();
    this->useRecoMu = useRecoMu;
    this->useRecoJets = useRecoJets;

    // init cut values
    cMassSplit = massSplit;
    cEtaCentralSplit = etaCentralSplit;
    cJetPtMin = jetPtMin;
    cJetEtaMax = jetEtaMax;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionFEWZ::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    float dimu_mass = vars.dimuCand->mass;
    float eta0 = vars.muons->at(vars.dimuCand->iMu1).eta;
    float eta1 = vars.muons->at(vars.dimuCand->iMu2).eta;
    unsigned int njets = vars.validJets.size();

    // Should cut out all events that don't fall into the wide mass window in earlier selection stage
    // All events that pass are in window of min to max
    categoryMap["c_Wide"].inCategory = true;

    // Narrow goes from min to cMassSplit
    if(dimu_mass < cMassSplit) categoryMap["c_Narrow"].inCategory = true;

    // Both central
    if(TMath::Abs(eta0) < 0.8 && TMath::Abs(eta1) < 0.8) categoryMap["c_Central_Central"].inCategory = true;

    // Not both, but at least one is central
    else if(TMath::Abs(eta0) < 0.8 || TMath::Abs(eta1) < 0.8) categoryMap["c_Central_Not_Central"].inCategory = true;

    // One category that passes basic selections and has exactly one jet
    if(njets == 1) categoryMap["c_1Jet"].inCategory = true; 

    // Final Categories ///////////////////////////////////////////////////////
    if(categoryMap["c_Wide"].inCategory && categoryMap["c_Central_Central"].inCategory) categoryMap["c_Central_Central_Wide"].inCategory = true;
    if(categoryMap["c_Narrow"].inCategory && categoryMap["c_Central_Central"].inCategory) categoryMap["c_Central_Central_Narrow"].inCategory = true;

    if(categoryMap["c_Wide"].inCategory && categoryMap["c_Central_Not_Central"].inCategory) categoryMap["c_Central_Not_Central_Wide"].inCategory = true;
    if(categoryMap["c_Narrow"].inCategory && categoryMap["c_Central_Not_Central"].inCategory) categoryMap["c_Central_Not_Central_Narrow"].inCategory = true;

    if(categoryMap["c_Wide"].inCategory && categoryMap["c_1Jet"].inCategory) categoryMap["c_1Jet_Wide"].inCategory = true;
    if(categoryMap["c_Narrow"].inCategory && categoryMap["c_1Jet"].inCategory) categoryMap["c_1Jet_Narrow"].inCategory = true;

    return;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________LotsOfCategoriesRun2__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void LotsOfCategoriesRun2::initCategoryMap()
{
// Initialize the categories

    geometricNames = {TString("BB"), TString("OO"), TString("EE"), TString("BO"), TString("BE"), TString("OE")};

    ///////////////// INCLUSIVE //////////////////////////////
    categoryMap["c_ALL"] = Category("c_ALL");

    ///////////////// GEOMETRY //////////////////////////////
    categoryMap["c_BB"] = Category("c_BB", true);
    categoryMap["c_BO"] = Category("c_BO", true);
    categoryMap["c_BE"] = Category("c_BE", true);
    categoryMap["c_OO"] = Category("c_OO", true);
    categoryMap["c_OE"] = Category("c_OE", true);
    categoryMap["c_EE"] = Category("c_EE", true);

    ///////////////// PRESELECTION //////////////////////////////
    categoryMap["c_Preselection_Pass"] = Category("c_Preselection_Pass");

      ///////////////// AT LEAST ONE B-JET //////////////////////////////
      categoryMap["c_1b"] = Category("c_1b");

        ///////////////// TTH 2 EXTRA LEPTONS //////////////////////////////
        categoryMap["c_1b_TTH"]        = Category("c_1b_TTH");
        categoryMap["c_1b_TTH_2e"]     = Category("c_1b_TTH_2e",     true);
        categoryMap["c_1b_TTH_1e_1mu"] = Category("c_1b_TTH_1e_1mu", true);
        categoryMap["c_1b_TTH_2mu"]    = Category("c_1b_TTH_2mu",    true);
  
        ///////////////// TTH/BBH 0 EXTRA LEPTONS (HADRONS)///////////////////
        categoryMap["c_1b_TTH_BBH"]              = Category("c_1b_TTH_BBH");
        categoryMap["c_1b_TTH_BBH_Tight"]        = Category("c_1b_TTH_BBH_Tight",        true);
        categoryMap["c_1b_TTH_BBH_V_Hadronic_H"] = Category("c_1b_TTH_BBH_V_Hadronic_H", true);
  
        ///////////////// 1 B-JET EVENTS THAT DONT FIT ELSEWHERE////////////////
        categoryMap["c_1b_Leftovers"] = Category("c_1b_Leftovers");
  
      ///////////////// NO B-JETS //////////////////////////////
      categoryMap["c_0b"] = Category("c_0b");
  
        ///////////////// NOT V(lept)H (VBF, gF, V(had)H, ZvvH) /////////////////////
        categoryMap["c_0b_nonVlH"] = Category("c_0b_nonVlH");

        ///////////////// 2-jet (VBF, V(had)H, gF) //////////////////////////////
        categoryMap["c_0b_nonVlH_2j"]              = Category("c_0b_nonVlH_2j");
        categoryMap["c_0b_nonVlH_2j_VBF_Tight"]    = Category("c_0b_nonVlH_2j_VBF_Tight");
        categoryMap["c_0b_nonVlH_2j_VBF_Loose"]    = Category("c_0b_nonVlH_2j_VBF_Loose");
        categoryMap["c_0b_nonVlH_2j_V_Hadronic_H"] = Category("c_0b_nonVlH_2j_V_Hadronic_H");
        categoryMap["c_0b_nonVlH_2j_gF"]           = Category("c_0b_nonVlH_2j_gF");

        ///////////////// 01-jet (gF Tight, gF Loose, ZvvH) //////////////////////////////
        categoryMap["c_0b_nonVlH_01j"]          = Category("c_0b_nonVlH_01j");

          // ZvvH
          categoryMap["c_0b_nonVlH_01j_ZvvH"]     = Category("c_0b_nonVlH_01j_ZvvH");

          // gF Tight
          categoryMap["c_0b_nonVlH_01j_gF_Tight"] = Category("c_0b_nonVlH_01j_gF_Tight");

            // gF Tight Geometrized
            categoryMap["c_0b_nonVlH_01j_gF_Tight_BB"] = Category("c_0b_nonVlH_01j_gF_Tight_BB");
            categoryMap["c_0b_nonVlH_01j_gF_Tight_BO"] = Category("c_0b_nonVlH_01j_gF_Tight_BO");
            categoryMap["c_0b_nonVlH_01j_gF_Tight_BE"] = Category("c_0b_nonVlH_01j_gF_Tight_BE");
            categoryMap["c_0b_nonVlH_01j_gF_Tight_OO"] = Category("c_0b_nonVlH_01j_gF_Tight_OO");
            categoryMap["c_0b_nonVlH_01j_gF_Tight_OE"] = Category("c_0b_nonVlH_01j_gF_Tight_OE");
            categoryMap["c_0b_nonVlH_01j_gF_Tight_EE"] = Category("c_0b_nonVlH_01j_gF_Tight_EE");

          // gF Loose
          categoryMap["c_0b_nonVlH_01j_gF_Loose"] = Category("c_0b_nonVlH_01j_gF_Loose");

            // gF Loose Geometrized
            categoryMap["c_0b_nonVlH_01j_gF_Loose_BB"] = Category("c_0b_nonVlH_01j_gF_Loose_BB");
            categoryMap["c_0b_nonVlH_01j_gF_Loose_BO"] = Category("c_0b_nonVlH_01j_gF_Loose_BO");
            categoryMap["c_0b_nonVlH_01j_gF_Loose_BE"] = Category("c_0b_nonVlH_01j_gF_Loose_BE");
            categoryMap["c_0b_nonVlH_01j_gF_Loose_OO"] = Category("c_0b_nonVlH_01j_gF_Loose_OO");
            categoryMap["c_0b_nonVlH_01j_gF_Loose_OE"] = Category("c_0b_nonVlH_01j_gF_Loose_OE");
            categoryMap["c_0b_nonVlH_01j_gF_Loose_EE"] = Category("c_0b_nonVlH_01j_gF_Loose_EE");

        ///////////////// V(lept)H (...) /////////////////////
        categoryMap["c_0b_VlH"]    = Category("c_0b_VlH");

            // VlH according to the V decays
            categoryMap["c_0b_VlH_We"]         = Category("c_0b_VlH_We");
            categoryMap["c_0b_VlH_Wmu"]        = Category("c_0b_VlH_Wmu");
            categoryMap["c_0b_VlH_Ztautau"]        = Category("c_0b_VlH_Ztautau");
            categoryMap["c_0b_VlH_Zmumu"]      = Category("c_0b_VlH_Zmumu");
            categoryMap["c_0b_VlH_Zee"]        = Category("c_0b_VlH_Zee");
            categoryMap["c_0b_VlH_Leftovers"]  = Category("c_0b_VlH_Leftovers");

    ///////////////// FAIL PRESELECTION //////////////////////////////
    categoryMap["c_Preselection_Fail"] = Category("c_Preselection_Fail");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

LotsOfCategoriesRun2::LotsOfCategoriesRun2()
{
// initialize the default cut values in the constructor

    initCategoryMap();

        // Preselection
        c_pre_numExtraLeptonsMax = 2;
        c_pre_useTau = false;
        c_pre_numBJetsMin = 1;

        // b-jet categories
        c_1b_numExtraLeptons_tth     = 2; 
        c_1b_numExtraLeptons_tth_bbh = 0; 

            // ttH categories
            c_tth_numExtraElectrons = -999;

            // tth-bbh categories
            c_tth_bbh_mbb      = -999;
            c_tth_bbh_mt_bMET  = -999;
            c_tth_bbh_MET      = -999;

        // right side of b-jet test, fewer b-jets
        // 0 extra leptons, 1 or 2 extra leptons
        c_0b_numExtraLeptonsMin = 1; 

            // nonVlH (0 extra leptons)
            c_0b_nonVlH_njetsMin = 2;
    
                // >= c_0b_nonVlH_njets
                c_0b_nonVlH_2j_mjj_min_vbfTight    = 500;
                c_0b_nonVlH_2j_dEtajj_min_vbfTight = 2.5;

                c_0b_nonVlH_2j_mjj_min_vbfLoose    = 250;
                c_0b_nonVlH_2j_dEtajj_min_vbfLoose = 2.5;

                c_0b_nonVlH_2j_mjj_min_VhH = 60;
                c_0b_nonVlH_2j_mjj_max_VhH = 110;
                c_0b_nonVlH_2j_dEtajjMuMu_max_VhH = 1.5;

                // < c_0b_nonVlH_njets
                c_0b_nonVlH_01j_MET_min_ZvvH = 40;
                c_0b_nonVlH_01j_dimuPt_min_gfTight = 25;

                // muon geometry categorization for gf categories
                c_geo_bmax = 0.8;
                c_geo_omax = 1.6;
                c_geo_emax = 2.4;

            // VlH (1 or 2 extra leptons)
            c_0b_VlH_MET_min = 40;

                // VlH according to the V decays
                c_0b_VlH_We_num_e   = 1;         
                c_0b_VlH_We_num_mu  = 0;         

                c_0b_VlH_Wmu_num_e  = 0;        
                c_0b_VlH_Wmu_num_mu = 1;        

                c_0b_VlH_Ztautau_num_e  = 1;        
                c_0b_VlH_Ztautau_num_mu = 1;        

                c_0b_VlH_Zmumu_num_e  = 0;
                c_0b_VlH_Zmumu_num_mu = 2;

                c_0b_VlH_Zee_num_e  = 2; 
                c_0b_VlH_Zee_num_mu = 0; 

}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void LotsOfCategoriesRun2::evaluateMuonGeometry(VarSet& vars)
{
// Determine which geometric category the event belongs to

    // Geometric Categories
    // Barrel Barrel
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < c_geo_bmax && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < c_geo_bmax) 
        categoryMap["c_BB"].inCategory = true;

    // Overlap Overlap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta)>=c_geo_bmax 
      && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta)<c_geo_omax 
      && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta)>=c_geo_bmax 
      && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta)<c_geo_omax) 
        categoryMap["c_OO"].inCategory = true;

    // Endcap Endcap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= c_geo_omax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= c_geo_omax) 
        categoryMap["c_EE"].inCategory = true;

    // Barrel Overlap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < c_geo_bmax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= c_geo_bmax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < c_geo_omax) 
        categoryMap["c_BO"].inCategory = true;

    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < c_geo_bmax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= c_geo_bmax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < c_geo_omax) 
        categoryMap["c_BO"].inCategory = true;

    // Barrel Endcap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < c_geo_bmax 
      && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= c_geo_omax) 
        categoryMap["c_BE"].inCategory = true;

    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < c_geo_bmax && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= c_geo_omax) 
        categoryMap["c_BE"].inCategory = true;

    // Overlap Endcap
    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= c_geo_bmax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) < c_geo_omax 
      && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= c_geo_omax) 
        categoryMap["c_OE"].inCategory = true;

    if(TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) >= c_geo_bmax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu2).eta) < c_geo_omax 
       && TMath::Abs(vars.muons->at(vars.dimuCand->iMu1).eta) >= c_geo_omax) 
        categoryMap["c_OE"].inCategory = true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void LotsOfCategoriesRun2::evaluate(VarSet& vars)
{
    ///////////////// INCLUSIVE //////////////////////////////
    categoryMap["c_ALL"].inCategory = true;

    // figure out bb,oo,ee,bo,be,oe
    ///////////////// MUON GEOMETRY //////////////////////////////
    evaluateMuonGeometry(vars);

    ///////////////// PRESELECTION //////////////////////////////
    if(vars.validExtraMuons.size() + vars.validElectrons.size() <= c_pre_numExtraLeptonsMax) 
        categoryMap["c_Preselection_Pass"].inCategory = true;
    else
        categoryMap["c_Preselection_Fail"].inCategory = true;

   // Determine whether we are in the at least 1b-jet categories or 0b-jet categories
   if(categoryMap["c_Preselection_Pass"].inCategory)
   {
       //std::cout << "    pass preselection..." << std::endl;
       if(vars.validBJets.size() >= c_pre_numBJetsMin) 
           categoryMap["c_1b"].inCategory = true;
       else
           categoryMap["c_0b"].inCategory = true;
   }

       ///////////////// 1b CATEGORIES //////////////////////////////
       if(categoryMap["c_1b"].inCategory)
       {
           //std::cout << "    pass 1b..." << std::endl;
           if(vars.validExtraMuons.size() + vars.validElectrons.size() == c_1b_numExtraLeptons_tth)
               categoryMap["c_1b_TTH"].inCategory = true;
           else if(vars.validExtraMuons.size() + vars.validElectrons.size() == c_1b_numExtraLeptons_tth_bbh)
               categoryMap["c_1b_TTH_BBH"].inCategory = true;
           else 
               categoryMap["c_1b_Leftovers"].inCategory = true;
       }
           ///////////////// 1b-TTH (2 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["c_1b_TTH"].inCategory)
           {
               //std::cout << "    pass 1b TTH..." << std::endl;
           }
           ///////////////// 1b_TTH_BBH (0 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["c_1b_TTH_BBH"].inCategory)
           {
               //std::cout << "    pass 1b TTH_BBH..." << std::endl;
           }

       ///////////////// 0b CATEGORIES //////////////////////////////
       if(categoryMap["c_0b"].inCategory)
       {
           //std::cout << "    pass 0b..." << std::endl;
           // output event information  here to debug... we have 70% in this category for N_valid_whatevers in data
           if(vars.validExtraMuons.size() + vars.validElectrons.size() >= c_0b_numExtraLeptonsMin)
               categoryMap["c_0b_VlH"].inCategory = true;
           else
               categoryMap["c_0b_nonVlH"].inCategory = true;
       }
           ///////////////// 0b-VlH    (1,2 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["c_0b_VlH"].inCategory)
           {
               //std::cout << "    pass 0b_VlH..." << std::endl;
               if(vars.met->pt >= c_0b_VlH_MET_min)
               {
                    if(vars.validElectrons.size() == c_0b_VlH_We_num_e && vars.validExtraMuons.size() == c_0b_VlH_We_num_mu)    
                        categoryMap["c_0b_VlH_We"].inCategory = true;

                    else if(vars.validElectrons.size() == c_0b_VlH_Wmu_num_e && vars.validExtraMuons.size() == c_0b_VlH_Wmu_num_mu)    
                        categoryMap["c_0b_VlH_Wmu"].inCategory = true;

                    else if(vars.validElectrons.size() == c_0b_VlH_Ztautau_num_e && vars.validExtraMuons.size() == c_0b_VlH_Ztautau_num_mu)    
                        categoryMap["c_0b_VlH_Ztautau"].inCategory = true;

                    else    
                        categoryMap["c_0b_VlH_Leftovers"].inCategory = true;
               }
               else
               {
                    if(vars.validElectrons.size() == c_0b_VlH_Zmumu_num_e && vars.validExtraMuons.size() == c_0b_VlH_Zmumu_num_mu)    
                        categoryMap["c_0b_VlH_Zmumu"].inCategory = true;

                    else if(vars.validElectrons.size() == c_0b_VlH_Zee_num_e && vars.validExtraMuons.size() == c_0b_VlH_Zee_num_mu)    
                        categoryMap["c_0b_VlH_Zee"].inCategory = true;

                    else    
                        categoryMap["c_0b_VlH_Leftovers"].inCategory = true;
               }
           }

           ///////////////// 0b-nonVlH (0 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["c_0b_nonVlH"].inCategory)
           {
               //std::cout << "    pass 0b_non_VlH..." << std::endl;
               if(vars.validJets.size() >= c_0b_nonVlH_njetsMin)
                   categoryMap["c_0b_nonVlH_2j"].inCategory = true;
               else
                   categoryMap["c_0b_nonVlH_01j"].inCategory = true;
           }
               ///////////////// 0b-nonVlH_2j (0 extra lept) CATEGORIES //////////////////////////////
               if(categoryMap["c_0b_nonVlH_2j"].inCategory)
               {
                   //std::cout << "    pass 0b_non_VlHi_2j..." << std::endl;
                   TLorentzVector leadJet    = vars.validJets[0];
                   TLorentzVector subleadJet = vars.validJets[1];
                   TLorentzVector dijet      = leadJet + subleadJet;
           
                   float dEta = TMath::Abs(leadJet.Eta() - subleadJet.Eta());
                   float dijetMass = dijet.M();
                   float dEtajjMuMu = TMath::Abs(dijet.Eta() - vars.dimuCand->eta); 

                   if(dijetMass > c_0b_nonVlH_2j_mjj_min_vbfTight && dEta > c_0b_nonVlH_2j_dEtajj_min_vbfTight)
                       categoryMap["c_0b_nonVlH_2j_VBF_Tight"].inCategory = true; 

                   else if(dijetMass > c_0b_nonVlH_2j_mjj_min_vbfLoose && dEta > c_0b_nonVlH_2j_dEtajj_min_vbfLoose)
                       categoryMap["c_0b_nonVlH_2j_VBF_Loose"].inCategory = true; 

                   else if(dijetMass > c_0b_nonVlH_2j_mjj_min_VhH && dijetMass < c_0b_nonVlH_2j_mjj_max_VhH && dEtajjMuMu < c_0b_nonVlH_2j_dEtajjMuMu_max_VhH)
                       categoryMap["c_0b_nonVlH_2j_V_Hadronic_H"].inCategory = true; 

                   else
                       categoryMap["c_0b_nonVlH_2j_gF"].inCategory = true; 
               }

               ///////////////// 0b-nonVlH_01j (0 extra lept) CATEGORIES //////////////////////////////
               if(categoryMap["c_0b_nonVlH_01j"].inCategory)
               {
                   //std::cout << "    pass 0b_nonVlH_01j..." << std::endl;
                   if(vars.met->pt > c_0b_nonVlH_01j_MET_min_ZvvH)
                       categoryMap["c_0b_nonVlH_01j_ZvvH"].inCategory = true; 

                   else if(vars.dimuCand->pt >= c_0b_nonVlH_01j_dimuPt_min_gfTight)
                       categoryMap["c_0b_nonVlH_01j_gF_Tight"].inCategory = true; 

                   else
                       categoryMap["c_0b_nonVlH_01j_gF_Loose"].inCategory = true; 
               }
                   ///////////////// Geometrized 0b-nonVlH_01j_gF CATEGORIES //////////////////////////////
                   if(categoryMap["c_0b_nonVlH_01j_gF_Tight"].inCategory)
                   {
                       for(auto& geoName: geometricNames)
                           if(categoryMap["c_"+geoName].inCategory) categoryMap["c_0b_nonVlH_01j_gF_Tight_"+geoName].inCategory = true;
                   }
                   if(categoryMap["c_0b_nonVlH_01j_gF_Loose"].inCategory)
                   {
                       for(auto& geoName: geometricNames)
                           if(categoryMap["c_"+geoName].inCategory) categoryMap["c_0b_nonVlH_01j_gF_Loose_"+geoName].inCategory = true;
                   }
}
