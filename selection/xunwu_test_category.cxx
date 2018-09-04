#ifndef add_category_selection
#define add_category_selection

#include "VarSet.h"
#include "TH1D.h"
#include "TList.h"
#include "TXMLEngine.h"
#include <map>
#include <utility>
#include <iostream>

class Category
{
    public:
	Category(){};
	~Category(){};

	Category(TString key)
	{	
	    this->key = key;
	    this->name = key;
	    this->hide =false;
	}
	
	std::map<TString, std::vector< std::pair<int,long long int> > >eventsMap;

	bool inCategory = false;
	bool isTerminal = false;
	bool hide = false;
	
	TString key;
	TString name = "";

	std::map<TString, TH1D*> histoMap;
	TList* histoList = new TList();
	TList* signalList = new TList();
	TList* bkgList = new TList();
	TList* dataList = new TList();
};

class Categorizer
{
    public:
	std::map<TString, Category> categoryMap;
	virtual void initCategoryMap() = 0 ;
	virtual void evaluate(VarSet& vars)=0;
	void reset()
	{
	    for (auto &entry:categoryMap)
		entry.second.inCategory = false;
	};

	void outputResults()
	{
	    for (auto &entry:categoryMap)
		if(!entry.second.hide && entry.second.name != entry.second.key)
		    std::cout<<  " ( "  <<entry.second.name<<  " ) "  <<entry.first<< " : " <<entry.second.inCategory<< std::endl;
		else if (!entry.second.hide)
		    std::cout<< "   " <<entry.first<<  " : "  <<entry.second.inCategory<<std::endl;
	    std::cout<<std::endl;
	};


};

class x1_Categorizer : public Categorizer
{
    public:
	x1_Categorizer();

	float LeadPtMin;
	float SubleadPtMin;
	float METMax;
	float DimuMassMin;

	void evaluate(VarSet& vars);
	void initCategoryMap();
};


#endif



#include "ParticleTools.h"
#include "EventTools.h"
#include <sstream>

x1_Categorizer::x1_Categorizer()
{
    initCategoryMap();

    LeadPtMin=60;
    SubleadPtMin = 50;
    METMax = 40;
    DimuMassMin = 80;
}

void x1_Categorizer::initCategoryMap()
{
    categoryMap["c_ALL"]=Category("c_ALL");

    categoryMap["x1"]=Category("x1");
}

void x1_Categorizer::evaluate(VarSet& vars)
{
    categoryMap["c_ALL"].inCategory =true;

    MuonInfo& mu1 = vars.muons->at(vars.dimuCand->iMu1);
    MuonInfo& mu2 = vars.muons->at(vars.dimuCand->iMu2);
    if(vars.dimuCand->mass)
	if( (mu1.pt>LeadPtMin && mu2.pt>SubleadPtMin) || (mu1.pt>SubleadPtMin && mu2.pt>LeadPtMin))
	    categoryMap["x1"].inCategory = true;


}


