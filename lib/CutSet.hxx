///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// CutSet.hxx                                                            //
// ======================================================================//
// CutInfo represents the information for a  single cut.                 //
// Can turn cuts on and off get info                                     //
// from this object to create a histogram of the cutvar for validation   //
// or N-1 studies.                                                       //
// CutSet is a collection of CutInfo. Can output the cuts concat them    //
// into a single string. See which ones in the set failed, useful for    //
// synchronization or debugging.                                         //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_CUTSET
#define ADD_CUTSET

#include "VarSet.h"
#include "TString.h"
#include <iostream>

class CutInfo
{
    public:

        CutInfo(){};
        CutInfo(TString tstring, TString name, int id)
        {
            this->tstring = tstring;
            this->name = name;
            this->id = id;
        }

        TString tstring = ""; // represent the set of cuts as a tstring
        TString name = "";    // the name of the cut, useful when plotting
        TString cutname = ""; // The name of the variable in the Selection Object
        bool passed = true;   // flag whether the event passed the cut or not
        int id = -1;          // the number of the cut in the cutset vector
        bool ismin = false;   // the cut is a minimum bound or maximum bound
        bool on = 1;          // whether to use the cut in the N-1 plots
        float value = -9999;  // the value of the variable that was cut on
        int bins = 50;        // the number of bins needed for the histogram
        float min = 0;        // the minimum value for the histogram
        float max = 200;      // the maximum value for the histogram
        float* cutvalue = 0;  // the address of the selection variable defining the cut
                              // i.e. float* cutvalue = &SelectionCuts.cDimuMassMax
                              // useful when optimizing the cuts cut[i].cutvalue = somevalue

        std::vector<float> cutvalues;
      
        CutInfo operator&&(const CutInfo& c)
        {
           CutInfo cutinfo;
           cutinfo.tstring = this->tstring.EqualTo("")?c.tstring:this->tstring + " && " + c.tstring; 
           cutinfo.name = this->name.EqualTo("")?c.tstring:this->name + " and " + c.name; 
           cutinfo.passed = this->passed && c.passed;
           return cutinfo;
        }

        
};

class CutSet
{
    public:

        CutInfo cutinfo;                         // represent the set of cuts as a tstring
        std::vector<CutInfo> cuts;               // the collection of cuts

        void resetCuts()
        {
           for(unsigned int i=0; i<cuts.size(); i++)
           {
               cuts[i].value = -9999;      
               cuts[i].passed = true;
           }
        }

        void outputFailedCuts()
        {
           for(unsigned int i=0; i<cuts.size(); i++)
           {
               if(!cuts[i].passed)
               {
                   std::cout << "        +++ " <<cuts[i].tstring << " FAILED." << std::endl;
                   std::cout << "        +++ " <<cuts[i].name << ": " << cuts[i].value << std::endl;
               }
           }

        }

        void concatCuts(std::vector<CutInfo>& cuts)
        {
            CutInfo c;
            // concatenate active cuts together
            for(unsigned int i=0; i<cuts.size(); i++)
            {
                if(cuts[i].on)
                    c = c && cuts[i];
            } 
            cutinfo = c;
        }

        void concatCuts()
        {
            CutInfo c;
            // concatenate active cuts together
            for(unsigned int i=0; i<cuts.size(); i++)
            {
                if(cuts[i].on)
                    c = c && cuts[i];
            } 
            cutinfo = c;
        }

        void turnOnCut(int i)
        {
            cuts[i].on = true;
        }
};

#endif
