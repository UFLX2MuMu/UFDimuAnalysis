///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// MassCalibration.cxx                                                   //
// ======================================================================//
// Plot Z Mass Resolution vs mu+/- eta, mu+/- phi, etc.                  //
// Creates mass histograms in different bins of some x variable.         //
// The object can then fit the histograms with a Voigtian and plot       //
// the voigt fit mean/resolution vs x.                                   //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "MassCalibration.h"

#include "TMinuit.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"

#include <sstream>
#include <cmath>
#include <iostream>


//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

MassCalibration::MassCalibration()
{
    massname = "mass";
    xname = "phi_plus";
    massmin = 86.2;
    massmax = 96.2;
    massbins = 50;
    xmin = -3.14;
    xmax = 3.14;
    xbins = 25;
    fitsig = 1;

    histos = std::vector<TH1D*>();
    binning = std::vector<Float_t>();
    vfis = std::vector<VoigtFitInfo>();
    init();
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

MassCalibration::MassCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                           Float_t massmax, Int_t massbins, Float_t xmin, Float_t xmax, Int_t xbins, Float_t voigt_gamma)
{
    this->xname = xname;
    this->massname = massname;
    this->massmin = massmin;
    this->massmax = massmax;
    this->massbins = massbins;
    this->xmin = xmin;
    this->xmax = xmax;
    this->xbins = xbins;
 
    this->fitsig = fitsig;
    this->voigt_gamma = voigt_gamma;

    histos = std::vector<TH1D*>();
    binning = std::vector<Float_t>();
    vfis = std::vector<VoigtFitInfo>();
    init();
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

MassCalibration::MassCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                           Float_t massmax, Int_t massbins, std::vector<Float_t> binning, Float_t voigt_gamma)
{
    this->xname = xname;
    this->massname = massname;
    this->massmin = massmin;
    this->massmax = massmax;
    this->massbins = massbins;
    this->binning = binning;
 
    this->fitsig = fitsig;
    this->voigt_gamma = voigt_gamma;

    histos = std::vector<TH1D*>();
    vfis = std::vector<VoigtFitInfo>();
    init(this->binning);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void MassCalibration::init()
{
    TString basename = "Z_Peak_Calibration_";
    basename+=massname+"_"+xname;
    variableBinning=false;

    binning.push_back(xmin);
    Float_t interval = (xmax - xmin)/xbins;
    Float_t upper_bin_edge = xmin;

    // Set up the bin edges 
    for(unsigned int i=0; i<xbins; i++)
    {
        upper_bin_edge += interval;
        binning.push_back(upper_bin_edge);
    }

    // print binning for debugging
    //std::cout << "printing binning..." << std::endl;
    //for(unsigned int i=0; i<binning.size(); i++)
    //{
    //    std::cout << binning[i] << std::endl;
    //}

    // Set up the TH1Ds
    for(unsigned int i=0; i<binning.size()-1; i++)
    {
        TString range = Form("_%5.2f_to_%5.2f", binning[i], binning[i+1]);
        TString title = basename+range;
        title.ReplaceAll(" ","");
        title.ReplaceAll(".","p");
        histos.push_back(new TH1D(title, title, massbins, massmin, massmax));
        histos[i]->GetXaxis()->SetTitle(massname);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void MassCalibration::init(std::vector<Float_t>& binning)
{
    TString basename = "Z_Peak_Calibration_";
    basename+=massname+"_"+xname;

    if(binning.size() < 2) std::cout << "!!!! ERROR in MassCalibration::init, binning.size() < 2 !!! \n";

    variableBinning = true;
    xmin = binning[0];
    xmax = binning[binning.size()-1];
    xbins = binning.size()-1;

    // Set up the bin edges
    // Set up the TH1Ds
    for(unsigned int i=0; i<binning.size()-1; i++)
    {
        TString range = Form("_%5.2f_to_%5.2f", binning[i], binning[i+1]);
        TString title = basename+range;
        title.ReplaceAll(" ","");
        title.ReplaceAll(".","p");
        histos.push_back(new TH1D(title, title, massbins, massmin, massmax));
        histos[i]->GetXaxis()->SetTitle(massname);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

int MassCalibration::whichTH1D(Float_t xvalue)
{
    // constant binning
    if(!variableBinning)
    {
        Float_t interval = (xmax - xmin)/xbins;
        int bin = (xvalue-xmin)/interval;

        if(xvalue == xmin) bin = 0;
        if(xvalue == xmax) bin = binning.size()-2;
        if(xvalue < xmax && bin > binning.size()-2) bin = binning.size()-2;

        if(xvalue < xmin || xvalue > xmax) bin = -1;

        return bin;
    }
    // variable binning
    else
    {
        if(xvalue < xmin || xvalue > xmax) return -1;

        int begin = 0;
        int end = binning.size()-1;
        int check = (begin+end)/2;

        int bin = -9999;
        while(bin < -1)
        {
            if(xvalue >= binning[check] && xvalue < binning[check+1]) // xvalue in check bin 
            {
                bin = check;
                //std::cout << " checkv: " << binning[check] << ", check+1v: " << binning[check+1] << ", xvalue: " << xvalue << std::endl;
            }
            else if(xvalue < binning[check]) // xvalue is in the left half, set end to half-way 
                end = check;

            else if(end == check+1 || check == begin) // xvalue not in check bin and nowhere else to go
            {
                std::cout << "begin: " << begin << ", end: " << end << ", check: " << check << ", xvalue: " << xvalue << std::endl;
                std::cout << "beginv: " << binning[begin] << ", endv: " << binning[end] << ", checkv: " << binning[check] << std::endl;
                bin = -1;
            }

            else if(xvalue >= binning[check]) // xvalue is in the right half, set begin to right half 
                begin = check;

            check = (begin+end)/2;
        }
        return bin;
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void MassCalibration::fill(Float_t xvalue, Float_t massvalue)
{
    int h = whichTH1D(xvalue);
    //if(h>=histos.size())
    //{
    //    std::cout << Form("whichHist: %d, histos.size(): %d, binning.size(): %d, xmin: %f, xmax: %12.10f, xvalue: %12.10f, massvalue: %f \n", 
    //                       h, histos.size(), binning.size(), xmin, xmax, xvalue, massvalue); 
    //}
    if(h>=0) 
        histos[h]->Fill(massvalue);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void MassCalibration::fill(Float_t xvalue, Float_t massvalue, double weight)
{
    int h = whichTH1D(xvalue);
    //if(h>=histos.size())
    //{
    //    std::cout << Form("whichHist: %d, histos.size(): %d, binning.size(): %d, xmin: %f, xmax: %12.10f, xvalue: %12.10f, massvalue: %f \n", 
    //                       h, histos.size(), binning.size(), xmin, xmax, xvalue, massvalue); 
    //}
    if(h>=0) 
        histos[h]->Fill(massvalue, weight);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

VoigtFitInfo MassCalibration::fit(TH1D* inhist, Float_t x, Float_t x_err)
{
    //std::cout << Form("==== Fitting %s ==== \n", inhist->GetTitle());
    
    VoigtFitInfo vfi;
    vfi.x = x;
    vfi.x_err = x_err;
    

    //if histo empty, set mean=91, error=91, resolution=0
    if(inhist->GetEntries() == 0)
    {
        vfi.vmean = 0;
        vfi.vmean_err = 0;
        vfi.vsigma = 0;
        vfi.vsigma_err = 0;

        return vfi;
    }

    // Grab the results we are interested in.
    vfi.hmean = inhist->GetMean();
    vfi.hmean_err = inhist->GetMeanError();
    
    vfi.hrms = inhist->GetRMS();
    vfi.hrms_err = inhist->GetRMSError();
    
    // Display fit info on canvas.
    gStyle->SetOptFit(0011);
    // Fit things.
    TString fitname = TString("fit_")+massname+"_"+xname+Form("_%5.2f", x);
    fitname.ReplaceAll(" ", "");
    fitname.ReplaceAll(".", "p");
    TF1* fit = new TF1(fitname, "[0]*TMath::Voigt(x - [1], [2], [3])", massmin, massmax);
    fit->SetParNames("Constant", "Mean", "Sigma", "Gamma");
    
    float initial_sigma = vfi.hrms - voigt_gamma;
    if(initial_sigma < 0) initial_sigma = vfi.hrms/100;

    // Reasonable initial guesses for the fit parameters
    fit->SetParameters(inhist->GetMaximum(), inhist->GetMean(), initial_sigma, voigt_gamma);
    fit->SetParLimits(2, 0, 10*initial_sigma);
    fit->SetParLimits(1, massmin, massmax);

    //std::cout << Form("sigma = %f, in (%f, %f)", initial_sigma, 0, 2*vfi.hrms);

    // Fix the intrinsic width to the theoretical value
    fit->FixParameter(3, voigt_gamma);

    // Sometimes we have to fit the data a few times before the fit converges
    bool converged = 1;
    int ntries = 0;
    
    // Make sure the fit converges.
   while( converged != 0 ) 
    {
      if(ntries >= 50) break;
      fit->SetParameter(2, initial_sigma);
      inhist->Fit(fitname, "q");

      // Fit to a width of fitsig sigmas
      TFitResultPtr ptr = inhist->Fit(fitname,"q","", fit->GetParameter(1) - fitsig*vfi.hrms, fit->GetParameter(1) + fitsig*vfi.hrms);
      //inhist->Fit(fitname,"q","", massmin, massmax);
      //TString sconverge = gMinuit->fCstatu.Data();
      converged = ptr;
      //converged = sconverge.Contains(TString("CONVERGED"));
      ntries++; 
    }

    vfi.fit = fit;
    vfi.vmean = fit->GetParameter(1);
    vfi.vmean_err = fit->GetParError(1);
    vfi.vsigma = fit->GetParameter(2);
    vfi.vsigma_err = fit->GetParError(2);
    vfi.vgamma = fit->GetParameter(3);
    vfi.vgamma_err = fit->GetParError(3);

    return vfi;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void MassCalibration::fit()
{
    for(unsigned int h=0; h<histos.size(); h++)
    {
        VoigtFitInfo vfi = fit(histos[h], (binning[h+1] + binning[h])/2, (binning[h+1] - binning[h])/2);
        vfis.push_back(vfi);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void MassCalibration::plot()
{
    fit();
    mean_vs_x = new TGraphErrors();
    resolution_vs_x = new TGraphErrors();
    for(unsigned int i=0; i<vfis.size(); i++)
    {
        mean_vs_x->SetPoint(i, vfis[i].x, vfis[i].vmean);
        mean_vs_x->SetPointError(i, vfis[i].x_err, vfis[i].vmean_err);

        resolution_vs_x->SetPoint(i, vfis[i].x, vfis[i].vsigma);
        resolution_vs_x->SetPointError(i, vfis[i].x_err, vfis[i].vsigma_err);
    }
    mean_vs_x->SetTitle(massname+"_fit_mean_vs_"+xname);
    mean_vs_x->SetName("mean_"+massname+"_"+xname);
    mean_vs_x->GetXaxis()->SetTitle(xname);
    mean_vs_x->GetYaxis()->SetTitle(massname+"_mean");

    resolution_vs_x->SetTitle(massname+"_fit_resolution_vs_"+xname);
    resolution_vs_x->SetName("resolution_"+massname+"_"+xname);
    resolution_vs_x->GetXaxis()->SetTitle(xname);
    resolution_vs_x->GetYaxis()->SetTitle(massname+"_resolution");
}
