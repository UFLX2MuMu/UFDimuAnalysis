##############################################
# BGSpectrumFitter.py                        #
##############################################
# fit the bkg spectrum using the histos      #
# created by categorize.cxx                  #
#                                            #
##############################################

#============================================
# import
#============================================

import PDFDatabase as pdfs
import prettytable
import string
import re
import argparse
import math as math
import numpy as np
from ROOT import *

import sys
sys.argv.append( '-b-' )

#============================================
# code
#============================================

class BGSpectrumFitter:
# Object to fit the background spectrum in data or MC. Takes in rootfiles from categorize.cxx.
# Can fit blinded or unblinded histograms. Could be used to fit the signal with a few tweaks.

    infilename = ''
    category = ''
    data_hist = 0
    bg_dy_hist = 0
    bg_ttbar_hist = 0
    bg_diboson_hist = 0
    bg_not_dy_hist = 0
    bg_all_hist = 0
    tfile = 0
    nuisance_params = []

    def __init__(self, infilename, category):
        self.infilename = infilename
        self.category = category
        self.tfile = TFile(infilename)
        self.setHists()
    
    def setHists(self):
    # use the naming convention from categorize.cxx to automatically grab the different histograms that might be fit
        self.data_hist      = self.tfile.Get('net_histos/'+self.category+"_Net_Data")
        self.bg_dy_hist      = self.tfile.Get('net_histos/'+self.category+"_Drell_Yan_")
        self.bg_ttbar_hist   = self.tfile.Get('net_histos/'+self.category+"_TTbar_Plus_SingleTop")
        self.bg_diboson_hist = self.tfile.Get('net_histos/'+self.category+"_Diboson_plus")
        self.bg_all_hist  = self.tfile.Get('net_histos/'+self.category+"_Net_Bkg")

        self.bg_not_dy_hist  = self.bg_ttbar_hist.Clone()
        self.bg_not_dy_hist.Add(self.bg_diboson_hist)
        self.bg_not_dy_hist.SetName(self.category+"_TTBar_Diboson_plus")
        self.bg_not_dy_hist.SetTitle(self.category+"_TTBar_Diboson_plus")
    
    
    def getX(self, histo):
    # need the roofit x variable to do the fit, so create it from the histogram
        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78
    
        # Get the dimu mass histogram to use for limit setting
        nbins   = histo.GetNbinsX()
        massmin = histo.GetBinLowEdge(1)
        massmax = massmin + nbins*histo.GetBinWidth(1)
    
        #----------------------------------------
        # create di-muon mass variable
        # syntax:
        # <name>[initial-val, min-val, max-val]
        #----------------------------------------
        x = RooRealVar('x','x',120, massmin, massmax)
        x.SetTitle('m_{#mu#mu}')
        x.setUnit('GeV')
  
        return x
    
    def fit(self, histo, pdfMmumu, x, xmin=100, xmax=150, blinded=True, save=True):
    # fit the spectrum and save the fit and the data/MC histogram to .png and .root

        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78

        # create binned dataset from histogram
        # needs to be named data_obs for higgs combine limit setting
        data = RooDataHist('data_obs', 'data_obs', RooArgList(x), histo)

        #----------------------------------------
        # fit and plot
        #----------------------------------------
        x.setRange("window",xmin, xmax)  # whole window for MC
        x.setRange("left",xmin, 120+0.1)    # exclude signal region for Data
        x.setRange("right",130-0.1, xmax)

      
        if "Data" in histo.GetName() and blinded:
            pdfMmumu.fitTo(data, RooFit.Range("left,right"))
        else: 
            pdfMmumu.fitTo(data, RooFit.Range("window"))

        xframe = x.frame(RooFit.Name(histo.GetName()+"_Fit"), RooFit.Title(histo.GetName()+"_Fit"))
        xframe.GetXaxis().SetNdivisions(505)
        data.plotOn(xframe)
        pdfMmumu.plotOn(xframe, RooFit.Name(pdfMmumu.GetName()))
        pdfMmumu.paramOn(xframe, RooFit.Format("NELU", RooFit.AutoPrecision(2)), RooFit.Layout(0.3, 0.95, 0.92) )
        chi2 = xframe.chiSquare()

        print "chi2    :     %7.3f"               % chi2
        print

        c1 = TCanvas(histo.GetName()+"_"+pdfMmumu.GetName()+"_c", histo.GetName()+"_"+pdfMmumu.GetName(), 10, 10, 600, 600)
        xframe.Draw()
        t = TLatex(.6,.6,"#chi^{2}/ndof = %7.3f" % chi2);  
        t.SetNDC(kTRUE);
        t.Draw();

        f = xframe.findObject(pdfMmumu.GetName());

        # save the fit to the data/MC via the TCanvas
        if(save):
            c1.SaveAs('img/'+c1.GetName()+'.png')
            c1.SaveAs('rootfiles/'+histo.GetName()+"_"+pdfMmumu.GetName()+'.root');

        return f;
