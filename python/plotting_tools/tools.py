##############################################
# tools.py                                   #
##############################################

#============================================
# import
#============================================

import string
import re
import argparse
import numpy as np
from ROOT import *
from array import *

import sys
import tdrstyle

#sys.argv.append( '-b-' )

#============================================
# code
#============================================

#colors = [58, 2, 8, 36, 91, 46, 50, 30, 9, 29, 3, 42, 98, 62, 74, 20, 29, 32, 49, 12, 3, 91]
#colors = [58, 2, 8, 30, 91, 42, 50, 30, 9, 29, 3, 42, 98, 62, 74, 20, 29, 32, 49, 12, 3, 91]
#colors = [225, 205, 12, 8, 220, 29, 30, 35, 28, 18]
#colors = [225, 205, 8, 29, 12, 220]
colors = [225, 8, 205, 29, 12, 220]

def overlay(plot_list, title="overlay", savename="", name="", xtitle="", ytitle="", yrange=(-999,-999), ldim=[0.71,0.82,1.0,1.0], 
            draw="fill-transparent", autocolor=True):
    """ Overlay histograms or tgraphs onto the same canvas. plot_list = [hist/graph, hist/graph, hist/graph, ...]"""

    if savename == "": savename = title
    legend = TLegend(ldim[0], ldim[1], ldim[2], ldim[3], "", "brNDC")
    c = TCanvas()
    for i,plot in enumerate(plot_list):
        if 'fill' in draw: 
            if 'transparent' in draw: 
                plot.SetFillStyle(3001)
                #plot.SetFillColorAlpha(colors[i], 0.55)
            else: 
                plot.SetFillColor(colors[i])

        if 'point' in draw:
            plot.SetMarkerSize(1)
        else:
            if(autocolor and plot.InheritsFrom("TH1")): 
                plot.SetLineColor(1)
                plot.SetLineWidth(1)
            elif(autocolor): 
                plot.SetLineColor(colors[0])
                plot.SetLineWidth(2)
            plot.SetMarkerSize(2)

        if(i==0):
            if plot.InheritsFrom("TH1"): 
                plot.SetStats(0)
                if 'transparent' in draw: 
                    plot.Draw("hist")
                elif 'point' in draw: 
                    plot.Draw("P")
                else: 
                    plot.Draw()
            else: 
                plot.Draw()

            if yrange != (-999,-999): 
                plot.GetYaxis().SetRangeUser(yrange[0], yrange[1])

            if xtitle!="": 
                plot.GetXaxis().SetTitle(xtitle)

            if ytitle!="": 
                plot.GetYaxis().SetTitle(ytitle)

            c.SetTitle(title)
            c.SetName(title)
            if name!="": c.SetName(name)

        else:
            if plot.InheritsFrom("TH1") and 'transparent' in draw: 
                plot.Draw("hist SAME")

            elif plot.InheritsFrom("TH1") and 'point' in draw: 
                plot.Draw("P SAME")

            else: 
                plot.Draw("SAME")

        if 'fill' in draw: 
            legend.AddEntry(plot, plot.GetTitle(), "f")
        elif 'point' not in draw: 
            legend.AddEntry(plot, plot.GetTitle(), "l")
        elif 'point' in draw: 
            legend.AddEntry(plot, plot.GetTitle(), "p")

        legend.Draw("SAME")
        plot_list[0].SetTitle(title)

    legend.SetFillStyle(0)
    c.cd()
    c.SetGridx(0)
    c.SetGridy(0)
    tdrstyle.cmsPrel(35900, 13, False, prelim=True, onLeft=True, textScale=1.5)
    plot_list[0].GetYaxis().SetTitleOffset(1.55)
    c.SaveAs('../out/imgs/%s.png' % savename)
    c.SaveAs('../out/rootfiles/%s.root' % savename)

##################################################################
#-----------------------------------------------------------------
##################################################################

def stack(stack_list, unc_map={}, title="stack", name="stack", xtitle="Mass (GeV)", ytitle="", yrange=(-999,-999), ldim=[0.71,0.82,1.0,1.0], autocolor=True):
    """ Stack some histograms to compare to another histogram. Take the ratio chist/stack and put it in the pad in the bottom. """
    stack = THStack() 
    stack.SetTitle(title)
    stack.SetName(name)
    legend = TLegend(ldim[0], ldim[1], ldim[2], ldim[3], "", "brNDC")

    # keep track of these to auto zoom the y-axis
    minmax = 999999999

    for i,hist in enumerate(stack_list):
        #print "%s: %f" % (hist.GetName(), hist.Integral())
        #  set the colors, etc
        if autocolor:
            hist.SetLineColor(1)
            hist.SetFillColor(colors[i])
        hist.SetStats(0)

        if(i==0):
            minmax = hist.GetMaximum();
        if(hist.GetMaximum() < minmax):  minmax = hist.GetMaximum();   # the lowest maximum of the histograms

        stack.Add(hist) 
        legend.AddEntry(hist, hist.GetTitle(), "f")

    #print  "stack integral: %f" % stack.GetStack().Last().Integral()
    return stack, legend, minmax;

##################################################################
#-----------------------------------------------------------------
##################################################################

def add(histo_list, name="Net_MC", title="Net MC"):
    """ Sum the histos given in the list of (file, name) pairs. """
    hsum = 0
    for i,hist in enumerate(histo_list):
        #print "%s: %f" % (hist.GetName(), hist.Integral())
        if(i==0):
            gROOT.cd()
            hist.Clone()
            hsum = hist.Clone(name) 
            hsum.SetTitle(title)
        
        else:
            hsum.Add(hist)

    #print "%s: %f" % (hsum.GetName(), hsum.Integral())
    return hsum

##################################################################
#-----------------------------------------------------------------
##################################################################

def get(item_list):
    """ Return a list of ROOT objects given a list of (file, name) pairs. """
    rlist = []
    for i,item_info in enumerate(item_list):
        # first get the histograms from the rootfiles
        filename, itemname = item_info
        print "&&& Grabbing %s, %s" % (filename, itemname)
        tfile = TFile(filename)
        gROOT.cd() # AHHHHH
        item = tfile.Get(itemname).Clone()
        #print item.GetName()
        rlist.append(item)
    return rlist

##################################################################
#-----------------------------------------------------------------
##################################################################

def ratio(histo_list, name="data_mc_ratio", title="data_mc_ratio"):
    """ Make a ratio histogram by dividing the Nth histogram by the sum of the 0 to N-1 histograms. """
    hsum = add(histo_list[0:-1])
    hlast = histo_list[-1].Clone(name)
    #print "%s: %f, %s: %f, ratio: %f" % (hlast.GetName(), hlast.Integral(), hsum.GetName(), hsum.Integral(), hlast.Integral()/hsum.Integral())
    hlast.Divide(hsum)
    return hlast

##################################################################
#-----------------------------------------------------------------
##################################################################

def subtractHistos(histo, up_or_down):
    """ Get the difference between two histograms. """
    diff = up_or_down.Clone(up_or_down.GetName()+"_minus_"+histo.GetName()) 
    diff.Add(histo, -1)
    print "%s: %f, %s: %f, diff: %f" % (up_or_down.GetName(), up_or_down.Integral(), histo.GetName(), histo.Integral(), up_or_down.Integral()-histo.Integral())
    return diff

##################################################################
#-----------------------------------------------------------------
##################################################################

def errorBandHisto(histo, up_or_down_list):
    """ Calculate the error band histograms given a histogram and a list of systematics (up/down histos). Assumes the list is ordered
        like [systematic1_up, systematic1_down, systematic2_up, systematic2_down, ...] or [systematic1_down, systematic1_up, 
        systematic2_down, systematic2_up, ...]"""

    err_up_hist = histo.Clone(histo.GetName()+"_err_up")
    err_down_hist = histo.Clone(histo.GetName()+"_err_down")

    # Clear the bin and the bin error for the up/down error histos
    # then set the inior_tial up/down error in each bin to the nominal error
    for b in range(1,histo.GetNbinsX()+1): 
        err_up_hist.SetBinContent(b, 0)
        err_up_hist.SetBinError(b, 0)
        err_up_hist.SetBinContent(b, histo.GetBinError(b)*histo.GetBinError(b))
        err_down_hist.SetBinContent(b, 0)
        err_down_hist.SetBinError(b, 0)
        err_down_hist.SetBinContent(b, histo.GetBinError(b)*histo.GetBinError(b))

    # Get the bin by bin differences
    for i in range(0, len(up_or_down_list)-1, 2):
        h1 = up_or_down_list[i] 
        h2 = up_or_down_list[i+1] 
        diff1 = subtractHistos(histo, h1) 
        diff2 = subtractHistos(histo, h2)
        for b in range(1,histo.GetNbinsX()+1): 
            diff1b = diff1.GetBinContent(b)
            diff2b = diff2.GetBinContent(b)
            maxdiff = max(diff1b, diff2b)
            mindiff = max(diff1b, diff2b)
            err_up_hist.SetBinContent(b, err_up_hist.GetBinContent(b)+maxdiff*maxdiff)
            err_down_hist.SetBinContent(b, err_down_hist.GetBinContent(b)+mindiff*mindiff)

    # Add +/- net_err = sqrt(err0^2 + err1^2 + ...) to the nominal value to get the upper and lower histos    
    up_hist = histo.Clone(histo.GetName()+"_up")
    down_hist = histo.Clone(histo.GetName()+"_down")
    band_hist = histo.Clone(histo.GetName()+"_err_band")

    for b in range(1,histo.GetNbinsX()+1): 
        up_hist.SetBinContent(b, up_hist.GetBinContent(b) + TMath.Sqrt(err_up_hist.GetBinContent(b)))
        down_hist.SetBinContent(b, down_hist.GetBinContent(b) - TMath.Sqrt(err_down_hist.GetBinContent(b)))
        band_hist.SetBinContent(b, (up_hist.GetBinContent(b) + down_hist.GetBinContent(b))/2)
        band_hist.SetBinError(b, (up_hist.GetBinContent(b) - down_hist.GetBinContent(b))/2)

    print "\n up_hist - down_hist = %f \n" % (up_hist.Integral() - down_hist.Integral())
    return band_hist, up_hist, down_hist

##################################################################
#-----------------------------------------------------------------
##################################################################

# Stack the histograms on the top part of the canvas and put a ratio plot on the bottom part of the canvas
def stackAndRatio(histo_list, rlist=0, unc_map={}, title="stack", log=True, name="stack", xtitle="Mass (GeV)", ytitle="", yrange=(-999,-999), ytitleratio="Data/MC", ldim=[0.71,0.82,1.0,1.0], yrange_ratio=(-999,-999), csize=(-999,-999), savename="", errband="MC Error Band"):

    if name=="stack" and title!=stack:
        name=title

    if(rlist==0):
        rlist = histo_list

    for h in histo_list:
        h.SetStats(False)

    for h in rlist:
        h.SetStats(False)

    if csize == (-999,-999):
        csize = (800,800)

    if savename == "":
        savename=name

    print "\n --- Nominal histograms --- "
    for h in histo_list:
        print h.GetName()

    print "\n --- Systematics histograms --- "
    for key,hlist in unc_map.iteritems():
        print "\n" + key
        for h in hlist:
            print h.GetName()

    print ""

    net_mc_hist = add(histo_list[:-1], "NET_MC_X") # the nominal net mc

    # For each systematic (JES_up, JES_down, PU_up, PU_down)
    # add up the inidivdual MC histos into a net histogram
    net_mc_unc_list = []
    net_mc_ratio_unc_list = []
    for key, hlist in unc_map.iteritems():
        net_mc_unc_list.append(add(hlist, "Net_MC_"+key))
        hlist.append(net_mc_hist)
        net_mc_ratio_unc_list.append(ratio(hlist, name="data_mc_ratio_"+key))

    # Get the error band histo given the individual up/down contributions in the net_mc_unc_list
    band_hist = 0
    up_hist = 0
    down_hist = 0

    # get the error band histo that accounts for the systematics
    if len(net_mc_unc_list) != 0: 
        band_hist, up_hist, down_hist = errorBandHisto(net_mc_hist, net_mc_unc_list)

    c = TCanvas()
    c.SetCanvasSize(csize[0],csize[1]);
    pad1 = TPad("pad1", "pad1", 0,0.3,1,1.0)
    pad1.SetBottomMargin(0.0125)
    pad1.Draw()
    pad1.cd()

    s, legend, minmax = stack(histo_list[0:-1], title=title, name=name, xtitle=xtitle, ytitle=ytitle, yrange=yrange, ldim=ldim)
    
    s.Draw("hist")
    s.GetXaxis().SetTitle("")
    s.GetYaxis().SetTitle(ytitle)
    s.GetYaxis().SetLabelFont(43); 
    s.GetYaxis().SetLabelSize(15);
    s.GetXaxis().SetLabelSize(0.000001);
    
    # overlay the mc error band
    bandhist = s.GetStack().Last().Clone() # set the default error band to the statistical uncertainty of the Net MC

    # if there are systematics provided then use that error band histo
    if band_hist != 0:
        bandhist = band_hist
        up_hist.SetLineColor(1)
        down_hist.SetLineColor(1)

    # set the style for the error band histogram
    bandhist.SetFillStyle(3001)
    bandhist.SetLineColor(12)
    bandhist.SetFillColor(12)
    bandhist.SetMarkerStyle(0)
    bandhist.Draw("E2SAME")
    legend.AddEntry(bandhist, errband, "f")
    gPad.Modified()
    
    if(yrange==(-999,-999) and log): yrange = (minmax*0.1, s.GetMaximum()*10e5)
    elif(yrange==(-999,-999)): yrange = (0, s.GetMaximum()*1.5)
    if(log and yrange[0] > 0): 
        gPad.SetLogy(1)
    s.SetMaximum(yrange[1])
    s.SetMinimum(yrange[0])
    
    hlast = histo_list[-1]
    hlast.SetMarkerStyle(20)
    hlast.SetLineColor(1)
    hlast.SetFillColor(0)
    legend.AddEntry(hlast, hlast.GetTitle(), "p")
    hlast.Draw("epsames")
    legend.Draw("SAME");

    # change to bottom pad and draw the ratio plot
    c.cd()
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetGridy()
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.25)
    pad2.Draw()
    pad2.cd()

    hratio = ratio(rlist) # the nominal ratio plot usually of data/mc
 
    if yrange_ratio[0] == -999: 
        hratio.SetMinimum(0.58)
    else: 
        hratio.SetMinimum(yrange_ratio[0])

    if yrange_ratio[1] == -999: 
        hratio.SetMaximum(1.42)
    else: 
        hratio.SetMaximum(yrange_ratio[1])

    # The histograms for the error band with systematics on the ratio
    rband_hist = 0
    rup_hist = 0
    rdown_hist = 0

    # Create the error band from the systematics if we have the histograms available
    if len(net_mc_ratio_unc_list) != 0:
        print "\n --- Ratio plots for error band ---"
        for h in net_mc_ratio_unc_list:
            print h.GetName()
        print ""

        rband_hist, rup_hist, rdown_hist = errorBandHisto(ratio([net_mc_hist, net_mc_hist], "net_mc_to_net_mc_ratio"), net_mc_ratio_unc_list)

    rbandhist = hratio # set the default error band to the statistical uncertainty
    if rband_hist != 0:
        rbandhist = rband_hist

    # set the style for the ratio error band
    rbandhist.SetFillStyle(3001)
    rbandhist.SetLineColor(12)
    rbandhist.SetFillColor(12)
    rbandhist.SetMarkerStyle(0)
    hratio.SetMarkerStyle(20)
    hratio.Draw("ep")
    rbandhist.Draw("E2SAME")
    #rup_hist.Draw("lsame")
    #rdown_hist.Draw("lsame")

    hratio.SetTitle("")
    hratio.GetYaxis().SetTitle(ytitleratio);
    hratio.GetYaxis().SetNdivisions(505);
    hratio.GetYaxis().SetTitleSize(25); # 20
    hratio.GetYaxis().SetTitleFont(43);
    hratio.GetYaxis().SetTitleOffset(0.6*1.55);
    hratio.GetYaxis().SetLabelFont(43); 
    hratio.GetYaxis().SetLabelSize(15);

    hratio.GetXaxis().SetTitle(xtitle);
    hratio.GetXaxis().SetTitleSize(25); # 20
    hratio.GetXaxis().SetTitleFont(43);
    hratio.GetXaxis().SetTitleOffset(0.8*4.);
    hratio.GetXaxis().SetLabelFont(43); 
    hratio.GetXaxis().SetLabelSize(15);

    pad1.cd()
    #tdrstyle.cmsPrel(35900, 13, True, onLeft=True, textScale=1.5)
    tdrstyle.cmsPrel(35900, 13, False, prelim=True, onLeft=True, textScale=1.5)

    c.SaveAs('../out/imgs/'+savename+".png")
    c.SaveAs('../out/rootfiles/'+savename+".root")

##################################################################
#-----------------------------------------------------------------
##################################################################

# Figure out the error squared for the ratio plot based upon the amount in the numerator and denominator
def ratioError2(numerator, numeratorError2, denominator, denominatorError2):
    """ Propogation of errors for N_numerator/N_denominator
        Returns error squared, takes in numerator error squared and denom err squared as well as N_num and N_denom"""
    if(numerator   <= 0): return 0
    if(denominator <= 0): return 0
    ratio = numerator/denominator;

    return  ratio*ratio*( numeratorError2/numerator/numerator + denominatorError2/denominator/denominator );

##################################################################
#-----------------------------------------------------------------
##################################################################

# Figure out the rebinning scheme based upon the max allowed error in each bin in the ratio plot
def getRebinEdges(numerator, denominator, max_err = 0.10):
    """The ratio plots are a bit crazy with huge errors sometimes, so we want to rebin with variable binning
       such that the error is always low in each of the ratio plot bins"""

    newBins = []    # will return this as the new binning scheme based upon the ratio plot errors
    errVec = []     # error in each bin
    numName = numerator.GetName()
    isMassData = False
    if("dimu_mass" in numerator.GetXaxis().GetTitle() and "Data" in numName): isMassData = True
    #print "numerator name: %s, numerator xtitle: %s, isMassData: %d" %(numName, numerator.GetXaxis().GetTitle(), isMassData)

    # check if the two histos are binned the same way
    compatible = numerator.GetNbinsX() == denominator.GetNbinsX() and numerator.GetBinLowEdge(1) == denominator.GetBinLowEdge(1) \
                 and numerator.GetBinLowEdge(numerator.GetNbinsX()) == denominator.GetBinLowEdge(denominator.GetNbinsX())
    if(not compatible):
        print "!!! ERROR: tools.getBinning >> numerator and denominator histograms do not have the same binning scheme!"
        return []

    if(numerator.Integral() <= 0 or denominator.Integral() <=0):
        print "!!! ERROR: tools.getBinning >>  numerator or denominator integral <= 0"
        return []

    # first low edge is the lowest edge by default
    newBins.append(numerator.GetBinLowEdge(1))     # edges for the new binning scheme
    sumNum = 0                                     
    sumErr2Num = 0
    sumDenom = 0
    sumErr2Denom = 0
    lastBinIsCutBin = 0

    # Collect bins together until the error is low enough for the corresponding ratio plot bin.
    # Once the error is low enough call that the new bin, move on and repeat.
    for i in range(1,numerator.GetNbinsX()+1):
        sumNum += numerator.GetBinContent(i)
        sumErr2Num += numerator.GetBinError(i)*numerator.GetBinError(i)

        sumDenom += denominator.GetBinContent(i)
        sumErr2Denom += denominator.GetBinError(i)*denominator.GetBinError(i)

        # squared error for the ratio, using propogation of errors
        # 0 if sumNum or sumDenom is 0
        ratioErr2 = ratioError2(sumNum, sumErr2Num, sumDenom, sumErr2Denom)

        percentError = 0
        if(ratioErr2>0): 
            percentError = TMath.Sqrt(ratioErr2)/(sumNum/sumDenom)

        # If the bins suddenly drop to zero or go up from zero then this is probably due to some cut that was applied
        # and we dont' want to combine the 0 bins on one side of the cut with the bins after the cut. It messes up the comparison plots.
        # it's better to see exactly where the cut caused everything to drop to zero.
        if(i!=numerator.GetNbinsX() and TMath.Abs(numerator.GetBinContent(i)-numerator.GetBinContent(i+1)) > 10
           and (numerator.GetBinContent(i) == 0 or numerator.GetBinContent(i+1) == 0)):
            newBins.append(numerator.GetBinLowEdge(i)+numerator.GetBinWidth(i))
            errVec.append(percentError)

            sumNum = 0
            sumErr2Num = 0

            sumDenom = 0
            sumErr2Denom = 0
            lastBinIsCutBin = True
        # we have the minimum error necessary create the bin in the vector
        # or just make a bin if we are in the blinded region of the mass spectrum
        elif(percentError !=0 and percentError <= max_err or (isMassData and numerator.GetBinLowEdge(i) >= 120 
             and numerator.GetBinLowEdge(i) < 130)):
            newBins.append(numerator.GetBinLowEdge(i)+numerator.GetBinWidth(i))
            errVec.append(percentError)

            sumNum = 0
            sumErr2Num = 0

            sumDenom = 0
            sumErr2Denom = 0

            lastBinIsCutBin = False
        # we got to the end of the histogram and the last bins can't add up
        # to get the error low enough, so we merge these last bins with the 
        # last bin set up in the new variable binning scheme.
        if(i==numerator.GetNbinsX() and (percentError > max_err or percentError == 0)):
            # get rid of the last bin added, but make sure there is at least one bin. Don't remove the low edge of the zero-th bin.
            # la la we need at least a value for the beginning and end, size 1 vector doesn't make sense and rebinning will fail.
            if(len(newBins) > 1 and not lastBinIsCutBin): 
                del newBins[-1]
            # replace with the end bin value in the numerator histo
            newBins.append(numerator.GetBinLowEdge(i)+numerator.GetBinWidth(i))
    return newBins

##################################################################
#-----------------------------------------------------------------
##################################################################

# Rebin histograms based upon the list of edges
def rebin(hlist, edges):
    """Use the a variable binning scheme to rebin the histograms in the list."""
    rebin_hlist = []
    for h in hlist:
        hrebin = h.Rebin(len(edges)-1, h.GetName()+"_", array('d',edges))
        rebin_hlist.append(hrebin)
        
    return rebin_hlist

##################################################################
#-----------------------------------------------------------------
##################################################################

def scaleByBinWidth(hlist, edges):
    """ Scale each bin by its width. The smallest width is the reference. For example, a bin 3x the smallest will be scaled by 1/3."""
    min_diff = 9999999999
    scale = []
    retlist = []

    # find minimum bin width
    for i in range(len(edges)-1):
        scale.append(edges[i+1]-edges[i])   # bin width for each bin
        if(edges[i+1] - edges[i] < min_diff): min_diff = edges[i+1]-edges[i]

    scale = np.array(scale)
    scale = scale/min_diff

    for h in hlist:
        hclone = h.Clone(h.GetName()+"_")
        # scale each bin in hrebin so that it is reduced by binwidth/minwidth
        for i,b in enumerate(scale):
            hclone.SetBinContent(i+1, 1/b*hclone.GetBinContent(i+1))
            hclone.SetBinError(i+1, 1/b*hclone.GetBinError(i+1))
        retlist.append(hclone)

    return retlist

##################################################################
#-----------------------------------------------------------------
##################################################################

def getMinMaxY(plot):
    """ Get the min and max of a plot. """
    miny = -999
    maxy = -999

    if plot.InheritsFrom("TH1"):
        miny = plot.GetMinimum()
        maxy = plot.GetMaximum()
    else:
        yvals = plot.GetY()
        miny = min(yvals)
        maxy = max(yvals)

    return (miny, maxy)

##################################################################
#-----------------------------------------------------------------
##################################################################

# Recently implemented smoothing by averaging w/ neighbors to get better results.
# Haven't tested this recent change to make sure it's not buggy.
def fwhm(h, start=1):
    """ Calculate the FWHM of a peaked histogram by scanning the bin contents. Not stable for low stats. """
    maxbin = -999
    leftbin = -999
    rightbin = -999
    max = -999
    half_max = -999 

    # find max
    for i in range(start, h.GetNbinsX()):
        val = 0
        if i-1 >= start and i+1 < h.GetNbinsX(): 
            val = (h.GetBinContent(i-1) + h.GetBinContent(i) + h.GetBinContent(i+1))/3

        elif i+1 < h.GetNbinsX(): 
            val = (h.GetBinContent(i) + h.GetBinContent(i+1))/2

        elif i-1 >= start:
            val = (h.GetBinContent(i) + h.GetBinContent(i+1))/2

        else:
            val = h.GetBinContent(i)

        if val > max:
            max = val
            maxbin = i

    half_max = max/2

    # find half-max bin on left of max
    for i in range(start, maxbin):
        val = 0
        if i-1 >= start and i+1 < h.GetNbinsX(): 
            val = (h.GetBinContent(i-1) + h.GetBinContent(i) + h.GetBinContent(i+1))/3

        elif i+1 < h.GetNbinsX(): 
            val = (h.GetBinContent(i) + h.GetBinContent(i+1))/2

        elif i-1 >= start:
            val = (h.GetBinContent(i) + h.GetBinContent(i+1))/2

        else:
            val = h.GetBinContent(i)

        if val >= half_max:
            leftbin = i
            break

    # find half-max bin on right of max
    for i in range(maxbin, h.GetNbinsX()):
        val = 0
        if i-1 >= start and i+1 < h.GetNbinsX(): 
            val = (h.GetBinContent(i-1) + h.GetBinContent(i) + h.GetBinContent(i+1))/3

        elif i+1 < h.GetNbinsX(): 
            val = (h.GetBinContent(i) + h.GetBinContent(i+1))/2

        elif i-1 >= start:
            val = (h.GetBinContent(i) + h.GetBinContent(i+1))/2

        else:
            val = h.GetBinContent(i)

        if val <= half_max:
            rightbin = i
            break

    left = h.GetBinCenter(leftbin);
    right = h.GetBinCenter(rightbin);
    full_width_half_max = right - left;

    return leftbin, rightbin, full_width_half_max

##################################################################
#-----------------------------------------------------------------
##################################################################

# Fit the signal with a gaussian and get the FWHM boundaries
#def fwhm_fit(h, fitmin=121, fitmax=129):
def fwhm_fit(h, fitmin=115, fitmax=135, double=True):
    """ Fit the signal with a Gaussian and get the FWHM and left and right bins. """
    from math import sqrt
    ntries = 0
    converged = False
    fit = 0

    if not double: 
        fit = TF1("Gaussian", "[0]*TMath::Gaus(x, [1], [2])", fitmin, fitmax) 
    else: 
        fit = TF1("Gaussian", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5])", fitmin, fitmax) 

    fit.SetParNames("Normalization1", "Mean1", "Sigma1", "Normalization2", "Mean2", "Sigma2")
    fit.SetParameters(0.7*h.GetMaximum(), h.GetMean(), 1.2, 0.3*h.GetMaximum(), h.GetMean(), 4);

    for i in range(0, 50):
        h.Fit("Gaussian", "q", "", fitmin, fitmax)
        sconverge = gMinuit.fCstatu
        if "CONVERGED" in sconverge:
            converge = True
            break

    hfit = fit.GetHistogram()

    norm1 = fit.GetParameter(0);
    gauss_mean1 = fit.GetParameter(1)
    gauss_sigma1 = fit.GetParameter(2)

    norm2 = fit.GetParameter(3);
    gauss_mean2 = fit.GetParameter(4)
    gauss_sigma2 = fit.GetParameter(5)

    gauss_mean = gauss_mean1
    sigma = gauss_sigma1

    full_width_half_max = sigma*2.355
    left = gauss_mean - full_width_half_max/2;
    right = gauss_mean + full_width_half_max/2;
    leftbin1 = h.GetXaxis().FindBin(left)
    rightbin1 = h.GetXaxis().FindBin(right)

    leftbin = hfit.FindFirstBinAbove(hfit.GetMaximum()/2);
    rightbin = hfit.FindLastBinAbove(hfit.GetMaximum()/2);
    left = hfit.GetBinCenter(leftbin)
    right = hfit.GetBinCenter(rightbin)
    
    leftbin = h.GetXaxis().FindBin(left);
    rightbin = h.GetXaxis().FindBin(right);
    full_width_half_max = right - left;

    print "\ndLeft %d, dRight %d" % (leftbin1-leftbin, rightbin1-rightbin)
    print "NBINS : %d vs %d \n" % (hfit.GetNbinsX(), h.GetNbinsX())
    return leftbin, rightbin, full_width_half_max

##################################################################
#-----------------------------------------------------------------
##################################################################

# Fit the bkg with bwz redux and return the fit 
def fit_bkg(h, fitmin=110, fitmax=160):
    """Fit the background with BWZRedux and return the fit."""
    ntries = 0
    converged = False
    fit = TF1("bwz_redux", "[0]*TMath::Exp([1]*x/100 + [2]*(x/100)^2)/(TMath::Power(x-91.2,[3]) + TMath::Power(2.5/2, [3]))", fitmin, fitmax) 
    fit.SetParNames("Normalization", "a1", "a2", "a3")
    fit.SetParameters(h.GetMaximum(), 1.39, 0.46, -0.26);

    for i in range(0, 50):
        h.Fit("bwz_redux", "q", "", fitmin, fitmax)
        sconverge = gMinuit.fCstatu
        if "CONVERGED" in sconverge:
            converge = True
            break

    a1 = fit.GetParameter(1)
    a2 = fit.GetParameter(2)
    a3 = fit.GetParameter(3)
    return fit
