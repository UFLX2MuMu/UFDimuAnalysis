##############################################
# make_tables.py                             #
##############################################

import tools as tools
import sys
from ROOT import *
from categories_and_samples import *
import re
import argparse
import numpy as np

gROOT.SetBatch()

parser = argparse.ArgumentParser(description='Parse cmd line arguments.')
parser.add_argument('f')
opts = parser.parse_args()
filename = opts.f

#[function(number) for number in numbers if condition(number)]

#============================================
# code
#============================================

print "\n /// Overlaying backgrounds for %s... \n" % filename

# List of histograms to make the overlay from -----------------
in_list = []
for c in bdt_categories:
    for s in net_samples_all:
        if s == 'Net_Data':
            in_list.append((filename, 'net_histos/%s_%s' % (c, s))),

print "\n =========== GET ================== \n"
hlist = tools.get(in_list)
edges = [x for x in range(110,161) if x%2 == 0]
for e in edges:
    print e 
hlist = tools.rebin(hlist, edges)

inclusive = 0

# Normalize all of the backgrounds to the same scale
for i,h in enumerate(hlist):
    h.Scale(1/h.Integral())
    #h.SetName(bdt_categories[i])
    h.SetTitle(bdt_categories[i])
    if 'cAll' in h.GetName():
        h.SetTitle("Inclusive")
        inclusive = h

#print "\n =========== OVERLAY ================== \n"
#tools.overlay(hlist, title="Background Shapes", savename="compare_bkg_shapes", xtitle="dimu_mass", ytitle="Events / 0.2 GeV", 
#              ldim=[0.85, 0.4, 0.98, 0.98], draw='fill-transparent')

print "\n =========== STACKS ================== \n"
list_to_stack = []

for h in hlist:
    if 'cAll' in h.GetName():
        continue

    l = [h,inclusive]
    e = tools.getRebinEdges(h, inclusive, max_err=0.05)
    lrebin = tools.scaleByBinWidth(tools.rebin(l, e), e)
    tools.stackAndRatio(lrebin, log=False, name=h.GetTitle(), savename="bkg_shapes_stack_%s" % h.GetTitle(), title=h.GetTitle(), 
                        ytitleratio="%s/%s" % ("Inclusive",h.GetTitle()), errband="Error Band", ldim=[0.6, 0.55, 0.92, 0.92])

#for l in list_to_stack:
#    h = l[0]
#    tools.stackAndRatio(l, log=False, name=h.GetName(), savename="bkg_shapes_stack_%s" % h.GetName(), title=h.GetName(), 
#                        ytitleratio="%s/%s" % (h.GetName(), 'inclusive'))
