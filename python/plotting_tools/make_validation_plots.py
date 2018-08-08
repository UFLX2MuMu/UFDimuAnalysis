##############################################
# make_validation_plots.py                   #
##############################################

import tools as tools
import sys
from ROOT import *
from categories_and_samples import *
import re
import argparse

gROOT.SetBatch()

parser = argparse.ArgumentParser(description='Parse cmd line arguments.')
parser.add_argument('f')
opts = parser.parse_args()
filename = opts.f

#str.find()
#m = re.search("\d", s1)
#if m:
#    print "Digit found at position %d" % m.start()
# i = [x.isdigit() for x in sequence].index(True)

# [(i,c) for i,c in enumerate('xdtwkeltjwlkejt7wthwk89lk') if c.isdigit()]
# {c:i for i,c in enumerate('xdtwkeltjwlkejt7wthwk89lk') if c.isdigit()}

#haystack = "xdtwkeltjwlkejt7wthwk89lk"
#digits   = "012345689"
#found    = [haystack.index(dig) for dig in digits if dig in haystack]
#firstdig = min(found) if found else None

#filename[validate.start()+len('validate_'):uscore.start()]

#============================================
# code
#============================================

categories = ['cAll', 'c14']

tfile = 0
print "\n /// Making validation plots for %s ... \n" % filename

# parse the file to get necessary naming information

blinded = re.search("blinded_", filename)
validate = re.search("validate_", filename)
dashnum = re.search("_-\d", filename)
num = re.search("_\d", filename)

if blinded and dashnum: 
    varname = filename[blinded.start()+len('blinded_'):dashnum.start()]

elif blinded and num: 
    varname = filename[blinded.start()+len('blinded_'):num.start()]

elif validate and dashnum: 
    varname = filename[validate.start()+len('validate_'):dashnum.start()]

elif validate and num: 
    varname = filename[validate.start()+len('validate_'):num.start()]

vartitle = varname.replace('_',' ')
vartitle = vartitle.replace('pt','p_{t}')
vartitle = vartitle.replace('dPhi','#Delta#phi')
vartitle = vartitle.replace('dEta','#Delta#eta')

vartitle = vartitle.replace('mass Roch','m')
vartitle = vartitle.replace('mass PF','m')
vartitle = vartitle.replace('mass KaMu','m')
vartitle = vartitle.replace('mass','m')


if 'dijet' not in varname:
    if 'jet1' in varname:
        vartitle = vartitle.replace('jet1', '')
        vartitle += '(j_{1})'
    if 'jet2' in varname:
        vartitle = vartitle.replace('jet2', '')
        vartitle += '(j_{2})'

if 'dimu' in vartitle:
    vartitle = vartitle.replace('dimu','')
    vartitle += '(#mu^{+} #mu^{-})'

if 'mu1' in vartitle:
    vartitle = vartitle.replace('mu1','')
    vartitle += '(#mu_{1})'

if 'mu2' in vartitle:
    vartitle = vartitle.replace('mu2','')
    vartitle += '(#mu_{2})'

if 'dijet1' in vartitle:
    vartitle = vartitle.replace('dijet1','')
    vartitle += '(jj_{1})'

if 'dijet2' in vartitle:
    vartitle = vartitle.replace('dijet2','')
    vartitle += '(jj_{2})'

if 'dEta' not in varname: vartitle = vartitle.replace('eta','#eta')
if 'dPhi' not in varname: vartitle = vartitle.replace('phi','#phi')

if 'abs' in vartitle: 
    vartitle = vartitle.replace('abs ', '|')
    vartitle = vartitle+'|'

if "mass" in varname or "pt" in varname: 
    units = ' (GeV)'

else: units = ''

if 'bdt' in vartitle:
    vartitle.replace('bdt', 'BDT')

# Drell Yan Naming bug ....
dy_bug = False

tfile = TFile(filename)
tfile_dir = tfile.GetDirectory('net_histos')
for key in tfile_dir.GetListOfKeys():
    name = key.GetName()
    #print "(%d) Looping over %s ..." % (dy_bug, key)
    if 'Drell_Yan_' in name and name[-1] == '_':
        dy_bug = True
        break

# Make the stacks
for c in categories:
    in_list = []
    hlist = []
    in_unc_map = {'JES_up': [], 'JES_down': [], 'PU_up': [], 'PU_down' : []}
    unc_hist_map = {}

    for s in net_samples:
        if s == 'Drell_Yan' and dy_bug:
            s+='_'
        #print "   /// (%d) Looking at sample %s ..." % (dy_bug, s)
        # Histos to get for this category
        in_list.append( (filename, 'net_histos/%s_%s' % (c, s)) )

        # Up and Down uncertainty histograms to get for this category, e.g. JES_up, JES_down, ...
        if s == 'Net_Data': 
            continue       # Data does not have up/down uncertainty histograms, only MC

        for key,value in in_unc_map.iteritems():
            value.append( (filename, 'net_histos/%s_%s_%s' % (c, s, key)) ) 
    
    # list of histograms to make the stack and ratio plots for this category
    hlist = tools.get(in_list)

    # map to lists of up/down uncertainty histograms for this category
    # 'JES_up' -> list of JES_up histos for this category
    for in_key, in_value in in_unc_map.iteritems():
        unc_hist_map[in_key] = tools.get(in_value) 

    #hmc = tools.add(hlist[0:-1])
    #hdata = hlist[-1]
    #newBins = tools.getRebinEdges(hdata, hmc, max_err=0.05)
    #print newBins
    #rebin_hlist = tools.rebin(hlist, newBins) #rebinned to var bin width
    #srebin_hlist = tools.scaleByBinWidth(rebin_hlist, newBins) # rebinned to var bin width and wide bins are scaled down
    
    print "\n =========== STACK AND RATIO ====== \n" 
    
    tools.stackAndRatio(hlist, unc_map=unc_hist_map, savename=varname+'_%s' % c, title=('%s '+ vartitle) % c, ldim=[0.6, 0.55, 0.92, 0.92], xtitle=vartitle+units)
    #tools.stackAndRatio(srebin_hlist, savename=varname+'_%s' % c, title=('%s '+ vartitle) % c, ldim=[0.6, 0.55, 0.92, 0.92], xtitle=vartitle+units)
    #tdrstyle.cmsPrel(35900, 13, False, onLeft=True)
    c = TCanvas();
