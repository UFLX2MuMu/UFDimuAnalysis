##############################################
# make_signal_validation_plots.py            #
##############################################

import tools as tools
import sys
from ROOT import *
from categories_and_samples import *
import re
import argparse
import tdrstyle 

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
    for s in net_samples:
        if 'H2Mu' not in s and 'VH' not in s:
            continue
        if s == 'Drell_Yan' and dy_bug:
            s+='_'
        #print "   /// (%d) Looking at sample %s ..." % (dy_bug, s)
        in_list.append( (filename, 'net_histos/%s_%s' % (c, s)) )
    
    hlist = tools.get(in_list)
    
    #hmc = tools.add(hlist[0:-1])
    #hdata = hlist[-1]
    #newBins = tools.getRebinEdges(hdata, hmc, max_err=0.1)
    #print newBins
    #rebin_hlist = tools.rebin(hlist, newBins) #rebinned to var bin width
    #srebin_hlist = tools.scaleByBinWidth(rebin_hlist, newBins) # rebinned to var bin width and wide bins are scaled down
    
    print "\n =========== STACK AND RATIO ====== \n" 
    xtitle=vartitle+units
    canvas = TCanvas()
    savename='sig_kinematics_%s_%s' % (varname, c)
    stack, legend, minmax = tools.stack(hlist, title=('%s '+ vartitle) % c, ldim=[0.6, 0.65, 0.92, 0.92])
    stack_sum = stack.GetStack().Last()
    maximum = -999
    for i in range(2, stack_sum.GetNbinsX()):
        if stack_sum.GetBinContent(i) > maximum:
            maximum = stack_sum.GetBinContent(i)

    if 'n' != varname[0]: 
        stack.SetMaximum(maximum*1.30)
    stack.Draw("hist")
    if stack_sum.GetBinContent(1) > 5* stack_sum.GetBinContent(2) and 'n' != varname[0]:
        stack.GetXaxis().SetRangeUser(stack_sum.GetBinLowEdge(2), stack.GetXaxis().GetXmax());


    legend.Draw("SAME")
    stack.GetXaxis().SetTitle(xtitle)
    tdrstyle.cmsPrel(35900, 13, True, onLeft=True, textScale=1.5)
    canvas.SaveAs('../out/imgs/'+savename+'.png')
    canvas.SaveAs('../out/rootfiles/'+savename+'.root')
