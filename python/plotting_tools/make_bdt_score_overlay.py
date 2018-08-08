##############################################
# make_bdt_score_overlay.py                  #
##############################################

import tools as tools
from categories_and_samples import *
from ROOT import *

#============================================
# code
#============================================

filename = '../../UFDimuAnalysis_v2/bin/rootfiles/validate_bdt_score_-1_1_categories2_36814_dyAMC-J_minpt20_b1_sig-xlumi1.root'
in_list = []
c = 'cAll'

for s in net_samples:
    if 'Data' in s or 'Diboson' in s or 'VH' in s:
        continue
    in_list.append( (filename, 'net_histos/%s_%s' % (c, s)) )

plot_list = tools.get(in_list)

i=2
for h in plot_list:
    h.Scale(1/h.Integral())
    i+=1

tools.overlay(plot_list, title="BDT Score", savename="bdt_score_overlay_%s" % c, 
              draw="fill-transparent", xtitle="BDT score", ytitle="", ldim=[0.5, 0.7, 0.88, 0.88], yrange =(0,0.15))
