##############################################
# mass_dep.py                                #
##############################################

import tools as tools
import sys
sys.argv.append( '-b-' )
from ROOT import *

#============================================
# code
#============================================

dirname   = '../validation/'
filename = 'validate_bdt_score_-1_1_categories2_36814_dyAMC-J_minpt20_b0_sig-xlumi1.root'
category = 'cAll'
vbf_samples = ['H2Mu_VBF', 'H2Mu_VBF_120', 'H2Mu_VBF_130']
gg_samples = ['H2Mu_gg', 'H2Mu_gg_120', 'H2Mu_gg_130']

filename = dirname+filename
gg_getlist  = []
vbf_getlist = []

for sample in gg_samples:
    gg_getlist.append((filename, 'signal_histos/%s_%s' % (category, sample)))

for sample in vbf_samples:
    vbf_getlist.append((filename, 'signal_histos/%s_%s' % (category, sample)))

print "\n =========== GET ================== \n" 
gg_hlist = tools.get(gg_getlist)
vbf_hlist = tools.get(vbf_getlist)

# normalize histograms
for h in gg_hlist:
    h.Scale(1/h.Integral())

for h in vbf_hlist:
    h.Scale(1/h.Integral())

# -- ggf -----------------------
gg_h125 = gg_hlist[0]
gg_h120 = gg_hlist[1]
gg_h130 = gg_hlist[2]

gg_h125.SetTitle("GGF_M125");
gg_h120.SetTitle("GGF_M120");
gg_h130.SetTitle("GGF_M130");

gg_l120 = [gg_h120, gg_h125]
gg_l130 = [gg_h130, gg_h125]

# -- vbf -----------------------
vbf_h125 = vbf_hlist[0]
vbf_h120 = vbf_hlist[1]
vbf_h130 = vbf_hlist[2]

vbf_h125.SetTitle("VBF_M125");
vbf_h120.SetTitle("VBF_M120");
vbf_h130.SetTitle("VBF_M130");

vbf_l120 = [vbf_h120, vbf_h125]
vbf_l130 = [vbf_h130, vbf_h125]

print "\n =========== STACK AND RATIO ====== \n" 
tools.stackAndRatio(gg_l120, title='BDT_Score_GGF_M125_vs_M120', ytitleratio='M125/M120', log=False, yrange=(0,0.06), xtitle="bdt_score", 
                    ldim=[0.23, 0.65, 0.55, 0.83])
tools.stackAndRatio(gg_l130, title='BDT_Score_GGF_M125_vs_M130', ytitleratio='M125/M130', log=False, yrange=(0,0.06), xtitle="bdt_score",
                    ldim=[0.23, 0.65, 0.55, 0.83])
tools.stackAndRatio(vbf_l120, title='BDT_Score_VBF_M125_vs_M120', ytitleratio='M125/M120', log=False, yrange=(0,0.08), xtitle="bdt_score", 
                    ldim=[0.23, 0.65, 0.55, 0.83])
tools.stackAndRatio(vbf_l130, title='BDT_Score_VBF_M125_vs_M130', ytitleratio='M125/M130', log=False, yrange=(0,0.08), xtitle="bdt_score",
                    ldim=[0.23, 0.65, 0.55, 0.83])

