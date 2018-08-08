##############################################
# make_muon_calibration_plots.py             #
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

tfile = 0
print "\n /// Making muon calibration plots for %s ... \n" % filename

# parse the file to get necessary naming information

masscal = re.search("masscalibration_Z_", filename)
medium = re.search("_medium", filename)

if masscal and medium: 
    varname = filename[masscal.start()+len('masscalibration_Z_'):medium.start()]

vartitle = varname.replace('_',' ')

vartitle = vartitle.replace('eta', '#eta')
vartitle = vartitle.replace('phi', '#phi')
vartitle = vartitle.replace('pt', 'p_{t}')

print varname
print '\"'+vartitle+'\"'

if 'plus' in vartitle:
    vartitle = vartitle.replace(' plus', '')
    vartitle += '(#mu^{+})'

if 'minus' in vartitle:
    vartitle = vartitle.replace(' minus', '')
    vartitle += '(#mu^{-})'

if 'dimu' in vartitle:
    vartitle = vartitle = vartitle.replace('dimu','')
    vartitle += '(#mu^{+} #mu^{-})'


if "pt" in vartitle: 
    units = ' (GeV)'
else: 
    units = ''

print varname
print '\"'+vartitle+'\"'

# put plots here
smean_data = []
smean_mc = []
sres_data = []
sres_mc=[]

calibs = ['PF', 'Roch', 'KaMu']

# Make 
for c in calibs:
    #print "   /// (%d) Looking at sample %s ..." % (dy_bug, s)
    smean_data.append( (filename, 'plots/mean_Net_Data_mass_%s_%s' % (c, varname)) )
    smean_mc.append(   (filename, 'plots/mean_ZJets_AMC_mass_%s_%s' % (c, varname)) )

    sres_data.append( (filename, 'plots/resolution_Net_Data_mass_%s_%s' % (c, varname)) )
    sres_mc.append(   (filename, 'plots/resolution_ZJets_AMC_mass_%s_%s' % (c, varname)) )

mean_data = tools.get(smean_data)
mean_mc = tools.get(smean_mc)
res_data = tools.get(sres_data)
res_mc = tools.get(sres_mc)
    
mean_all = []
mean_all.extend(mean_data)
mean_all.extend(mean_mc)

res_all = []
res_all.extend(res_data)
res_all.extend(res_mc)

mean_minmax = (999999,-999999)
res_minmax = (999999,-999999)

# get y min and max of all mean plots
for i in mean_all:
    if "PF" in i.GetTitle():
        i.SetTitle("Uncalibrated (Particle Flow)")
    if "KaMu" in i.GetTitle():
        i.SetTitle("Kalman Filter Muon Corrections")
    if "Roch" in i.GetTitle():
        i.SetTitle("Rochester Corrections")

    minmax = tools.getMinMaxY(i)

    if minmax[0] < mean_minmax[0]:
        mean_minmax = (minmax[0], mean_minmax[1])
    if minmax[1] > mean_minmax[1]:
        mean_minmax = (mean_minmax[0], minmax[1])

miny, maxy = mean_minmax
mean_minmax = (miny-0.2*(maxy-miny), 1.5*(maxy-miny)+maxy)

for i in res_all:
    if "PF" in i.GetTitle():
        i.SetTitle("Uncalibrated (Particle Flow)")
    if "KaMu" in i.GetTitle():
        i.SetTitle("Kalman Filter Corrections")
    if "Roch" in i.GetTitle():
        i.SetTitle("Rochester Corrections")

    minmax = tools.getMinMaxY(i)

    if minmax[0] < res_minmax[0]:
        res_minmax = (minmax[0], res_minmax[1])
    if minmax[1] > res_minmax[1]:
        res_minmax = (res_minmax[0], minmax[1])

miny, maxy = res_minmax
res_minmax = (miny-0.2*(maxy-miny), 1.5*(maxy-miny)+maxy)

print mean_minmax
print res_minmax
    
mean_datak = [mean_data[0], mean_data[2]]
res_datak = [res_data[0], res_data[2]]

mean_datak[0].SetLineColor(1)
mean_datak[1].SetLineColor(8)
res_datak[0].SetLineColor(1)
res_datak[1].SetLineColor(8)

print "\n =========== OVERLAY DATA ====== \n" 
tools.overlay(mean_datak, savename='zcal_data_mean_'+varname, title=('Z Calibration Mean in Data vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=mean_minmax, ytitle='Z-Peak Mean (GeV)', autocolor=False)
tools.overlay(res_datak, savename='zcal_data_res_'+varname, title=('Z Calibration Resolution in Data vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=res_minmax, ytitle='Z-Peak Resolution (GeV)', autocolor=False)

print "\n =========== OVERLAY DATA/MC PF ====== \n" 

mean_data[0].SetTitle('Data (Particle Flow)')
mean_mc[0].SetTitle('Drell Yan MC (Particle Flow)')
res_data[0].SetTitle('Data (Particle Flow)')
res_mc[0].SetTitle('Drell Yan MC (Particle Flow)')

tools.overlay([mean_data[0], mean_mc[0]], savename='zcal_pf_mc-data_mean_'+varname, title=('Z Calibration Mean vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=mean_minmax, ytitle='Voigt Fit Mean (GeV)')
tools.overlay([res_data[0], res_mc[0]], savename='zcal_pf_mc-data_res_'+varname, title=('Z Calibration Resolution vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=res_minmax, ytitle='Voigt Fit Resolution (GeV)')

print "\n =========== OVERLAY DATA/MC Roch ====== \n" 

mean_data[1].SetTitle('Data (Rochester)')
mean_mc[1].SetTitle('Drell Yan MC (Rochester)')
res_data[1].SetTitle('Data (Rochester)')
res_mc[1].SetTitle('Drell Yan MC (Rochester)')

mean_data[1].SetLineColor(1)
mean_mc[1].SetLineColor(8)
res_data[1].SetLineColor(1)
res_mc[1].SetLineColor(8)

tools.overlay([mean_data[1], mean_mc[1]], savename='zcal_roch_mc-data_mean_'+varname, title=('Z Calibration Mean vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=mean_minmax, ytitle='Z-Peak Mean (GeV)', autocolor=False)
tools.overlay([res_data[1], res_mc[1]], savename='zcal_roch_mc-data_res_'+varname, title=('Z Calibration Resolution vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=res_minmax, ytitle='Z-Peak Resolution (GeV)', autocolor=False)

print "\n =========== OVERLAY DATA/MC KaMu ====== \n" 

mean_data[2].SetTitle('Data (Kalman Filter)')
mean_mc[2].SetTitle('Drell Yan MC (Kalman Filter)')
res_data[2].SetTitle('Data (Kalman Filter)')
res_mc[2].SetTitle('Drell Yan MC (Kalman Filter)')

mean_data[2].SetLineColor(1)
mean_mc[2].SetLineColor(8)
res_data[2].SetLineColor(1)
res_mc[2].SetLineColor(8)

tools.overlay([mean_data[2], mean_mc[2]], savename='zcal_kamu_mc-data_mean_'+varname, title=('Z Calibration Mean vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=mean_minmax, ytitle='Z-Peak Mean (GeV)', autocolor=False)
tools.overlay([res_data[2], res_mc[2]], savename='zcal_kamu_mc-data_res_'+varname, title=('Z Calibration Resolution vs '+ vartitle), ldim=[0.3,0.63,0.88,0.88], 
             xtitle=vartitle+units, yrange=res_minmax, ytitle='Z-Peak Resolution (GeV)', autocolor=False)
