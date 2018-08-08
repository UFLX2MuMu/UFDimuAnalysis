#!/bin/bash
#./categorize --binning=-1 --categories=tree_categorization_final.xml --systematics="JES_up JES_down PU_up PU_down"
#./masscalibration --reductionFactor=1 --peak=JPsi --x=phi_plus
#./masscalibration --reductionFactor=1 --peak=JPsi --x=phi_minus
#./masscalibration --reductionFactor=1 --peak=JPsi --x=eta_plus
#./masscalibration --reductionFactor=1 --peak=JPsi --x=eta_minus
#./masscalibration --reductionFactor=1 --peak=JPsi --x=pt_plus
#./masscalibration --reductionFactor=1 --peak=JPsi --x=pt_minus
#./masscalibration --reductionFactor=1 --peak=JPsi --x=dimu_pt
#
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=phi_plus
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=phi_minus
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=eta_plus
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=eta_minus
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=pt_plus
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=pt_minus
#./masscalibration --reductionFactor=1 --peak=Upsilon --x=dimu_pt
#
#./masscalibration --reductionFactor=1 --peak=JPsi --x=phi_plus --leadPtMin=40  
#./masscalibration --reductionFactor=1 --peak=JPsi --x=phi_minus --leadPtMin=40
#./masscalibration --reductionFactor=1 --peak=JPsi --x=eta_plus --leadPtMin=40
#./masscalibration --reductionFactor=1 --peak=JPsi --x=eta_minus --leadPtMin=40
#./masscalibration --reductionFactor=1 --peak=JPsi --x=pt_pluss --leadPtMin=40
#./masscalibration --reductionFactor=1 --peak=JPsi --x=pt_minus --leadPtMin=40
#./masscalibration --reductionFactor=1 --peak=JPsi --x=dimu_pt --leadPtMin=40

./masscalibration --reductionFactor=1 --peak=Z --x=phi_plus
./masscalibration --reductionFactor=1 --peak=Z --x=phi_minus
./masscalibration --reductionFactor=1 --peak=Z --x=eta_plus
./masscalibration --reductionFactor=1 --peak=Z --x=eta_minus
./masscalibration --reductionFactor=1 --peak=Z --x=pt_plus
./masscalibration --reductionFactor=1 --peak=Z --x=pt_minus
./masscalibration --reductionFactor=1 --peak=Z --x=dimu_pt

