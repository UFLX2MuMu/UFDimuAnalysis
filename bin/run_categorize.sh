#!/bin/bash
./categorize --binning=-4 --var=dimu_mass_Roch --categories=tree_categorization_final.xml --systematics="JES_up JES_down PU_up PU_down" --sig_xlumi=0
./categorize --binning=-4 --var=dimu_mass_KaMu --categories=tree_categorization_final.xml --systematics="JES_up JES_down PU_up PU_down" --sig_xlumi=0
