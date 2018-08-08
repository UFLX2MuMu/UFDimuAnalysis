#!/bin/bash

declare -a vars=("bdt_score" 
                 "dimu_mass_Roch" 
                 "dimu_mass_KaMu" 
                 "mu1_pt" 
                 "mu1_eta" 
                 "mu2_pt" 
                 "mu2_eta" 
                 "dimu_pt" 
                 "dimu_eta" 
                 "dimu_abs_dEta" 
                 "dimu_abs_dPhi" 
                 "jet1_eta" 
                 "jet1_pt" 
                 "jet2_eta" 
                 "jet2_pt" 
                 "dijet1_mass" 
                 "dijet2_mass" 
                 "dijet1_abs_dEta" 
                 "dijet2_abs_dEta" 
                 "nJets" 
                 "nJetsCent" 
                 "nJetsFwd" 
                 "nBMed" 
                 "nEle" 
                 "nExtraMu" 
                 "nVertices" 
                 "MET")

for var in "${vars[@]}"
do
   ./categorize --categories=tree_categorization_final.xml --var=$var  
   ./categorize --categories=tree_categorization_final.xml --var=$var --binning=1
done

#./grid.sh
