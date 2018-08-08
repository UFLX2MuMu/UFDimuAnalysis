#!/bin/bash

declare -a vars=("bdt_score" 
                 "dimu_mass_Roch" 
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
    montage imgs/*${var}*_T_*b0.png imgs/*${var}*_root_*.png -tile 4x4 -geometry +2+2 grids/${var}_b0_grid.png
    montage imgs/*${var}*_T_*b1.png imgs/*${var}*_root_*.png -tile 4x4 -geometry +2+2 grids/${var}_b1_grid.png
done

declare -a categories=("root"
                       "c0"
                       "c1"
                       "c2" 
                       "c3"
                       "c4"
                       "c5"
                       "c6"
                       "c7" 
                       "c8"
                       "c9"
                       "c10"
                       "c11"
                       "c12"
                      )

for c in "${categories[@]}"
do
    montage imgs/*${c}*b0.png -tile 4x4 -geometry +2+2 grids/${c}_b0_grid.png
    montage imgs/*${c}*b1.png -tile 4x4 -geometry +2+2 grids/${c}_b1_grid.png
done
