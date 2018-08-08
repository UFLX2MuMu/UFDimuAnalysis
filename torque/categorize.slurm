#!/bin/bash
#SBATCH --job-name=categorize
#SBATCH -o categorize_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=puno@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=8:00:00
#SBATCH --array=0-27

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

var=${vars[$SLURM_ARRAY_TASK_ID]}

#Set the compiler architecture and CMSSW
export SCRAM_ARCH=slc6_amd64_gcc620
export CMS_PATH=/cvmfs/cms.cern.ch
source ${CMS_PATH}/cmsset_default.sh
export CVSROOT=:pserver:anonymous@cmssw.cvs.cern.ch:/local/reps/CMSSW

# source CMSSW libs to get root
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc620/cms/cmssw/CMSSW_9_0_0_pre5/src
eval `scram runtime -sh`
cd /home/puno/h2mumu/UFDimuAnalysis_v2/bin/

# Add analysis objects to library path (can't find this shared object library otherwise)
LD_LIBRARY_PATH=/home/puno/h2mumu/UFDimuAnalysis_v2/bin:${LD_LIBRARY_PATH}

# run the executable
date
hostname
pwd
echo "JOB ID: ${SLURM_ARRAY_JOB_ID}"
echo "ARRAY ID: ${SLURM_ARRAY_TASK_ID}"
echo "VAR TO PLOT: $var"
echo ""
./categorize --categories=2 --var=$var --nthreads=10 --binning=0 --systematics="JES_up JES_down PU_up PU_down" 
./categorize --categories=2 --var=$var --nthreads=10 --binning=1 --systematics="JES_up JES_down PU_up PU_down"
date
