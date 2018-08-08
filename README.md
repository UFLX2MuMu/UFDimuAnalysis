# UFDiMuAnalysis
```
* Code used to process the Analyzed TTrees for the h2mu analysis at UF
    - Grabs info from the TTree to ...
    - make plots, categorize, synchronize with other groups, perform significance studies, resolution studies, etc
    - there are a variety of objects to help with
        % getting objects out of the TTree, plotting, bookkeeping, weighting, event selection, collection cleaning, parallelization, etc

* Create executables in the bin directory. Edit the makefile and set MAIN to the file you want to compile, e.g. MAIN=categorize.
    - This will compile categorize.cxx into categorize. Then run categorize with ./categorize.
    - bin/categorize.cxx is the most important executable
    - check it out to see how to use the code

* the makefile creates a shared library, so you need to update your LD_LIBRARY_PATH so that it can find it in the bin directory 
    - in bash: LD_LIBRARY_PATH=/path/to/UFDimuAnalysis/bin:${LD_LIBRARY_PATH}
    - in tcsh: setenv LD_LIBRARY_PATH /path/to/UFDimuAnalysis/bin:${LD_LIBRARY_PATH}

* bin/categorize.cxx produces the histograms and data vs MC stacks for the run1, run2, or autocategorizer categories
    - ./categorize --categories=1
    - ./categorize --categories=2
    - ./categorize --categories=some_autocategorizer_xml_output.xml
    - see categorize.cxx for the other options
    - limit setting requires --sig_xlumi=0 --binning=-4 

* create c++ executables in bin, python scripts in python
* the rest of the directories have objects that help with the h2mu analysis
    - Check out the ABOUT file in other folders to see how they help 
    - most of the .cxx/.h files should have short descriptions at the top and lots of comments
```

# UF H2Mu Analysis Code
https://github.com/acarnes/UFDimuAnalysis

## 1 - Installing the code

Clone or fork the repo from the above github above into some directory on lxplus or the ufhpc. 

Then you need to get access to some gcc compiler newer than gcc620 and you need to add to your LD_LIBRARY_PATH so that the system knows about some shared objects the library makes. Here's a snippet I took from my bash profile showing how to do that on the ufhpc.

```
# You need access to the slc6_amd64_gcc620 c++ compiler
# I get that by cmsenv'ing something higher than CMSSW_9_X_X
# on the ufhpc I just cd to a preinstalled one
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc620/cms/cmssw/CMSSW_9_0_0_pre5/src
eval `scram runtime -sh`             # robust way to cmsenv
cd /home/puno/

# Let the system know about the shared objects the makefile creates
LD_LIBRARY_PATH=/path/to/your/UFDimuAnalysis/bin:${LD_LIBRARY_PATH}
```

## 2 - Running the code

The bin folder (/path/to/your/UFDimuAnalysis/bin) has all of the cxx files for the various plotting and studies. All of the other directories have objects the executables call to help perform the studies (the objects help with cuts, categorization, keeping a database of the variables we use, cleaning object collections, etc). If you need to create some other cxx files to run other studies make them in the bin directory.

The code runs on samples stored at the ufhpc or CERN depending on your location. In the cxx scripts there is a GetSamples(..., "UF/CERN") function. You need to choose UF or CERN appropriately in the function to access our data and mc samples. The detailed sample information is stored in SampleDatabase.cxx. The samples we run on in this library are created using the code from https://github.com/acarnes/UfHMuMuCode. The UfHMuMuCode runs over the MINIAOD using CMSSW and CRAB and creates the root files listed in SampleDatabase.cxx. 

To run anything in the bin directory you need to edit the makefile and change MAIN=whatever to MAIN=study, where study.cxx is the cxx file you want to compile and run.  Then you will create the executable by compiling with

```
	make
```

Then you run the study with 

```
./study option1 option2
```

Where study is the executable compiled from study.cxx via make and option1 and option2 are some inputs to the executable. See some of the sh files in the bin directory that show how to run some of the studies.

Try this with example.cxx. This file is a very basic example showing how to use the different objects of the library to make selections, categorize, and plot. It makes plots of the ggf dimuon mass in the run1 categories. After making, run it with ./example

Here is a summary of the most used cxx files.

### categorize.cxx

Plot histograms of any variable in all of the categories. It outputs a root file with stacks comparing Data/MC and the individual histograms for the data and the signal and bkg MC. It is already set up to run on the run1 2012 categories, the run2 2016 BDT categories, or any autocategorizer categories. If you want to run it on another set of categories you just need to create a "categorizer" (see selection/CategorySelection.cxx and the .h) and add an option in categorize.cxx to use that. 

You may also customize the cuts by swapping out the cut objects the script references (see selection/MuonSelection.cxx and selection/EventSelection.cxx and their .h).  This script also makes the histograms needed for limit setting by plotting the dimuon mass in each category.

After compiling with the makefile run the script like so to plot the dimuon mass with the data blinded in the signal region.

```
./categorize --categories=2 --var=dimu_mass_KaMu --nthreads=10 --binning=0
```

Binning controls whether the plots are blinded in the signal region, the range of the plots, and the number of bins. Negative binning unblinds the mass in the signal region 120-130 GeV. The unblinding is used for the mass plots for limit setting. The option --sig_xlumi=0 is also required for the limit setting mass histos. It makes sure the signal normalization is given by the efficiecy times acceptance rather than scaling by the xsec times lumi to get the total number of expected events. Higgs combine will automatically scale to the number of events using the official yellow report values rather than the xsecs from MC.  

To debug or run over fewer events you can use --reductionFactor=10. This would reduce the events you run over by 10x, then scale the plots so that you get the correct numbers for the luminosity. There are other options you can check out as well in the script. It is well commented. If you want histograms for the up and down systematics you can use --systematics="JES_up JES_down PU_up PU_down". The systematics can be used to get error bands on the Data/MC comparison plots that account for these uncertainties. 

There is a slurm script at UFDimuAnalysis/torque/categorize.slurm that submits a batch job on the ufhpc to plot all of the variables used in the analysis in the different run2categories. The categorize.cxx script is parallelized assigning threads per sample. 

The output files from categorize.cxx are set to go to bin/rootfiles. So make sure you have that directory or the output files won't be saved. This is true for a lot of the scripts. Check where they save the output and make sure you have the appropriate directory structure. Or you can change the save location in the script to save them wherever you want. 

### masscalibration.cxx

This study plots the mean and the resolution for a dimuon mass peak (Z, J/Psi, Upsilon) vs some variable x ( x = mu+/- phi, eta, pt, or dimu_pt). It makes histograms of the mass peak for x_i < x < x_j and x_j < x < x_k and x_k < x < x_l and so on until there are mass histograms for the whole range of x. You set the binning and can even set a variable binning. 

The mass histograms for each x range are then fit with Voigtian and the mean and experimental resolution are extracted for each histogram. The mean and resolution are then plotted vs x. 

Here are some examples running the script to plot the Z-peak mean and resolution vs the different variables.

./masscalibration --reductionFactor=1 --peak=Z --x=phi_plus
./masscalibration --reductionFactor=1 --peak=Z --x=phi_minus
./masscalibration --reductionFactor=1 --peak=Z --x=eta_plus
./masscalibration --reductionFactor=1 --peak=Z --x=eta_minus
./masscalibration --reductionFactor=1 --peak=Z --x=pt_plus
./masscalibration --reductionFactor=1 --peak=Z --x=pt_minus
./masscalibration --reductionFactor=1 --peak=Z --x=dimu_pt

### classify.cxx

This script runs ROOT's TMVA to classify signal vs background using MC in a certain mass window passing given selections. This is a simplified version of TMVAClassification_H2Mu.cxx, which was used to train the BDTs for the 2016 run2 categories. The script outputs an xml file with all of the information defining the trained classifier. The resulting xml file can be used in TMVA to apply the classifier to new data/MC. People call this file the "weights" file.
 
TMVAClassificationApplication_H2Mu.cxx is an example showing how to apply the trained classifier to new data/MC. The categorize.cxx script does this when it applies an xml file in TMVA to get the bdt_score for an event. The categorizer then uses the bdt_score to place the event into the correct 2016 run2 BDT/mu_eta category. tools/TMVATools help with this.  

classify.cxx will also output a root file that shows the training and testing results, which can be analyzed by looking at it in a TBrowser or by using the TMVAGui.  You can check out the ROC to see how well the classifier did in the TMVAGui. 

Run using

```
./classify
```

There are no options. Edit the code to change the classifier and its settings or whatever else. Add options if you would like to.

### outputToDataframe.cxx

This script outputs a flat ntuple and a csv of all the features available to the analysis (see lib/VarSet.cxx and the .h).  It does this for every MC/Data sample in a given mass window. These files are used as input to the autocategorizer to get automatically create optimum categories based upon a set of features. 

You can also use this to work with our data and mc samples outside of ROOT. You could, for example, use the .csv files to classify signal vs background with sci-kit learn or keras (neural nets).

Run using

```
./outputToDataframe
```

There aren't any options for this one. Edit the code to change the behavior. You can train a classifier on these flat ntuples using classifyDataFrame.cxx. With classifyDataFrame the selections are already made and the computation time is spent only on classifying. The classify script tests each event from the larger ttrees of SampleDatabase.cxx and adds them to the training and testing sets one at a time if they pass the selections.    

### outputCounts.cxx

This script outputs the # signal,  # background, S/B, and significance in each category to csv. The script works with the different significance metrics in bin/SignificanceMetrics.hxx. This is useful to evaluate the behavior of the categorization.

Run using

```
./outputCounts
```

There are no options for this one. Edit the code to change the behavior.

### synchronize.cxx

Synchronize with other groups in the analysis. The script outputs event info for a list of certain events (or all of the events). Generally, it's best to direct the output to a LOG file of some sort since there is a ton of output. The list of events to investigate can be imported from csv if there are many of them. The script can also output lists of events in each category to csv.  You can run the script for specific samples or for all of the samples. 

Run using

```
./synchronize
```

There are no options for this one. Edit the code to modify the behavior.

## 2 - The Supporting Objects of the Library

There are many objects in the UFDimuAnalysis library that help with the studies run in bin. These are covered by directory.

### lib/
These objects are the backbone of the library. 

Sample.cxx is an object describing one of our samples. It keeps track of the name of the sample, the cross section, whether it is MC, Data, or Background, the rootfile(s) with the data for the sample. bin/SampleDatabase.cxx creates Sample objects and links the rootfile(s) for the datasets to the objects, setting the names, the sample type, etc. The object also keeps track of the number of events, the way to weight the event based upon the xsec and the lumi, and some other things.

VarSet.cxx/.h keeps track of all of our objects and the collections of objects. It also provides the ability to access any variable used in the analysis using the name. The muon pt for the 0th muon in the vector of muons in the event can be accessed like so sample.vars.muons->at(0).pt. The VarSet also has collections for the dimuon pairs, electrons, jets, etc. To access a variable by name it must be added to the map in the object that maps a TString to a Function. Once in the map the value of a feature l(ike the pt of the first muon of the dimuon pair) can be accessed like so sample.vars.getValue("mu1_pt"). 

DiMuPlottingSystem.cxx has a bunch of functions to help with plotting stuff.

MassCalibration.cxx has the functionality to do the bin/masscalibration.cxx studies. It makes mass histograms for different bins of some x variable. The object also fits the histograms and extracts the fit info plotting a TGraph of the mean and resolution of a given peak vs x. 

### lib/analyzer_objects
These are the objects stored into the rootfiles used in this library (the ones in SampleDatabase.cxx. The makefile creates a shared object library that tells root about our object structure so that it can pull the objects out of the rootfiles correctly.


### selection/
The objects here help with various selections. We have the MuonSelection and EventSelection objects used by various studies to cut events from the analysis in order to reduce the number to a feasible amount for plots and limits without losing much signal. You can implement the interfaces to create your own selection objects. Just add to the cxx/h files following the examples already there.

The CategorySelection files define our categorizer objects. The categorizer objects place events into categories based upon their kinematics and object counts. They also keep track of the names of each category and the plots and lists of events for each category. The CategorySelection objects do not cut any events from the analysis, they simply place the events passing the muon and event selections into categories. You can implement the interface to create your own categorizer and use it in bin/categorize.cxx. Just add to the CategorySelection.cxx/h files.

The collection cleaning objects also do not cut any events from our analysis. These are used to remove untrustworthy objects from certain collections. For instance we might remove fake jets from the jets collection or we might remove electrons that aren't isolated from the collection of electrons in the event.

### threadpool/
This provides us with a threadpool object that lets us parallelize different scripts. Categorize.cxx and some other scripts use this to parallelize by assigning one thread to each sample. When there are 10 threads in the pool the code will run over 10 samples in parallel, and whenever one of the 10 threads finishes it will automatically be assigned to the next sample in the queue. In this way the threads are always doing productive work grabbing new samples as soon as they are free.

### tools/
The objects in here help output event information, read TMVA xml files, and do some other random things. This is a collection of miscellaneous tools. 

### torque/
I put my batch submission scripts in here. The ufhpc used to use torque for batch submission but now it uses slurm. So the directory name is outdated. 

## 3 - Python Utilities
python/ has utilities to read in FEWZ files and make histograms from them (python/fewz). There are also utilities to fit the background (python/fit_bkg). See h2mu_poly.py for a background fit example. There is also an object that shows how to make a basic workspace+datacard for limit setting in higgs combine (python/limit_setting). This is the prototype for the actual limit setting code we used. There are no systematics or interpolations or anything like that in the workspace it's a very basic prototype. There are also a bunch of plotting utilities to make the charts and the plots for the Analysis Note. 

## 4 - The Autocategorizer
https://github.com/acarnes/bdt/ has the code for the autocategorization. But you need to use the binned_categorizer branch. The master is a boosted decision tree package I made to do some trigger work way earlier. The code was repurposed to make the autocategorizer. 

The autocategorizer makes optimum categories for events binned along the pdf variable -- the pdf variable is the one the signal and background histograms/fits are along for higgs combine. The autocategorizer minimizes the expected p-value on the background onlyi. It minimizes the expected p-value by categorizing to maximize the difference between the S+B and B-only hypotheses. Different sensitivity metrics (S/sqrt(B), Asimov, Punzi, etc) are used to quantify the difference between the hypothesis rating how the difference compares to an expected fluctuation.  

 The algorithm needs to know the class (signal/bkg), the weight (xsec, lumi,i scale_factors, mc_weight/sum_weights), the bin of each event, and the features. It takes in the csv files or flat ntuples with the sig/bkg, weight, bin, and feature info. It outputs an xml file of the categories with the cuts for each variable in the form of a decision tree. 

You can see how the code runs in bdt/studies/h2mumu/BasicTrainAndTest.cxx. There is a makefile to compile the executable. Just run `make` then run with `./BasicTrainAndTest`. See run.sh for an example of how to run the code with all of the options and their descriptions. LoadEvents.hxx has code to load the events from the csv/ntuples. I was running on the ntuples on the ufhpc at /home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt ,  so you can copy those over to wherever you are working or set up your code on the uf hpc and run on them directly. 

The code is pretty well commented so check it out. The probability density function we put into higgs combine for H->mumu is along the dimuon mass spectrum. So the bin values used in the autocategorizer were along that variable. I had twenty bins (if I remember right) representing which 0.5 GeV mass bin in the 120-130 GeV signal region the event fell into. A bin value of -1 means an event fell outside the signal window. These bin=-1 events can optionally be used used to scale the mc in the window based upon the data/mc ratio in the control region outside the window and for some other optional corrections. If you are using the bdt_score as the pdf for higgs combine then you would bin the signal/bkg the same way along that variable.

The bin/outputToDataframe.cxx script in the UFDimuAnalysis repo creates the needed csv/ntuples. The output xml file can be automatically used by UFDimuAnalysis/bin/categorize.cxx through the XMLCategorizer class. You just need to set the appropriate option like ./categorize --categories=output_xml_file.xml. Because of string->function map in VarSet.h/.cxx, our library can automatically calculate the value of a feature based upon the name of the feature. The xml file has the cuts based upon the feature names, which along with our varset map allows us to automate the autocategorizer->plotting process.


## 5 - Limits
The limit setting code takes in the root files from categorize.cxx with the options --binning=-4 --var=dimu_mass_KaMu --sig_xlumi=0. Binning -4 sets an unblinded singal region, fine binning, and a wide range, while sig_xlumi=0 makes sure the signal histograms are not scaled by the xsec times lumi since higgs combine will automatically do this. You can also run bias studies. Here is the info needed to set up the code. 

```
cd ~
mkdir h2mu_limit_setting # directory for your limit setting code
# install latest higgs combine via https://cms-hcomb.gitbooks.io/combine/content/part1/
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/ # CMSSW_X_X_X is the version installed for higgs combine 

## Set up the viktor analysis code (see https://github.com/uiowahep/Analysis)
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/
git clone https://github.com/uiowahep/Analysis ViktorAnalysis
mkdir build
mkdir build/ViktorAnalysis
cd build/ViktorAnalysis
cmake ../../ViktorAnalysis
cd ../../ViktorAnalysis

## Source the correct CMSSW
exec /bin/bash
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/
eval `scramv1 runtime -sh`

## source viktor's env.sh script so that python knows where to find all of the scripts
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/ViktorAnalysis
source $PWD/config/env.sh
```

For the latest combine you may have to edit Modeling/combine/generate_precombine.py and change setPhysicsModelParameters to setParameters and change nuisancesToFreeze to freezeParameters.

Before running the code always source the setup.sh file. Make sure to replace my directory paths with your paths for your limit setting library path and your higgs combine cmssw src path. 

To run the limits use run_limits.sh. To get plots from those files use get_limits.sh. Similarly for the bias studies use run_bias.sh and get_bias.sh. The label and jobLabel in the run and get scripts need to match. These use the information from one of the settings files (e.g. https://github.com/uiowahep/Analysis/blob/master/Configuration/higgs/UF_AMC_settings.py )

You need to make your own settings.py file there and call it using --mode UF_AMC in the run/get.sh scripts. You will need to edit the code to get it to read in your own settings file. Do a grep for "UF_AMC" to see where in the code it looks for this input and add yours appropriately.
 
You also need to make sure the directories from your Configuration/higgs/X_settings.py are all linked correctly. Most notably, 

```
cmsswDir        = '/afs/cern.ch/work/a/acarnes/public/h2mumu/limit_setting/combine/CMSSW_7_4_7/src/'
projectDirToUse = '/afs/cern.ch/work/a/acarnes/public/h2mumu/limit_setting/out/limits_and_bias/'
jobLabel        = 'AMC'
```

cmsswDir needs to match your install location for your higgs combine as per the instructions above. projectDirToUse is where the limit setting code will automatically create a bunch of directories and save output, log files, etc. projectDirToUse needs to match that in the run_limts.sh and run_bias.sh and get_limits.sh and get_bias.sh scripts. 

For example you will see this pattern in the run/get.sh files "cmsswDir/somesavedir/jobLabel/$label" e.g. 

```
cd /afs/cern.ch/work/a/acarnes/public/h2mumu/limit_setting/out/limits_and_bias/ftest/AMC/$label/
```
$label is from the .sh script itself, ftest is a directory created by the limit setting library automatically to store some output from the ftest, AMC is my jobLabel, and the long path before all of that is my projectDirToUse. Replace yours appropriately. You can edit the files so that all of these are bash variables at the top and people don't have to do find replace in the future.

You also need to make sure you tell the system where to find the output from categorize.cxx using your settings.py file. I had a file here. 

```
inputFIleUF = '/afs/cern.ch/work/a/acarnes/public/h2mumu/rfiles/validate_UNBLINDED_dimu_mass_Roch_90_200_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b-4_sig-xlumi0.root'
```

The settings.py also controls which background functions are used, ftested, etc. Moreover, you can set the signal model and the sig/bkg model parameters there. The code needs to be updated so that you can use different models in each category. 

Viktor wrote this code so I'm not intimately familiar with all the peculiarities.

## 6 - Creating the Samples
https://github.com/acarnes/UfHMuMuCode
Andrew Brinkerhoff has edited this code a bunch since I last used it. So maybe he can add to this to make it clearer and more detailed. 

This code creates CRAB jobs for use on the grid (I ran all of my jobs from lxplus). It will grab samples from DAS (MC and Data) and create the .root files needed by UFDimuAnalysis. It creates our objects from raw data, makes basic selections, and then saves .root files that are fed into our plotting, mass calibration studies, autocategorizer, etc. 

I was last creating the crab submission jobs using https://github.com/acarnes/UfHMuMuCode/blob/master/UFDiMuonsAnalyzer/test/make_crab_script.py , but I think that Andrew Brinkerhoff is using a newer version  in the crab directory of the repo at https://github.com/acarnes/UfHMuMuCode/blob/master/UFDiMuonsAnalyzer/crab/make_crab_script.py .
