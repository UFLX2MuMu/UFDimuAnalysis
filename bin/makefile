# MakeFile

#ROOTLIBS = $(shell root-config --libdir) 
#ROOTINCS = $(shell root-config --incdir)  

CC = g++ 
LIBFLAGS = `root-config --libs` -O3 -lXMLIO -lMLP -lMinuit -lTMVA -lTMVAGui

LIBDIR = ../lib/
SDIR = ../selection/
CDIR = ../selection/collection_cleaning/
TDIR = ../tools/
THREADDIR = ../threadpool/
OBJDIR = ../lib/analyzer_objects/

CFLAGS = -I ${CDIR} -I${OBJDIR} -I${LIBDIR} -I${SDIR} -I${TDIR} -I${THREADDIR} -c -g `root-config --cflags`

# The name of the .cxx file you want to compile into an executable
MAIN = classifyDataFrame
#MAIN = TMVAClassification_H2Mu
#MAIN = TMVAClassificationApplication_H2Mu
#MAIN = categorize
#MAIN = dyAndDataForFEWZ
#MAIN = outputCounts
#MAIN = synchronize
#MAIN = example
#MAIN = dycompare
#MAIN = masscalibration
#MAIN = fakes
#MAIN = outputToDataframe
#MAIN = listXMLNodes

MAINRULES1 = ${LIBDIR}Sample.o ${LIBDIR}VarSet.o ${LIBDIR}MassCalibration.o ${SDIR}EventSelection.o ${SDIR}MuonSelection.o ${SDIR}CategorySelection.o  
MAINRULES2 = ${CDIR}EleCollectionCleaner.o ${CDIR}JetCollectionCleaner.o ${CDIR}MuonCollectionCleaner.o ${TDIR}TMVATools.o
MAINRULES3 = ${LIBDIR}DiMuPlottingSystem.o ${TDIR}EventTools.o ${TDIR}PUTools.o ${TDIR}ParticleTools.o libAnalysisObjects.so ${MAIN}.oo 
MAINDEPS   = ${THREADDIR}ThreadPool.hxx ${LIBDIR}BranchSet.h SampleDatabase.cxx ${CDIR}CollectionCleaner.hxx ${LIBDIR}VarSet.h
DEPS       = ${LIBDIR}Cut.h ${LIBDIR}CutSet.hxx SignificanceMetrics.hxx ${CDIR}CollectionCleaner.hxx ${LIBDIR}VarSet.h

# ANALYZER objects that were stored into the ttree need to be linked to ROOT
# by creating a dictionary
OBJH := $(shell find $(OBJDIR) -maxdepth 1 -name '*.h') 
OBJCC := $(shell find $(OBJDIR) -maxdepth 1 -name '*.cc') 

# execute all of the rules to make the main program
# $@ is the left of the :, and $^ is everything on the right
# so main requires execution of all MAINRULES, run these rules then compile by linking all MAINRULES
${MAIN}: $(MAINRULES1) $(MAINRULES2) $(MAINRULES3)
	$(CC) $(LIBFLAGS) -o $@ $^

# define a rule for each .o rule
# %.o  all files ending in .o, % is the filename without the suffix
# $< is the first item on the right side of the :, and $@ is the left side of the :
%.o: %.cxx %.h $(DEPS)
	$(CC) $(CFLAGS) $< -o $@ 

#MAIN has no MAIN.h and needs to be compiled with its own dependencies, so we have a different rule
${MAIN}.oo: ${MAIN}.cxx $(MAINDEPS)
	$(CC) $(CFLAGS) ${MAIN}.cxx -o ${MAIN}.oo 

# Creates a dictionary out of the analyzer objects in OBJDIR and Linkdef.h
# Then creates a shared library object from the dictionary, which we link to our executable
# need this so that ROOT TTrees know how to pull out our analysis objects
MyDict.cxx: ${OBJH} ${OBJDIR}linkdef/Linkdef.h
	rootcint -f $@ -c `root-config --cflags` -p $^

libAnalysisObjects.so: MyDict.cxx ${OBJCC}
	g++ -shared -o$@ `root-config --cflags` -fPIC -I$(ROOTSYS)/include $^

clean:
	rm MyDict*
	rm *.oo
	rm *.so

cleanall:
	rm ${LIBDIR}*.o
	rm ${SDIR}*.o
	rm ${TDIR}*.o
	rm ${CDIR}*.o
	rm MyDict*
	rm *.oo
	rm *.so
