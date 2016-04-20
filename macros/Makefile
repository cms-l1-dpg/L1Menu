######### GENERAL STUFF: DON NOT CONSIDER
ODIR     =  objs
CXX      =  g++
CXXFLAGS =  -std=c++11
CXXFLAGS += -Wall -fPIC -O3
#CXXFLAGS  += -g -pg
#CXXFLAGS  += -Wextra
#CXXFLAGS += -Wextra -Weffc++ #Asking for troubles
#CXXFLAGS += -Wno-reorder #Dirty fix
CXXFLAGS  += $(subst -I, -isystem , $(shell root-config --cflags))

LD         = g++
LDFLAGS    =
LIBS       = $(shell root-config --libs)   -lTable

ifneq ($(shell echo $$CMSSW_BASE), )
  CXXFLAGS  += -isystem $(CMSSW_BASE)/src/
  CXXFLAGS  += -isystem $(CMSSW_RELEASE_BASE)/src
  CXXFLAGS  += -isystem $(shell scram tool info boost | awk -F"=" '/INCLUDE=(.*)/{print $$NF}')
  LIBS += -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ \
		  -lFWCoreFWLite -lDataFormatsL1TGlobal
  LIBS += -L$(shell scram tool info boost | awk -F"=" '/LIBDIR=(.*)/{print $$NF}') \
		  -lboost_program_options -lboost_filesystem
endif

TestOBJS = L1Ntuple L1AlgoFactory L1Menu2016 L1Plot L1TnP L1uGT
OBJS := $(TestOBJS:%=$(ODIR)/%.o)

PROGRAM = testMenu2016 RunL1Menu2016

MKBIN = $(CXX) $(CXXFLAGS) `root-config  --libs --cflags` -lMinuit -lGenVector

$(ODIR)/%.o : %.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : %.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

all : $(PROGRAM)

comparePlots : comparePlots.cc 
	$(MKBIN) $< -o $@ 

RunL1Menu2016 : $(ODIR)/RunL1Menu2016.o $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) $^ $(LIBS)  -o $@
	@echo "done"

testMenu2016 : $(ODIR)/testMenu2016.o $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) $^ $(LIBS)  -o $@
	@echo "done"

clean:
	rm -f *exe *.o $(PROGRAM) $(OBJS) $(ODIR)/*.o


