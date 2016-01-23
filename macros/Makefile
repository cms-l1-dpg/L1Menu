######### GENERAL STUFF: DON NOT CONSIDER
CXX = g++
CXXFLAGS  = -std=c++11
CXXFLAGS  += -Wall -fPIC -O2 -g -pg
#CXXFLAGS  += -Wextra -Weffc++ #Asking for troubles
CXXFLAGS  += $(subst -I, -isystem , $(shell root-config --cflags))

LD         = g++
LDFLAGS    =
LIBS       = $(shell root-config --libs)   -lTable

ifneq ($(shell echo $$CMSSW_BASE), )
  CXXFLAGS  += -isystem $(CMSSW_RELEASE_BASE)/src
  CXXFLAGS  += -isystem $(CMSSW_BASE)/src/
  CXXFLAGS  += -isystem $(shell echo $$CMSSW_FWLITE_INCLUDE_PATH | cut -d':' -f4)
  LIBS += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -lFWCoreFWLite
  #LIBS +=  -L$(subst include, lib, $(shell echo $$CMSSW_FWLITE_INCLUDE_PATH | cut -d':' -f4)) -lboost_program_options
endif

OBJS = L1Ntuple.o L1AlgoFactory.o

PROGRAM = testMenu2016
#PROGRAM = RunL1Menu2016

MKBIN = $(CXX) $(CXXFLAGS) `root-config  --libs --cflags` -lMinuit -lGenVector

%.o : %.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

%.o : %.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

all : $(PROGRAM)

comparePlots.exe : comparePlots.cc 
	$(MKBIN) $< -o $@ 

RunL1Menu2016 : RunL1Menu2016.o $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) $^ $(LIBS)  -o $@
	@echo "done"

testMenu2016 : testMenu2016.o L1Menu2016.o  $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) $^ $(LIBS)  -o $@
	@echo "done"

clean:
	rm -f *exe *.o $(PROGRAM) $(OBJS)


