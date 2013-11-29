


######### GENERAL STUFF: DON NOT CONSIDER
CXX = g++

CXXFLAGS = # -Wall -fPIC -O2 -g -pg

MKBIN = $(CXX) $(CXXFLAGS) `root-config  --libs --cflags` -lMinuit -lGenVector

all : comparePlots.exe

comparePlots.exe : comparePlots.cc 
	$(MKBIN) $< -o $@ 


clean:
	rm -f *exe


