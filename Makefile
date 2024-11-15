# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2018 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc


# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -ldl -g -std=c++20
CXX_ALL:=$(CXX_COMMON)
CXX_ALL+= -lpythia8

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "Usage: make mainXX"

# The Makefile configuration.
Makefile.inc:
	$(error Error: PYTHIA must be configured, please run "./configure"\
                in the top PYTHIA directory)

# PYTHIA libraries.
$(PREFIX_LIB)/libpythia8.so :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)


pythia2root.so: pythia2rootDct.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ -c -w -I$(ROOT_INCLUDE) $(CXX_SHARED) $(CXX_ALL) -pthread -std=c++17 -m64 -L$(ROOT_LIB) -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -Wl,-rpath,$(ROOT_LIB) -lm -ldl -rdynamic
#`$(ROOT_BIN)root-config --cflags` `$(ROOT_BIN)root-config --glibs` 
pythia2rootDct.cc: pythia2root.h pythia2rootLinkDef.h
	export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT_LIB);\
	 $(ROOT_BIN)rootcint -f $@ -c -I$(PREFIX_INCLUDE) $^

pythia2fastjet: $$@.cc $(PREFIX_LIB)/libpythia8.a 
ifeq ($(FASTJET3_USE),true)
	$(CXX) $< -o $@ -w -I$(FASTJET3_INCLUDE) $(CXX_ALL)\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet -lRecursiveTools -lNsubjettiness -lfastjettools \
	 -Wl,-rpath,./
else
	@echo "Error: $@ requires fastjet"
endif


fastjet_example_bare: $$@.cc 
ifeq ($(FASTJET3_USE),true)
	$(CXX) $< -o $@ -w -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet -lRecursiveTools -lNsubjettiness -lfastjettools \
	 -Wl,-rpath,./
else
	@echo "Error: $@ requires fastjet"
endif


clustering: $$@.cc
	$(CXX) $< -o $@ -w -I$(FASTJET3_INCLUDE) -I$(ROOT_INCLUDE) $(CXX_SHARED) $(CXX_ALL) -pthread -std=c++17 -m64 -L$(ROOT_LIB) -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -Wl,-rpath,$(ROOT_LIB) -lm -ldl -rdynamic -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet -lRecursiveTools -lNsubjettiness -lfastjettools -Wl,-rpath,./



# Internally used tests, without external dependencies.
test% : test%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ $(CXX_ALL) $(GZIP_INC) $(GZIP_FLAGS)

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat;\
	rm -f weakbosons.lhe; rm -f Pythia8.promc; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;\
	rm -f pythia2root pythia2fastjet
