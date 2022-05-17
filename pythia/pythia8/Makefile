CXXFLAGS:=$(shell root-config --cflags) -I$(PYTHIA8)/include
    # $(shell root-config --libs) -L$(PYTHIA8)/lib -lpythia8

# Define directory paths
FASTJETCONFIG := fastjet-config

# Compiler options
FASTJETLIBS := $(shell $(FASTJETCONFIG) --libs)
FASTJETLIBS := -Wl,-rpath,$(FASTJET_ROOT)/lib -lm -L$(CGAL_ROOT)/lib -lCGAL -lgmp -L$(FASTJET_ROOT)/lib -lfastjettools -lfastjet
FASTJETFLAGS := $(shell $(FASTJETCONFIG) --cxxflags)

ROOTLIBS := $(shell root-config  --libs)
ROOTFLAGS := $(shell root-config --cflags)

CXXFLAGS += $(FASTJETFLAGS) $(ROOTFLAGS)
FLAGS += $(FASTJETLIBS) $(ROOTLIBS)
LDFLAGS += $(FASTJETLIBS) -lfastjetcontribfragile $(ROOTLIBS) $(HepMClib)

# .SUFFIXES:      .o .cxx

# pythia8_simple: pythia8_simple.o

