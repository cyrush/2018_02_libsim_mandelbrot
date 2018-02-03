###############################################################################
# Makefile for building mandelbrot simulation using libsim. 
###############################################################################

############################## USER EDIT SECTION ################################

# Fill in the Path to the VisIt installation. This is the directory that contains 
# VisIt's 2.13.0 directory. Mac users may want to
# use /path/to/VisIt.app/Contents/Resources to locate the "visit" directory
# within an app bundle.
VISITHOME=path/to/VisIt.app/Contents/Resources/

# Set this to the version of VisIt that you use
VISITVERSION=2.13.0

# Choose one, depending on your system architecture
VISITARCH=darwin-x86_64
#VISITARCH=linux-x86_64

# Edit your compiler and its settings
CXX=clang++
CPPFLAGS=
CXXFLAGS=-O3
LDFLAGS=
LIBS=

#################################################################################
SIMDIR=$(VISITHOME)/$(VISITVERSION)/$(VISITARCH)/libsim/V2

SIM_CXXFLAGS=-I$(SIMDIR)/include
SIM_LDFLAGS=-L$(SIMDIR)/lib
SIM_LIBS=-lsimV2 -ldl

SRC=mandelbrot.C patch.C
OBJ=$(SRC:.C=.o)


SRC_BATCH=mandelbrot_batch.C patch.C
OBJ_BATCH=$(SRC_BATCH:.C=.o)

all: mandelbrot mandelbrot_batch

clean:
	rm -f mandelbrot mandelbrot_batch $(OBJ)

mandelbrot: $(OBJ)
	$(CXX) -o mandelbrot $(OBJ) $(LDFLAGS) $(SIM_LDFLAGS) $(SIM_LIBS) $(LIBS)

mandelbrot_batch: $(OBJ_BATCH)
	$(CXX) -o mandelbrot_batch $(OBJ_BATCH) $(LDFLAGS) $(SIM_LDFLAGS) $(SIM_LIBS) $(LIBS)


.C.o:
	$(CXX) $(CXXFLAGS) $(SIM_CXXFLAGS) $(CPPFLAGS) -c $<