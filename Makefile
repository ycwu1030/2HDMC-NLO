# Makefile for 2HDMC 
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
# Choose your C++ compiler here (in general g++ on Linux systems):
CC = $(LT)/../bin/f++
LDFLAGS=-lgsl -lgslcblas

#Optimisation level, eg: -O3
OPT=-O3
#OR debug level: -g(n=1,2,3)
DEBUG= 

VPATH=src 

CFLAGS= -Wall $(DEBUG) $(OPT) -I$(LT)/../include
CFLAGS+=-DNLOCOUPLINGS

OBJDIR=lib
LIBDIR=$(OBJDIR)
SOURCES=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp CouplingFunctionTypeI.cpp CouplingFunctionTypeII.cpp CouplingFunctionTypeLS.cpp CouplingFunctionTypeFL.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
LIB=lib2HDMC.a
LDFLAGS+=-L$(LIBDIR) -l2HDMC -lgsl -lgslcblas -lm -L$(LT) -looptools -L/usr/local/gfortran/lib -lgfortran
#LDFLAGS+=-L$(LIBDIR) -l2HDMC
LIBS=
PROG=CalcPhys.x CalcGen.x CalcHiggs.x CalcHybrid.x CalcHMSSM.x CalcMSSM.x CalcInert.x CalcLH.x Demo.x

# To use HiggsBounds/HiggsSignals for Higgs constraints, put the corresponding libraries
# in the "lib" subdirectory and uncomment the following lines.
# Note that both programs: HiggsBounds (> 4.2.0) and HiggsSignals (> 1.4.0) must be present
#CFLAGS+=-DHiggsBounds
#LDFLAGS+=-L$(LIBDIR) -lHS -lHB -lgfortran
#SOURCES+=HBHS.cpp

.PHONY: directory lib clean distclean

all: lib $(PROG)

directory: $(OBJDIR)
	@ mkdir -p $(OBJDIR)

$(OBJDIR)/%.o : %.cpp %.h directory
	$(CC) $(CFLAGS) -c $< -o $@

lib: $(addprefix $(LIBDIR)/, $(LIB))

$(addprefix $(LIBDIR)/, $(LIB)): $(addprefix $(OBJDIR)/, $(OBJECTS))
	@ echo "Making library $@"
	@ ar rcs $@ $(addprefix $(OBJDIR)/, $(OBJECTS))

%.x: %.cpp $(addprefix $(LIBDIR)/, $(LIB))
	@ $(CC) $< -Isrc $(CFLAGS) $(LDFLAGS) $(addprefix $(LIBDIR)/, $(LIBS)) -o $@

clean:
	@ echo "Cleaning library"
	@ rm -f $(addprefix $(OBJDIR)/, *.o)
	@ rm -f $(addprefix $(LIBDIR)/, $(LIB))

distclean:
	@ make -s clean
	@ echo "Cleaning programs"
	@ rm -f $(PROG)
