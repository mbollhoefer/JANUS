#!/bin/csh
#

STARTDIR=$(PWD)/..
MYSTARTDIR=$(STARTDIR)

include ../user.mk
MYPLATFORM=$(PLATFORM)
include ../makefiles/platforms.mk
include ../makefiles/default.mk

# where is your "mex" script located?
# even when matlab is installed, the mex command is not necessarily
# installed as well
# usually "mex" is in the same place as the "matlab" command
# if your ar lucky, the command is installed
CMEX=mex

MEXDIR=../matlab

DOBJECTS=$(MEXDIR)/DGNLbilu.$(EXT)\
         $(MEXDIR)/ZGNLbilu.$(EXT)\
         $(MEXDIR)/DSPDbilu.$(EXT)\
         $(MEXDIR)/ZHPDbilu.$(EXT)\
         $(MEXDIR)/bilusol.$(EXT) \
         $(MEXDIR)/bilusolt.$(EXT) \
         $(MEXDIR)/bilusolh.$(EXT) \
         $(MEXDIR)/bildlsol.$(EXT) \
         $(MEXDIR)/bildlsols.$(EXT) \
         $(MEXDIR)/bilunnz.$(EXT)\
         $(MEXDIR)/DSYMbilu.$(EXT)\
         $(MEXDIR)/ZHERbilu.$(EXT)\
         $(MEXDIR)/ZSYMbilu.$(EXT)\
         $(MEXDIR)/DSYMbfspai.$(EXT)\
         $(MEXDIR)/ZSYMbfspai.$(EXT)\
         $(MEXDIR)/ZHERbfspai.$(EXT)\
         $(MEXDIR)/DSYMselbinv.$(EXT)\
         $(MEXDIR)/ZSYMselbinv.$(EXT)\
         $(MEXDIR)/ZHERselbinv.$(EXT)\
         $(MEXDIR)/DSPDcgbilu.$(EXT)\
         $(MEXDIR)/ZHPDcgbilu.$(EXT)\
         $(MEXDIR)/DSYMqmrbilu.$(EXT)\
         $(MEXDIR)/ZSYMqmrbilu.$(EXT)\
         $(MEXDIR)/ZHERqmrbilu.$(EXT)\
         $(MEXDIR)/Dgmresbilu.$(EXT)\
         $(MEXDIR)/Zgmresbilu.$(EXT)
	 
  
ZOBJECTS=

OBJECTS=$(DOBJECTS) $(ZOBJECTS)

all: $(OBJECTS)

clean:
	rm -f $(OBJECTS)

.SUFFIXES: .c .$(EXT)

$(DOBJECTS): ../lib/$(PLATFORM)/libjanus.a
$(ZOBJECTS): ../lib/$(PLATFORM)/libjanus.a 
$(MEXDIR)/%.$(EXT): %.c ../lib/$(PLATFORM)/libjanus.a






# 64 BIT ------------------------------
# gcc/gfortran with mexopts77_gfortran.sh 
#	$(CMEX) -O CFLAGS='$(CFLAGS) -fPIC -fopenmp' LDFLAGS='$(LDFLAGS) -fopenmp' -largeArrayDims -D_LONG_INTEGER_  $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR) $<  -L../lib/$(PLATFORM) -ljanus -lsparspak -lamd -lsuitesparseconfig -lmwlapack -lmwblas -lblaslike -lm -lc -lgfortran ../lib/$(PLATFORM)/MC*.o
	$(CMEX) -O CFLAGS='$(CFLAGS) -fPIC -fopenmp' LDFLAGS='$(LDFLAGS) -fopenmp' -largeArrayDims -D_LONG_INTEGER_  $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR) $< -L../lib/$(PLATFORM) -ljanus -lsparspak -lmtmetis -lamd -lsuitesparseconfig -lmwlapack -lmwblas -lblaslike -lm -lc -lgfortran ../lib/$(PLATFORM)/MC*.o
#	$(CMEX) -O CFLAGS='$(CFLAGS) -fPIC -fopenmp' LDFLAGS='$(LDFLAGS) -fopenmp' -largeArrayDims -D_LONG_INTEGER_ $(MATCHING) $(FORTRANNAMES)  -outdir $(MEXDIR) $< -L../lib/$(PLATFORM) -ljanus -lsparspak -lmtmetis  -lamd -lsuitesparseconfig -lmwlapack -lmwblas -lblaslike -lm -lc -lgfortran $(MKL_LIBPATH) -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5 ../lib/$(PLATFORM)/MC*.o
# 64 BIT ------------------------------


$(MEXDIR)/D*.$(EXT): ../lib/$(PLATFORM)/libjanus.a  
$(MEXDIR)/Z*.$(EXT): ../lib/$(PLATFORM)/libjanus.a 
