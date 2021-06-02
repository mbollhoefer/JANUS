#!/bin/csh
#
# makefile for ILUPACK
# 
STARTDIR=$(PWD)
MYSTARTDIR=$(STARTDIR)

include user.mk
MYPLATFORM=$(PLATFORM)
include makefiles/platforms.mk
include makefiles/default.mk

# ----------------------   program specific settings   -----------------------

DIRSPARSPAK=sparspak/

DIRPRECONDITIONERS=preconditioners/
DIRKRYLOV=krylov/
DIRORDERINGS=orderings/
DIRBLASLIKE=blaslike/


MYSTARTDIR=$(PWD)

.PHONY:  rpreconditioners rkrylov rorderings  rtools  rlibs \
         cpreconditioners ckrylov corderings  ctools  clibs \
         clean clear blaslike sparspak matlab


all: rpreconditioners rkrylov rorderings  rtools rlibs \
     cpreconditioners ckrylov corderings ctools clibs \
     blaslike sparspak matlab

real: rpreconditioners blaslike rkrylov rorderings rtools rlibs 

complex: cpreconditioners blaslike ckrylov corderings ctools clibs

rpreconditioners: 
	@ echo build double real preconditioners
	@ cd preconditioners; \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=D"

cpreconditioners: 
	@ echo build double complex preconditioners
	@ cd preconditioners; \
        $(MAKE) "ARITHMETIC="    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=Z"


rkrylov: 
	@ echo build double real krylov subspace methods
	@ cd krylov; \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=D"

ckrylov: 
	@ echo build double complex krylov subspace methods
	@ cd krylov; \
        $(MAKE) "ARITHMETIC="    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=Z"


rorderings: 
	@ echo build double real orderings
	@ cd orderings; \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=D"

corderings: 
	@ echo build double complex orderings
	@ cd orderings; \
        $(MAKE) "ARITHMETIC="    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=Z"
	

blaslike: 
	@ echo build blaslike library
	@ cd blaslike;\
        $(MAKE) "MYPLATFORM=$(PLATFORM)" "MYSTARTDIR=$(STARTDIR)"

sparspak:
	@ echo build sparspak library
	@ cd sparspak;\
        $(MAKE) "MYPLATFORM=$(PLATFORM)" "MYSTARTDIR=$(STARTDIR)"


rtools:  
	@ echo build double real tools
	@ cd tools;\
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=D"

ctools:  
	@ echo build double complex tools
	@ cd tools;\
        $(MAKE) "ARITHMETIC="    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=Z"






rlibs:   rpreconditioners rkrylov rorderings  rtools  sparspak blaslike
	@ echo build double real libraries
	@ cd lib;\
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=D"

clibs:   cpreconditioners ckrylov corderings ctools   sparspak blaslike
	@ echo build double complex libraries
	@ cd lib;\
        $(MAKE) "ARITHMETIC="    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=Z"


matlab:  
	@ echo build matlab executables
	@ cd matlabsrc;\
        $(MAKE) 



clean:
	rm -rf  */D/*.o */Z/*.o */*.o */*.mod 


clear:
	rm -rf */*.o */*/*.o */D/* */Z/* */*.mod;\
        rm -rf lib/$(PLATFORM)/lib*.a\
        rm -rf ../routines/*.mex*


