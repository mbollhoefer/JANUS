# makefile using 64 bit long integer


# --------------------------------
# C sample drivers
# main driver that reads in a matrix in matrix market format,
# the user defines which of the blocking strategies to turn
# on/off and select drop tolerance
# 1. version that uses the iterative solver without scaling
#    the system in advance (sometimes not recommended)
 MAIN=janusdriver
# 2. version that applies the iterative solver to the diagonally
#    scaled system (sometimes superior)
# MAIN=janusdriver_scaled
#
# C++ variant of janusdriver
# MAIN=janusdriver_cpp
#
#
# simple real general driver, double precision 
# MAIN=dgnljanus
#
# simple real SPD driver, double precision
# MAIN=dspdjanus
#
# simple real symmetric indefinite driver, double precision
# MAIN=dsymjanus
#
# simple complex general driver, double precision
# MAIN=zgnljanus
#
# simple complex HPD driver, double precision
# MAIN=zhpdjanus
#
# simple complex Hermitian indefinite driver, double precision
# MAIN=zherjanus
#
# simple complex symmetric indefinite driver, double precision
# MAIN=zsymjanus
#
#
#
# --------------------------------
# -------- further drivers -------
# --------------------------------
# factored approximate inverse via Neumann series based on block incomplete Cholesky/LDL
# real symmetric driver
# MAIN=bfspaidsym
# complex Hermitian driver
# MAIN=bfspaizher
# complex symmetric driver
# MAIN=bfspaizsym
# --------------------------------
# selected inverse based on block incomplete Cholesky/LDL
# real symmetric driver
# MAIN=dsymselbinv
# complex Hermitian driver
# MAIN=zherselbinv
# complex symmetric driver
# MAIN=zsymselbinv

#MKL=-D_USE_MKL_
#USE_METIS4=-D_USE_METIS4_

ifeq ($(MKL),-D_USE_MKL_)
MKL_INCLUDE=-I$(MKLROOT)/include 
MKL_LIBPATH=-L$(MKLROOT)/lib 
MKL_LIBS=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core  -lpthread -liomp5
else
MKL_LIBS=-llapack -lblas -lblaslike
endif

ifeq ($(USE_METIS4),-D_USE_METIS4_)
METIS4_LIB=-lmetis -lcamd
else
METIS4_LIB=
endif




$(MAIN).out: $(MAIN).o mmio.o
#	ifort -nofor_main  -O3 -fPIC -qopenmp -i8 -o$(MAIN).out $(MAIN).o mmio.o -L ../lib/INTEL64_long/ -ljanus  -lmtmetis   $(METIS4_LIB) -lamd -lsuitesparseconfig -lsparspak $(MKL_LIBPATH) $(MKL_LIBS) 
	gfortran -O -fPIC -fopenmp -fdefault-integer-8 -fimplicit-none -m64 -mcmodel=medium -o$(MAIN).out $(MAIN).o mmio.o -L ../lib/GNU64_long/ -ljanus  -lmtmetis $(METIS4_LIB) -lamd -lsuitesparseconfig -lsparspak $(MKL_LIBPATH) $(MKL_LIBS) -lstdc++
#	gfortran -g -fsanitize=address,signed-integer-overflow -fPIC -fopenmp -fdefault-integer-8 -fimplicit-none -m64 -mcmodel=medium -o$(MAIN).out $(MAIN).o -L ../lib/GNU64_long/ -ljanus  -lmtmetis  $(METIS4_LIB) -lamd -lsuitesparseconfig -lsparspak -llapack -lblas  -lstdc++ -lasan

# $(MAIN).o: ../lib/INTEL64_long/libjanus.a
$(MAIN).o: ../lib/GNU64_long/libjanus.a
# $(MAIN).o: ../lib/GNU64_long/libjanus.a -lasan


.c.o:
#	icc -c -O3 -fPIC -qopenmp -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL)                               -I ../include $(MKL_INCLUDE) -I. $*.c
#	icc -c -O3 -fPIC -qopenmp -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL) -ipo -prof-use -prof-dir=/tmp -I ../include $(MKL_INCLUDE) -I. $*.c
#	icc -c -g  -fPIC -qopenmp -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL)                               -I ../include $(MKL_INCLUDE) -I. $*.c
	gcc -c -O -fPIC -fopenmp -m64 -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL) -mcmodel=medium -I ../include $(MKL_INCLUDE) -I. $*.c
#	gcc -c -g -fsanitize=address -fsanitize=signed-integer-overflow -fPIC -fopenmp -m64 -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL) -mcmodel=medium -I ../include $(MKL_INCLUDE) -I. $*.c
#	gcc -c -g                                                       -fPIC -fopenmp -m64 -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL) -mcmodel=medium -I ../include $(MKL_INCLUDE) -I. $*.c

.cpp.o:
#	icc -c -O3 -fPIC -qopenmp -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL)                                -I ../include $(MKL_INCLUDE) -I. $*.cpp
#	icc -c -O3 -fPIC -qopenmp -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64 $(MKL) -ipo  -prof-use -prof-dir=/tmp -I ../include $(MKL_INCLUDE) -I. $*.cpp
	g++ -c -O -fPIC -fopenmp -m64 -D__UNDERSCORE__ -D_LONG_INTEGER_ -DMKL_ILP64  $(MKL) -mcmodel=medium -I ../include $(MKL_INCLUDE) -I. $*.cpp

clean:
	rm -f *.o *.out
