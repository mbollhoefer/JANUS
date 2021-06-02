# decide which type of matching should be used
MATCHING=-D_MC64_MATCHING_
# decide whether to use Intel MKL
#MKL=-D_USE_MKL_

ifeq ($(MKL),-D_USE_MKL_)
MKL_INCLUDE=-I$(MKLROOT)/include 
MKL_LIBPATH=-L$(MKLROOT)/lib 
else
endif

METIS4_LIB=



LIBPRECONDITIONERS=$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)janus.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bilu.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlspd.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsym.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsyms.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bfspaisym.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bfspaisyms.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symselbinv.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symselbinvs.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)spdapp_supernode.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symcompute_supernode.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symcompute_supernodes.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bilusol.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsol.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsols.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bilusolt.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bilusolh.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilu1t_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilupt_blocks_omp.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlpt_blocks_omp.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlpt_blocks_omp.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlpt_blockss.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilu_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bilu1t_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildl1t_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildl1t_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildl1t_blockss.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildl_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildl_blockss.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)cosine_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)mcosine_blocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)mcosine_sblocks.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)biludriver.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildldriver.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildldrivers.o\
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bicholdriver.o

OBJECTPRECONDITIONERS=

PRECONDITIONERS=      $(LIBPRECONDITIONERS) $(OBJECTPRECONDITIONERS)

LIBKRYLOV=   $(DIRKRYLOV)$(FLOAT)/$(FLOAT)gmres.o        \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)pcg.o          \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)bisinit.o      \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)brkdn.o        \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)mgsro.o        \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)stopbis.o      \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)tidycg.o       \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)sqmr.o         \
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)sqmrs.o\
             $(DIRKRYLOV)$(FLOAT)/$(FLOAT)zdotc2.o        

OBJECTKRYLOV=

KRYLOV=      $(LIBKRYLOV) $(OBJECTKRYLOV)


LIBORDERINGS=$(DIRORDERINGS)$(FLOAT)/$(FLOAT)permamd.o\
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64amd.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64null.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64rcm.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64amd.o \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64amds.o \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64rcm.o  \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64rcms.o \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_null.o  \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_nulls.o \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_mtmetis.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64_mtmetis.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_mtmetis.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_mtmetiss.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symwm.o \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symwms.o \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)spdpermPP.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)spdpermrcm.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)indPPF3.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)qsortR2I.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permnull.o	      \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permrcm.o	           \
             $(DIRORDERINGS)$(FLOAT)/$(FLOAT)swm.o


OBJECTORDERINGS=

ORDERINGS=      $(LIBORDERINGS) $(OBJECTORDERINGS)






LIBTOOLS=   $(DIRTOOLS)$(FLOAT)/$(FLOAT)qsortgnl.o                \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsortgnl.o                 \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)swapj.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)swapm.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsort.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qsort.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janussolver.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janussol.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janussolt.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janussolh.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janusnnz.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)etree.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)compute_supernode.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)generic_matvec.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)matvec.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)mattvec.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)mathvec.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)symmatvec.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)symmatvecs.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)generic_printmatrix.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)printmatrix.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)AMGinit.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)spdAMGinit.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)symAMGinit.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)symAMGinits.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janusdelete.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janusinit.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)janusdefault.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)blocksparsedelete.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)sparsedelete.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)blocksparseinit.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)sparseinit.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSetupGraph.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSetupGraph_epsilon.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)Free.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)geteps.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)Malloc.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)Realloc.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)scale.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)spdscale.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)symscale.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)clear.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)evaluatetime2.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)cc_etimes.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)ddot2.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsorti.o

OBJECTTOOLS=

TOOLS=      $(LIBTOOLS)  $(OBJECTTOOLS)






SPARSPAK=$(DIRSPARSPAK)genrcm.o \
         $(DIRSPARSPAK)rcm.o    \
         $(DIRSPARSPAK)rootls.o \
         $(DIRSPARSPAK)degree.o \
         $(DIRSPARSPAK)revrse.o \
         $(DIRSPARSPAK)fnroot.o







STARTDIR=$(PWD)

MYSTARTDIR=$(STARTDIR)







# where are the headers
INCDIR=$(MYSTARTDIR)/include

# where are the libraries
LIBDIR=$(MYSTARTDIR)/lib/$(MYPLATFORM)


.SUFFIXES: .c .f .F .mod .o .a
.DEFAULT: main

$(FLOAT)/$(FLOAT)%.o: %.c
	$(CC)  $(CCFLAGS)  -I$(INCDIR) $(MATCHING) $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $@ $<

$(FLOAT)/%.o: %.c
	$(CC)  $(CCFLAGS)  -I$(INCDIR) $(MATCHING) $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $@ $<


$(FLOAT)/$(FLOAT)%.o: %.F
	$(FCOMPILE)

$(FLOAT)/%.o: %.F
	$(FCOMPILE)

$(FLOAT)/%.o: %.f
	$(FCOMPILE)





