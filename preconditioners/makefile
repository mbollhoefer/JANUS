include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk


main: $(PRECONDITIONERS)




$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsols.o:                $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsol.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildl1t_blockss.o:          $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildl1t_blocks.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlpt_blockss.o:          $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlpt_blocks_omp.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsyms.o:                $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildlsym.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bfspaisyms.o:              $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bfspaisym.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symselbinvs.o:              $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symselbinv.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symcompute_supernodes.o:    $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symcompute_supernode.o
$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildldrivers.o:             $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)bildldriver.o
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/ilupackmacros.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/ilupack*.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/namesilupack.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/long_integer.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/f2c.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/blas.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/namesblas.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/lapack.h
$(DIRPRECONDITIONERS)$(FLOAT)/*.o: ../include/nameslapack.h
