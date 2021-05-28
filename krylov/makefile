include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk


main: $(KRYLOV)


$(DIRKRYLOV)$(FLOAT)/*.o: ../include/ilupackmacros.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/ilupack*.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/namesilupack.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/long_integer.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/f2c.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/blas.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/namesblas.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/lapack.h
$(DIRKRYLOV)$(FLOAT)/*.o: ../include/nameslapack.h
$(DIRKRYLOV)$(FLOAT)/$(FLOAT)sqmrs.o:          $(DIRKRYLOV)$(FLOAT)/$(FLOAT)sqmr.o
