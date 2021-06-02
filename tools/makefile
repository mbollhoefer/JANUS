include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk


main: $(TOOLS)


$(DIRTOOLS)$(FLOAT)/*.o: ../include/ilupackmacros.h
$(DIRTOOLS)$(FLOAT)/*.o: ../include/ilupack*.h
$(DIRTOOLS)$(FLOAT)/*.o: ../include/namesilupack.h

$(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSYMmergematricess.o: $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSYMmergematrices.o
