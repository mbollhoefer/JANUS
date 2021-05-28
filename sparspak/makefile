include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk


.c.o:
	$(CC)  $(CCFLAGS)  -I$(INCDIR)  $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $@ $<

.f.o:
	echo $(FF);$(FF)  $(FFFLAGS)  -I$(INCDIR) -c -o $@ $<

.F.o:
	echo $(FF);$(FF)  $(FFFLAGS)  -I$(INCDIR)  $(LONGINTEGER) -c -o $@ $<


main : $(SPARSPAK)

$(DIRSPARSPAK)*.o: ../include/ilupack_fortran.h
