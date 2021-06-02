include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk

DIRPRECONDITIONERS=../preconditioners/
DIRKRYLOV=../krylov/
DIRORDERINGS=../orderings/
DIRTOOLS=../tools/
DIRBLASLIKE=../blaslike/
DIRSPARSPAK=../sparspak/

MULTIPLEOBJECTS=$(LIBPRECONDITIONERS) \
                $(LIBKRYLOV) \
                $(LIBORDERINGS) \
                $(LIBTOOLS) 

main: $(MYPLATFORM)/libjanus.a $(MYPLATFORM)/libsparspak.a \
      $(MYPLATFORM)/libblaslike.a 


$(MYPLATFORM)/libjanus.a: $(MULTIPLEOBJECTS) $(NOTDISTRIBUTED)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libblaslike.a:  $(BLASLIKE)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libsparspak.a:  $(SPARSPAK)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@




