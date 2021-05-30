include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk

DIRPRECONDITIONERS=../preconditioners/
DIRKRYLOV=../krylov/
DIRORDERINGS=../orderings/
DIRTOOLS=../tools/
DIRNOTDISTRIBUTED=../notdistributed/
DIRSPARSPAK=../sparspak/
DIRMETIS=../metis/
DIRMTMETIS=../mtmetis/
DIRSUITESPARSECONFIG=../SuiteSparse_config/
DIRAMD=../AMD/
DIRCAMD=../CAMD/

MULTIPLEOBJECTS=$(LIBPRECONDITIONERS) \
                $(LIBKRYLOV) \
                $(LIBORDERINGS) \
                $(LIBTOOLS) 

ifeq ($(USE_METIS4),-D_USE_METIS4_)
main: $(MYPLATFORM)/libjanus.a $(MYPLATFORM)/libsparspak.a \
      $(MYPLATFORM)/libmetis.a  \
      $(MYPLATFORM)/libsuitesparseconfig.a  \
      $(MYPLATFORM)/libamd.a $(MYPLATFORM)/libcamd.a 
else
main: $(MYPLATFORM)/libjanus.a $(MYPLATFORM)/libsparspak.a \
      $(MYPLATFORM)/libsuitesparseconfig.a  \
      $(MYPLATFORM)/libamd.a 
endif
#      $(MYPLATFORM)/libmtmetis.a  \


$(MYPLATFORM)/libjanus.a: $(MULTIPLEOBJECTS) $(NOTDISTRIBUTED)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libsparspak.a:  $(SPARSPAK)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libmetis.a:  $(METIS)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libmtmetis.a:  $(MTMETIS)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libsuitesparseconfig.a:  $(ALLSUITESPARSECONFIG)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libamd.a:  $(ALLAMD)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(MYPLATFORM)/libcamd.a:  $(ALLCAMD)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@




