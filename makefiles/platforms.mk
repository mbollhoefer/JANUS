ifeq ($(MYPLATFORM),GNU32)
include $(MYSTARTDIR)/makefiles/makefile.includeGNU32

else
ifeq ($(MYPLATFORM),GNU64)
include $(MYSTARTDIR)/makefiles/makefile.includeGNU64

else
ifeq ($(MYPLATFORM),GNU64_long)
include $(MYSTARTDIR)/makefiles/makefile.includeGNU64_long

else
ifeq ($(MYPLATFORM),INTEL32)
include $(MYSTARTDIR)/makefiles/makefile.includeINTEL32

else
ifeq ($(MYPLATFORM),INTEL64)
include $(MYSTARTDIR)/makefiles/makefile.includeINTEL64

else
ifeq ($(MYPLATFORM),INTEL64_long)
include $(MYSTARTDIR)/makefiles/makefile.includeINTEL64_long

else
ifeq ($(MYPLATFORM),PGI32)
include $(MYSTARTDIR)/makefiles/makefile.includePGI32

else
ifeq ($(MYPLATFORM),PGI64)
include $(MYSTARTDIR)/makefiles/makefile.includePGI64

else
ifeq ($(MYPLATFORM),PGI64_long)
include $(MYSTARTDIR)/makefiles/makefile.includePGI64_long

else
ifeq ($(MYPLATFORM),MACOSX32)
include $(MYSTARTDIR)/makefiles/makefile.includemacosx32

else
ifeq ($(MYPLATFORM),MACOSX64)
include $(MYSTARTDIR)/makefiles/makefile.includemacosx64

else
ifeq ($(MYPLATFORM),MACOSX64_long)
include $(MYSTARTDIR)/makefiles/makefile.includemacosx64_long

else
ifeq ($(MYPLATFORM),aix)
include $(MYSTARTDIR)/makefiles/makefile.includeaix

else
ifeq ($(MYPLATFORM),alpha)
include $(MYSTARTDIR)/makefiles/makefile.includealpha

else
ifeq ($(MYPLATFORM),sun)
include $(MYSTARTDIR)/makefiles/makefile.includesun

else
ifeq ($(MYPLATFORM),hpux)
include $(MYSTARTDIR)/makefiles/makefile.includehpux

else
ifeq ($(MYPLATFORM),irix)
include $(MYSTARTDIR)/makefiles/makefile.includeirix

else
ifeq ($(MYPLATFORM),sgi)
include $(MYSTARTDIR)/makefiles/makefile.includesgi

endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

