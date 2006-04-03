#
# Rules.mk
# written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 11:41:04 EET 2001
#

.PHONY: all
all: $(LIBRARIES) $(PROGRAMS)

# TAGS and tags are for people working with vi(m) or emacs who wants to
# use the tag feature of the editors (VERY useful...)
#.PHONY: TAGS tags
TAGS: $(tag_src)
	$(ETAGS) $(tag_src)

tags: $(tag_src)
	$(CTAGS) $(tag_src)

.PHONY: clean realclean distclean
clean:
	-@set -e; \
	echo "rm -f $(all_objs) *.mod *.d"; \
	rm -f $(all_objs) *.mod *.d; 
ifneq ($(PROGRAMS),)
	@set -e; for i in $(PROGRAMS); do \
	echo "rm -f $(bindir)/$$i"; rm -f $(bindir)/$$i; done
endif
ifneq ($(LIBRARIES),)
	@set -e; for i in $(LIBRARIES); do \
	echo "rm -f $(libdir)/$$i"; rm -f $(libdir)/$$i; done
endif

realclean: clean
	@set -e; \
	rm -f tags TAGS y.tab.c y.tab.h lex.yy.c; 

distclean: realclean
	-rm -f config.h Makefile Config.mk config.log config.status

.PHONY: install install_bin install_lib install_all 
install: install_bin

install_all: install_bin install_lib

install_bin: $(PROGRAMS) $(inst_bindir) 
ifneq ($(PROGRAMS),)
	@set -e; \
	for i in $(PROGRAMS); do \
		$(INSTALL_PROGRAM) $(bindir)/$$i $(inst_bindir); \
		echo "Installed $$i in $(inst_bindir)"; \
	done
endif

install_lib: $(LIBRARIES) $(inst_libdir)
ifneq ($(LIBRARIES),)
	@set -e; \
	for i in $(LIBRARIES); do \
		$(INSTALL_DATA) $(libdir)/$$i $(inst_libdir); \
		echo "Installed $$i in $(inst_libdir)"; \
	done
endif

$(inst_bindir):
	@$(mkinstalldirs) $(inst_bindir)

$(inst_libdir):
	@$(mkinstalldirs) $(inst_libdir)

%.o: %.c
	$(CC) $(CFLAGS) -c $< 

%.o: %.f
	$(F90) $(FFLAGS) -c $< 

%.o: %.f90
	$(F90) $(FFLAGS) -c $< 

(%.o): %.o
	$(AR) rc $@ $<

