#
#  Compiler options.
#
# jeff's flags
#FFLAGS= +O3 +es +U77 +E3
#
FFLAGS= $(MYFFLAGS)
#
# Library path name
#
LIBNAM= $(MYLIB)/libharm.a
#
#
#  Compile, archive and clean.
#
.f.a:
	f77 $(FFLAGS) -c  $<
	ar rv $@ $*.o
	rm -f $*.o
#
#  List all the target objects.
#
$(LIBNAM): \
	$(LIBNAM)(angles.subs.o) \
	$(LIBNAM)(avgspher.subs.o)
#
#  Set index.
#
$(LIBNAM): ; ranlib $(LIBNAM)
