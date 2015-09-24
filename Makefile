##############################################################################
# (GNU) Makefile for spectral element solvers.
#
# $Id: Makefile,v 8.1 2015/04/20 11:14:13 hmb Exp $
##############################################################################

# -- $MAKE supplies the path to GNU make. Since it's typical that make
#    *is* gmake on modern UNIX machines, we make no distinction (you may
#    have to modify).

MAKE = make

tar:
	tar cvf semtex.tar *
	gzip semtex.tar

# ----------------------------------------------------------------------------
dist:
	tar cvf semtex.tar			        \
	README Makefile include src lib mesh sm doc	\
	dns elliptic utility
	gzip semtex.tar

# ----------------------------------------------------------------------------
srcdist:
	tar cvf srcdist.tar				  \
	README Makefile include src veclib femlib mesh sm \
	dns elliptic utility
	gzip srcdist.tar

# ----------------------------------------------------------------------------
# Run this to compile libraries required by programs.

libs:
	cd veclib;		\
	$(MAKE) -s headers;

	cd femlib;		\
	$(MAKE) -s headers;

	cd src;			\
	$(MAKE) -s

	cd veclib;		\
	$(MAKE) -s clean;	\
	$(MAKE) -s;		\
	$(MAKE) -s install

	cd femlib;		\
	$(MAKE) -s clean;	\
	$(MAKE) -s;		\
	$(MAKE) -s install

# ----------------------------------------------------------------------------
# Make version of femlib with MPI.

parlib:
	cd femlib;		\
	$(MAKE) -s headers;

	cd femlib;		\
	$(MAKE) -s clean;	\
	$(MAKE) -s MPI=1;	\
	$(MAKE) -s install MPI=1


# ----------------------------------------------------------------------------
# Run this to compile all serial executables.

all: libs
	cd utility;  $(MAKE) clean; $(MAKE) -s all
	cd elliptic; $(MAKE) clean; $(MAKE) -s
	cd dns;      $(MAKE) clean; $(MAKE) -s; $(MAKE) -s tbcs

# ----------------------------------------------------------------------------
# Compile parallel executables.

parallel: parlib
	cd dns; $(MAKE) MPI=1; $(MAKE) tbcs MPI=1
	cd elliptic; $(MAKE) MPI=1

# ----------------------------------------------------------------------------
# Run a regression test of the (serial) DNS solver.

test:  libs
	cd utility; $(MAKE) -s clean; $(MAKE) -s essential
	cd dns; $(MAKE) -s clean; $(MAKE) -s;
	cd test ; \
	rm -f compare;   ln -s ../utility/compare   . ;	\
	rm -f enumerate; ln -s ../utility/enumerate . ; \
	rm -f dns;       ln -s ../dns/dns . ; \
	./testregress dns

# ----------------------------------------------------------------------------
# Run test of parallel version of DNS solver: do "make test" first.
# Also may need to edit the file test/testregress_mp to set MPI command.

partest: parlib
	cd dns ; $(MAKE) -s clean ; $(MAKE) -s ALIAS=1 MPI=1
	cd test ; rm -f dns_mp; ln -s ../dns/dns_mp . ; testregress_mp dns_mp

# ----------------------------------------------------------------------------
# Clean up.

clean:
	cd veclib;   $(MAKE) clean
	cd femlib;   $(MAKE) clean
	cd utility;  $(MAKE) distclean
	cd elliptic; $(MAKE) distclean
	cd dns;      $(MAKE) distclean
