##############################################################################
# Makefile for dns_ff, spectral element solver for unsteady
# Navier--Stokes problems, with generalised body forcing added.
#
# $Id: Makefile,v 8.1 2015/04/20 11:14:13 hmb Exp $
##############################################################################

PROG = dns

SEM   = ..
VPATH = $(SEM)/src

ifdef MPI		# -- Make a parallel executable.
  PROG := $(PROG)_mp
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis auxfield bcmgr boundary bsys condition domain   \
	   edge element family feml field geometry history  \
           integration matrix mesh misc particle statistics \
           fieldforce data2df svv
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = $(addsuffix .h,$(SEMFILES)) sem.h

# ----------------------------------------------------------------------------
# Default build for Navier--Stokes solver.
#
NS_FILES = drive integrate nonlinear dnsanalysis
NS_OBJ   = $(addsuffix .o,$(NS_FILES))

$(PROG): sem $(SEMOBJ) $(NS_OBJ)
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

# ----------------------------------------------------------------------------
# Clean up executables and anything else
#
distclean:
	$(MAKE) clean
	rm -f dns dns_mp

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.cpp
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.cpp
endif
