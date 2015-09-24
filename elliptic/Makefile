##############################################################################
# Makefile for spectral element Helmholtz solver.
#
# $Id: Makefile,v 8.1 2015/04/20 11:14:13 hmb Exp $
##############################################################################

SEM   = ..
VPATH = $(SEM)/src

ifdef MPI
  PROG = elliptic_mp
else
  PROG = elliptic
endif

# ----------------------------------------------------------------------------
#

SEMFILES  = auxfield bcmgr boundary bsys condition data2df domain edge \
	    element family feml field geometry integration matrix mesh \
	    misc svv
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = sem.h

# ----------------------------------------------------------------------------
# Default build.
#
ELLOBJ = drive.o helmholtz.o

$(PROG): sem $(SEMOBJ) $(ELLOBJ)
	$(LD) -o $@ $(ELLOBJ) $(SEMOBJ) $(LDFLAGS)

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
	rm -f elliptic elliptic_mp


# ----------------------------------------------------------------------------
# Special machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
