##############################################################################
# Remake finite element routine library (libfem.a).  Gnu make required.
#
# $Id: Makefile,v 8.1 2015/04/20 11:14:14 hmb Exp $
##############################################################################

SEM    = ..
SEMSRC = $(SEM)/src
ARCH   = $(shell uname -s)

ifdef MPI
  LIB = libfem_mp.a
else
  LIB = libfem.a
endif

ifdef DEBUG_FFT
  DEFINES += -DDEBUG_FFT
endif

# ----------------------------------------------------------------------------
# Create dependency list.
#

INITIAL    = initial

OPERATORS  = polyops  operators  matops  polylib  filter

FFT	   = temfftd  canfft  fourier

MEMORY     = mapping  family  message

IMPORT     = netlib

ALL        = $(INITIAL) $(OPERATORS) $(FFT) $(MEMORY) $(IMPORT)

LMD        = $(foreach routine, $(ALL), $(routine).o)

# ----------------------------------------------------------------------------
# Make library, default action.
#
$(LIB): $(LMD)
	$(AR) $(ARFLAGS) $(LIB) $(LMD)

ifeq ($(ARCH),Linux-alpha)
	ranlib $(LIB)
endif
ifeq ($(ARCH),Linux-iX86)
	ranlib $(LIB)
endif

# ----------------------------------------------------------------------------
# Make and install parallel library version.
parallel:
	rm -f message.o
	make MPI=1
	make install MPI=1

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(SEMSRC)/Makefile

# ----------------------------------------------------------------------------
# Install header files.
#
headers:
	if test ! -d  $(SEM)/include	;	then    \
		mkdir $(SEM)/include	;	fi

	rm -f $(SEM)/include/cfemdef.h	;	\
	rm -f $(SEM)/include/cfemlib.h  ;	\
	rm -f $(SEM)/include/polylib.h  ;	\
	rm -f $(SEM)/include/femlib.h	;			

	cp cfemdef.h	$(SEM)/include	;	\
	cp cfemlib.h	$(SEM)/include	;	\
	cp polylib.h	$(SEM)/include  ;	\
	cp femlib.h	$(SEM)/include

# ----------------------------------------------------------------------------
# Install library & associated header files.
#
install:
	$(MAKE) headers

	if test ! -d $(SEM)/lib 			;	then	\
		mkdir $(SEM)/lib 			; 	fi;	\
	if test ! -d $(SEM)/lib/$(ARCH) 		;	then	\
		mkdir $(SEM)/lib/$(ARCH)		;	fi;	\
	if test -r $(LIB) 				;	then	\
		rm -f $(SEM)/lib/$(ARCH)/$(LIB)		;		\
		cp $(LIB) $(SEM)/lib/$(ARCH)		;	fi
	if test $(RANLIB)				;	then	\
		$(RANLIB) $(SEM)/lib/$(ARCH)/$(LIB)	;	fi


# ----------------------------------------------------------------------------
# Ensure any new default parser variables are picked up.
#
initial.o: defaults.h

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.
#
ifeq ($(ARCH),OSF1)
temfftd.o: temfftd.F
	f90 -c -fast -inline speed -tune host temfftd.F -O5 -unroll 2
temffts.o: temffts.F
	f90 -c -fast -inline speed -tune host temffts.F -O5 -unroll 2
endif

# ----------------------------------------------------------------------------
# Rules to help with .F routines.
#
matops.o:  matops.F
temfftd.o: temfftd.F
