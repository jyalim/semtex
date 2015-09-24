##############################################################################
# Remake linear algebra primitive library (libvec.a), and maintain header
# files for Veclib, BLAS, & LAPACK calls.
#
# $Id: Makefile,v 8.1 2015/04/20 11:14:19 hmb Exp $
##############################################################################

SEM    = ..
SEMSRC = $(SEM)/src
LIB    = libvec.a
ARCH   = $(shell uname -s)

# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

UTIL       = util

MEMORY     = memory

MATH       = xcopy   xfill   xneg    xvneg   xsadd   xvadd   xssub  \
             xvsub   xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero  \
	     xspow   xvsgn   xznan

OTHER      = xvabs   xvamax  xvpow   xvexp   xvlg10  xvlog   xvatan \
             xvatn2  xvcos   xvsin   xvsqrt  xrand   xvhypot xvtanh

TRIAD      = xsvmvt  xsvpvt  xsvtsp  xsvtvm  xsvtvp  xsvvmt  xsvvpt \
	     xsvvtm  xsvvtp  xvvmvt  xvvpvt  xvvtvm  xvvtvp  xvvvtm \
	     xsvvtt  xvvvtt

RELATIONAL = xseq    xsge    xsle    xslt    xsne

REDUCTION  = icount  ifirst  ixmax   ixmin   xsum    lany    lxsame

CONVERSION = vdble   xvfloa  vsngl   xbrev

MISC       = xscatr      xgathr      	xgathr_scatr		\
	     xscatr_sum  xgathr_sum  	xgathr_scatr_sum	\
	     xramp       xcndst      	xmask       		\
	     xvpoly      xpolint    	xclip			\
	     xspline     xsplquad	xiclip			\
	     xmxv        xmxva	  				\
	     xvvtvvtp    xsvvttvp	 xvvtvvtm

FORTRAN    = xmxm							

ALL        = $(UTIL)       $(MEMORY)    $(MATH)       $(OTHER)   $(TRIAD) \
             $(RELATIONAL) $(REDUCTION) $(CONVERSION) $(MISC)    $(FORTRAN)    
#LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))
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
# Standard rules and definitions.
#
include $(SEMSRC)/Makefile

# ----------------------------------------------------------------------------
# Install header files.
#
headers:
	if test ! -d  $(SEM)/include	;	then    \
		mkdir $(SEM)/include	;	fi

	rm -f $(SEM)/include/cveclib.h;	\
	rm -f $(SEM)/include/veclib.h;	\
	rm -f $(SEM)/include/utility.h;	\
	rm -f $(SEM)/include/blas.h;	\
	rm -f $(SEM)/include/lapack.h

	cp cveclib.h	$(SEM)/include	;	\
	cp veclib.h    	$(SEM)/include	;	\
	cp utility.h	$(SEM)/include	;	\
	cp blas.h	$(SEM)/include	;	\
	cp lapack.h	$(SEM)/include

# ----------------------------------------------------------------------------
# Install library & associated header files.
#
install:
	$(MAKE) headers

	if test ! -d    $(SEM)/lib			;	then    \
		mkdir   $(SEM)/lib			;	fi;     \
	if test ! -d    $(SEM)/lib/$(ARCH)		;	then    \
		mkdir   $(SEM)/lib/$(ARCH)		;	fi;     \
	if test -r      $(LIB)			;	then    \
		rm -f   $(SEM)/lib/$(ARCH)/$(LIB)	;               \
		cp $(LIB) $(SEM)/lib/$(ARCH)		;	fi
	if test $(RANLIB)			;	then	\
		$(RANLIB) $(SEM)/lib/$(ARCH)/$(LIB) 	;	fi
