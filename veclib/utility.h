#ifndef UTILITY_H
#define UTILITY_H
///////////////////////////////////////////////////////////////////////////////
// Utility.h:  C++ wrapper for veclib utility and memory routines.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// $Id: utility.h,v 8.1 2015/04/20 11:14:19 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef.h>

template<class T> inline T sqr(T x)             { return x * x;            }
template<class T> inline T sgn(T x)             { return (x < 0) ? -1 : 1; }
template<class T> inline T clamp(T t, T a, T b) { return max(min(t,b),a);  }

// Max & Min are part of STL now, so have been removed from the above.

#ifndef M_PI
const double M_PI   = 3.14159265358979323846;
#endif

#ifndef TWOPI
const double TWOPI  = 6.28318530717958647692;
#endif

#ifndef PI_180
const double PI_180 = 0.01745329251994329576;
#endif

const double EPSm3  = 1.0e-3;
const double EPSm4  = 1.0e-4;
const double EPSm5  = 1.0e-5;
const double EPSm6  = 1.0e-6;
const double EPSSP  = 6.0e-7;
const double EPSm7  = 1.0e-7;
const double EPSm8  = 1.0e-8;
const double EPSm12 = 1.0e-12;
const double EPSDP  = 6.0e-14;
const double EPSm14 = 1.0e-14;
const double EPSm20 = 1.0e-20;
const double EPSm30 = 1.0e-30;

const int StrMax    = STR_MAX;

enum lev {WARNING, ERROR, REMARK};

extern "C" {
  void       message (const char *routine, const char *txt, int level);

  double     dclock  ();
  float      sclock  ();

  double   *dvector  (int_t nl, int_t nh);
  double  **dmatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
  double ***d3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		      int_t dl, int_t dh);

  float    *svector  (int_t nl, int_t nh);
  float   **smatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
  float  ***s3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		      int_t dl, int_t dh);

  int_t   *ivector  (int_t nl, int_t nh);
  int_t  **imatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
  int_t ***i3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		     int_t dl, int_t dh);
  
  void freeDvector  (double    *v, int_t nl);
  void freeDmatrix  (double   **m, int_t nrl, int_t ncl);
  void freeD3matrix (double  ***t, int_t nrl, int_t ncl, int_t ndl);

  void freeSvector  (float     *v, int_t nl);
  void freeSmatrix  (float    **m, int_t nrl, int_t ncl);
  void freeS3matrix (float   ***t, int_t nrl, int_t ncl, int_t ndl);

  void freeIvector  (int_t   *v, int_t nl);
  void freeImatrix  (int_t  **m, int_t nrl, int_t ncl);
  void freeI3matrix (int_t ***t, int_t nrl, int_t ncl, int_t ndl);
}

#endif

