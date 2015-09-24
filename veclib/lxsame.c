/*****************************************************************************
 * lxsame: return 1 if all elements of x any y match, else zero.
 *
 * Floating point versions use absolute tolerances on the allowable difference.
 *
 * $Id: lxsame.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#define EPSSP   6.0e-7
#define EPSDP   6.0e-14
#define EPSm30  1.0e-30

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif


int_t lisame (int_t n,
	      const int_t* x, int_t incx,
	      const int_t* y, int_t incy)
{ 
  register int_t i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) if (x[i*incx] != y[i*incy]) return 0;
  
  return 1;
}


int_t ldsame (int_t n,
	      const double* x, int_t incx,
	      const double* y, int_t incy)
{ 
  register int_t i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) if (fabs (x[i*incx] - y[i*incy]) > EPSDP) return 0;

  return 1;
}


int_t lssame (int_t n,
	      const float* x, int_t incx,
	      const float* y, int_t incy)
{ 
  register int_t i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

#if  defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) if (fabs  (x[i*incx] - y[i*incy]) > EPSSP) return 0;
#else
  for (i = 0; i < n; i++) if (fabsf (x[i*incx] - y[i*incy]) > EPSSP) return 0;
#endif
  
  return 1;
}
