/*****************************************************************************
 * xmask:  conditional assignment:  if (y[i]) z[i] = w[i]; else z[i] = x[i].
 *
 * $Id: xmask.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dmask (int_t n,
	    const double*  w, int_t incw,
	    const double*  x, int_t incx,
	    const int_t* y, int_t incy,
	          double*  z, int_t incz)
{
  register int_t i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = (y[i*incy]) ? w[i*incw] : x[i*incx];
}


void imask (int_t n,
	    const int_t* w, int_t incw,
	    const int_t* x, int_t incx,
	    const int_t* y, int_t incy,
	          int_t* z, int_t incz)
{
  register int_t i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = (y[i*incy]) ? w[i*incw] : x[i*incx];
}


void smask (int_t n,
	    const float*   w, int_t incw,
	    const float*   x, int_t incx,
	    const int_t* y, int_t incy,
	          float*   z, int_t incz)
{
  register int_t i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = (y[i*incy]) ? w[i*incw] : x[i*incx];
}
