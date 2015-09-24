/*****************************************************************************
 * xcndst:  conditional assignment:  if (y[i]) z[i] = x[i].
 *
 * $Id: xcndst.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dcndst (int_t n,
	     const double*  x, int_t incx,
	     const int_t* y, int_t incy,
	           double*  z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) if (y[i*incy]) z[i*incz] = x[i*incx];
}


void icndst (int_t n,
	     const int_t* x, int_t incx,
	     const int_t* y, int_t incy,
	           int_t* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) if (y[i*incy]) z[i*incz] = x[i*incx];
}


void scndst (int_t n,
	     const float*   x, int_t incx,
	     const int_t* y, int_t incy,
	           float*   z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) if (y[i*incy]) z[i*incz] = x[i*incx];
}
