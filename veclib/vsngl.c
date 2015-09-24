/*****************************************************************************
 * y[i] = (float) x[i]
 *
 * $Id: vsngl.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

void vsngl (int_t n, const double *x, int_t incx, float *y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (float) x[i*incx];
}
