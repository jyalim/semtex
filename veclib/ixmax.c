/*****************************************************************************
 * ixmax: index of maximum value in x.
 *
 * $Id: ixmax.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>


int_t idmax (int_t n, const double* x, int_t incx)
{
  register int_t i, imax;
  register double  xmax;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = x[0];
  imax = 0;

  for (i = 1; i < n; i++)
   if (x[i*incx] > xmax) {
      xmax = x[i*incx];
      imax = i;
    }

  return imax;
}


int_t iimax (int_t n, const int_t* x, int_t incx)
{
  register int_t i, xmax, imax;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = x[0];
  imax = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] > xmax) {
      xmax = x[i*incx];
      imax = i;
    }

  return imax;
}


int_t ismax (int_t n, const float* x, int_t incx)
{
  register int_t i, imax;
  register float   xmax;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = x[0];
  imax = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] > xmax) {
      xmax = x[i*incx];
      imax = i;
    }

  return imax;
}
