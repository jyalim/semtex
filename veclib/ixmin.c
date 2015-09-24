/*****************************************************************************
 * ixmin: index of minimum value in x.
 *
 * $Id: ixmin.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>


int_t idmin (int_t n, const double* x, int_t incx)
{
  register int_t i, imin;
  register double  xmin;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = x[0];
  imin = 0;

  for (i = 1; i < n; i++)
   if (x[i*incx] < xmin) {
      xmin = x[i*incx];
      imin = i;
    }

  return imin;
}


int_t iimin (int_t n, const int_t* x, int_t incx)
{
  register int_t i, xmin, imin;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = x[0];
  imin = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] < xmin) {
      xmin = x[i*incx];
      imin = i;
    }

  return imin;
}


int_t ismin (int_t n, const float* x, int_t incx)
{
  register int_t i, imin;
  register float   xmin;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = x[0];
  imin = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] < xmin) {
      xmin = x[i*incx];
      imin = i;
    }

  return imin;
}
