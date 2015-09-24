/*****************************************************************************
 * xcopy:  y[i] = x[i].
 *
 * Use memcpy for cases where both skips are unity.
 *
 * $Id: xcopy.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <string.h>
#include <cfemdef.h>


#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dcopy (int_t n, const double* x, int_t incx,
                             double* y, int_t incy)
{
  register int_t i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (double));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}


void icopy (int_t n, const int_t* x, int_t incx,
                             int_t* y, int_t incy)
{
  register int_t i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (int_t));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}


void scopy (int_t n, const float* x, int_t incx,
                             float* y, int_t incy)
{
  register int_t i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (float));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}
