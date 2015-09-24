/*****************************************************************************
 * xseq:  y[i] = alpha == x[i].
 *
 * $Id: xseq.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void iseq (int_t n, int_t alpha,
	   const int_t* x, int_t incx,
	         int_t* y, int_t incy)
{
  register int_t  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] == alpha) ? 1 : 0;
}
