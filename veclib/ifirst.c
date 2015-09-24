/*****************************************************************************
 * ifirst:  index of first non-zero value in x.
 *
 * $Id: ifirst.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>


int_t ifirst (int_t n, const int_t* x, int_t incx)
{ 
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) if (x[i*incx]) return i;
  
  return -1;
}
