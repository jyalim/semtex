/*****************************************************************************
 * xbrev:  byte-reversal routines.
 *
 * $Id: xbrev.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

int_t iformat (void)
/* ------------------------------------------------------------------------- *
 * Return 1 if machine floating-point format is IEEE little-endian,
 * 0 if IEEE big-endian, -1 for unrecognized format.
 * ------------------------------------------------------------------------- */
{
  union { float  f; int i;    unsigned char c[4]; } v;
  union { double d; int i[2]; unsigned char c[8]; } u;
  int_t reverse = (-1);
  u.d = 3;
  v.i = 3;
  if      (u.c[0] == 64 && u.c[1] == 8 && v.c[3] == 3) reverse = 0;
  else if (u.c[7] == 64 && u.c[6] == 8 && v.c[0] == 3) reverse = 1;

  return (reverse);
}


void format (char* s)
/* ------------------------------------------------------------------------- *
 * Fill s with a string describing machine's floating-point storage format.
 * ------------------------------------------------------------------------- */
{
  switch (iformat ()) {
  case -1:
    sprintf (s, "unknown");
    break;
  case 1:
    sprintf (s, "IEEE little-endian");
    break;
  case 0: default:
    sprintf (s, "IEEE big-endian");
    break;
  }
}


void dbrev (int_t n,
	    const double* x, int_t incx,
	          double* y, int_t incy)
{
  register char    *cx, *cy, d;
  register int_t i, j;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    cx = (char*) (x + i * incx);
    cy = (char*) (y + i * incy);
    for (j = 0; j < 4; j++) { 
      d       = cx[j];
      cy[j]   = cx[7-j]; 
      cy[7-j] = d;
    }
  }
}


void ibrev (int_t n,
	    const int_t* x, int_t incx,
	          int_t* y, int_t incy)
{
  register char    *cx, *cy, d;
  register int_t i, j;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    cx = (char*) (x + i * incx);
    cy = (char*) (y + i * incy);
    for (j = 0; j < 2; j++) { 
      d       = cx[j];
      cy[j]   = cx[3-j]; 
      cy[3-j] = d;
    }
  }
}


void sbrev (int_t n,
	    const float* x, int_t incx,
	          float* y, int_t incy)
{
  register char    *cx, *cy, d;
  register int_t i, j;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    cx = (char*) (x + i * incx);
    cy = (char*) (y + i * incy);
    for (j = 0; j < 2; j++) { 
      d       = cx[j];
      cy[j]   = cx[3-j]; 
      cy[3-j] = d;
    }
  }
}

