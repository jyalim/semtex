/*****************************************************************************
 * xpolint - Polynomial interpolation.
 *
 * The following function performs polyomial interpolation through a given
 * set of points using Neville's algorithm.  Given arrays xa[1..n] and
 * ya[1..n], each of length n, this routine returns a value y and error
 * estimate dy.
 *
 * Reference: Numerical Recipes, 2nd edn, pp. 108-110.
 * xpoly() provide simplified interfaces which do array offsets.
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cfemdef.h>
#include <cveclib.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

double dpoly (int_t n, double x, const double* xp, const double* yp)
{
  double value, err;

  dpolint (xp - 1, yp - 1, n, x, &value, &err);
  return value;
}


float spoly (int_t n, float  x, const float* xp, const float* yp)
{
  float value, err;

  spolint (xp - 1, yp - 1, n, x, &value, &err);
  return value;
}


void dpolint (const double *xa, const double *ya, int_t n,
	            double x,         double *y,  double *dy)
{
  register int_t ns  = 1;
  double           dif = fabs (x - xa[1]),
                   *c  = dvector (1, n),
                   *d  = dvector (1, n);
  register int_t i, m;
  double           den, dift, ho, hp, w;

  for (i = 1; i <= n; i++) {
    if ((dift = fabs (x-xa[i])) < dif) {
      ns  = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
    
  *y = ya[ns--];
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      ho = xa[i] - x;
      hp = xa[i+m] -x;
      w  = c[i+1] - d[i];
      if ((den = ho-hp) == 0.0)
 	message("dpolint()", "two x values the same (within roundoff)", ERROR);
      den  = w  / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    *y += (*dy=((ns<<1) < (n-m) ? c[ns+1] : d[ns--]));
  }

  freeDvector (c, 1);
  freeDvector (d, 1);
}


void spolint (const float *xa, const float *ya, int_t n,
	            float x,         float *y,  float *dy)
{
  register int_t ns  = 1;
#if  defined(__uxp__) || defined(_SX)
  float            dif = (float) fabs (x - xa[1]),
#else
  float            dif = fabsf(x - xa[1]),
#endif
                   *c  = svector(1, n),
                   *d  = svector(1, n);
  register int_t i, m;
  float            den, dift, ho, hp, w;

  for (i = 1; i <= n; i++) {
#if  defined(__uxp__) || defined(_SX)
    if ((dift = (float) fabs  (x - xa[1])) < dif) {
#else
    if ((dift =         fabsf (x - xa[i])) < dif) {
#endif
      ns  = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
    
  *y = ya[ns--];
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      ho = xa[i] - x;
      hp = xa[i+m] -x;
      w  = c[i+1] - d[i];
      if( (den = ho-hp) == 0.0F)
	message("spolint()", "two x values the same (within roundoff)", ERROR);
      den  = w  / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    *y += (*dy=((ns<<1) < (n-m) ? c[ns+1] : d[ns--]));
  }

  freeSvector (c, 1);
  freeSvector (d, 1);
}
