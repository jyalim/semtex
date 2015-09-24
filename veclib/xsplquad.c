/*****************************************************************************
 * xsplquad.c: use cubic spline coefficients to compute definite integral
 * of a piecewise polynomial representation.
 *
 * $Id: xsplquad.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

double dsplquad (const double* xa ,
		 const double* ya ,
		 const double* ya2,
		 const int_t n  ,
		 const double  xlo,
		 const double  xhi)
/* ------------------------------------------------------------------------- *
 * Return the definite integral of a cubic spline interpolant from xlo to xhi.
 *
 * Input vectors xa & ya (length n) give the knot points for the spline,
 * while ya2 supplies the cubic spline coefficients, which are the second
 * derivative of the interpolant at the knots.  Vector xa must be strictly
 * increasing, and xlo and xhi must fall within the interval xa[0] and
 * xa[n - 1].  If xlo > xhi, or if xlo or xhi are out of range, return 0.0.
 * ------------------------------------------------------------------------- */
{
  register int_t k, kk, kp, nint_t;
  register double  h, sum = 0.0;
  int_t          jlo, jhi, klo, khi;
  double           a1, a2, a3, a4, c1, c2, c3, c4, aa1, aa2, aa3, aa4;

  if (xlo >= xhi || xlo < xa[0] || xhi > xa[n - 1]) return sum;

  /* -- Bisection to locate intervals for xlo & xhi. */

  jlo = klo = 0;
  jhi = khi = n - 1;

  while (jhi - jlo > 1) {
    k = (jlo + jhi) >> 1;
    if   (xa[k] > xlo) jhi = k;
    else               jlo = k;
  }

  while (khi - klo > 1) {
    k = (klo + khi) >> 1;
    if   (xa[k] > xhi) khi = k;
    else               klo = k;
  }

  if (nint_t = klo - jlo) {	/* -- xlo & xhi are in different intervals. */

    /* -- End intervals. */

    h  = xa[jhi] - xa[jlo];
    c1 = ya2[jlo] * h / 12.0 - ya[jlo] / (h + h);
    c2 = ya2[jhi] * h / 12.0 - ya[jhi] / (h + h);
    c3 = ya2[jhi] / (24.0 * h);
    c4 = ya2[jlo] / (24.0 * h);
    a2 = xlo     - xa[jhi];
    a3 = xa[jhi] - xa[jlo];
    a4 = xlo     - xa[jlo];
    
    aa2 = a2 * a2;
    aa3 = a3 * a3;
    aa4 = a4 * a4;

    sum -= c1 * aa2 +  c2 * (aa3 - aa4);
    sum += c3 * (aa3 * aa3 - aa4 * aa4) + c4 * aa2 * aa2;

    h  = xa[khi] - xa[klo];
    c1 = ya2[klo] * h / 12.0 - ya[klo] / (h + h);
    c2 = ya2[khi] * h / 12.0 - ya[khi] / (h + h);
    c3 = ya2[khi] / (24.0 * h);
    c4 = ya2[klo] / (24.0 * h);
    a1 = xhi     - xa[khi];
    a2 = xa[klo] - xa[khi];
    a3 = xhi     - xa[klo];
    
    aa1 = a1 * a1;
    aa2 = a2 * a2;
    aa3 = a3 * a3;

    sum += c1 * (aa1 - aa2) - c2 * aa3;
    sum += c3 *  aa3 * aa3 - c4 * (aa1 * aa1 - aa2 * aa2);

    /* -- Internal intervals. */

    for (k = 1; k < nint_t; k++) {
      kk   = jlo + k;
      kp   = kk  + 1;
      h    = xa[kp] - xa[kk];
      a1   = h * h;
      c1   = ya2[kk] * h / 12.0 - ya[kk] / (h + h);
      c2   = ya2[kp] * h / 12.0 - ya[kp] / (h + h);
      c3   = ya2[kp] / (24.0 * h);
      c4   = ya2[kk] / (24.0 * h);
      sum += a1 * ((c3 + c4) * a1 - (c1 + c2));
    }

  } else {			/* -- xlo & xhi are in same interval. */

    h  = xa[jhi] - xa[jlo];
    c1 = ya2[jlo] * h / 12.0 - ya[jlo] / (h + h);
    c2 = ya2[jhi] * h / 12.0 - ya[jhi] / (h + h);
    c3 = ya2[jhi] / (24.0 * h);
    c4 = ya2[jlo] / (24.0 * h);
    a1 = xhi - xa[jhi];
    a2 = xlo - xa[jhi];
    a3 = xhi - xa[jlo];
    a4 = xlo - xa[jlo];
    
    aa1 = a1 * a1;
    aa2 = a2 * a2;
    aa3 = a3 * a3;
    aa4 = a4 * a4;

    sum += c1 * (aa1       -       aa2) - c2 * (aa3       -       aa4);
    sum += c3 * (aa3 * aa3 - aa4 * aa4) - c4 * (aa1 * aa1 - aa2 * aa2);
  }

  return sum;
}


float ssplquad (const float*  xa ,
		const float*  ya ,
		const float*  ya2,
		const int_t n  ,
		const float   xlo,
		const float   xhi)
/* ------------------------------------------------------------------------- *
 * Return the definite integral of a cubic spline interpolant from xlo to xhi.
 *
 * Input vectors xa & ya (length n) give the knot points for the spline,
 * while ya2 supplies the cubic spline coefficients, which are the second
 * derivative of the interpolant at the knots.  Vector xa must be strictly
 * increasing, and xlo and xhi must fall within the interval xa[0] and
 * xa[n - 1].  If xlo > xhi, or if xlo or xhi are out of range, return 0.0.
 * ------------------------------------------------------------------------- */
{
  register int_t k, kk, kp, nint_t;
  register float   h, sum = 0.0;
  int_t          jlo, jhi, klo, khi;
  float            a1, a2, a3, a4, c1, c2, c3, c4, aa1, aa2, aa3, aa4;

  if (xlo >= xhi || xlo < xa[0] || xhi > xa[n - 1]) return sum;

  /* -- Bisection to locate intervals for xlo & xhi. */

  jlo = klo = 0;
  jhi = khi = n - 1;

  while (jhi - jlo > 1) {
    k = (jlo + jhi) >> 1;
    if   (xa[k] > xlo) jhi = k;
    else               jlo = k;
  }

  while (khi - klo > 1) {
    k = (klo + khi) >> 1;
    if   (xa[k] > xhi) khi = k;
    else               klo = k;
  }

  if (nint_t = klo - jlo) {	/* -- xlo & xhi are in different intervals. */

    /* -- End intervals. */

    h  = xa[jhi] - xa[jlo];
    c1 = ya2[jlo] * h / 12.0 - ya[jlo] / (h + h);
    c2 = ya2[jhi] * h / 12.0 - ya[jhi] / (h + h);
    c3 = ya2[jhi] / (24.0 * h);
    c4 = ya2[jlo] / (24.0 * h);
    a2 = xlo     - xa[jhi];
    a3 = xa[jhi] - xa[jlo];
    a4 = xlo     - xa[jlo];
    
    aa2 = a2 * a2;
    aa3 = a3 * a3;
    aa4 = a4 * a4;

    sum -= c1 * aa2 +  c2 * (aa3 - aa4);
    sum += c3 * (aa3 * aa3 - aa4 * aa4) + c4 * aa2 * aa2;

    h  = xa[khi] - xa[klo];
    c1 = ya2[klo] * h / 12.0 - ya[klo] / (h + h);
    c2 = ya2[khi] * h / 12.0 - ya[khi] / (h + h);
    c3 = ya2[khi] / (24.0 * h);
    c4 = ya2[klo] / (24.0 * h);
    a1 = xhi     - xa[khi];
    a2 = xa[klo] - xa[khi];
    a3 = xhi     - xa[klo];
    
    aa1 = a1 * a1;
    aa2 = a2 * a2;
    aa3 = a3 * a3;

    sum += c1 * (aa1 - aa2) - c2 * aa3;
    sum += c3 *  aa3 * aa3 - c4 * (aa1 * aa1 - aa2 * aa2);

    /* -- Internal int_tervals. */

    for (k = 1; k < nint_t; k++) {
      kk   = jlo + k;
      kp   = kk  + 1;
      h    = xa[kp] - xa[kk];
      a1   = h * h;
      c1   = ya2[kk] * h / 12.0 - ya[kk] / (h + h);
      c2   = ya2[kp] * h / 12.0 - ya[kp] / (h + h);
      c3   = ya2[kp] / (24.0 * h);
      c4   = ya2[kk] / (24.0 * h);
      sum += a1 * ((c3 + c4) * a1 - (c1 + c2));
    }

  } else {			/* -- xlo & xhi are in same interval. */

    h  = xa[jhi] - xa[jlo];
    c1 = ya2[jlo] * h / 12.0 - ya[jlo] / (h + h);
    c2 = ya2[jhi] * h / 12.0 - ya[jhi] / (h + h);
    c3 = ya2[jhi] / (24.0 * h);
    c4 = ya2[jlo] / (24.0 * h);
    a1 = xhi - xa[jhi];
    a2 = xlo - xa[jhi];
    a3 = xhi - xa[jlo];
    a4 = xlo - xa[jlo];
    
    aa1 = a1 * a1;
    aa2 = a2 * a2;
    aa3 = a3 * a3;
    aa4 = a4 * a4;

    sum += c1 * (aa1       -       aa2) - c2 * (aa3       -       aa4);
    sum += c3 * (aa3 * aa3 - aa4 * aa4) - c4 * (aa1 * aa1 - aa2 * aa2);
  }

  return sum;
}
