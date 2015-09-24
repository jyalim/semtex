//////////////////////////////////////////////////////////////////////////////
// testDLT.C: test run 1D Discrete Legendre Transforms.
//
// $Id: testDLT.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>

#include <iostream>
#include <iomanip>

using namespace std;

#include <Array.h>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>

static double pnleg (const double, const int);


int main ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  int            i, j, k;
  const int      np = 9, nm = np - 1;
  vector<double> work (3*np+np*(np+1));
  double         ck, *x, *y, *z, *tab;
  const double   *g, *w;

  x   = work();
  y   = x + np;
  z   = y + np;
  tab = z + np;

  // -- Get GLL grid points & weights, prepare table of Legendre polys.

  Femlib::quad (LL, np, np, &g, 0, &w, 0, 0, 0, 0);

  for (i = 0; i < np; i++) {
    tab[Veclib::row_major (np, i, np)] = (i < nm) ?  0.5*(i+i+1) : 0.5*nm;
    for (j = 0; j < np; j++) 
      tab[Veclib::row_major (i, j, np)] = pnleg (g[j], i);
  }

  // -- Create original data.

  Veclib::vrandom (np, x, 1);

  // -- Forward DLT.

  for (k = 0; k < np; k++) {
    ck = (k < nm) ? 0.5 * (k + k + 1) : 0.5 * nm;
    y[k] = 0.0;
    for (j = 0; j < np; j++)
//      y[k] += x[j] * pnleg (g[j], k) * w[j];
      y[k] += x[j] * tab[Veclib::row_major (k, j, np)] * w[j];
//    y[k] *= ck;
    y[k] *= tab[Veclib::row_major (np, k, np)];
  }

  // -- Inverse DLT.

  for (j = 0; j < np; j++) {
    z[j] = 0.0;
    for (k = 0; k < np; k++)
//      z[j] += y[k] * pnleg (g[j], k);
      z[j] += y[k] * tab[Veclib::row_major (k, j, np)];
  }

  // -- Print everything up.

  for (i = 0; i < np; i++)
    cout << setw(14) << x[i] << setw(14) << y[i] << setw(14) << z[i] << endl;

  return (EXIT_SUCCESS);
}


static double pnleg (const double z,
		     const int    n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the nth order Legendre polynomial at z, based on the
 * recursion formula for Legendre polynomials.
 * ------------------------------------------------------------------------- */
{
  register integer k;
  register double  dk, p1, p2, p3;

  if (n == 0) return 1.0;

  p1 = 1.0;
  p3 = p2 = z;

  for (k = 1; k < n; k++) {
    dk = (double) k;
    p3 = ((2.0*dk + 1.0)*z*p2 - dk*p1) / (dk + 1.0);
    p1 = p2;
    p2 = p3;
  }
 
  return p3;
}
