//////////////////////////////////////////////////////////////////////////////
// testDLT3.C: 2D Discrete Legendre Transforms with filtering.
//
// $Id: testDLT3.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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

#define SIZE 6


int main ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  int            i, j, k, l, p, q;
  const int      np = SIZE, np2 = np * np, np4 = np2 * np2;
  double         ci, *x, *y, *z, *t, *FW, *FT, *BW, *BT, *mas, *sam;
  double         ts, tf;
  const double   *w, *tab;

  x   = new double [8*np2+np*(np+1)+2*np];
  y   = x   + np2;
  z   = y   + np2;
  t   = z   + np2;
  FW  = t   + np2;
  FT  = FW  + np2;
  BW  = FT  + np2;
  BT  = BW  + np2;
  mas = BT  + np2;
  sam = mas + np;
  tab = sam + np;

  // -- Get GLL grid points & weights, prepare table of Legendre polys.

  Femlib::legCoef (np, &tab);
  Femlib::quad    (LL, np, np, 0, 0, &w, 0, 0, 0, 0);

  // -- Create filter coefficients and inverse.

  Femlib::erfcFilter (np-1, 4, (integer) (0.5 * (np - 1)), 0.75, mas);
  for (i = 0; i < np; i++) sam[i] = 1.0 / mas[i];

  // -- Create forward & inverse tensor-product transform matrices.

  for (i = 0; i < np; i++) {
    ci = tab[Veclib::row_major (np, i, np)];
    for (j = 0; j < np; j++) {
      FW[Veclib::row_major (i, j, np)] =
	mas[i] * ci * w[j] * tab[Veclib::row_major (i, j, np)];
      BW[Veclib::row_major (i, j, np)] =
	sam[j] * tab[Veclib::row_major (j, i, np)];
    }
  }

  // -- And their transposes.

  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      FT[Veclib::row_major (i, j, np)] = FW[Veclib::row_major (j, i, np)];
      BT[Veclib::row_major (i, j, np)] = BW[Veclib::row_major (j, i, np)];
    }
  }

  // -- Create original data.

  Veclib::vrandom (np2, x, 1);

  // -- Forward DLT.

  Blas::mxm (FW, np, x,  np, t, np);
  Blas::mxm (t,  np, FT, np, y, np);


  // -- Inverse DLT.

  Blas::mxm (BW, np, y,  np, t, np);
  Blas::mxm (t,  np, BT, np, z, np);

  // -- Print everything up.

  cout << "Problem size: " << np << " x " << np << endl;

  for (i = 0; i < np2; i++)
    cout << setw(14) << x[i] << setw(14) << y[i] << setw(14) << z[i] << endl;

  return (EXIT_SUCCESS);
}
