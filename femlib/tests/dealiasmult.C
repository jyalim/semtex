//////////////////////////////////////////////////////////////////////////////
// dealiasmult.C: given a table of N (u, v) pairs, assumed to lie at
// the GLL knot points in [-1, 1], compute dealaised (band-limited)
// product of the two Lagrange interpolants of the points, w. Print up
// the collocation (pointwise) multiple, and the dealiased product.
//
// Given a vectors {u}_{N x 1}, {v}_{N x 1}, we first project these
// onto the longer vectors {u^+}_{3N/2 x 1}, {v^+}_{3N/2 x 1} by
// interpolation of the Lagrange interpolants, via matrix
// multiplication.
//
//   {u^+}_{3N/2 x 1} = I^+_{3N/2 x N} {u}_{N x 1}
//   {v^+}_{3N/2 x 1} = I^+_{3N/2 x N} {v}_{N x 1}
//
// Then we form the product
//
//   {w^+} = {u^+} * {w^+} where * represents vector multiply.
//
// This product is not polluted by aliasing; the task is to project
// back to the smaller space only those polynomial coefficients which
// are appropriate to the space of {N-1}-degree polynomials. This is
// done by projecting to a space of orthogonal polynomials, filtering
// and projecting back. This is done via the discrete polynomial
// transform, with a (diagonal) filter matrix L, which has 1 as the
// diagonal entry for the first N rows, 0 in all other locations.
//
//   {w^+_f} = B L B^{-1} {w^+},
//
// where all matrices are 3N/2 x 3N/2, and the vectors are length
// 3N/2.  The Bs are matrices of basis functions (Legendre
// polynomials, given the knots and weights we are using).
//
// Finally we project back to the space of GLL Lagrange polynomials:
//
//   {w}_{N x 1} = I^-_{N x 3N/2} {w^+_f}_{3N/2 x 1}
//
// The matrices I^- B L B^{-1} may be premultiplied.
//
// Usage:  dealiasmult [-v] [file]
//
// $Id: dealiasmult.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

#include <cfemdef>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>

static char prog[] = "dealiasmult";

void getargs (int       argc,
	      char**    argv,
	      int&      verb,
	      istream*& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: dealiasmult [file]\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'v':
      verb = 1;
      break;
    default:
      cerr << usage; exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    file = new ifstream (*argv);
    if (file -> bad()) message (prog, "unable to open input file", ERROR);
  } else file = &cin;
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  istream*       input;
  int            i, j, N, Np, verb = 0;
  vector<double> u, v, w, up, vp, wp, tp;
  double         ud, vd;
  const double   **Ip, **Im, *B, *Bm;
  vector<double> L;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, verb, input);

  while (*input >> ud >> vd) {
    u.insert (u.end(), ud);
    v.insert (v.end(), vd);
  }

  w .resize (N = u.size());
  up.resize (Np = (3*N + 1) / 2);
  vp.resize (Np);
  wp.resize (Np);
  tp.resize (Np);
  
  if (verb) cerr << "N = " << N << ", 3N/2 (integer)= " << Np << endl;

  // -- Retrieve the interpolation and basis function matrices.

  Femlib::mesh    (GLL, GLL, N, Np, 0, &Ip, 0, 0, 0);
  Femlib::mesh    (GLL, GLL, Np, N, 0, &Im, 0, 0, 0);
  Femlib::legTran (Np, &B, 0, &Bm, 0, 0, 0);

  // -- Project u & v up to vp & vp.

  Blas::mxv (Ip[0], Np, &u[0], N, &up[0]);
  Blas::mxv (Ip[0], Np, &v[0], N, &vp[0]);

  // -- Multiply together pointwise.

  Veclib::vmul (Np, &up[0], 1, &vp[0], 1, &wp[0], 1);

#if 0
  // -- Set the lowpass filter matrix L.

  L.resize (sqr(Np));
  Veclib::zero (sqr(Np), &L[0], 1);
  Veclib::fill (N, 1.0, &L[0], (Np + 1));

  if (verb) {
    cout << "-- [L]:" << endl;
    for (i = 0; i < Np; i++) {
      for (j = 0; j < Np; j++) {
	if   (j == 0) cerr         << L[Veclib::row_major (i, j, Np)];
	else          cerr << '\t' << L[Veclib::row_major (i, j, Np)];
      }
      cout << endl;
    }
  }

  // -- Lowpass.

  Blas::mxv (B,     Np, &wp[0], Np, &tp[0]);
  Blas::mxv (&L[0], Np, &tp[0], Np, &wp[0]);
  Blas::mxv (Bm,    Np, &wp[0], Np, &tp[0]);
#else
  vector<double> F1(sqr(Np)), F2(sqr(Np));

  L.resize(Np);
  Veclib::fill (N, 1.0, &L[0], 1);
  Veclib::zero (Np-N, &L[N], 1);

  if (verb) for (i = 0; i < Np; i++) cerr << L[i] << endl;

  // -- Create and premultiply all the matrices.

  for (i = 0; i < Np; i++)
    Veclib::smul (Np, L[i], B + i*Np, 1, &F1[0] + i*Np, 1);
  
  if (verb) {
    for (i = 0; i < Np; i++) {
      for (j = 0; j < Np; j++) {
	if   (j == 0) cerr         << F1[Veclib::row_major (i, j, Np)];
	else          cerr << '\t' << F1[Veclib::row_major (i, j, Np)];
      }
      cout << endl;
    }
  }
  
  Blas::mxm (Bm, Np, &F1[0], Np, &F2[0], Np);
  
  // -- Apply.

  Blas::mxv (&F2[0], Np, &wp[0], Np, &tp[0]);
  
#endif
 
  // -- Project back down to original subspace.

  Blas::mxv (Im[0], N, &tp[0], Np, &w[0]);

  // -- Print up aliased and dealiased product.

  for (i = 0; i < N; i++) 
    cout << u[i]*v[i] << '\t' << w[i] << endl;

  return EXIT_SUCCESS;
}
