//////////////////////////////////////////////////////////////////////////////
// testsine.C: test FFT routines on simple cos and sin waves.
//
// $Id: testsine.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>

using namespace std;

#include <Array.h>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>


int main ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int      nz   = 32;
  const int      np   = 2;
  const int      ntot = nz * np;
  int            i, j, nfax, ifax[64];
  int            ip, iq, ir, ipqr2;
  double         trig[512];
  double         *x, *y, *z, *w;
  vector<double> work  (4*nz*np);
  vector<double> work2 (3 * nz + 15);
  double*        tmp  = work2();
  double*        Wtab = tmp + nz;
  double*        ptr;

  x = work();
  y = x + ntot;
  z = y + ntot;
  w = z + ntot;

// Veclib::vrandom (ntot, x, 1);

  for (i = 0; i < nz; i++) {
    x[i*np  ] = cos(TWOPI*i/nz);
    x[i*np+1] = sin(TWOPI*i/nz);
  }
  Veclib::copy (ntot, x, 1, y, 1);
  Veclib::copy (ntot, x, 1, z, 1);

  // -- Echo original data.

  cout << "-- Original data:" << endl;

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }

  cout << endl << "-- FFTPACK:   " << endl;

  Femlib::rffti (nz, Wtab);
  ptr = x;

  if (nz & 1)			// -- nz is odd.
    for (i = 0; i < np; i++, ptr++) {
      Veclib::copy  (nz, ptr, np, tmp, 1);
      Femlib::rfftf (nz, tmp, Wtab);
      Veclib::copy  (nz, tmp, 1, ptr, np);
    }
  else
    for (i = 0; i < np; i++, ptr++) {
      Veclib::copy  (nz, ptr, np, tmp, 1);
      Femlib::rfftf (nz, tmp, Wtab);
      Veclib::copy  (nz - 2, tmp + 1, 1, ptr + 2 * np, np);
      ptr[0]  = tmp[0];
      ptr[np] = tmp[nz - 1];
    }
  Blas::scal (ntot, 1.0 / nz, x, 1);

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }

  cout << endl << "-- Temperton:  " << endl;

  Femlib::primes235 (i = nz, ip, iq, ir, ipqr2);
  Femlib::setpf     (trig, nz, ip, iq, ir);

  Femlib::mpfft (y, w, np, nz, ip, iq, ir, trig, +1);
  Blas::scal    (ntot, 1.0/nz, y, 1);

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << y[i*np+j];
    cout << endl;
  }

  cout << endl << "-- DFTr front end:" << endl;

  Femlib::DFTr  (z, nz, np, +1);

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << z[i*np+j];
    cout << endl;
  }

  return EXIT_SUCCESS;
}
