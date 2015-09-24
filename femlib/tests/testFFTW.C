//////////////////////////////////////////////////////////////////////////////
// testFFTW.C: test run FFTW routines against FFTPACK.
//
// $Id: testFFTW.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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

#include "/usr/local/include/rfftw.h"

int main ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int      nz   = 64;
  const int      np   = 20000;
  const int      ntot = nz * np;
  int            i, j, nfax, ifax[64];
  int            ip, iq, ir;
  double         trig[512];
  double         *x, *y, *w;
  vector<double> work  (3*nz*np);
  vector<double> work2 (3 * nz + 15);
  double*        tmp  = work2();
  double*        Wtab = tmp + nz;
  double*        ptr;
  double         stime, ftime;

  x = work();
  y = x + ntot;
  w = y + ntot;

  rfftw_plan fplan = rfftw_create_plan_specific 
    (nz, FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_OUT_OF_PLACE, x, np, y, np);

  Veclib::vrandom (ntot, x, 1);
  Veclib::copy    (ntot, x, 1, y, 1);

  // -- Echo original data.
  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }
  */
  // -- Compare results on forward transform.
  cout << "-- Forward transform" << endl;

  Femlib::rffti (nz, Wtab);
  ptr = x;

  stime = dclock();
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
  ftime = dclock();
  Blas::scal (ntot, 1.0 / nz, x, 1);
  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }
  */
  cout << "FFTPACK: " << ftime - stime << " seconds" << endl;

  stime = dclock();
  rfftw (fplan, np, y, np, 1, w, np, 1);
  Veclib::copy (np, w, 1, y, 1);
  Veclib::copy (np, w+np*(nz>>1), 1, y+np, 1);
  for (i = 1; i < (nz>>1); i++) {
    Veclib::copy (np, w+np*i,      1, y+np*(i+i),   1);
    Veclib::copy (np, w+np*(nz-i), 1, y+np*(i+i+1), 1);
  }
  ftime = dclock();
  Blas::scal    (ntot, 1.0/nz, y, 1);
  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << y[i*np+j];
    cout << endl;
  }
  */
  cout << "FFTW:    " << ftime - stime << " seconds" << endl;

  Veclib::vsub (ntot, x, 1, y, 1, w, 1);
  cout << "   Maximum difference: " << w[Blas::iamax (ntot, w, 1)] << endl;

#if 0 // -- The remainder of this hasn't been converted to use FFTW yet.
  cout << "-- Repeat:" << endl;

  Veclib::vrandom (ntot, x, 1);
  Veclib::copy    (ntot, x, 1, y, 1);

  stime = dclock();
  Femlib::mpfft (y, w, np, nz, ip, iq, ir, trig, +1);
  ftime = dclock();
  Blas::scal    (ntot, 1.0/nz, y, 1);
  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << y[i*np+j];
    cout << endl;
  }
  */
  cout << "Temperton:  " << ftime - stime << " seconds" << endl;

  ptr = x;  

  stime = dclock();
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
  ftime = dclock();
  Blas::scal (ntot, 1.0 / nz, x, 1);

  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }
  */

  Veclib::vsub (ntot, x, 1, y, 1, w, 1);
  cout << "FFTPACK:   " << ftime - stime << " seconds" << endl;
  cout << "   Maximum difference: " << w[Blas::iamax (ntot, w, 1)] << endl;

  // ------------------------------------------------------------------------
  // -- Compare results on inverse transform.

  cout << "-- Back transform:" << endl;

  ptr = x;
  if (nz & 1)
    for (i = 0; i < np; i++, ptr++) {
      Veclib::copy  (nz, ptr, np, tmp, 1);
      Femlib::rfftb (nz, tmp, Wtab);
      Veclib::copy  (nz, tmp, 1, ptr, np);
    }
  else
    for (i = 0; i < np; i++, ptr++) {
      tmp[nz - 1] = ptr[np];
      tmp[0]      = ptr[0];
      Veclib::copy  (nz - 2, ptr + 2 * np, np, tmp + 1, 1);
      Femlib::rfftb (nz, tmp, Wtab);
      Veclib::copy  (nz, tmp, 1, ptr, np);
    }

  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }
  cout << endl;
  */

  Femlib::mpfft (y, w, np, nz, ip, iq, ir, trig, -1);

  Veclib::vsub (ntot, x, 1, y, 1, w, 1);
  cout << "   Maximum difference: " << w[Blas::iamax (ntot, w, 1)] << endl;

  /*
  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << y[i*np+j];
    cout << endl;
  }
  */

  // ------------------------------------------------------------------------
  // -- Multiple calls.

  cout << "-- Multiple calls" << endl;

  for (j = 0; j < 16; j++) {
    Femlib::mpfft (y, w, np, nz, ip, iq, ir, trig, +1);
    Blas::scal    (ntot, 1.0 / nz, y, 1);
    Femlib::mpfft (y, w, np, nz, ip, iq, ir, trig, -1);
  }

  Veclib::vsub (ntot, x, 1, y, 1, w, 1);
  cout << "   Maximum difference: " << w[Blas::iamax (ntot, w, 1)] << endl;
#endif
#if 0
  // ------------------------------------------------------------------------
  // -- Test front end routines in fourier.c

  cout << "-- DFTr front end" << endl;

  Femlib::mrcft (y, np, nz, w, nfax, ifax, trig, +1);
  Blas::scal    (ntot, 1.0 / nz, y, 1);
  Femlib::DFTr  (x, nz, np, +1);

  Veclib::vsub (ntot, x, 1, y, 1, w, 1);
  cout << "   Maximum difference: " << w[Blas::iamax (ntot, w, 1)] << endl;
#endif
  return EXIT_SUCCESS;
}






