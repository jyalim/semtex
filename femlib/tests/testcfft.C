//////////////////////////////////////////////////////////////////////////////
// testFFT.C: test run multiple FFT routines against versions known to work.
//
// $Id: testcfft.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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


int main ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int      nz   = 16;
  const int      np   = 2;
  const int      ntot = nz * np;
  int            i, j, ip, iq, ir, nfax, ifax[64];
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

  Veclib::vrandom (ntot, x, 1);
  Veclib::copy    (ntot, x, 1, y, 1);

  // -- Echo original data.

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << x[i*np+j];
    cout << endl;
  }
  cout << endl;

  // -- Compare results on forward transform.
  cout << "-- Forward transform" << endl;

  Femlib::preft (nz>>1, nfax, ifax, trig);
  Femlib::fft1 (x, w, nz>>1, nfax, ifax, +1, trig, np);

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << w[i*np+j];
    cout << endl;
  }
  
  cout << endl;

  Femlib::setpf (trig, nz>>1, ip, iq, ir);
  Femlib::gpfa  (y, y + np, trig, np + np, 1, nz>>1, ip, iq, ir, np, -1);

  for (i = 0; i < nz; i++) {
    for (j = 0; j < np; j++)
      cout << setw(15) << y[i*np+j];
    cout << endl;
  }

  Veclib::vsub (ntot, w, 1, y, 1, w, 1);
  cout << "   Maximum difference: " << w[Blas::iamax (ntot, w, 1)] << endl;

  return EXIT_SUCCESS;
}






