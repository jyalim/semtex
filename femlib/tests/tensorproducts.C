//////////////////////////////////////////////////////////////////////////////
// tensorproducts.C: test utility for Femlib::tpr2d, and Femlib:grad2.
//
// $Id: tensorproducts.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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

#define TPR2D 0
#define GRAD2 1


int main()
{
  const int NR  = 3;
  const int NS  = 2;
  const int NEL = 1;

#if TPR2D

  int i, j, k;
  vector<real> work (NR*NR+NS*NS+3*NR*NS);
  real*        y  = &work[0];
  real*        x  = y  + NR*NR;
  real*        dv = x  + NS*NS;
  real*        dt = dv + NR*NS;
  real*        t  = dt + NR*NS;

  Veclib::vrandom (NR*NS, dv, 1);
  Veclib::vrandom (NS*NS,  x, 1);
  Veclib::vrandom (NS*NR, dt, 1);

  Femlib::tpr2d (x, y, t, dv, dt, NR, NS, NEL);

  cout << "--DV" << endl;

  for (i = 0; i < NR; i++) {
    for (j = 0; j < NS; j++) {
      if (j == 0) cout << dv[Veclib::row_major (i, j, NS)];
      else        cout << '\t' <<  dv[Veclib::row_major (i, j, NS)];
    }
    cout << endl;
  }

  cout << "-- X" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NS; j++) {
      if (j == 0) cout << x[Veclib::row_major (i, j, NS)];
      else        cout << '\t' <<  x[Veclib::row_major (i, j, NS)];
    }
    cout << endl;
  }

  cout << "--DT" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << dt[Veclib::row_major (i, j, NR)];
      else        cout << '\t' << dt[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  cout << endl;

  cout << "-- Y = DV X DT" << endl;

  for (i = 0; i < NR; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << y[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  y[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  Blas::mxm (x,  NS, dt, NS, t, NR);
  Blas::mxm (dv, NR, t,  NS, y, NR);

  cout << endl;

  cout << "-- Y (check using mxm)" << endl;

  for (i = 0; i < NR; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << y[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  y[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

#elif GRAD2

  int i, j, k;
  vector<real> work (NR*NR+NS*NS+4*NR*NS);
  real*        x  = &work[0];
  real*        y  = x  + NS*NR;
  real*        ur = y  + NS*NR;
  real*        us = ur + NS*NR;
  real*        dv = us + NS*NR;
  real*        dt = dv + NS*NS;

  Veclib::vrandom (NS*NR,  x, 1);
  Veclib::vrandom (NS*NR,  y, 1);
  Veclib::vrandom (NS*NS, dv, 1);
  Veclib::vrandom (NR*NR, dt, 1);
  Veclib::zero    (NS*NR, ur, 1);
  Veclib::zero    (NS*NR, us, 1);

  cout << "-- X" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << x[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  x[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  cout << "-- Y" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << y[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  x[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  cout << "--DV" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NS; j++) {
      if (j == 0) cout << dv[Veclib::row_major (i, j, NS)];
      else        cout << '\t' <<  dv[Veclib::row_major (i, j, NS)];
    }
    cout << endl;
  }

  cout << "--DT" << endl;

  for (i = 0; i < NR; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << dt[Veclib::row_major (i, j, NR)];
      else        cout << '\t' << dt[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  Femlib::grad2 (x, y, ur, us, dv, dt, NR, NS, NEL);

  cout << endl;

  cout << "-- UR = X DT" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << ur[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  ur[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  cout << "-- US = DV Y" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << us[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  us[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  Blas::mxm (x,  NS, dt, NR, ur, NR);
  Blas::mxm (dv, NS, y,  NS, us, NR);

  cout << endl;

  cout << "-- UR (check using mxm)" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << ur[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  ur[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }

  cout << "-- US (check using mxm)" << endl;

  for (i = 0; i < NS; i++) {
    for (j = 0; j < NR; j++) {
      if (j == 0) cout << us[Veclib::row_major (i, j, NR)];
      else        cout << '\t' <<  us[Veclib::row_major (i, j, NR)];
    }
    cout << endl;
  }
#endif

  return EXIT_SUCCESS;
}
