//////////////////////////////////////////////////////////////////////////////
// dpt.C: carry out 1D forward/inverse discrete polynomial transforms.
//
// Usage: dpt -<type> [-i] [file]
// where
// -i       ... specifies inverse transform
// -<type>  ... shape function kind
//              type == l ==> Legendre polynomials
//              type == m ==> 'Modal'  polynomials
//
// Forward transform: N values are input from file; these are assumed
// to be nodal values at the GLL points on [-1, 1].  Compute the 1D DPT
// and print up the spectral coefficients.
// 
// Inverse transform: the N values are assumed to be spectral
// coefficients; compute the nodal values.
//
// $Id: dpt.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#include <cfemdef>
#include <Array.h>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>
#include <Stack.h>

static char prog[] = "dpt";

static void getargs (int       argc,
		     char**    argv,
		     int&      invt,
		     char&     ptyp,
		     istream*& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: dpt -<type> [-i] [file]\n"
    "where\n"
    "-i       ... specifies inverse transform\n"
    "-<type>  ... shape function kind\n"
    "             type == l ==> Legendre polynomials\n"
    "             type == m ==> 'Modal'  polynomials\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'i':
      invt = INVERSE;
      break;
    default:
      ptyp = *argv[0];
      break;
    }
 
  if (!(ptyp == 'l' || ptyp == 'm')) { cerr << usage; exit (EXIT_FAILURE); }

  if (argc == 1) {
    file = new ifstream (*argv);
    if (file -> bad()) message (prog, "unable to open input file", ERROR);
  } else file = &cin;
}


int loadVals (istream&        file,
	      vector<double>& val )
// ---------------------------------------------------------------------------
// Get nodal values from file, return number of values.
// ---------------------------------------------------------------------------
{
  int           ntot, num = 0;
  double        datum;
  Stack<double> data;

  while (file >> datum) { data.push (datum); num++; }

  val.setSize (ntot = num);

  while (num--) val[num] = data.pop();

  return ntot;
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  istream*       file;
  int            i, ngll, dir = FORWARD;
  char           polytype;
  vector<double> u, v;
  const double   *F, *I;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, dir, polytype, file);

  v.setSize (ngll = loadVals (*file, u));

  switch (polytype) {
  case 'l': Femlib::legTran (ngll, &F, 0, &I, 0, 0, 0); break;
  case 'm': Femlib::modTran (ngll, &F, 0, &I, 0, 0, 0); break;
  }

  switch (dir) {
  case FORWARD: Blas::mxv (F, ngll, u(), ngll, v()); break;
  case INVERSE: Blas::mxv (I, ngll, u(), ngll, v()); break;
  }

  for (i = 0; i < ngll; i++) cout << i << '\t' << v[i] << endl;

  return EXIT_SUCCESS;
}
