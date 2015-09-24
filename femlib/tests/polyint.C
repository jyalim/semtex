//////////////////////////////////////////////////////////////////////////////
// polyint.C: given a table of N (x, y) pairs, compute the
// coefficients of the interpolating polynomial of order N-1, and
// print up the values of the interpolation polynomial at I locations,
// equispaced on the interval spanned by the x values. The x values do
// not have to be in sorted order.
//
// Usage: polyint -i <num> [file]
// where
// -i <num> ... number of interpolant points on [xmin, xmax]
//
// $Id: polyint.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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
#include <cnr77>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>

static char prog[] = "polyint";

void getargs (int       argc,
	      char**    argv,
	      int&      nint,
	      istream*& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: polyint -i <num> [file]\n"
    "where\n"
    "-i <num> ... number of interpolant points on [xmin, xmax]\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'i':
      if (*++argv[0]) nint = atoi (*argv);
      else { --argc;  nint = atoi (*++argv); }
      break;
    default:
      cerr << usage; exit (EXIT_FAILURE);
      break;
    }

  if (!nint) { cerr << usage; exit (EXIT_FAILURE); }

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
  int            i, j, nknot, nint = 0;
  vector<double> x, y, cof;
  double         xc, yc, xmin = FLT_MAX, xmax = -FLT_MAX, dy;

  cerr.precision (8);
  cerr.setf (ios::fixed, ios::floatfield);
  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, nint, input);

  while (*input >> xc >> yc) {
    x.insert (x.end(), xc);
    y.insert (y.end(), yc);
    xmin = min (xc, xmin);
    xmax = max (xc, xmax);
  }

  cof.resize (nknot = x.size());
  
  Recipes::polcoe (&x[0], &y[0], nknot, &cof[0]);

  for (i = 0; i < nknot; i++)
    cerr << cof[i] << endl;

  for (i = 0; i < nint; i++) {
    xc = i * (xmax - xmin)/(nint - 1) + xmin;
    Recipes::polint (&x[0], &y[0], nknot, xc, yc, dy);
    cout << xc << '\t' << yc << endl;
  }

  return EXIT_SUCCESS;
}
