//////////////////////////////////////////////////////////////////////////////
// lagint.C: given a table of N values, and their associated knots (2N
// values in total), return I values of the N-1 polynomial interpolant
// for locations uniformly distributed on the interval, computed using
// the Lagrange polynomials.
//
// NB: this program can be used to compute Lagrange interpolant shape
// functions, by setting the value at one of the knot points equal to
// 1, the others to zero.
//
// Usage: lagint -I <num> [file]
// where
// -I <num> ... number of interpolant points on [-1, 1]
//
// File contains <num> 2-tuples, each supplying a knot location and
// the corresponding value there.
//
// $Id: lagint.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stack>

using namespace std;

#include <cfemdef.h>
#include <veclib.h>
#include <femlib.h>
#include <blas.h>
#include <utility.h>

static char prog[] = "lagint";

void getargs (int       argc,
	      char**    argv,
	      int_t&    nint,
	      istream*& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: lagint -I <num> [file]\n"
    "where\n"
    "-I <num> ... number of interpolant points on [-1, 1]\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'I':
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


int_t loadVals (istream&        file,
		vector<real_t>& knot,
		vector<real_t>& val )
// ---------------------------------------------------------------------------
// Get knot locations and nodal values from file, return number of values.
// ---------------------------------------------------------------------------
{
  int_t         ntot, num = 0;
  real_t        datum;
  stack<real_t> data;

  while (file >> datum) { data.push (datum); num++; }

  if (num % 2) message (prog, "need even number of (knots & values)", ERROR);

  ntot = num >>= 1;

  knot.resize (ntot);
  val. resize (ntot);

  while (num--) {
    val [num] = data.top(); data.pop();
    knot[num] = data.top(); data.pop();
  }

  return ntot;
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  istream*       input;
  int_t          i, j, ngll, nint = 0;
  vector<real_t> u, v, x, z, II, IT;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, nint, input);

  ngll = loadVals (*input, z, v);

  x.resize (nint);
  u.resize (nint);

  II.resize (nint*ngll);	// -- Lagrange interpolant operator (flat).
  IT.resize (ngll*nint); 	// -- Not used but required by LagrangeInt.

  for (i = 0; i < nint; i++) x[i] = z[0] + i * (z[ngll-1]-z[0])/(nint - 1);

  Femlib::LagrangeInt (ngll, &z[0], nint, &x[0], &II[0], &IT[0]);

  Blas::mxv (&II[0], nint, &v[0], ngll, &u[0]);

  for (i = 0; i < nint; i++) cout << x[i] << '\t' << u[i] << endl;

  return EXIT_SUCCESS;
}
