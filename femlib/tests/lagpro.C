//////////////////////////////////////////////////////////////////////////////
// lagpro.C: given a table of N values, assumed located at the GLL
// points in [-1, 1], return values of the original interpolating
// polynomial, order N-1, at I GLL points; also return the values of
// the order I-1 interpolating polynomial for the new I GLL points
// evaluated at the original N GLL points.
//
// Usage: lagpro -i <num> [file]
// where
// -i <num> ... number of new GLL points on [-1, 1]
//
// $Id: lagpro.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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

static char prog[] = "lagpro";

void getargs (int       argc,
	      char**    argv,
	      int&      nint,
	      istream*& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: lagpro -i <num> [file]\n"
    "where\n"
    "-i <num> ... number of new GLL points on [-1, 1]\n";

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
  istream*       input;
  int            i, nold, nnew = 0;
  vector<double> u, v;
  const double   **IF, **IB;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, nnew, input);

  nold = loadVals (*input, v);

  u.setSize (nnew);

  if (nold == nnew)
    Veclib::copy (nold, v(), 1, u(), 1);
  else {
    Femlib::mesh (GLL, GLL, nold, nnew, 0, &IF, 0, 0, 0);
    Femlib::mesh (GLL, GLL, nnew, nold, 0, &IB, 0, 0, 0);
    Blas::mxv    (*IF, nnew, v(), nold, u());
    Blas::mxv    (*IB, nold, u(), nnew, v());
  }

  for (i = 0; i < nnew; i++) cout << i << '\t' << u[i] << endl;
  cout << endl;
  for (i = 0; i < nold; i++) cout << i << '\t' << v[i] << endl;

  return EXIT_SUCCESS;
}
