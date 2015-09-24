//////////////////////////////////////////////////////////////////////////////
// lagder.C: given a table of N values, assumed located at the GLL
// points in [-1, 1], return values of the derivative at the same points.
//
// Usage: lagder [file]
//
// $Id: lagder.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
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

static char prog[] = "lagder";

void getargs (int       argc,
	      char**    argv,
	      istream*& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: lagder [file]\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    default:
      cerr << usage; exit (EXIT_FAILURE);
      break;
    }

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
  int            i, np;
  vector<double> u, v;
  const double   **DV;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, input);

  np = loadVals (*input, v);

  u.setSize (np);

  Femlib::mesh (GLL, GLL, np, np, 0, 0, 0, &DV, 0);
  Blas::mxv    (*DV, np, v(), np, u());

  for (i = 0; i < np; i++) cout << i << '\t' << u[i] << endl;

  return EXIT_SUCCESS;
}
