//////////////////////////////////////////////////////////////////////////////
// filter.C: Generate/test filtering functions.
//
// Usage: filter -a <num> -N <num> -p <num> -r <num>
// where
// -N <num> ... supplies the number of points in the filter
// -r <num> ... supplies the filter rollof point     [0, 1]
// -p <num> ... supplies the filter order
// -a <num> ... supplies attenuation factor at high frequencies
//
// $Id: testBVD.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>

using namespace std;

#include <cfemdef>
#include <Array.h>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>

static void getargs (int     argc ,
		     char**  argv ,
		     real&   atten,
		     real&   roll ,
		     real&   order,
		     int&    N    )
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: filter -a <num> -p <num> -r <num> -N <num>\n"
    "where\n"
    "-a <num> ... supplies attenuation factor at high frequencies [0, 1]\n"
    "-r <num> ... supplies the filter rolloff                     [0, 1]\n"
    "-p <num> ... supplies the filter order\n"
    "-N <num> ... supplies the number of points in the filter\n";

  if (argc != 9) {
    cerr << usage;
    exit (EXIT_FAILURE);
  }
  
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'a':
      --argc;
      atten = atof (*++argv);
      break;
    case 'N':
      --argc;
      N = atoi (*++argv);
      break;
    case 'p':
      --argc;
      order = atof (*++argv);
      break;
    case 'r':
      --argc;
      roll = atof (*++argv);
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int  i, N;
  real atten, roll, order;

  getargs (argc, argv, atten, roll, order, N);

  vector<real> work (N);
  real         *filter = work();

  Femlib::erfcFilter (N-1, order, roll, atten, filter);

  for (i = 0; i < N; i++) 
    cout << i << '\t' << filter[i] << endl;
  
  return EXIT_SUCCESS;
}
