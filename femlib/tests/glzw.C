//////////////////////////////////////////////////////////////////////////////
// glzw.C: print out Gauss--Lobatto nodes and weights on [-1, 1].
//
// Copyright (c) 1999 <--> $Date: 2015/04/20 11:14:14 $, Hugh Blackburn
//
// Usage: glzw -N <num> [-w <a,b>]
// where
// -N <num> ... supplies the number of GL nodes.
// -w <a,b> ... a and b are the exponents in (1+x)^a*(1-x)^b (default: 0,0)
//
// The default weights 0,0 generate quadrature points for Legendre polys.
//
// $Id: glzw.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <cstring>

#include <iostream>
#include <iomanip>

using namespace std;

#include <cfemdef.h>
#include <femlib.h>
#include <utility.h>

static char prog[] = "glzw";


static void getargs (int     argc,
		     char**  argv,
		     int_t&  N   ,
		     real_t& a   ,
		     real_t& b   )
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: glzw -N <num> [-w <a,b>]\n"
    "where\n"
    "-N <num> ... supplies the number of GL nodes\n"
    "-w <a,b> ... a and b are exponents in (1+x)^a*(1-x)^b (default: 0,0)\n";
  char *tok, *pspec;

  if (argc < 3) {
    cerr << usage;
    exit (EXIT_FAILURE);
  }
  
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'N':
      if (*++argv[0]) N = atoi (*argv);
      else { --argc;  N = atoi (*++argv); }
      break;
    case 'w':
      if (*++argv[0])
	pspec = *argv;
      else {
	--argc;
	pspec = *++argv;
      }
      if (tok = strtok (pspec, ",")) {
	a = atof (tok);
      } else {
	message (prog, "couldn't parse number a from string", ERROR);
      }
      if (tok = strtok (0, ",")) {
	b = atof (tok);
      } else {
	message (prog, "couldn't parse number b from string", ERROR);
      }
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
  int_t        i, N;
  real_t       alpha = 0.0, beta = 0.0;
  const real_t *z, *w;

  getargs (argc, argv, N, alpha, beta);

  Femlib::quadrature (&z, &w, 0, 0, N, LL, alpha, beta);

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);
  for (i = 0; i < N; i++) 
    cout << i << '\t' << z[i] << '\t' << w[i] << endl;
  
  return EXIT_SUCCESS;
}
