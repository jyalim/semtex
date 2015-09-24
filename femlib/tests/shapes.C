//////////////////////////////////////////////////////////////////////////////
// shapes.C: print out shape functions evaluated at points equispaced
// on [-1, 1].
//
// Copyright (c) 1999 <--> $Date: 2015/04/20 11:14:14 $, Hugh Blackburn
//
// Usage: shapes -n <num> -i <num> -t <type>
// where
// -n <num>  ... number of shape functions (columns)
// -i <num>  ... number of interpolant points on [-1, 1] (rows)
// -<type>   ... shape function kind
//               type == l ==> Legendre polynomials
//               type == m ==> 'Modal'  polynomials
//               type == L ==> Lagrange interpolants through GL nodes
//
// NB: the equivalent polynomial order is ONE LESS THAN the number of
// shape functions in the case of Lagrange interpolants.
//
// $Id: shapes.C,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>

using namespace std;

#include <cfemdef.h>
#include <Array.h>
#include <veclib.h>
#include <femlib.h>
#include <blas.h>
#include <utility.h>


static void getargs (int     argc ,
		     char**  argv ,
		     int&    n    ,
		     int&    i    ,
		     char&   ptype)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: shapes -n <num> -i <num> -<type>\n"
    "where\n"
    "-n <num> ... number of shape functions (columns)\n"
    "-i <num> ... number of interpolant points on [-1, 1] (rows)\n"
    "-<type>  ... shape function kind:\n"
    "             type == l ==> Legendre polynomials\n"
    "             type == m ==> 'Modal'  polynomials\n"
    "             type == L ==> Lagrange polynomials through GL nodes\n\n"
    "NB: the equivalent polynomial order is ONE LESS THAN the number of\n"
    "shape functions in the case of Lagrange interpolants.\n";

  if (argc != 6) {
    cerr << usage;
    exit (EXIT_FAILURE);
  }
  
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'i':
      if (*++argv[0]) i = atoi (*argv);
      else { --argc;  i = atoi (*++argv); }
      break;
    case 'n':
      if (*++argv[0]) n = atoi (*argv);
      else { --argc;  n = atoi (*++argv); }
      break;
    default:
      ptype = *argv[0];
      if (!(ptype == 'l' || ptype == 'm' || ptype == 'L')) {
	cerr << usage;
	exit (EXIT_FAILURE);
      }
      break;
    }
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int  i, j, ns, ni;
  char polytype;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, ns, ni, polytype);

  double* x = dvector (0, ni-1);
  for (i = 0; i < ni; i++)
    x[i] = -1.0 + i*2.0/(ni-1);
  double** s = dmatrix (0, ni-1, 0, ns-1);


  switch (polytype) {

  case 'L': {			// -- Lagrange polynomials.

    double*  z  = dvector (0, ns-1);
    double*  w  = dvector (0, ns-1);
    double** II = dmatrix (0, ni-1, 0, ns-1);
    double** IT = dmatrix (0, ns-1, 0, ni-1);

    Femlib::GLLzw       (ns, z, w);
    Femlib::LagrangeInt (ns, z, ni, x, II, IT);

    for (j = 0; j < ns; j++)
      for (i = 0; i < ni; i++) s[i][j] = II[i][j];

  } break;
  
  case 'l':			// -- Legendre polynomials.

    for (j = 0; j < ns; j++)
      for (i = 0; i < ni; i++)
	s[i][j] = Femlib::LegendreVal (j, x[i]);
    break;
  
  case 'm':			// -- Modal polynomials.

    for (j = 0; j < ns; j++)
      for (i = 0; i < ni; i++)
	s[i][j] = Femlib::ModalVal (j, x[i]);
    break;
  }

  for (i = 0; i < ni; i++) {
    cout << x[i];
    for (j = 0; j < ns; j++)
      cout << '\t' << s[i][j];
    cout << endl;
  }

  return EXIT_SUCCESS;
}
