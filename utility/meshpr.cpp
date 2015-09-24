///////////////////////////////////////////////////////////////////////////////
// meshpr.C: utility to generate mesh nodes from mesh description file.
//
// Copyright (c) 1995 <--> $Date: 2015/09/08 21:43:00 $, Hugh Blackburn
//
// Usage: meshpr [options] file
//   options:
//   -h       ... display this message
//   -c       ... disable checking of mesh connectivity
//   -s       ... list surfaces not determined by mesh connectivity (only)
//   -v       ... set verbose output
//   -u       ... set uniform spacing [Default: GLL]
//   -3       ... produce 3D mesh output: Np*Np*Nz*Nel*(x y z)
//   -n <num> ... override element order to be num
//   -z <num> ... override number of planes to be num
//   -b <num> ... override wavenumber beta to be <num> (3D)
//
// Prism-compatible output.
//
// Note that option 's' does not print mesh node locations but instead
// lists element sides that are free from internal element
// connectivity. This option could be used to provide a default list
// of surfaces as a starting point for editing if this information is
// not yet determined. -s ==> -c.
//
// --
// This file is part of Semtex.
//
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: meshpr.cpp,v 8.3 2015/09/08 21:43:00 hmb Exp $";

#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

#include "cfemdef.h"
#include "femlib.h"
#include "utility.h"
#include "mesh.h"

static char prog[] = "meshpr";
static void getargs (int_t, char**, char*&, int_t&, bool&, bool&,
		     int_t&, int_t&, bool&, int_t&, real_t&);

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// From FEML file named on command line, generate mesh knot
// information and print up on standard output.
// ---------------------------------------------------------------------------
{
  // -- Set defaults & parse command line.

  char*  session = 0;
  int_t  verb    = 0,
         np      = 0,
         nz      = 0,
         basis   = GLJ;
  real_t beta    = -1.;
  bool   check = true, surf = false, threed = false;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, verb, check, surf, np, nz, threed, basis, beta);

  // -- Set up to read from file, initialize Femlib parsing.

  FEML feml (session);

  if   (verb) Femlib::ivalue ("VERBOSE", verb);
  if   (np)   Femlib::ivalue ("N_P",     np  );
  else  np =  Femlib::ivalue ("N_P");
  if   (nz)   Femlib::ivalue ("N_Z",     nz  );
  else  nz =  Femlib::ivalue ("N_Z");

  if (nz > 1 && beta > 0.0) Femlib::value ("BETA", beta);

  // -- Build mesh from session file information.

  Mesh M (&feml, check);

  if (surf) {		       // -- Generate listing of mesh-egde valency.
    M . assemble (surf);
  } else {		       // -- Standard functionality.
    // -- Generate mesh knots and print up.

    const int_t    NEL  = M.nEl();
    const int_t    NTOT = np * np;
    const real_t   dz   = Femlib::value ("TWOPI/BETA") / nz;
    register int_t ID, j, k;
    vector<real_t> x (np*np), y (np*np), unimesh (np);
    real_t         *mesh_r, *mesh_s;
    const real_t   *zero_r, *zero_s;
    real_t         z;

    if (!threed) cout
		   << np  << " "
		   << np  << " "
		   << nz  << " "
		   << NEL << " NR NS NZ NEL"<< endl;

    if (basis == TRZ) {
      Femlib::equispacedMesh (np, &unimesh[0]);
      zero_r = zero_s = &unimesh[0];
    } else {
      Femlib::quadrature (&zero_r, 0, 0, 0, np, GLJ, 0.0, 0.0);
      Femlib::quadrature (&zero_s, 0, 0, 0, np, GLJ, 0.0, 0.0);
    }

    if (threed) {

      // -- Print_t out x, y, z for every mesh location, in planes.

      nz = (nz > 1) ? nz : 0;
      for (k = 0; k <= nz; k++) {
	z = k * dz;
	for (ID = 0; ID < NEL; ID++) {
	  M.meshElmt (ID, np, zero_r, zero_r, &x[0], &y[0]);
	  for (j = 0; j < NTOT; j++)
	    cout << x[j] << '\t' << y[j] << '\t' << z << endl;
	}
      }

    } else {

      // -- Print_t out x-y mesh.

      std::cout.precision(16);
      for (ID = 0; ID < NEL; ID++) {
	M.meshElmt (ID, np, zero_r, zero_r, &x[0], &y[0]);
	for (j = 0; j < NTOT; j++)
	  cout << setw(20) << x[j] << setw(24) << y[j] << endl;
      }

      // -- Print_t out z-mesh.

      if (nz > 1) for (j = 0; j <= nz; j++) cout << setw(20) << j * dz << endl;
    }
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int     argc   ,
		     char**  argv   ,
		     char*&  session,
		     int_t&  verb   ,
		     bool&   check  ,
		     bool&   surf   ,
		     int_t&  np     ,
		     int_t&  nz     ,
		     bool&   threed ,
		     int_t&  basis  ,
		     real_t& beta   )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: meshpr [options] session\n"
    "options:\n"
    "  -h       ... display this message\n"
    "  -c       ... disable checking of mesh connectivity\n"
    "  -s       ... list surfaces not determined by mesh connectivity (only)\n"
    "  -v       ... set verbose output\n"
    "  -u       ... set uniform spacing [Default: GLL]\n"
    "  -3       ... produce 3D mesh output: Np*Np*Nz*Nel*(x y z)\n"
    "  -n <num> ... override number of element knots to be num\n"
    "  -z <num> ... override number of planes to be num\n"
    "  -b <num> ... override wavenumber beta to be <num> (3D)\n";
  char err[StrMax], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      for (verb = 1; *++argv[0] == 'v'; verb++);
      break;
    case 'b':
      if (*++argv[0]) beta = atof (*argv);
      else { --argc;  beta = atof (*++argv); }
      break;
    case 'c':
      check = false;
      break;
    case 's':
      surf  = true;
      check = false;
      break;
    case 'u':
      basis = TRZ;
      break;
    case '3':
      threed = true;
      break;
    case 'n':
      if (*++argv[0]) np = atoi (*argv);
      else { --argc;  np = atoi (*++argv); }
      break;
    case 'z':
      if (*++argv[0]) nz = atoi (*argv);
      else { --argc;  nz = atoi (*++argv); }
      break;
    default:
      sprintf (err, "illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}
