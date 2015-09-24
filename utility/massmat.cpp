///////////////////////////////////////////////////////////////////////////////
// massmat.C: print out the system mass matrix value for all (2D)
// points in the mesh, element-by-element.  Ordering is the same as
// the mesh locations printed out by meshpr. If CYLINDRICAL = 1 in
// session file, mass matrix values include weighting by radius
// (y). Printout in ASCII.
//
// Copyright (c) 2013 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
// Synopsis:
// --------
// massmat [-h] session
//
// Description: 
// ----------- 
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

static char RCS[] = "$Id: massmat.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>

static char  prog[]  = "massmat";
static int_t verbose = 0;
static void  getargs  (int, char**, char*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char*              session = 0;
  int_t              NP, NEL, NPNP, i, j;
  vector<real_t>     mass;
  FEML*              F;
  Mesh*              M;
  Element*           E;
  Geometry::CoordSys space;

  // -- Initialize.

  Femlib::initialize (&argc, &argv);
  getargs            (argc, argv, session);

  // -- Set output format and precesion. Edit these if you like.

  cout << scientific;
  cout << setprecision (12);

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);

  NEL = M -> nEl();  
  NP  = Femlib::ivalue ("N_P");

  NPNP = NP * NP;

  space = (Femlib::ivalue ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (NP, 1, NEL, space);
  mass.resize   (NPNP);

  for (i = 0; i < NEL; i++) {
    E = new Element (i, NP, M);
    Veclib::fill (NPNP, 1.0, &mass[0], 1);
    E -> weight (&mass[0]);
    if (Geometry::cylindrical()) E -> mulY (&mass[0]);

    for (j = 0; j < NPNP; j++) cout << mass[j] << endl;
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: massmat [options] session\n"
    "options:\n"
    "-h ... print this message\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if      (argc == 1)   session = argv[0];
  else                  message (prog, usage, ERROR);
}
