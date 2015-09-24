///////////////////////////////////////////////////////////////////////////////
// mapmesh.C: a filter to map the NODES in a session file according to
// formulae supplied on the command line. Can be used e.g. to map a
// rectangular domain to a boundary-fitted shape.  It is the user's
// responsiblity to ensure that the mapping does not produce Jacobians
// that are not positive-definite (for example a reflection would
// violate this requirement).
//
// Copyright (c) 2010 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
// Usage:
// -----
// mapmesh [-x <string>] [-y <string>] session
//   -x <string> ... x <-- f(x, y), f is supplied by string
//   -y <string> ... y <-- g(x, y), g is supplied by string
//   session     ... name of a semtex session file
//
// Strings are STR_MAX (see cfemdef.h) maximum length, and acceptable to
// the semtex parser (or calc).
//
// Examples:
// --------
// mapmesh -y acos(PI*y/2.0) session
// maps the y locations of nodes: y = cos^(-1)(PI/2*y)
//
// mapmesh -x x*cos(PI/4)-y*sin(PI/4) -y x*sin(PI/4)+y*cos(PI/4) session
// rotates the mesh by 45 degrees.
//
// Files:
// -----
// Input file "session" is assumed to be a valid session file. Only
// the x, y locations of the NODES section are operated on, the
// remainder is repeated verbatim (so that any plain text file passes
// through unmodified).
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

static char RCS[] = "$Id: mapmesh.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>

static char prog[] = "mapmesh";
static void getargs (int, char**, char*&, char*&, char*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char                *session = NULL, *xmap = NULL, *ymap = NULL, err[64];
  vector<const char*> present;
  int_t               i, nKey, nNode;
  vector<int_t>       id;
  vector<real_t>      x, y, mapx, mapy;
  real_t              z;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, xmap, ymap, session);

  FEML F (session);

  nKey = F.sections (present);

  for (i = 0; i < nKey; i++)
    if (strcmp (present[i], "NODES")) {
      F.echo (cout, present[i]);
      cout << endl;
    }

  // -- Now all that is left is the NODES section, which is what we want.
  //    Read it all in.

  nNode = F.attribute ("NODES", "NUMBER");
  istream& input = F.stream ();

  id.resize   (nNode);
  x.resize    (nNode);
  y.resize    (nNode);
  mapx.resize (nNode);
  mapy.resize (nNode);

  if (nNode < 4) {
    sprintf (err, "At least 4 Nodes are needed, found %1d declared", nNode);
    message (prog, err, ERROR);
  }

  for (i = 0; i < nNode; i++) {

    while (input.peek() == '#') // -- Skip comments.
      input.ignore (StrMax, '\n');

    input >> id[i] >> x[i] >> y[i] >> z;

    if (id[i] > nNode) {
      sprintf (err, "Node ID %1d exceeds attribution (%1d)", id[i], nNode);
      message (prog, err, ERROR);
    }
  }

  // -- Now we do the mappings requested on the command line.

  if (xmap) {
    Femlib::prepVec ("x y", xmap);
    Femlib__parseVec (nNode, &x[0], &y[0], &mapx[0]);
  } else mapx = x;
  if (ymap) {
    Femlib::prepVec ("x y", ymap);
    Femlib__parseVec (nNode, &x[0], &y[0], &mapy[0]);
  } else mapy = y;

  // -- And print up the section with the mapped NODES.

  cout << "<NODES NUMBER=" << nNode << '>' << endl;

  for (i = 0; i < nNode; i++)
    cout << id[i] << '\t' << mapx[i] << '\t' << mapy[i] << "\t0" << endl;

  cout << "</NODES>" << endl;

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& xmap   ,
		     char*& ymap   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: mapmesh [options] session\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -x <string> ... x <-- f(x, y), f is supplied by string\n"
    "  -y <string> ... y <-- g(x, y), g is supplied by string\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'x':
      if (*++argv[0]) xmap = *argv;
      else { xmap = *++argv; argc--;}
      break;
    case 'y':
      if (*++argv[0]) ymap = *argv;
      else { ymap = *++argv; argc--;}
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}
