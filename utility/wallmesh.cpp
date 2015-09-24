///////////////////////////////////////////////////////////////////////////////
// wallmesh.C
//
// Copyright (c) 2004 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
// SYNOPSIS
// --------
// This is a filter for the mesh output of meshpr: reproduce on cout
// only the points that are fall on "wall" surfaces, as defined in the
// session file.  Read from stdin or optionally from a named file.
// "Obviously" the session file used should match that used by meshpr.
//
// USAGE
// -----
// wallmesh [options] session [mesh.file]
// options:
// -h ... print this message
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

static char RCS[] = "$Id: wallmesh.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>

static char prog[] = "wallmesh";

static void getargs    (int, char**, char*&, istream*&);
static void readMesh   (istream&,
			vector<real_t>&, vector<real_t>&, vector<real_t>&);
static void printWalls (int_t, int_t, int_t, BCmgr*, vector<Element*>&,
			vector<real_t>&, vector<real_t>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char*            session;
  istream*         meshfile;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  vector<Element*> elmt;
  int_t            i, nel, np, nz;
  vector<real_t>   x, y, z;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, meshfile);

  file = new FEML (session);
  mesh = new Mesh (file);
  nel  = mesh -> nEl();
  np   = Femlib::ivalue ("N_P");
  nz   = Femlib::ivalue ("N_Z");
  Femlib::ivalue ("NEL", nel);

  Geometry::set (np, nz, nel, Geometry::Cylindrical);

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  bman = new BCmgr (file, elmt);

  if (bman -> nWall() == 0) message (prog, "no walls found", ERROR);

  readMesh (*meshfile, x, y, z);
  
  printWalls (np, nz, nel, bman, elmt, x, y);

  if (nz > 1) for (i = 0; i <= nz; i++) cout << z[i] << endl;

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc,
		     char**    argv,
		     char*&    sess,
		     istream*& mesh)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: wallmesh [options] session [mesh.file]\n"
                 "options:\n"
                 "  -h ... display this message\n";
  char err[StrMax], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h': cerr << usage; exit (EXIT_SUCCESS); break;
    default:
      sprintf (err, "illegal option: %c\n", c);
      message (prog, err, ERROR); break;
    }

  switch (argc) {
  case 1:
    sess = argv[0];
    mesh = &cin;
    break;
  case 2:
    sess = argv[0];
    mesh = new ifstream (argv[1]);
    if (mesh -> bad()) message (prog, "couldn't open mesh file", ERROR);
    break;
  default:
    cerr << usage;
    exit (EXIT_FAILURE);
    break;
  }  
}


static void readMesh (istream&        file,
		      vector<real_t>& x   ,
		      vector<real_t>& y   ,
		      vector<real_t>& z   )
// ---------------------------------------------------------------------------
// Here we read in meshpr output and do a few rudimentary checks that
// it conforms with the current session file's information.
// ---------------------------------------------------------------------------
{
  char  buf[StrMax];
  int_t i, nr, ns, nz, nel, ntot;

  file.getline (buf, StrMax);
  if (!strstr (buf, "NR NS NZ NEL"))
    message (prog, "input not a mesh file", ERROR);

  string s (buf);
  istringstream ss (buf);
  ss  >> nr >> ns >> nz >> nel;

  if (nr  != Femlib::ivalue ("N_P") ||
      nz  != Femlib::ivalue ("N_Z") ||
      nel != Femlib::ivalue ("NEL") ||
      nr  != ns)
    message (prog, "input mesh does not match session file", ERROR);

  ntot = nr * ns * nel;

  x.resize (ntot);
  y.resize (ntot);
  z.resize ((nz > 1) ? nz+1 : 0);

  for (i = 0; i < ntot; i++) file >> x[i] >> y[i];

  if (nz > 1) for (i = 0; i <= nz; i++) file >> z[i];

  if (!file) message (prog, "reached end of mesh file prematurely", ERROR);
}


static void printWalls (int_t             np  ,
			int_t             nz  ,
			int_t             nel ,
			BCmgr*            bman,
			vector<Element*>& elmt,
			vector<real_t>&   x   ,
			vector<real_t>&   y   )
// ---------------------------------------------------------------------------
// This is the guts of the program.  Go through the list of surfaces,
// and when we find one that's a wall, we output the appropriate data
// from internal mesh storage.
// ---------------------------------------------------------------------------
{
  const int_t        np2   = np * np;
  const int_t        Nedge = bman -> nBCedges();
  const int_t        Nwall = bman -> nWall();
  vector<BCtriple*>& edge  = bman -> getBCedges();
  int_t              i, j, k, s;
  vector<real_t>     xs (np), ys (np);

  cout << np << " 1 " << nz << " " << Nwall << " NR NS NZ NEL" << endl;
  
  for (i = 0; i < Nedge; i++)
    if (strstr (bman -> groupInfo (edge[i] -> group), "wall")) {
      k = edge[i] -> elmt;
      s = edge[i] -> side;
      elmt[k] -> sideGet (s, &x[k*np2], &xs[0]);
      elmt[k] -> sideGet (s, &y[k*np2], &ys[0]);
      for (j = 0; j < np; j++) cout << xs[j] << "\t" << ys[j] << endl;
    }
}
