//////////////////////////////////////////////////////////////////////////////
// drive.C: compute solution to elliptic problem, optionally compare to
// exact solution (see getoptions(), below).
//
// Copyright (c) 1994<-->$Date: 2015/04/20 11:14:13 $, Hugh Blackburn
//
// USAGE:
// -----
// elliptic [options] session
//   options:
//   -h       ... print this message
//   -i       ... use iterative solver
//   -v[v...] ... increase verbosity level
//
// If session.frc is found, use this field file as forcing for the
// elliptic problem, otherwise use the 'forcing' string in the USER
// section; failing that, set forcing to zero.
//
// Author
// ------
// Hugh Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@eng.monash.edu.au
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
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: drive.cpp,v 8.1 2015/04/20 11:14:13 hmb Exp $";

#include <sem.h>

static char prog[] = "elliptic";
static void getargs    (int, char**, char*&);
static void getoptions (FEML*, char*&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&, AuxField*&);
static void getforcing (const char*, const char*, AuxField*);

void Helmholtz (Domain*, AuxField*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char             *session, *forcefunc = 0, *exact = 0;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  AuxField*        forcefld;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain, forcefld);

  getoptions (file, forcefunc, exact);
  getforcing (session, forcefunc, forcefld);

  domain -> restart();

  Helmholtz (domain, forcefld);

  ROOTONLY if (exact) domain -> u[0] -> errors (mesh, exact);

  domain -> dump();

  Femlib::finalize();

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  const char routine[] = "getargs";
  const char usage[] =
    "Usage: %s [options] session\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -i       ... use iterative solver\n"
    "  -v[v...] ... increase verbosity level\n";
  char buf[StrMax];
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      ROOTONLY {
	sprintf (buf, usage, prog);
	cout << buf;
      }
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      Femlib::ivalue ("ITERATIVE", static_cast<int_t>(1));
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE", Femlib::ivalue("VERBOSE")+1);
      while (*++argv[0] == 'v');
      break;
    default:
      ROOTONLY { sprintf (buf, usage, prog); cout << buf; }
      exit (EXIT_FAILURE);
      break;
    }
  
  if (argc != 1) message (routine, "no session definition file", ERROR);

  session = *argv;
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain ,
			AuxField*&        forcing)
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.

  VERBOSE cout << "Building boundary condition manager ..." << endl;

  bman = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, elmt, bman);

  VERBOSE cout << "done" << endl;

  // -- Build the forcing field.

  VERBOSE cout << "Building forcing ...";

  forcing = new AuxField(new real_t[(size_t)Geometry::nTotProc()],nz,elmt,'f');

  VERBOSE cout << "done" << endl;
}


static void getoptions (FEML*  feml ,
			char*& forcf,
			char*& exact)
// ---------------------------------------------------------------------------
// Try to load forcing function string and exact solution string from USER
// section of FEML file.  The section is not required to be present.
// 
// Expect something in the form:
// <USER>
// forcing 0
// exact   sin(TWOPI*x)*sinh(TWOPI*y)/sinh(TWOPI)
// </USER>
//
// Either or both of the two strings may be absent.
// ---------------------------------------------------------------------------
{
  char routine[] = "options";
  char s[StrMax];

  if (feml -> seek ("USER")) {
    feml -> stream().ignore (StrMax, '\n');

    while (feml -> stream() >> s) {
      if (strcmp (s, "</USER>") == 0) break;

      upperCase (s);
      if (strcmp (s, "FORCING") == 0)
	feml -> stream() >> (forcf = new char [StrMax]);
      else if (strcmp (s, "EXACT") == 0)
	feml -> stream() >> (exact = new char [StrMax]);
    }

    if (strcmp (s, "</USER>") != 0)
      message (routine, "couldn't sucessfully close <USER> section", ERROR);
  }
}


static void getforcing (const char* session  , 
			const char* forcefunc,
			AuxField*   forcefld )
// ---------------------------------------------------------------------------
// If file session.frc is found, use the contents (of the first field
// variable) to initialise forcefld, failing that use the string
// forcefunc, otherwise initialise to zero.  Finally, Fourier
// transform.
// ---------------------------------------------------------------------------
{
  const char routine[] = "getforcing";
  char       restartfile[StrMax];
  
  ROOTONLY cout << "-- Forcing          : ";
  ifstream file (strcat (strcpy (restartfile, session), ".frc"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << restartfile;
      cout.flush();
    }

    // -- Strip header and check the data conforms.

    int_t         np, nz, nel, ntot, nfields;
    int_t         npchk,  nzchk, nelchk, swab = 0;
    char          s[StrMax], f[StrMax];

    if (file.getline(s, StrMax).eof())
      message (routine, "forcing file is empty", ERROR);

    file.getline(s,StrMax).getline(s,StrMax);

    string        ss(s);
    istringstream sss (ss);

    sss >> np    >> np    >> nz    >> nel;

    forcefld -> describe (f);

    sss.clear();
    sss.str (ss = f);
    sss >> npchk >> npchk >> nzchk >> nelchk;

    if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
    if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
    if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
    ntot = np * np * nz * nel;
    if (ntot != Geometry::nTot())
      message (routine, "declared sizes mismatch", ERROR);

    file.getline(s,StrMax).getline(s,StrMax);
    file.getline(s,StrMax).getline(s,StrMax);
    file.getline(s,StrMax).getline(s,StrMax);

    nfields = 0; while (isalpha (s[nfields])) nfields++;
    if (!nfields) message (routine, "no fields declared in forcing", ERROR);

    file.getline (s, StrMax);
    Veclib::describeFormat (f);

    if (!strstr (s, "binary"))
      message (routine, "input field file not in binary format", ERROR);
  
    if (!strstr (s, "endian"))
      message (routine, "input field file in unknown binary format", WARNING);
    else {
      swab = ((strstr (s, "big") && strstr (f, "little")) ||
	      (strstr (f, "big") && strstr (s, "little")) );
      ROOTONLY {
	if (swab) cout << " (byte-swapping)";
	cout.flush();
      }
    }

    // -- Read data, byteswap if required. 

    file >> *forcefld;
    if (swab) forcefld -> reverse();
    file.close();

    forcefld -> transform (FORWARD);
    ROOTONLY forcefld -> zeroNyquist();

  } else if (forcefunc) {
    ROOTONLY cout << "set to string: " << forcefunc;
    (*forcefld = forcefunc).transform (FORWARD);
    ROOTONLY forcefld -> zeroNyquist();

  } else {
    ROOTONLY cout << "set to zero";
    *forcefld = 0.0;
  }

  ROOTONLY cout << endl;
}
