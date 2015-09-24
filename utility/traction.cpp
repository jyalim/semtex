//////////////////////////////////////////////////////////////////////////////
// traction.C: Compute tractions on wall boundaries from field file.
//
// Copyright (c) 2006 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
// USAGE:
// -----
// traction session [file]
//
// Essentially this carries out the same computation as is done during
// execution of dns, but as a standalone utility.
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

static char RCS[] = "$Id: traction.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>

static char prog[] = "traction";
static void getargs    (int, char**, char*&, istream*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&);
static bool getDump    (Domain*, istream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys         system;
  char*                      session;
  istream*                   file;
  FEML*                      F;
  Mesh*                      M;
  BCmgr*                     B;
  Domain*                    D;
  vector<Element*>           E;
  int_t                      _nwall, _nline, _npad, DIM;
  vector<real_t>             _work;

  Femlib::initialize (&argc, &argv);

  // -- Read command line.

  getargs (argc, argv, session, file);

  // -- Set up domain.

  preprocess (session, F, M, E, B, D);

  DIM = Geometry::nDim();

  // -- Loop over all dumps in field file, compute and print traction.

  while (getDump (D, *file)) {

    // -- Input was in physical space but we need Fourier.

    D -> transform (FORWARD);

    // -- Set up to compute wall shear stresses.    
    
    const int_t npr = Geometry::nProc();
    const int_t np  = Geometry::nP();
    const int_t nz  = Geometry::nZProc();

    // -- Allocate storage area: 3 = 1 normal component + 2 tangential.
    
    _nwall = B -> nWall();
    _nline = np * _nwall;
    _npad  = 3  * _nline;

    // -- Round up length for Fourier transform/exchange.

    if   (npr > 1) _npad += 2 * npr - _npad % (2 * npr);
    else           _npad += _npad % 2;

    _work.resize (_npad * nz);

    // --------------------------------------------------------------------

    const int_t    nP  = Geometry::nP();
    const int_t    nZ  = Geometry::nZ();
    int_t          i, j, k;
    real_t*        plane;
    vector<real_t> buffer (_nline);

    // -- Load the local storage area.

    Veclib::zero (_work.size(), &_work[0], 1);

    if (DIM == 3 || D -> nField() == 4)
      Field::traction (&_work[0], &_work[_nline], &_work[2*_nline], _nwall,
		       _npad, D->u[3],D->u[0],D->u[1],D->u[2]);
    else
      Field::traction (&_work[0], &_work[_nline], &_work[2*_nline], _nwall,
		       _npad, D->u[2],D->u[0],D->u[1]);

    // -- Inverse Fourier transform (like Field::bTransform).

    if (nZ > 1)
      if (nZ == 2)
	Veclib::copy (_npad, &_work[0], 1, &_work[_npad], 1);
      else
	Femlib::DFTr (&_work[0], nZ, _npad, INVERSE);

    // -- Write to file.
    
    // -- Header: this will be a lot like a standard header.
    //    Output normal and tangential tractions, 'n', 't', 's'.

    const char *Hdr_Fmt[] = { 
      "%-25s "                "Session\n",
      "%-25s "                "Created\n",
      "%-5d1    %-5d %-10d"   "Nr, Ns, Nz, Elements\n",
      "%-25d "                "Step\n",
      "%-25.6g "              "Time\n",
      "%-25.6g "              "Time step\n",
      "%-25.6g "              "Kinvis\n",
      "%-25.6g "              "Beta\n",
      "%-25s "                "Fields written\n",
      "%-25s "                "Format\n" };
      
    char   s1[StrMax], s2[StrMax];
    time_t tp (time (0));

    strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
	
    sprintf (s1, Hdr_Fmt[0], D->name);                  cout << s1;
    sprintf (s1, Hdr_Fmt[1], s2);                       cout << s1;
    sprintf (s1, Hdr_Fmt[2], nP, nZ, _nwall);           cout << s1;
    sprintf (s1, Hdr_Fmt[3], D->step);                  cout << s1;
    sprintf (s1, Hdr_Fmt[4], D->time);                  cout << s1;
    sprintf (s1, Hdr_Fmt[5], Femlib::value ("D_T"));    cout << s1;
    sprintf (s1, Hdr_Fmt[6], Femlib::value ("KINVIS")); cout << s1;
    sprintf (s1, Hdr_Fmt[7], Femlib::value ("BETA"));   cout << s1;
    sprintf (s1, Hdr_Fmt[8], "nts");                    cout << s1;
    sprintf (s2, "binary "); Veclib::describeFormat  (s2 + strlen (s2));
    sprintf (s1, Hdr_Fmt[9], s2);                       cout << s1;

    // -- Data.
  
    for (j = 0; j < 3; j++)
      for (i = 0; i < nZ; i++) {
	plane = &_work[i*_npad + j*_nline];
	cout.write (reinterpret_cast<char*>(plane),
		    static_cast<int_t>(_nline * sizeof (real_t))); 
      }
  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int        argc   ,
		     char**     argv   ,
		     char*&     session,
		     istream*&  file   )
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session [file]"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -v[vv..] ... increase verbosity level\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE", Femlib::ivalue ("VERBOSE") + 1);
      while (*++argv[0] == 'v');
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  switch (argc) {
  case 1:
    session = argv[0];
    file = &cin;
    break;
  case 2:
    session = argv[0];
    file = new ifstream (argv[1]);
    if (file -> bad()) {
      cerr << usage;
      sprintf (buf, "unable to open field file: %s", argv[1]);
      message (prog, buf, ERROR);
    }
    break;
  default:
    cerr << usage;
    message (prog, "session file not supplied", ERROR);
    break;
  }  
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain )
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
}


static bool getDump (Domain*    D   ,
		     istream&   dump)
// ---------------------------------------------------------------------------
// Read next set of field dumps from file.
// ---------------------------------------------------------------------------
{
  dump >> *D;
  return dump.good ();
}
