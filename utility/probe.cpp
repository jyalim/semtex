///////////////////////////////////////////////////////////////////////////////
// probe.C: extract results from a field file at a set of 3D points.
//
// Copyright (c) 1997 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
// Synopsis
// --------
// Probe has three different user interfaces --- the internal
// mechanics of data extraction are the same in all cases, but the way
// the points are defined and data are output differ for each interface.
//
// Common information
// ------------------
// If the x--y values of a point cannot be located in the (2D) mesh,
// it is ignored, but if the z values falls outside [0, TWOPI/BETA],
// Fourier interpolation will be applied on assumption of periodicity.
//
// The field file must be in binary format, and the value of BETA in the
// session file will override the value given in the field file header.
//
// Interface 1: Extract data at set of points
// -----------
//
// Usage: probe [-h] [-p file] -s session dump [file]
//
// This is the most general form.  Extract data at a specified set
// of 3D points, given either on standard input or in a named file.
// Points must have three coordinates, and any number of them can be given,
// e.g.
//         1.32461       0.514135      1.00
//         1.31102       0.509459     -0.25
//            ..             ..         ..
//
// Points can either be supplied on standard input or in a named file.
//
// Output is always ASCII format.  Each line of output contains the values
// for the fields in the file in columns, in the order they were written
// to the field file.
//
// Interface 2: Extract data along a straight line
// -----------
//
// Usage: probeline [-h] -p "[n:]x0,y0,z0,dx,dy,dz" -s session dump
//
// Extract data along the straight line defined by the parameters.
// The number of points extracted is specified as <num>.  Output
// format is the same as for probe.
//
// Interface 3: Extract data at 2D array of points on x-y, x-z or y-z planes
// -----------
//
// Usage: probeplane [-h] [options] -s session dump
// options:
// -xy "x0,y0,dx,dy" ... xy-cutting plane
// -xz "x0,z0,dx,dz" ... xz-cutting plane
// -yz "y0,z0,dy,dz" ... yz-cutting plane
// -orig #           ... origin of the cutting plane along orthogonal axis
// -nx #             ... resolution along the x-axis
// -ny #             ... resolution along the y-axis
// -swap             ... swap x <--> y in output file (rotate)
// -tec              ... write TECPLOT-formatted ASCII output
// -sm               ... write SM-formatted binary output
// -0		     ... output zero if point is outside mesh (no warning)
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

static char RCS[] = "$Id: probe.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>

static const int_t NPTS = 64;	// -- Default number of points for line/plane.

typedef enum {			// -- Flags for coordinate axes.
  None = 0,
  X    = 'x',
  Y    = 'y',
  Z    = 'z'
} AXIS;

static char  *prog;
static void  getargs     (int, char**, char*&, char*&, int_t&,
			  char*&, char*&, char*&);
static int_t loadPoints  (istream&, vector<Point*>&);
static int_t linePoints  (vector<Point*>&);
static int_t planePoints (vector<Point*>&, Mesh*);
static void  findPoints  (vector<Point*>&, vector<Element*>&,
			  vector<Element*>&, vector<real_t>&, vector<real_t>&);
static int_t getDump     (ifstream&, vector<AuxField*>&, vector<Element*>&,
			  const int_t, const int_t, const int_t);
static void  putData     (const char*, const char*, const char*, int_t,
			  int_t, vector<AuxField*>&, vector<Element*>&,
			  vector<Point*>&, vector<vector<real_t> >&);
static void  Finterp     (vector<AuxField*>&, const Point*, const Element*,
			  const real_t, const real_t, const int_t,
			  real_t*, int_t*, real_t*);
static bool  doSwap      (const char*);
static char* root        (char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char                    *session, *dump, *format;
  char                    *interface = 0, *points = 0;
  int_t                   NP, NZ,  NEL;
  int_t                   i, j, k, nf, ntot = 0, rotswap = 0;
  ifstream                fldfile;
  istream*                pntfile;
  FEML*                   F;
  Mesh*                   M;
  const real_t*           knot;
  vector<int_t>           iwork(15);
  vector<real_t>          rwork, r, s, datum;
  vector<Point*>          point;
  vector<Element*>        elmt;
  vector<Element*>        Esys;
  vector<AuxField*>       u;
  vector<vector<real_t> > data;

  // -- Initialize.

  prog = *argv;
  Femlib::initialize (&argc, &argv);

  // -- Set defaults for probeplane interface.

  Femlib::ivalue ("SIZED"  , None);
  Femlib::ivalue ("NX"     , NPTS);
  Femlib::ivalue ("NY"     , NPTS);
  Femlib::ivalue ("ORTHO"  , Z   );
  Femlib::ivalue ("PRINT_OUTSIDE", 0 );
  Femlib:: value ("OFFSET" , 0.0 );
  Femlib:: value ("X_MIN"  , 0.0 );
  Femlib:: value ("Y_MIN"  , 0.0 );
  Femlib:: value ("Z_MIN"  , 0.0 );
  Femlib:: value ("X_DELTA", 0.0 );
  Femlib:: value ("Y_DELTA", 0.0 );
  Femlib:: value ("Z_DELTA", 0.0 );

  // -- Parse command line.

  getargs (argc, argv, interface, format, rotswap, session, dump, points);

  // -- Check presence of field file before proceeding.

  fldfile.open (dump, ios::in);
  if (!fldfile) message (prog, "no field file", ERROR);
  
  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);

  NEL = M -> nEl();  
  NP  = Femlib::ivalue ("N_P");
  NZ  = Femlib::ivalue ("N_Z");
  
  Geometry::set (NP, NZ, NEL, Geometry::Cartesian);
  Esys.resize   (NEL);

  for (k = 0; k < NEL; k++) Esys[k] = new Element (k, NP, M);
  
  // -- Set up FFT work areas.

  rwork.resize  (3*NZ);
  Femlib::rffti (NZ, &rwork[0], &iwork[0]);

  // -- Construct list of points.

  if (strcmp (interface, "probe") == 0) {
    if (points) {
      pntfile = new ifstream (points);
      if (pntfile -> bad()) message (prog, "unable to open point file", ERROR);
    } else pntfile = &cin;

    ntot = loadPoints (*pntfile, point);
  } 
  else if (strcmp (interface, "probeline")  == 0) ntot = linePoints  (point);
  else if (strcmp (interface, "probeplane") == 0) ntot = planePoints (point,M);

  // -- Locate points in the Mesh.

  findPoints (point, Esys, elmt, r, s);

  // -- Load field file.

  if (!(getDump (fldfile, u, Esys, NP, NZ, NEL)))
    message (prog, "no data extracted", ERROR);

  datum.resize (nf = u.size());
  data.resize  (ntot);
  for (i = 0; i < ntot; i++) data[i].resize (nf);

  // -- Interpolate within it.

  for (i = 0; i < ntot; i++)
    if (elmt[i]) {
      Finterp (u,point[i],elmt[i],r[i],s[i],NZ,&rwork[0],&iwork[0],&datum[0]); 
      for (j = 0; j < nf; j++) data[i][j] = datum[j];
    } else
      for (j = 0; j < nf; j++) data[i][j] = 0.0;

  // -- Output collected data.

  putData (dump, interface, format, ntot, rotswap, u, elmt, point, data);
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc     ,
		     char** argv     ,
		     char*& interface,
		     char*& format   ,
		     int_t& swap     ,
		     char*& session  ,
		     char*& dump     ,
		     char*& points   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  format = new char [16];
  strcpy (format, "free");	// -- Default output format.
  interface = *argv;
  char err[StrMax];

  if (strcmp (interface, "probe") == 0) {

    char usage[] =
      "Usage: probe [options] -s session dump\n"
      "  options:\n"
      "  -h      ... print this message\n"
      "  -p file ... name file of point data [Default: stdin]\n";

    while (--argc  && **++argv == '-')
      switch (*++argv[0]) {
      case 'h':
	cout << usage;
	exit (EXIT_SUCCESS);
	break;
      case 'p':
	if (*++argv[0])
	  points = *argv;
	else {
	  --argc;
	  points = *++argv;
	}
	break;
      case 's':
	if (*++argv[0])
	  session = *argv;
	else {
	  --argc;
	  session = *++argv;
	}
	break;
      default:
	cerr << usage;
	exit (EXIT_FAILURE);
	break;
      }

  } else if (strcmp (interface, "probeline") == 0) {

    char usage[] =
      "Usage: probeline [-h] -p \"[n:]x0,y0,z0,dx,dy,dz\" -s session dump\n";
    char *tok, *pspec;
    int set = 0;

    while (--argc && **++argv == '-')
      switch (*++argv[0]) {
      case 'h':
	cout << usage;
	exit (EXIT_SUCCESS);
	break;
      case 'p':
	if (*++argv[0])
	  pspec = *argv;
	else {
	  --argc;
	  pspec = *++argv;
	}
	if (strchr (pspec, ':')) {
	  if   (tok = strtok (pspec, ":")) Femlib::value ("NPTS", atof (tok));
	  else                             Femlib::value ("NPTS", NPTS);
	  if (tok = strtok (0, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    set = 1;
	  } else {
	    message (prog, "couldn't x0 parse number from string", ERROR);
	  }
	} else {
	  Femlib::value ("NPTS", NPTS);
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    set = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	}
	while (tok = strtok (0, ","))
	  switch (++set) {
	  case 2:
	    Femlib::value ("Y_MIN",   atof (tok));
	    break;
	  case 3:
	    Femlib::value ("Z_MIN",   atof (tok));
	    break;
	  case 4:
	    Femlib::value ("X_DELTA", atof (tok));
	    break;
	  case 5:
	    Femlib::value ("Y_DELTA", atof (tok));
	    break;
	  case 6:
	    Femlib::value ("Z_DELTA", atof (tok));
	    break;
	  default:
	    message (prog, "too many numbers in point string", ERROR);
	    break;
	  }
	if (set != 6) {
	  sprintf (err, "wrong number of parameters to line string (%1d)",set);
	  message (prog, err, ERROR);
	}
	break;
      case 's':
	if (*++argv[0])
	  session = *argv;
	else {
	  --argc;
	  session = *++argv;
	}
	break;
      default:
	cerr << usage;
	exit (EXIT_FAILURE);
	break;
      }
    if (!set)
      message (prog, "no points set", ERROR);

  } else if (strcmp (interface, "probeplane") == 0) {

    char *tok, *pspec, *usage =
      "Usage: probeplane [options] -s session dump\n"
      "options:\n"
      "-h                ... print this message\n"
      "-xy \"x0,y0,dx,dy\" ... xy-cutting plane\n"
      "-xz \"x0,z0,dx,dz\" ... xz-cutting plane\n"
      "-yz \"y0,z0,dy,dz\" ... yz-cutting plane\n"
      "-orig #           ... origin of the cutting plane along ortho axis\n"
      "-nx #             ... resolution along the x-axis\n"
      "-ny #             ... resolution along the y-axis\n"
      "-swap             ... swap output x <--> y (rotate)\n"
      "-tec              ... write TECPLOT-formatted ASCII output\n"
      "-sm               ... write SM-formatted binary output\n"
      "-0	         ... output zero if point is outside mesh (instead of warning)\n";

    int_t nset = 0;

    while (--argc && **++argv == '-')
      switch (*++argv[0]) {
      case 'h':
	cout << usage;
	exit (EXIT_SUCCESS);
	break;
      case 'n':
	switch (argv[0][1]) {
	case 'x':
	  --argc; ++argv;
	  Femlib::value ("NX", atof (*argv));
	  break;
	case 'y':
	  --argc; ++argv;
	  Femlib::value ("NY", atof (*argv));
	  break;
	default:
	  message (prog, "can only specify nx or ny", ERROR);
	  break;
	}
	break;
      case 'o':
	--argc;
	Femlib::value ("OFFSET", atof (*++argv));
	break;
      case 's':
	if      (argv[0][1] == 'm') strcpy (format, "sm");
	else if (argv[0][1] == 'w') swap = 1;
	else { --argc; session = *++argv; }
	break;
      case 't':
	strcpy (format, "tecplot");
	break;
      case '0':
	Femlib::ivalue ("PRINT_OUTSIDE", 1);
	break;
      case 'x':
	switch (argv[0][1]) {
	case 'y':
	  Femlib::value ("ORTHO", Z);
	  --argc; pspec = *++argv;
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    nset = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	  while (tok = strtok (0, ","))
	    switch (++nset) {
	    case 2:
	      Femlib::value ("Y_MIN",   atof (tok));
	      break;
	    case 3:
	      Femlib::value ("X_DELTA", atof (tok));
	      break;
	    case 4:
	      Femlib::value ("Y_DELTA", atof (tok));
	      break;
	    default:
	      message (prog, "too many numbers in point string", ERROR);
	      break;
	    }
	  if (nset != 4)
	    message (prog, "need 4 parameters for cutting plane", ERROR);
	  Femlib::value ("SIZED", 1.0);
	  break;
	case 'z':
	  Femlib::value ("ORTHO", Y);
	  --argc; pspec = *++argv;
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    nset = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	  while (tok = strtok (0, ","))
	    switch (++nset) {
	    case 2:
	      Femlib::value ("Z_MIN",   atof (tok));
	      break;
	    case 3:
	      Femlib::value ("X_DELTA", atof (tok));
	      break;
	    case 4:
	      Femlib::value ("Z_DELTA", atof (tok));
	      break;
	    default:
	      message (prog, "too many numbers in point string", ERROR);
	      break;
	    }
	  if (nset != 4)
	    message (prog, "need 4 parameters for cutting plane", ERROR);
	  Femlib::value ("SIZED", 1.0);
	  break;
	default:
	  message (prog, "extents can be xy, xz, or yz", ERROR);
	  break;
	}
	break;
      case 'y':
	switch (argv[0][1]) {
	case 'z':
	  Femlib::value ("ORTHO", X);
	  --argc; pspec = *++argv;
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("Y_MIN", atof (tok));
	    nset = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	  while (tok = strtok (0, ","))
	    switch (++nset) {
	    case 2:
	      Femlib::value ("Z_MIN",   atof (tok));
	      break;
	    case 3:
	      Femlib::value ("Y_DELTA", atof (tok));
	      break;
	    case 4:
	      Femlib::value ("Z_DELTA", atof (tok));
	      break;
	    default:
	      message (prog, "too many numbers in point string", ERROR);
	      break;
	    }
	  if (nset != 4)
	    message (prog, "need 4 parameters for cutting plane", ERROR);
	  Femlib::value ("SIZED", 1.0);
	  break;
	default:
	  message (prog, "extents can be xy, xz, or yz", ERROR);
	  break;
	}
	break;
      default:
	cerr << usage;
	exit (EXIT_FAILURE);
	break;
      }

  }

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static int_t loadPoints (istream&        pfile,
			 vector<Point*>& point)
// ---------------------------------------------------------------------------
// Probe point input for the "probe" interface, from file.
// ---------------------------------------------------------------------------
{
  int_t         ntot, num = 0;
  real_t        x, y, z;
  Point*        datum;
  stack<Point*> data;

  while (pfile >> x >> y >> z) {
    datum = new Point;
    datum -> x = x;
    datum -> y = y;
    datum -> z = z;
    data.push (datum);
    num++;
  }

  ntot = num;
  point.resize (ntot);

  while (num--) { point[num] = data.top(); data.pop(); }

  return ntot;
}


static int_t linePoints (vector<Point*>& point)
// ---------------------------------------------------------------------------
// Probe point generation for the "probeline" interface.
// ---------------------------------------------------------------------------
{
  int_t  i, ntot = Femlib::ivalue ("NPTS");
  real_t xmin = Femlib::value ("X_MIN");
  real_t ymin = Femlib::value ("Y_MIN");
  real_t zmin = Femlib::value ("Z_MIN");
  real_t dx   = (ntot == 1) ? 0.0 : Femlib::value ("X_DELTA") / (ntot - 1.0);
  real_t dy   = (ntot == 1) ? 0.0 : Femlib::value ("Y_DELTA") / (ntot - 1.0);
  real_t dz   = (ntot == 1) ? 0.0 : Femlib::value ("Z_DELTA") / (ntot - 1.0);

  point.resize (ntot);

  for (i = 0; i < ntot; i++) {
    point[i] = new Point;
    point[i] -> x = xmin + i * dx;
    point[i] -> y = ymin + i * dy;
    point[i] -> z = zmin + i * dz;
  }

  return ntot;
}


static int_t planePoints (vector<Point*>& point,
			  Mesh*           mesh )
// ---------------------------------------------------------------------------
// Probe point generation for the "probeplane" interface.
// ---------------------------------------------------------------------------
{
  int_t        i, j, k;
  const int_t  nx     = Femlib::ivalue ("NX");
  const int_t  ny     = Femlib::ivalue ("NY");
  const int_t  ortho  = Femlib::ivalue ("ORTHO");
  const real_t offset = Femlib:: value ("OFFSET");
  const int_t  ntot   = nx * ny;
  real_t       x0, y0, z0, dx, dy, dz;
  Point*       p;

  point.resize (ntot);

  switch (ortho) {
  case X:
    x0 = offset;
    y0 = Femlib::value ("Y_MIN");
    z0 = Femlib::value ("Z_MIN");
    dy = Femlib::value ("Y_DELTA") / (nx - 1.0);
    dz = Femlib::value ("Z_DELTA") / (ny - 1.0);
    for (k = 0, j = 0; j < ny; j++)
      for (i = 0; i < nx; i++, k++) {
	point[k] = p = new Point;
	p -> x = x0;
	p -> y = y0 + i * dy;
	p -> z = z0 + j * dz;
      }
    break;
  case Y:
    x0 = Femlib::value ("X_MIN");
    y0 = offset;
    z0 = Femlib::value ("Z_MIN");
    dx = Femlib::value ("X_DELTA") / (nx - 1.0);
    dz = Femlib::value ("Z_DELTA") / (ny - 1.0);
    for (k = 0, j = 0; j < ny; j++)
      for (i = 0; i < nx; i++, k++) {
	point[k] = p = new Point;
	p -> x = x0 + i * dx;
	p -> y = y0;
	p -> z = z0 + j * dz;
      }
    break;
  case Z:
    if (!(Femlib::ivalue ("SIZED"))) {
      Point lo, hi;
      mesh -> extent (lo, hi);
      Femlib::value ("X_MIN", lo.x);
      Femlib::value ("Y_MIN", lo.y);
      Femlib::value ("X_DELTA", hi.x - lo.x);
      Femlib::value ("Y_DELTA", hi.y - lo.y);
    }
    x0 = Femlib::value ("X_MIN");
    y0 = Femlib::value ("Y_MIN");
    dx = Femlib::value ("X_DELTA") / (nx - 1.0);
    dy = Femlib::value ("Y_DELTA") / (ny - 1.0);
    z0 = offset;
    for (k = 0, j = 0; j < ny; j++)
      for (i = 0; i < nx; i++, k++) {
	point[k] = p = new Point;
	p -> x = x0 + i * dx;
	p -> y = y0 + j * dy;
	p -> z = z0;
      }
    break;
  default:
    break;
  }
  
  return ntot;
}


static void findPoints (vector<Point*>&   point,
			vector<Element*>& Esys ,
			vector<Element*>& elmt ,
			vector<real_t>&   rloc ,
			vector<real_t>&   sloc )
// ---------------------------------------------------------------------------
// Locate points within elements, set Element pointer & r--s locations.
// ---------------------------------------------------------------------------
{
  int_t          i, k;
  real_t         x, y, z, r, s;
  const int_t    NEL   = Esys .size();
  const int_t    NPT   = point.size();
  const bool     guess = true;
  vector<real_t> work(static_cast<size_t>
		      (max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));

  elmt.resize (NPT);
  rloc.resize (NPT);
  sloc.resize (NPT);

  for (i = 0; i < elmt.size(); i++) elmt[i] = 0;

  cerr.precision (8);

  for (i = 0; i < NPT; i++) {
    x = point[i] -> x;
    y = point[i] -> y;
    z = point[i] -> z;
    for (k = 0; k < NEL; k++) {
      r = s = 0.0;
      if (Esys[k] -> locate (x, y, r, s, &work[0], guess)) {
	elmt[i] = Esys[k];
	rloc[i] = r;
	sloc[i] = s;
	break;
      }
    }

    if (!elmt[i] && Femlib::ivalue ("PRINT_OUTSIDE") == 0){
      cerr << "point ("
	   << setw(15) << x << ","
	   << setw(15) << y << ","
	   << setw(15) << z << ") is not in the mesh" << endl;
    }
  }
}


static int_t getDump (ifstream&          file,
		      vector<AuxField*>& u   ,
		      vector<Element*>&  Esys,
		      const int_t        np  ,
		      const int_t        nz  ,
		      const int_t        nel )
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  char  buf[StrMax], fields[StrMax];
  int_t i, nf, npnew, nznew, nelnew;
  bool  swab;

  if (file.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> npnew >> nznew >> nznew >> nelnew;
  file.getline (buf, StrMax);
  
  if (np != npnew || nz != nznew || nel != nelnew)
    message (prog, "size of dump mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields on first pass.

  if (u.size() == 0) {
    u.resize (nf);
    for (i = 0; i < nf; i++)
      u[i] = new AuxField (new real_t[Geometry::nTotal()], nz, Esys,fields[i]);
  } else if (u.size() != nf) 
    message (prog, "number of fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return file.good();
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping is required to make sense of binary input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}


static void Finterp (vector<AuxField*>& u    ,
		     const Point*       P    ,
		     const Element*     E    ,
		     const real_t       r    , 
		     const real_t       s    ,
		     const int_t        NZ   ,
		     real_t*            rwork,
		     int_t*             iwork, 
		     real_t*            data )
// ---------------------------------------------------------------------------
// Carry out 2DxFourier interpolation.
// ---------------------------------------------------------------------------
{
  register int_t  i, k, Re, Im;
  register real_t phase;
  const int_t     NF    = u.size();
  const int_t     NZH   = NZ >> 1;
  const int_t     NHM   = NZH - 1;
  const real_t    betaZ = P -> z * Femlib::value("BETA");
  const real_t*   Wtab  = rwork;
  real_t*         work  = rwork + NZ;
  real_t*         temp  = rwork + NZ + NZ;

  for (i = 0; i < NF; i++)	// -- For each field.

    if (NZ == 1)

      // -- Just 2D.

      data[i] = u[i] -> probe (E, r, s, 0);

    else {

      // -- 2D interpolation.
      
      for (k = 0; k < NZ; k++) work[k] = u[i] -> probe (E, r, s, k);

      // -- Fourier interpolation.

      Femlib::rfftf (NZ, work, temp, Wtab, iwork);
      Blas::scal    (NZ, 2.0/NZ, work, 1);

      if (NZ & 1) {

	data[i] = 0.5 * work[0];
	for (k = 1; k <= NZH; k++) {
	  Im       = k  + k;
	  Re       = Im - 1;
	  phase    = k * betaZ;
	  data[i] += work[Re] * cos (phase) - work[Im] * sin (phase);
	}

      } else {

	data[i] = 0.5 * work[0];
	for (k = 1; k <= NHM; k++) {
	  Im       = k  + k;
	  Re       = Im - 1;
	  phase    = k * betaZ;
	  data[i] += work[Re] * cos (phase) - work[Im] * sin (phase);
	}
	data[i] += 0.5 * work[NZ - 1] * cos (NZH * betaZ);

      }
    }
}


static char *root (char *s) {
  char *p = s + strlen(s)-1;
  while (p > s && *p != '.') p--;
  if (p != s) *p = '\0';
  return s;
}


static void putData (const char*              dump     ,
		     const char*              interface,
		     const char*              format   ,
		     int_t                    ntot     ,
		     int_t                    swap     ,
		     vector<AuxField*>&       u        ,
		     vector<Element*>&        elmt     ,
		     vector<Point*>&          point    ,
		     vector<vector<real_t> >& data     )
// ---------------------------------------------------------------------------
// Handle all the different output formats, according to the chosen interface.
// ---------------------------------------------------------------------------
{
  int_t i, j, k, n, nf = u.size();

  if (strstr (format, "free")) {
    cout.precision (6);
    for (i = 0; i < ntot; i++) {
      if (elmt[i] || Femlib::value ("PRINT_OUTSIDE") == 1) {
	cout << setw (5) << i + 1 << " " 
	     << setw(12) << point[i] -> x << " " 
	     << setw(12) << point[i] -> y << " " 
	     << setw(12) << point[i] -> z;
	for (j = 0; j < nf; j++)
	  cout << setw(15) << data[i][j];
	cout << endl;
      }
    }
    return;

  } else if (strstr (format, "sm") && strstr (interface, "probeplane")) {

    char      fname[StrMax], *base = root (strdup (dump));
    const int_t nx = Femlib::ivalue ("NX");
    const int_t ny = Femlib::ivalue ("NY");
    ofstream  out;

    for (n = 0; n < nf; n++) {
      sprintf (fname, "%s.%c", base, u[n] -> name());
      out.open  (fname);

      if (swap) {
	out.write (reinterpret_cast<const char*>(&ny), sizeof (int_t));
	out.write (reinterpret_cast<const char*>(&nx), sizeof (int_t));
	
	for (i = 0; i < nx; i++) {
	  for (j = 0; j < ny; j++) {
	    float tmp = static_cast<float>(data[i + j*nx][n]);
	    out.write (reinterpret_cast<char*>(&tmp), sizeof (float));
	  }
	}
      } else {
	out.write (reinterpret_cast<const char*>(&nx), sizeof (int_t));
	out.write (reinterpret_cast<const char*>(&ny), sizeof (int_t));

	for (k = 0, j = 0; j < ny; j++) {
	  for (i = 0; i < nx; i++, k++) {
	    float tmp = static_cast<float>(data[k][n]);
	    out.write (reinterpret_cast<char*>(&tmp), sizeof (float));
	  }
	}
      }
      out.close();
    }
  } else if (strstr (format, "tecplot") && strstr (interface, "probeplane")) {

    char      fname[StrMax], *base = root (strdup (dump));
    const int_t nx = Femlib::ivalue ("NX");
    const int_t ny = Femlib::ivalue ("NY");

    if (swap) {
        cerr << "Option -swap not yet implemented for tecplot output." << endl;
	return;
    }

    cout << "TITLE = \"" << base << "\"" << endl;
    cout << "VARIABLES = \"x\", \"y\", \"z\"";
    for (n = 0; n < nf; n++) cout << ", \"" << u[n] -> name() << "\"";
    cout << endl;
    cout << "ZONE I=" << nx << ", J=" << ny << ", F=Point, T=\"" << base << "\""<< endl;

    for (k = 0, j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++, k++) {
	cout << setw(12) << point[k] -> x << " "
	     << setw(12) << point[k] -> y << " "
	     << setw(12) << point[k] -> z;
        for (n = 0; n < nf; n++) {
	    cout << " " << setw(12) << data[k][n];
	}
	cout << endl;
      }
    }
    return;

  } else if (strstr (format, "tecplot")) {
    cerr << "tecplot output only implemented for probeplane interface." << endl; 
  }
}
