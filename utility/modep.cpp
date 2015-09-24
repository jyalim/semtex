///////////////////////////////////////////////////////////////////////////////
// modep.C: project velocity data in a nominated Fourier mode onto a
// spatial mode shape supplied in a file, return scalar value of integral.
//
// Copyright (c) 1999 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
// USAGE
// -----
// modep [options] -m shapefile session [file]
// options:
// -h         ... print this message.
// -f         ... input file is already in Fourier-transformed state.
// -n <num>   ... Fourier mode number (starts at 0, default value 1).
// 
// If file is not present, read from standard input.  Write to
// standard output.  All IO assumes binary field files.
//
// If mode number = 0, then the mode shape must be real, planar, and
// have a number of components matching the number of velocity
// components in the input field file.  For mode numbers 1 or greater,
// the requirements are the same but there should be two data planes
// in the mode shape for each velocity component, i.e. the shape is of
// a 2D complex velocity mode.  The number of velocity components to
// be dealt with corresponds to the number supplied in the shapefile
// ("uv" or "uvw"), and it is assumed that the ordering of the first
// few data fields of the input file matches.
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

static char RCS[] = "$Id: modep.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>


static char prog[] = "modep";
static void  getargs  (int, char**, int_t&, bool&, char*&, char*&, char*&);
static bool  getDump  (istream&, vector<AuxField*>&, vector<Element*>&,
		       const int_t, const int_t, const int_t);
static bool  doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char               *session = 0, *shape= 0, *dump = 0;
  istream            *shapefile, *fldfile;
  int_t              NP, NZ,  NEL;
  int_t              np, nel, ntot, i, nzMode, modeNum = 1;
  real_t             integral;
  real_t*            alloc;
  FEML*              F;
  Mesh*              M;
  Geometry::CoordSys space;
  bool               transform = true;
  vector<Element*>   Esys;
  vector<AuxField*>  u;			 // -- Velocity field.
  vector<AuxField*>  modeShape, extract; // -- Extract is a modal subset of u.
  AuxField*          energy;		 // -- Workspace for projection.

  // -- Initialize.

  Femlib::initialize (&argc, &argv);
  getargs            (argc, argv, modeNum, transform, session, shape, dump);
  cout.precision     (8);
  nzMode = (modeNum > 0) ? 2 : 1;

  shapefile = new ifstream (shape);
  if (shapefile -> bad()) message (prog, "no mode shape file", ERROR);
 
  if (dump) {
    fldfile = new ifstream (dump);
    if (fldfile -> bad()) message (prog, "no field file", ERROR);
  } else fldfile = &cin;

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);

  NEL = M -> nEl();  
  NP  = Femlib::ivalue ("N_P");
  NZ  = Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (NP, NZ, NEL, space);
  Esys.resize   (NEL);

  for (i = 0; i < NEL; i++) Esys[i] = new Element (i, NP, M);

  // -- Load mode shape file.

  getDump (*shapefile, modeShape, Esys, NP, nzMode, NEL);

  // -- Allocate storage for extraction and projection.

  extract.resize (modeShape.size());
  for (i = 0; i < modeShape.size(); i++) {
    alloc = new real_t [Geometry::planeSize() * nzMode];
    extract[i]  = new AuxField (alloc, nzMode, Esys, 'U' + i);
  }
  
  alloc  = new real_t [Geometry::planeSize() * nzMode];
  energy = new AuxField (alloc, nzMode, Esys);
  
  // -- Load field file and project onto supplied mode shape.

  getDump (*fldfile, u, Esys, NP, NZ, NEL);

  if (u.size() != modeShape.size())
    message (prog, "number of velocity components mismatch", ERROR);

  for (i = 0; i < u.size(); i++) {
    if (transform) u[i] -> transform (FORWARD); // -- Go to Fourier space.
    extract[i] -> extractMode (*u[i], modeNum);
  }

  energy -> innerProductMode (modeShape, extract);

  // -- Print out the KE-norm projection of the data onto the supplied
  //    mode shape followed by its projection onto the 1/4-wavelength
  //    rotation of the mode shape (each real numbers).

  cout << energy -> integral(0) << "  " << energy -> integral(1) << endl;

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     int_t& mode   ,
		     bool&  xform  ,
		     char*& session,
		     char*& shape  ,
		     char*& dump   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "modep [options] -m shapefile session [file]\n"
    " options:\n"
    " -h         ... print this message.\n"
    " -f         ... input file is already in Fourier-transformed state.\n"
    " -n <num>   ... Fourier mode number (starts at 0, default value 1).\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'f':
      xform = false;
      break;
    case 'n':
      if   (*++argv[0]) mode = atoi (*argv);
      else { --argc;    mode = atoi (*++argv); }
      break;
    case 'm':
      if   (*++argv[0]) shape = *argv;
      else { --argc;    shape = *++argv; }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if      (argc == 1)   session = argv[0];
  else if (argc == 2) { session = argv[0]; dump = argv[1]; }
  else                  message (prog, usage, ERROR);
}


static bool getDump (istream&           file,
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
  char    buf[StrMax], fields[StrMax];
  int_t   i, swab, nf, npnew, nznew, nelnew;
  real_t* alloc;

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
  nf = (nf > 3) ? 3 : nf; 	// -- At most 3 velocity components dealt with.

  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields on first pass.

  if (u.size() == 0) {
    u.resize (nf);
    for (i = 0; i < nf; i++) {
      alloc = new real_t [Geometry::nTotProc()];
      u[i]  = new AuxField (alloc, nz, Esys, fields[i]);
    }
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
