///////////////////////////////////////////////////////////////////////////////
// project.C:  Project solution files to different interpolation orders.
//
// Copyright (c) 1996 <--> $Date: 2015/05/28 07:49:14 $, Hugh Blackburn
//
// SYNOPSIS
// --------
// Process sem field file, project to new interpolation order on the
// same mesh.  Each dump in file is expected to be the same size.
// Also it is assumed that the field file represents a vector field
// dump, so that if the input file has N space dimensions, the first N
// fields represent vector components.  Input file must be binary
// format.
//
// USAGE
// -----
// project [options] [file]
// options:
// -h       ... print this message.
// -n <num> ... project elements to num x num.
// -z <num> ... project to <num> planes in the homogeneous direction.
// -w       ... Retain w components in 3D-->2D proj'n [Default: delete]
// -u       ... project elements to uniform internal grid [Default: GLL].
// -U       ... project from uniform grid to GLL.
// 
// If file is not present, read from standard input.  Write to
// standard output.
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

static char RCS[] = "$Id: project.cpp,v 8.2 2015/05/28 07:49:14 hmb Exp $";

#include <sem.h>

static int_t uniform = 0;


static int_t _index (const char* s, char c)
/* ------------------------------------------------------------------------- *
 * Return index of c in s, -1 if not found.
 * ------------------------------------------------------------------------- */
{
  int_t       i;
  const int_t len = strlen (s);

  for (i = 0; i < len; i++) if (s[i] == c) return i;

  return -1;
}


class Field2DF
// ============================================================================
// Canonical field class, each nr * ns element is defined on [-1,1] X [-1, 1].
// Data are arranged element-ordered in 2D planes to create a 3D scalar field.
// ============================================================================
{
friend istream& operator >> (istream&, Field2DF&);
friend ostream& operator << (ostream&, Field2DF&);

public:
  Field2DF  (const int_t nr, const int_t ns, const int_t nZ, const int_t nEl,
	     const char Name='\0');
  ~Field2DF () { delete data; delete plane; }

  char getName () { return name; }

  Field2DF& operator = (const Field2DF&);
  Field2DF& operator = (const real_t);

  Field2DF& transform (const int_t);
  Field2DF& reverse   ();
  
private:
  const char  name;
  const int_t nr, ns, nz, nel, nrns;
  int_t       nplane, ntot;
  real_t*     data;
  real_t**    plane;
};


Field2DF::Field2DF (const int_t nR  ,
		    const int_t nS  ,
		    const int_t nZ  ,
		    const int_t nEl ,
		    const char  Name) :

		    name       (Name),
                    nr         (nR  ),
                    ns         (nS  ),
		    nz         (nZ  ),
		    nel        (nEl ),
		    nrns       (nr * ns)
// ---------------------------------------------------------------------------
// Field2DF constructor. 
// ---------------------------------------------------------------------------
{
  register int_t i;
  
  nplane = nrns * nel;
  if (nplane > 1 && nplane & 1) nplane++;
  ntot   = nplane * nz;

  data  = new real_t  [ntot];
  plane = new real_t* [nz];

  for (i = 0; i < nz; i++) plane[i] = data + i * nplane;
  Veclib::zero (ntot, data, 1);
}


Field2DF& Field2DF::transform (const int_t sign)
// ---------------------------------------------------------------------------
// Carry out Fourier transformation in z direction.
// ---------------------------------------------------------------------------
{
  if (nz > 2) Femlib::DFTr (data, nz, nplane, sign);

  return *this;
}


Field2DF& Field2DF::operator = (const Field2DF& rhs)
// ---------------------------------------------------------------------------
// If the two fields conform, copy rhs's data storage to lhs.
//
// Otherwise perform projection/interpolation of rhs's data area to lhs.
// Interpolation ASSUMES THAT FOURIER TRANSFORMATION HAS ALREADY OCCURRED
// in z direction if rhs is 3D.  Truncation of Fourier modes occurs if this
// Field2DF has less modes than rhs (to avoid aliasing).
// ---------------------------------------------------------------------------
{
  if (rhs.nel != nel)
    message ("Field2DF::operator =", "fields can't conform", ERROR);

  if (rhs.nr == nr && rhs.ns == ns && rhs.nz == nz) // -- No project, just copy.
    Veclib::copy (ntot, rhs.data, 1, data, 1);

  else {			// -- Perform projection.

    register int_t  i, k;
    register real_t *LHS, *RHS;
    const real_t    *IN,  *IT;
    const int_t     nzm = min (rhs.nz, nz);
    vector<real_t>  work (rhs.nr * nr);
    real_t*         tmp = &work[0];

    if      (uniform == +1)
      Femlib::projection (&IN, &IT, rhs.nr, GLJ, 0.0, 0.0, nr, TRZ, 0.0, 0.0);
    else if (uniform == -1) 
      Femlib::projection (&IN, &IT, rhs.nr, TRZ, 0.0, 0.0, nr, GLJ, 0.0, 0.0);
    else
      Femlib::projection (&IN, &IT, rhs.nr, GLJ, 0.0, 0.0, nr, GLJ, 0.0, 0.0);

    for (k = 0; k < nzm; k++) {	// -- 2D planar projections.
      LHS = plane[k];
      RHS = rhs.plane[k];

      if (rhs.nr == nr && rhs.ns == ns)
	Veclib::copy (nplane, RHS, 1, LHS, 1);
      else
	for (i = 0; i < nel; i++, LHS += nrns, RHS += rhs.nrns) {
	  Blas::mxm (IN, nr, RHS, rhs.nr, tmp, rhs.nr);
	  Blas::mxm (tmp, nr, IT, rhs.nr, LHS,     nr);
	}
    }

    if ((i = nz - rhs.nz) > 0) // -- Zero pad for Fourier projections.
      Veclib::zero (i * nplane, data + rhs.ntot, 1);
  }

  return *this;
}


Field2DF& Field2DF::operator = (const real_t val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (ntot,      data, 1);
  else              Veclib::fill (ntot, val, data, 1);

  return *this;
}


Field2DF& Field2DF::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (ntot, data, 1, data, 1);

  return *this;
}


ostream& operator << (ostream&  strm,
		      Field2DF& F   )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  for (i = 0; i < F.nz; i++)
    strm.write ((char*) F.plane[i], F.nrns * F.nel * sizeof (real_t));

  return strm;
}


istream& operator >> (istream&  strm,
		      Field2DF& F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  for (i = 0; i < F.nz; i++)
    strm.read ((char*) F.plane[i], F.nrns * F.nel * sizeof (real_t));

  return strm;
}


static char prog[] = "project";
static void getargs  (int, char**, int_t&, int_t&, int_t&, istream*&);
static bool getDump  (istream&, ostream&, vector<Field2DF*>&,
		      int_t&, int_t&, int_t&, int_t&, int_t&, int_t&, char*);
static bool doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char              fields[StrMax];
  int_t             i, j, nEl, fInc;
  int_t             nRnew = 0, nSnew = 0, nZnew = 0, keepW = 0;
  istream*          input;
  vector<Field2DF*> Uold, Unew;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, nRnew, nZnew, keepW, input);

  // -- If a 2D projection is requested, enforce equal order in each direction.
  //    Otherwise, elemental structure (perhaps non-square) is left unchanged.
  //    This allows z-projection of non-square elements.

  nSnew = nRnew;
  
  while (getDump (*input,cout,Uold,nRnew,nSnew,nZnew,nEl,keepW,fInc,fields)) {

    Unew.resize (Uold.size() + fInc);

    switch (fInc) {

    case 0:			// -- Same number of fields out as in.
      for (i = 0; i < Uold.size(); i++) {
	Unew[i] = new Field2DF (nRnew, nSnew, nZnew, nEl, Uold[i] -> getName());
       *Unew[i] = *Uold[i];
      }
      break;

    case 1:			// -- Add a new blank field called 'w'.
      for (i = 0; i < Uold.size(); i++) {
	j = _index (fields, Uold[i] -> getName());
	Unew[j] = new Field2DF (nRnew, nSnew, nZnew, nEl, Uold[i] -> getName());
	*Unew[j] = *Uold[i];
      }
      j = _index (fields, 'w');
      Unew[j] = new Field2DF (nRnew, nSnew, nZnew, nEl, 'w');
      *Unew[j] = 0.0;
      break;

    case -1:			// -- Delete field called 'w'.
      for (j = 0, i = 0; i < Uold.size(); i++) {
	if (Uold[i] -> getName() == 'w') continue;
	Unew[j] = new Field2DF (nRnew, nSnew, nZnew, nEl, Uold[i] -> getName());
	*Unew[j] = *Uold[i];
	j++;
      }
      break;

    default:			// -- Whoops.
      message (prog, "unrecognized conversion", ERROR);
      break;
    }

    for (i = 0; i < Unew.size(); i++) {
      Unew[i] -> transform (INVERSE);
      cout << *Unew[i];
    }
  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     int_t&    np   ,
		     int_t&    nz   ,
		     int_t&    keepW,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: project [options] [file]\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -n <num> ... 2D projection onto num X num\n"
    "  -z <num> ... 3D projection onto num planes\n"
    "  -w       ... Retain w components in 3D-->2D proj'n [Default: delete]\n"
    "  -u       ... project to uniform grid from GLL\n"
    "  -U       ... project from uniform grid to GLL\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0]) np = atoi (*argv);
      else { --argc;  np = atoi (*++argv); }
      break;
    case 'z':
      if (*++argv[0]) nz = atoi (*argv);
      else { --argc; nz = atoi (*++argv); }
      break;
    case 'u':
      uniform =  1;
      break;
    case 'U':
      uniform = -1;
      break;
    case 'w':
      keepW   = 1;
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping of input is required to make sense of input.
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


static bool getDump (istream&           ifile ,
		     ostream&           ofile ,
		     vector<Field2DF*>& u     ,
		     int_t&             nrnew ,
		     int_t&             nsnew ,
		     int_t&             nznew ,
		     int_t&             nel   ,
		     int_t&             keepW ,
		     int_t&             finc  ,
		     char*              fields)
// ---------------------------------------------------------------------------
// Read next set of field dumps from ifile, put headers on ofile.
// Note that the output header matches the requested values of N_P and
// N_Z, but the input is done according to the values found in the
// original field dump.
//
// Convert to Fourier space.
// ---------------------------------------------------------------------------
{
  static const char* hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields written\n",
    "%-25s "    "Format\n"
  };
  char  buf[StrMax], fmt[StrMax], oldfields[StrMax];
  int_t i, j, swab, nf, nr, ns, nz;

  if (ifile.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  ofile << buf << endl;
  ifile.getline (buf, StrMax);
  ofile << buf << endl;

  // -- I/O numerical description of field sizes.

  ifile >> nr >> ns >> nz >> nel;
  ifile.getline (buf, StrMax);
  
  if (!nrnew) { nrnew = nr; nsnew = ns; }
  if (!nznew) nznew = nz;

  sprintf (fmt, "%1d %1d %1d %1d", nrnew, nsnew, nznew, nel);
  sprintf (buf, hdr_fmt[2], fmt);
  cout << buf;

#if 0
  if      (nz >  1 && nznew == 1 && !keepW) finc = -1;
  else if (nz == 1 && nznew >  1)           finc = +1;
  else                                      finc =  0;
#endif

  for (i = 0; i < 5; i++) {
   ifile.getline (buf, StrMax);
   ofile << buf << endl;
  }

  // -- I/O field names.

  ifile >> oldfields;
  nf = strlen (oldfields);

#if 1
  finc = 0; 			  // -- Default: fields out = fields in.
  if (strstr (oldfields, "uv")) { // -- Velocity vectors present.
    // 2D2C --> 3D3C.
    if (nz == 1 && nznew > 1 && !strchr (oldfields, 'w')) finc =  1;
    // 3D3C --> 2D2C.
    else if (nz > 1 && nznew == 1 && !keepW)              finc = -1;
  }
#else
  if (finc == +1 &&  strchr (oldfields, 'w')) finc = 0;
  if (finc == -1 && !strchr (oldfields, 'w')) finc = 0;
#endif

  for (j = 0, i = 0; i < nf; i++) {
    if (finc == -1 && oldfields[i] == 'w') continue;
    if (finc == +1 && oldfields[i] == 'v') {
      fields[j++] = 'v';
      fields[j++] = 'w';
    } else
      fields[j++] = oldfields[i];
  }
  fields[j] = '\0';
  sprintf (buf, hdr_fmt[8], fields);
  cout << buf;
  ifile.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  ifile.getline (buf, StrMax);

  swab = doSwap (buf);

  sprintf (buf, "binary ");
  Veclib::describeFormat (buf + strlen (buf));
  sprintf (fmt, hdr_fmt[9], buf);
  cout << fmt;

  if (u.size() != nf) {
    u.resize (nf);
    for (i = 0; i < nf; i++) u[i] = new Field2DF (nr,ns,nz,nel,oldfields[i]);
  }

  for (i = 0; i < nf; i++) {
    ifile >> *u[i];
    if (swab) u[i] -> reverse();
    u[i] -> transform (FORWARD);
  }

  return ifile.good();
}


