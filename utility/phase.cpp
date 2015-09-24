///////////////////////////////////////////////////////////////////////////////
// phase.C: do operations on 3D data file in phase/Fourier space.
//
// Copyright (c) 2002 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
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
//
// USAGE
// -----
// phase [options] [file]
// options:
// -h       ... print this message.
// -c       ... perform complex conjugation.
// -f       ... data are already Fourier transformed (do not transform).
// -z       ... take mode zero as complex (e.g. it is an eigenmode).
// -r       ... enforce reflection symmetry of velocity & pressure data.
// -s <num> ... shift data a fraction <num> of the fundamental wavelength.
// 
// If file is not present, read from standard input.  Write to
// standard output.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: phase.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>
#include <data2df.h>


static char  prog[] = "phase";
static void  getargs  (int,char**,bool&,bool&,bool&,bool&,real_t&,istream*&);
static int_t getDump  (istream&,ostream&,vector<Data2DF*>&);
static void  loadName (const vector<Data2DF*>&,char*);
static int_t doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int_t            i;
  bool             conj = false, cmplx = false, zero = false, symm = false;
  real_t           alpha;
  istream*         input;
  vector<Data2DF*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, conj, cmplx, zero, symm, alpha, input);
  
  while (getDump (*input, cout, u))
    for (i = 0; i < u.size(); i++) {
      if (!cmplx|zero) u[i] -> DFT1D (FORWARD);

      // -- Here are the various possible actions:

      if (fabs (alpha) > EPSDP) u[i] -> F_shift (alpha, zero);
      if (conj)                 u[i] -> F_conjugate    (zero);
      if (symm)                 u[i] -> F_symmetrize   (zero);

      if (!cmplx|zero) u[i] -> DFT1D (INVERSE);
      cout << *u[i];
    }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     bool&     conj ,
		     bool&     cmplx,
		     bool&     zero ,
		     bool&     symm ,
		     real_t&   shift,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: phase [options] [file]\n"
  "options:\n"
  "-h       ... print this message\n"
  "-c       ... perform complex conjugation (reflection).\n"
  "-f       ... data are already complex (do not Fourier transform).\n"
  "-z       ... take mode zero as complex (e.g. file is an eigenmode).\n"
  "-r       ... enforce reflection symmetry of velocity & pressure data.\n"
  "-s <num> ... shift data a fraction <num> of the fundamental wavelength.\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'c':
      conj = true;
      break;
    case 'f':
      cmplx = true;
      break;
    case 'z':
      zero = true;
      break;
    case 'r':
      symm = true;
      break;
    case 's':
      if   (*++argv[0]) shift = atof (*argv);
      else { --argc;    shift = atof (*++argv); }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (fabs(shift) < EPSDP && !(conj | cmplx | symm))
    message (prog, "no action signalled", ERROR);

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;
}


static void loadName (const vector<Data2DF*>& u,
		      char*                    s)
// --------------------------------------------------------------------------
// Load a string containing the names of fields.
// ---------------------------------------------------------------------------
{
  int_t i, N = u.size();

  for (i = 0; i < N; i++) s[i] = u[i] -> getName();
  s[N] = '\0';
}


static int_t doSwap (const char* ffmt)
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


static int_t getDump (istream&           ifile,
		      ostream&           ofile,
		      vector<Data2DF*>& u    )
// ---------------------------------------------------------------------------
// Read next set of field dumps from ifile, put headers on ofile.
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
  char  buf[StrMax], fmt[StrMax], fields[StrMax];
  int_t i, j, swab, nf, np, nz, nel;

  if (ifile.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  ofile << buf << endl;
  ifile.getline (buf, StrMax);
  ofile << buf << endl;

  // -- I/O Numerical description of field sizes.

  ifile >> np >> nz >> nz >> nel;
  ifile.getline (buf, StrMax);
  
  sprintf (fmt, "%1d %1d %1d %1d", np, np, nz, nel);
  sprintf (buf, hdr_fmt[2], fmt);
  ofile << buf;

  for (i = 0; i < 5; i++) {
   ifile.getline (buf, StrMax);
   ofile << buf << endl;
  }

  // -- I/O field names.

  ifile >> fields;
  nf = strlen (fields);
  for (j = 0, i = 0; i < nf; i++) fmt[j++] = fields[i];
  fmt[j] = '\0';
  sprintf (buf, hdr_fmt[8], fmt);
  ofile << buf;
  ifile.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  ifile.getline (buf, StrMax);

  swab = doSwap (buf);

  sprintf (buf, "binary ");
  Veclib::describeFormat (buf + strlen (buf));
  sprintf (fmt, hdr_fmt[9], buf);
  ofile << fmt;

  if (u.size() != nf) {
    u.resize (nf);
    for (i = 0; i < nf; i++) u[i] = new Data2DF (np, nz, nel, fields[i]);
  }

  for (i = 0; i < nf; i++) {
    ifile >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return ifile.good();
}
