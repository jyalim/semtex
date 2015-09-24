///////////////////////////////////////////////////////////////////////////////
// compare.C
//
// Copyright (c) 1996 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// SYNOPSIS
// --------
// Compute exact solution given in USER section of FEML file, subtract
// numerical solution (if present), print up infinity norm (largest
// error) and write field file of error field.  If no numerical
// solution is given, output is exact solution.  Only designed for
// serial execution.  When evaluating the exact solution, use the time
// value found in the field file.
//
// USAGE
// -----
// compare [options] session [field.file]
// options:
// -h ... print this message
// -n ... print "noise" if result is less than noise level
// -f ... forward Fourier transform output
//
// EXAMPLE
// -------
// Kovasznay flow in x--y plane.
// For a 3D solution, the USER section would read:
//
// <USER>
//   u = 1-exp(LAMBDA*x)*cos(2*PI*y)
//   v = LAMBDA/(2*PI)*exp(LAMBDA*x)*sin(2*PI*y)
//   w = 0.0
//   p = (1.0-exp(lambda*x))/2.0
// </USER>
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

static char RCS[] = "$Id: compare.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <ctime>
#include <sem.h>

static char  prog[]    = "compare";
const  int_t EXACT_MAX = 32;

static void getargs (int, char**, bool&, bool&, char*&, ifstream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Everything is in here, no file-scope subroutines.
// ---------------------------------------------------------------------------
{
  char               *session, *tok;
  ifstream           fieldfl;
  char               buf[StrMax], err[StrMax], fmt[StrMax], fields[StrMax];
  char               function[EXACT_MAX][StrMax];
  int_t              i, j, np, nz, nel, found;
  int_t              nexact = 0, nfields = 0;
  bool               swab = false, tran = false, noise = false;
  real_t             t;
  Geometry::CoordSys system;
  vector<Element*>   Esys;
  AuxField           *exact, *computed;
  const real_t       NOISE = 1e-12;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, tran, noise, session, fieldfl);

  Veclib::describeFormat (fmt);

  // -- Set up from session file.

  FEML* F = new FEML (session);
  Mesh* M = new Mesh (F);

  nel    =  M -> nEl();  
  np     =  Femlib::ivalue ("N_P");
  nz     =  Femlib::ivalue ("N_Z");
  system = (Femlib::ivalue ("CYLINDRICAL"))
    ? Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (np, nz, nel, system);

  if (F -> seek ("USER")) {	// -- Search for lines of form "char = string".

    F -> stream().ignore (StrMax, '\n');

    while (F -> stream().getline (buf, StrMax)) {
      if (strstr (buf, "=")) {
	strcpy (function[nexact], buf);
	nexact++;
      }
      tok = buf; while (*tok = toupper (*tok)) tok++;
      if (strstr (buf, "USER")) break;
    }
  }

  if (fieldfl) {

    // -- A field file exists.

    fieldfl.getline (buf, StrMax);
    if (!strstr (buf, "Session")) {
      sprintf (err, "%s is not a field file", argv[2]);
      message (prog, err, ERROR);
    }
    
    fieldfl .getline(buf, StrMax) .getline(buf, StrMax);

    // -- Check it conforms with the description in session file.

    string ss (buf);
    istringstream sss (ss);
    sss  >> np >> np >> nz >> nel;

    if (np != Geometry::nP()) {
      sprintf (err, "polynomial order mismatch (%1d <--> %1d)",
	       np, Geometry::nP());
      message (prog, err, ERROR);
    }
    if (nz != Geometry::nZ()) {
      sprintf (err, "number of z-planes mismatch (%1d <--> %1d)",
	       nz, Geometry::nZ());
      message (prog, err, ERROR);
    }
    if (nel != M -> nEl()) {
      sprintf (err, "number of elements mismatch (%1d <--> %1d)",
	       nel, Geometry::nElmt());
      message (prog, err, ERROR);
    }

    // -- Find the fields it contains, and set the time to found value.

    fieldfl .getline(buf, StrMax) .getline(buf, StrMax);

    sss.str (ss = buf);
    sss >> t;

    Femlib::value ("t", t);
    
    fieldfl .getline(buf, StrMax) .getline(buf, StrMax);
    fieldfl .getline(buf, StrMax) .getline(buf, StrMax);

    while (isalnum (buf[nfields])) {
      fields[nfields] = buf[nfields];
      nfields++;
    }
    fields[nfields] = '\0';
     
    // -- Is byte-swapping required?
    
    fieldfl.getline (buf, StrMax);
    if (!strstr (buf, "binary"))
      message (prog, "input field file not in binary format", ERROR);
    else if (!strstr (buf, "endian"))
      message (prog, "input field file in unknown binary format", WARNING);
    else {
      swab = (   (strstr (buf, "big") && strstr (fmt, "little"))
	      || (strstr (fmt, "big") && strstr (buf, "little")) );
    }
  }

  // -- Build Element information, then comparison AuxFields.

  Esys.resize (nel);
  for (i = 0; i < nel; i++) Esys[i] = new Element (i, np, M);

  exact    = new AuxField (new real_t[Geometry::nTotal()], nz, Esys);
  computed = new AuxField (new real_t[Geometry::nTotal()], nz, Esys);

  // -- Perform comparisons, output.

  if (nfields) {		// -- Have a set of data to check.

    fieldfl.seekg (0);
    for (i = 0; i < 10; i++) {	// -- Copy header information.
      fieldfl.getline (buf, StrMax);
      cout << buf << endl;
    }

    for (i = 0; i < nfields; i++) {

      fieldfl >> *computed;
      if (!fieldfl) {
	sprintf (err, "error reading input Field `%c'", fields[i]);
	message (prog, err, ERROR);
      }
      if (swab) computed -> reverse();

      for (found = 0, j = 0; j < nexact; j++) {
	tok = function[j];
	while (!isalpha (*tok)) tok++;
	if (found = tok[0] == fields[i]) {
	  tok    = strtok (function[j], "=");
	  tok    = strtok (0, "\0");
	  *exact = (const char*) tok;
	  break;
	}
      }

      if (!found) {
	sprintf (err, "no function description for field '%c'", fields[i]);
	message (prog, err, WARNING);
	*exact = 0.0;
      }

      *exact -= *computed;
      t = exact -> norm_inf();

      cerr.setf (ios::scientific, ios::floatfield);
      cerr << setprecision (4);
      cerr << "Field '" << fields[i] << "': norm_inf: ";
      if   (t < NOISE && noise ) cerr << "noise-level" << endl;
      else                       cerr << t             << endl;

      if (tran) exact -> transform (FORWARD);
      if (swab) exact -> reverse();
      cout << *exact;
    }

  } else {			// -- Create new data.

    const char *hdr_fmt[] = { 
      "%-25s "             "Session\n",
      "%-25s "             "Created\n",
      "%-5d%-5d%-5d%-10d " "Nr, Ns, Nz, Elements\n",
      "%-25d "             "Step\n",
      "%-25.6g "           "Time\n",
      "%-25.6g "           "Time step\n",
      "%-25.6g "           "Kinvis\n",
      "%-25.6g "           "Beta\n",
      "%-25s "             "Fields written\n",
      "%-25s "             "Format\n"
    };

    char   s1[StrMax], s2[StrMax];
    time_t tp (time (0));

    // -- Write the header.

    sprintf (s1, hdr_fmt[0], session);
    cout << s1;

    strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
    sprintf  (s1, hdr_fmt[1], s2);
    cout << s1;

    sprintf (s1, hdr_fmt[2], np, np, nz, nel);
    cout << s1;

    sprintf (s1, hdr_fmt[3], 0);
    cout << s1;
    
    sprintf (s1, hdr_fmt[4], Femlib::value ("t"));
    cout << s1;

    sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
    cout << s1;

    sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
    cout << s1;

    sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
    cout << s1;

    for (j = 0; j < nexact; j++) {
      tok = function[j];
      while (!isalpha (*tok)) tok++;
      fields[j] = *tok;
    }
    fields[j] = '\0';
    sprintf (s1, hdr_fmt[8], fields);
    cout << s1;

    sprintf (s2, "binary ");
    Veclib::describeFormat (s2 + strlen (s2));
    sprintf (s1, hdr_fmt[9], s2);
    cout << s1;
    
    // -- Compute and output fields.

    for (j = 0; j < nexact; j++) {
      tok    = strtok (function[j], "=");
      tok    = strtok (0, "\0");
      *exact = (const char*) tok;
      if (tran) exact -> transform (FORWARD);
      cout << *exact;
    }
  }
  
  cout.flush();
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     bool&     tran ,
		     bool&     noise,
		     char*&    sess ,
		     ifstream& fldf )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: compare [options] session [field.file]\n"
                 "options:\n"
                 "  -h ... display this message\n"
                 "  -n ... print 'noise-level' for small errors\n"
                 "  -t ... forward Fourier transform output\n";
  char err[StrMax], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h': cerr << usage; exit (EXIT_SUCCESS); break;
    case 'n': noise = true; break;
    case 't': tran  = true; break;
    default:
      sprintf (err, "illegal option: %c\n", c);
      message (prog, err, ERROR); break;
    }

  switch (argc) {
  case 1:
    sess = argv[0];
    fldf.clear (ios::badbit | fldf.rdstate());
    break;
  case 2:
    sess = argv[0];
    fldf.open (argv[1], ios::in);
    if (!fldf) message (prog, "couldn't open field file", ERROR);
    break;
  default:
    cerr << usage;
    exit (EXIT_FAILURE);
    break;
  }  
}
