///////////////////////////////////////////////////////////////////////////////
// misc.C: miscellaneous routines for I/O, memory management, service
// routines that don't fit class structures.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $ Hugh Blackburn
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
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: misc.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>
#include <data2df.h>


ostream& printVector (ostream&    strm,
		      const char* fmt , 
		      const int_t ntot,
		                  ... )
// ---------------------------------------------------------------------------
// Print up a variable number of numeric vectors on strm, in columns.
//
// The format specifier is gives the number and type of the vectors, with
// type specified by the first character.  The vectors must all be of the
// same type & length.
//
// Allowed types are: int ("i"), real ("r").
// Examples: four int_t vectors ==> fmt is "iiii".  Two reals ==> "rr".
// 
// Vectors are printed in a fixed field width of 15, regardless of type.
// ---------------------------------------------------------------------------
{
  char    routine[] = "printVector";
  int_t   nvect;
  va_list ap;

  nvect = strlen (fmt);
  if (! nvect   ) message (routine, "empty format string",   ERROR  );
  if (nvect > 10) message (routine, "more than 10 vectors?", WARNING);
  if (ntot  <  0) message (routine, "ntot < 0",              ERROR  );

  switch (fmt[0]) {

  case 'i': {
    int_t** u = new int_t* [nvect];
    va_start (ap, ntot);
    for (int_t k = 0; k < nvect; k++) u[k] = va_arg (ap, int_t*);
    va_end (ap);
    for (register int_t l = 0; l < ntot; l++) {
      for (register int_t j = 0; j < nvect; j++)
	strm << setw(15) << u[j][l];
      strm << endl;
    }
    delete [] u;
    break;
  }
  case 'r': {
    real_t** u = new real_t* [nvect];
    va_start (ap, ntot);
    for (int_t k = 0; k < nvect; k++) u[k] = va_arg (ap, real_t*);
    va_end (ap);
    for (register int_t l = 0; l < ntot; l++) {
      for (register int_t j = 0; j < nvect; j++)
	strm << setw(15) << u[j][l];
      strm << endl;
    }
    delete [] u;
    break;
  }
  default:
    message (routine, fmt, ERROR);
    break;
  }

  if (!strm) message (routine, "output failed", ERROR);

  return strm;
}


char* upperCase (char *s)
// ---------------------------------------------------------------------------
// Uppercase characters in string.
// ---------------------------------------------------------------------------
{
  char *z(s); while (*z = toupper (*z)) z++; return s;
}


void writeField (ostream&           file   ,
		 const char*        session,
		 const int_t        runstep,
		 const real_t       runtime,
		 vector<AuxField*>& field  )
// ---------------------------------------------------------------------------
// Write fields out to an opened file, binary nekton format.  Output is
// only done by the root processor.
// ---------------------------------------------------------------------------
{
  const char routine [] = "writeField";
  const char *hdr_fmt[] = { 
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

  char        s1[StrMax], s2[StrMax];
  time_t      tp (time (0));
  int_t       i;
  const int_t N = field.size();

  if (N < 1) return;

  ROOTONLY {
    sprintf (s1, hdr_fmt[0], session);
    file << s1;
    strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
    sprintf  (s1, hdr_fmt[1], s2);
    file << s1;

    field[0] -> describe (s2);
    sprintf (s1, hdr_fmt[2], s2);
    file << s1;

    sprintf (s1, hdr_fmt[3], runstep);
    file << s1;

    sprintf (s1, hdr_fmt[4], runtime);
    file << s1;

    sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
    file << s1;

    sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
    file << s1;

    sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
    file << s1;

    for (i = 0; i < N; i++) s2[i] = field[i] -> name();
    s2[i] = '\0';
    sprintf (s1, hdr_fmt[8], s2);
    file << s1;

    sprintf (s2, "binary ");
    Veclib::describeFormat (s2 + strlen (s2));
    sprintf (s1, hdr_fmt[9], s2);
    file << s1;
  }

  for (i = 0; i < N; i++) file << *field[i];

  ROOTONLY {
    if (!file) message (routine, "failed writing field file", ERROR);
    file << flush;
  }
}

void readField (istream&           file ,
                vector<AuxField*>& field)
// ---------------------------------------------------------------------------
// Read fields from an opened file, binary nekton format.
// ---------------------------------------------------------------------------
{
  const char routine [] = "readField";
  const int_t N = field.size();
  int_t i;

  if (N < 1) return;

  // -- read header, check sizes
  Header *hdr = new Header;
  file >> *hdr;

  ROOTONLY {
    if (hdr->nr != Geometry::nP() || hdr->ns != Geometry::nP())
      message(routine, "element size mismatch", ERROR);
    if (hdr->nz != Geometry::nZ())
      message (routine, "number of z planes mismatch", ERROR);
    if (hdr->nel != Geometry::nElmt())
      message (routine, "number of elements mismatch", ERROR);
  }

  // -- walk through fields, read appropriate
  char *type = hdr->flds;
  while (*type != 0)
  {
    bool skip = true;
    ROOTONLY cout << " type: " << *type;
    for (i = 0; i < N; i++) if (*type == field[i]->name())
    {
      file >> *field[i];
      ROOTONLY cout << "(reading)" << endl;
      skip = false;
    }
    if (skip) file.seekg (Geometry::nTot() * sizeof (real_t), ios::cur);
    type++;
  }
}
