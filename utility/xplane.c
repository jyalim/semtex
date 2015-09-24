/*****************************************************************************
 * xplane.c: from a 3D data file, extract a 2D plane of data on named
 * plane and output 2D data file.  Field must be binary format.  By
 * default, output plane 1. Optionally (-2) output two consecutive
 * planes, i.e. a Fourier mode (presuming the input has first been
 * Fourier transformed); in this case the requested plane number is
 * the first of the two desired plane indices.  Plane indexing is 1-based.
 *
 * Copyright (c) 2000 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
 *
 * USAGE
 * -----
 * xplane [-h] [-n plane] [-2] [input[.fld]
 *

 *
 * --
 * This file is part of Semtex.
 * 
 * Semtex is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * Semtex is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Semtex (see the file COPYING); if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 *****************************************************************************/

static char RCS[] = "$Id: xplane.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, FILE**, int*, int*);

static char prog[] = "xplane";
static const char *hdr_fmt[] = { 
  "%-25s "              "Session\n",
  "%-25s "              "Created\n",
  "%-5d%-5d%-5d%-10d "  "Nr, Ns, Nz, Elements\n",
  "%-25d "              "Step\n",
  "%-25.6g "            "Time\n",
  "%-25.6g "            "Time step\n",
  "%-25.6g "            "Kinvis\n",
  "%-25.6g "            "Beta\n",
  "%-25s "              "Fields written\n",
  "%-25s "              "Format\n"
};


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Wrapper.
 *
 * NB: internal plane index is one less than set on command line.
 * ------------------------------------------------------------------------- */
{
  char   buf[STR_MAX], fields[STR_MAX], fmt[STR_MAX];
  int    i, j, n, np, nz, nel, pindex = 0, twoplanes = 0;
  int    nfields, nplane, npts, ntot, nreal;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double **data;

  getargs (argc, argv, &fp_in, &pindex, &twoplanes);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    
    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    if (pindex >= nz) {
      sprintf (fields, "requested plane (%1d) not present (nz = %1d)",
	       pindex + 1, nz);
      message (prog, fields, ERROR);
    }
    
    fprintf (fp_out, hdr_fmt[2], np, np, (twoplanes) ? 2 : 1, nel);

    n = 7;
    while (--n) {
      fgets (buf, STR_MAX, fp_in);
      fputs (buf, fp_out);   
    }

    for (nfields = 0, i = 0; i < 25; i++) if (isalnum(buf[i])) nfields++;

    fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out);

    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);

    /* -- Set sizes, allocate storage. */

    nplane = np * np * nel;
    nreal  = (twoplanes) ? nplane + nplane : nplane;
    npts   = nz * nplane;
    ntot   = nfields * npts;
    data   = dmatrix (0, nfields - 1, 0, npts - 1);
    
    /* -- Read in all data fields. */

    for (i = 0; i < nfields; i++) {
      for (j = 0; j < nz; j++) {
	if (fread(data[i]+j*nplane,sizeof(double),nplane,fp_in) != nplane)
	  message (prog, "an error occured while reading", ERROR);
      }
    }

    /* -- Write out nominated plane(s). Note planes are consecutive. */

    for (i = 0; i < nfields; i++)
      if (fwrite(data[i]+pindex*nplane,sizeof(double),nreal,fp_out) != nreal)
	message (prog, "an error occured while writing", ERROR);
    
    freeDmatrix (data,  0, 0);
  } 
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc     ,
		     char** argv     ,
		     FILE** fp_in    ,
		     int*   plane    ,
		     int*   twoplanes)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "xplane [-h] [-n plane] [-2] [input[.fld]\n";

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0]) *plane = atoi (*argv);
      else {*plane = atoi (*++argv); argc--;}
      (*plane)--;
      break;
    case '2':
      *twoplanes = 1;
      break;
    default:
      fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
      break;
    }
  
  if (argc == 1)
    if ((*fp_in = fopen(*argv, "r")) == (FILE*) NULL) {
      sprintf(fname, "%s.fld", *argv);
      if ((*fp_in = fopen(fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open input file -- %s or %s\n",
		prog, *argv, fname);
	exit (EXIT_FAILURE);
      }
    }
}
