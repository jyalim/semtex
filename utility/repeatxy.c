/*****************************************************************************
 * repeatxy.c: read a field file, and output the number of repetitions
 * in x-y planes indicated on the command line.  Field must be
 * binary format.  By default, output the original file, no repeats.
 * Adust Nel in header as appropriate.
 *
 * This is the planar equivalent of repeatz.c
 *
 * Copyright (c) 2002 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn.
 *
 * USAGE
 * -----
 * repeatxy [-h] [-n <rep>] [input[.fld]
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

static char RCS[] = "$Id: repeatxy.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, FILE**, int*);
static char prog[] = "repeatxy";
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
 * ------------------------------------------------------------------------- */
{
  char   buf[STR_MAX], fmt[STR_MAX];
  int    i, j, k, np, nz, nel, nrep = 1;
  int    nfields, nplane, nptin, nptout, ntot, swab;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double *data;

  getargs (argc, argv, &fp_in, &nrep);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    
    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    fprintf (fp_out, hdr_fmt[2], np, np, nz, nel*nrep);

    i = 5;
    while (i--) { fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out); }

    fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out);
    for (nfields = 0, i = 0; i < 25; i++) if (isalnum(buf[i])) nfields++;

    fgets (buf, STR_MAX, fp_in);
    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);
    swab = (strstr (buf, "big") && strstr (fmt, "little")) || 
           (strstr (fmt, "big") && strstr (buf, "little"));
    strcat (strcpy (buf, "binary "), fmt);
    fprintf (fp_out, hdr_fmt[9], buf);
    
    /* -- Set sizes, allocate storage. */

    nplane = np * np * nel;
    ntot   = nz * nplane;
    data   = dvector (0, ntot - 1);
    
    /* -- Read and write all data fields. */

    for (i = 0; i < nfields; i++) {
      if (fread (data, sizeof (double), ntot, fp_in) != ntot)
	message (prog, "an error occured while reading", ERROR);
      if (swab) dbrev (ntot, data, 1, data, 1);

      for (j = 0; j < nz; j++)
	for (k = 0; k < nrep; k++)
	  if (fwrite (data+j*nplane, sizeof(double), nplane, fp_out) != nplane)
	    message (prog, "an error occured while writing", ERROR);
    }
  }

  freeDvector (data, 0);
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     FILE** fp_in,
		     int*   nrep )
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "repeatxy [-h] [-n <rep>] [input[.fld]\n";

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0]) *nrep = atoi (*argv);
      else           {*nrep = atoi (*++argv); argc--;}
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
