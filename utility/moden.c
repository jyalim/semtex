/*****************************************************************************
 * moden.c: from a 3D data file, compute 2D distribution of kinetic
 * energy in named mode, output 2D data file.  Field must be binary
 * format.  By default, output energy in mode 0.
 *
 * Copyright (c) 1999 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
 *
 * USAGE
 * -----
 * moden [-h] [-m <mode>] [-z] [input[.fld]
 *
 * -m nominates mode to select [Default: 0]
 * -z forces mode zero to be dealt with as complex 
 *    (to be used e.g. with a complex eigenmode).  In which case N_Z=2.
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
 * 02110-1301 USA.
 *****************************************************************************/

static char RCS[] = "$Id: moden.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, FILE**, int*, int*);
static int  _index  (const char*, char);

static char prog[] = "moden";
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
  char   buf[STR_MAX], fields[STR_MAX], fmt[STR_MAX];
  int    i, j, n, np, nz, nel, mode = 0, swab = 0, cmplx = 0;
  int    nfields, nplane, nplaneEven, nptsEven, ntot;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double **data, *plane, *vcmpt;

  getargs (argc, argv, &fp_in, &mode, &cmplx);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    
    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    if (2 * mode > nz) {
      sprintf (fields, "too many modes (%1d) for input (nz = %1d)", mode, nz);
      message (prog, fields, ERROR);
    }
    if (cmplx && nz != 2) {
      sprintf (fields, "need nz = 2 with full-complex single mode (%1d)", nz);
      message (prog, fields, ERROR);
    } 

    fprintf (fp_out, hdr_fmt[2], np, np, 1, nel);

    n = 6;
    while (--n) {
      fgets (buf, STR_MAX, fp_in);
      fputs (buf, fp_out);   
    }

    fgets  (fields, STR_MAX, fp_in);
    memset (fields+25, '\0', STR_MAX-25);
    for (nfields = 0, i = 0; i < 25; i++) if (isalpha(fields[i])) nfields++;
    if (nfields < 4) {
      if (!(strchr(fields, 'u') && strchr(fields, 'v')))
	message (prog, "need fields u, v to compute K.E.", ERROR);
    } else {
      if (!(strchr(fields, 'u') && strchr(fields, 'v') && strchr(fields, 'w')))
	message (prog, "need fields u, v, w to compute K.E.", ERROR);
    }
    fprintf (fp_out, hdr_fmt[8], "q");

    fgets (buf, STR_MAX, fp_in);
    for (i = 0; i < strlen (buf); i++) buf[i] = tolower (buf[i]);

    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);
    if (!strstr (buf, "endian"))
      message (prog, "input field file in unknown binary format", WARNING);
    else
      swab = ((strstr (buf, "big") && strstr (fmt, "little")) ||
	      (strstr (fmt, "big") && strstr (buf, "little")) );
    sprintf (buf, "%s %s", "binary", fmt);
    fprintf (fp_out, hdr_fmt[9], buf);

    /* -- Set sizes, allocate storage. */

    nplane     = np * np * nel;
    nplaneEven = (nplane & 1) ? nplane + 1 : nplane;
    nptsEven   = nz * nplaneEven;
    ntot       = nfields * nptsEven;

    data  = dmatrix (0, nfields - 1, 0, nptsEven - 1);
    plane = dvector (0, nplane  - 1);
    
    /* -- Read in all data fields. */

    dzero (ntot, data[0], 1);
    dzero (nplane, plane, 1);

    for (i = 0; i < nfields; i++) {
      for (j = 0; j < nz; j++) {
	if (fread (data[i] + j*nplaneEven, sizeof (double), nplane, fp_in)
	    != nplane)
	  message (prog, "an error occured while reading", ERROR);
      }
      if (swab)   dbrev (nptsEven, data[i], 1, data[i], 1);
      if (!cmplx) dDFTr (data[i], nz, nplaneEven, +1);
    }
    
    /* -- Compute K.E.: start by adding in real part. */

    for (i = 0; i < nfields - 1; i++) {
      vcmpt = data[_index (fields, 'u' + i)] + 2 * mode * nplaneEven;
      dvvtvp (nplane, vcmpt, 1, vcmpt, 1, plane, 1, plane, 1);
    }

    /* -- Add in imaginary part if not mode zero. */

    if (mode || cmplx) {
      for (i = 0; i < nfields - 1; i++) {
	vcmpt = data[_index (fields, 'u' + i)] + (2 * mode + 1) * nplaneEven;
	dvvtvp (nplane, vcmpt, 1, vcmpt, 1, plane, 1, plane, 1);
      }
    }

    /*  -- Normalize to make q = 0.5*UiUi. */

    dsmul (nplane, 0.5, plane, 1, plane, 1);
    
    if (fwrite (plane, sizeof (double), nplane, fp_out) != nplane)
      message (prog, "an error occured while writing", ERROR);

    freeDmatrix (data,  0, 0);
    freeDvector (plane, 0);
  } 
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     FILE** fp_in,
		     int*   mode ,
		     int*   cmplx)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "moden [-h] [-m mode] [-z] [input[.fld]\n";

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'm':
      if (*++argv[0]) *mode = atoi (*argv);
      else { *mode = atoi (*++argv); argc--; }
      break;
    case 'z':
      *cmplx = 1;
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


static int _index (const char* s,
		   char        c)
/* ------------------------------------------------------------------------- *
 * Return index of c in s, -1 if not found.
 * ------------------------------------------------------------------------- */
{
  int       i;
  const int len = strlen (s);

  for (i = 0; i < len; i++) if (s[i] == c) return i;

  return -1;
}
